#include <iostream>
#include <array>
#include "AMRSim.h"

void AmrSim::SwapElements(double * const arr, const int offset_a,
const int offset_b) {
  double tmp = arr[offset_a];
  arr[offset_a] = arr[offset_b];
  arr[offset_b] = tmp;
  return;
}

void AmrSim::UpdateBoundaries(int level) {
  dist_fn.at(level).FillBoundary(geom[level].periodicity());
  return;
}

void AmrSim::Collide(int level) {
  amrex::IntVect pos(0);
  std::array<double, NMODES> mode;
  double usq, TrS;
  double stress[NMODES][NMODES];
  int a, b, m, p;

  for (amrex::MFIter mfi(dist_fn.at(level)); mfi.isValid(); ++mfi) {
    amrex::FArrayBox& fab_dist_fn = dist_fn.at(level)[mfi];
    amrex::FArrayBox& fab_density = density.at(level)[mfi];
    amrex::FArrayBox& fab_velocity = velocity.at(level)[mfi];

    // will need to replace NX/NY/NZ with local domain sizes
    for (int k = 0; k < NZ; ++k) {
      pos.setVal(2, k);
      for (int j = 0; j < NY; ++j) {
        pos.setVal(1, j);
        for (int i = 0; i < NX; ++i) {
          pos.setVal(0, i);

          for (m = 0; m < NMODES; ++m) {
            mode[m] = 0.0;
            for (p = 0; p < NMODES; ++p) {
              mode[m] += fab_dist_fn(pos,p) * MODE_MATRIX[m][p];
            }
          }

          fab_density(pos) = mode[0];

          // no forcing is currently present in the model,
          // so we disregard uDOTf for now
          usq = 0.0;
          for (a = 0; a < NDIMS; ++a) {
            fab_velocity(pos, a) = mode[a+1] / fab_density(pos);
            usq += fab_velocity(pos, a) * fab_velocity(pos, a);
          }

          stress[0][0] = mode[4];
          stress[0][1] = mode[5];
          stress[0][2] = mode[6];

          stress[1][0] = mode[5];
          stress[1][1] = mode[7];
          stress[1][2] = mode[8];

          stress[2][0] = mode[6];
          stress[2][1] = mode[8];
          stress[2][2] = mode[9];

          // Form the trace
          TrS = 0.0;
          for (a = 0; a < NDIMS; ++a) {
            TrS += stress[a][a];
          }
          // Form the traceless part
          for (a = 0; a < NDIMS; ++a) {
            stress[a][a] -= (TrS / NDIMS);
          }

          // Relax the trace
          TrS -= OMEGA_B * (TrS - fab_density(pos)*usq);
          // Relax the traceless part
          for (a = 0; a < NDIMS; ++a) {
            for (b = 0; b < NDIMS; ++b) {
              stress[a][b] -= OMEGA_S * (stress[a][b] - fab_density(pos)
                                      * ( fab_velocity(pos,a)
                                          * fab_velocity(pos,b)
                                          - usq * DELTA[a][b]) );
            }
            stress[a][a] += (TrS / NDIMS);
          }

          // copy stress back into mode
          mode[4] = stress[0][0];
          mode[5] = stress[0][1];
          mode[6] = stress[0][2];

          mode[7] = stress[1][1];
          mode[8] = stress[1][2];

          mode[9] = stress[2][2];

          // Ghosts are relaxed to zero immediately
          mode[10] = 0.0;
          mode[11] = 0.0;
          mode[12] = 0.0;
          mode[13] = 0.0;
          mode[14] = 0.0;

          // project back to the velocity basis
          for (p = 0; p < NMODES; ++p) {
            fab_dist_fn(pos, p) = 0.0;
            for (m = 0; m < NMODES; ++m) {
              fab_dist_fn(pos, p) += mode[m] * MODE_MATRIX_INVERSE[p][m];
            }
          }

        } // i
      } // j
    } // k
  } // MFIter

  return;
}

void AmrSim::Propagate(int level) {
  // Explanation from subgrid source
  /* The propagation step of the algorithm.
  *
  * Copies (f[x][y][z][i] + R[i]) into f[x+dx][y+dy][z+dz][i]
  *
  * For each site, only want to access cells which are further
  * along in memory, so use: (expecting 7 directions, as
  * zero velocity doesn't move and we're swapping)
  * [1,0,0] [0,1,0] [0,0,1] i.e. 1, 3, 5
  * [1,1,1] [1,1,-1] [1,-1,1] [1,-1,-1] i.e. 7,8,9,10
  *
  * with, respectively:
  * 2,4,6
  * 14, 13, 12, 11
  *
  * We swap the value with the opposite direction velocity in
  * the target Lattice site. By starting in the halo, we ensure
  * all the real cells get fully updated. Note that the first
  * row must be treated a little differently, as it has no
  * neighbour in the [1,-1] position.
  * After this is complete, we have the distribution functions
  * updated, but with the fluid propagating AGAINST the velocity
  * vector. So need to reorder the values once we've done the
  * propagation step.
  */
  double * data;
  amrex::IntVect dims(0);

  for (amrex::MFIter mfi(dist_fn.at(level)); mfi.isValid(); ++mfi) {
    data = dist_fn.at(level)[mfi].dataPtr();
    dims = dist_fn.at(level)[mfi].length();

    // this looping is bad for FAB, but will need to reorder swaps too...
    for (int i = 0; i < dims[0]; ++i) {
      for (int j = 0; j < dims[1]; ++j) {
        for (int k = 0; k < dims[2]; ++k) {

          // [1,0,0]
          if (i < dims[0]-1)
          SwapElements(data, Lindex(i,j,k,1,dims), Lindex(i+1,j,k,2,dims));
          // [0,1,0]
          if (j < dims[1]-1)
          SwapElements(data, Lindex(i,j,k,3,dims), Lindex(i,j+1,k,4,dims));
          // [0,0,1]
          if (k < dims[2]-1)
          SwapElements(data, Lindex(i,j,k,5,dims), Lindex(i,j,k+1,6,dims));

          // [1,1,1]
          if (i < dims[0]-1 && j < dims[1]-1 && k < dims[2]-1)
          SwapElements(data, Lindex(i,j,k,7,dims), Lindex(i+1,j+1,k+1,14,dims));
          // [1,1,-1]
          if (i < dims[0]-1 && j < dims[1]-1 && k > 0)
          SwapElements(data, Lindex(i,j,k,8,dims), Lindex(i+1,j+1,k-1,13,dims));
          // [1,-1,1]
          if (i < dims[0]-1 && j > 0 && k < dims[2]-1)
          SwapElements(data, Lindex(i,j,k,9,dims), Lindex(i+1,j-1,k+1,12,dims));
          // [1,-1,-1]
          if (i < dims[0]-1 && j > 0 && k > 0)
          SwapElements(data, Lindex(i,j,k,10,dims), Lindex(i+1,j-1,k-1,11,dims));

          // reorder
          SwapElements(data, Lindex(i,j,k,1,dims), Lindex(i,j,k,2,dims));
          SwapElements(data, Lindex(i,j,k,3,dims), Lindex(i,j,k,4,dims));
          SwapElements(data, Lindex(i,j,k,5,dims), Lindex(i,j,k,6,dims));

          SwapElements(data, Lindex(i,j,k,7,dims), Lindex(i,j,k,14,dims));
          SwapElements(data, Lindex(i,j,k,8,dims), Lindex(i,j,k,13,dims));
          SwapElements(data, Lindex(i,j,k,9,dims), Lindex(i,j,k,12,dims));
          SwapElements(data, Lindex(i,j,k,10,dims), Lindex(i,j,k,11,dims));
        }
      }
    }

  }

  return;
}

void AmrSim::ErrorEst(int level, amrex::TagBoxArray& tba, double time, int ngrow) { return; }

void AmrSim::MakeNewLevelFromScratch(int level, double time, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) {
  // define MultiFabs
  density[level].define(ba, dm, 1, 0);
  velocity[level].define(ba, dm, NDIMS, 0);
  dist_fn[level].define(ba, dm, NMODES, HALO_DEPTH);

  return;
}

void AmrSim::MakeNewLevelFromCoarse(int level, double time, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) { return; }

void AmrSim::RemakeLevel(int level, double time, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) { return; }

void AmrSim::ClearLevel(int level) { return; }

AmrSim::AmrSim(int const nx, int const ny, int const nz,
double const tau_s, double const tau_b, int (&periodicity)[NDIMS],
amrex::RealBox& domain)
: AmrCore(&domain, 0, amrex::Vector<int>({AMREX_D_DECL(8,8,8)}), 0),
NX(nx), NY(ny), NZ(nz), NUMEL(nx*ny*nz), COORD_SYS(0),
PERIODICITY{ periodicity[0], periodicity[1], periodicity[2] },
TAU_S(tau_s), TAU_B(tau_b), OMEGA_S(1.0/(tau_s+0.5)),
OMEGA_B(1.0/(tau_b+0.5)), time_step(1, 0.0),
idx_domain( amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
            amrex::IntVect(AMREX_D_DECL(nx-1, ny-1, nz-1)),
            amrex::IndexType( {AMREX_D_DECL(0, 0, 0)} ) ),
phys_domain( domain ), ba_domain(idx_domain), dm(ba_domain)
{
  // As we don't use ParmParse the base class is initialised with dummy inputs.
  // We set the max grid size, blocking factor distribution map, BoxArray, and
  // geometry here...
  verbose = 1;
  SetMaxGridSize(amrex::IntVect(AMREX_D_DECL(nx-1,ny-1,nz-1)));
  SetBlockingFactor(amrex::IntVect(AMREX_D_DECL(nx-1,ny-1,nz-1)));
  SetDistributionMap(0, dm);
  SetBoxArray(0, ba_domain);
  geom.clear();
  // This is done to reset static variables in Geometry
  amrex::Geometry::Finalize();
  // note: there is a new constructor with RealBox& available in 19.03!
  geom.emplace_back(ba_domain[0], &phys_domain, COORD_SYS, PERIODICITY);

  // resize MF vectors
  int num_levels = max_level + 1;
  density.resize(num_levels);
  velocity.resize(num_levels);
  dist_fn.resize(num_levels);
};

void AmrSim::SetInitialDensity(const double rho_init) {
  initial_density.assign(NUMEL, rho_init);
  return;
}

void AmrSim::SetInitialDensity(const std::vector<double> rho_init) {
  initial_density = rho_init;
  return;
}

void AmrSim::SetInitialVelocity(const double u_init) {
  initial_velocity.assign(NUMEL, u_init);
  return;
}

void AmrSim::SetInitialVelocity(const std::vector<double> u_init) {
  initial_velocity = u_init;
  return;
}

void AmrSim::InitDistFunc() {

  return;
}

double AmrSim::GetDensity(const int i, const int j, const int k,
const int level) const {
  amrex::IntVect pos(i,j,k);
  for (amrex::MFIter mfi(density.at(level)); mfi.isValid(); ++mfi) {
    if (density.at(level)[mfi].box().contains(pos))
      return density.at(level)[mfi](pos);
  }
  amrex::Abort("Invalid index to density MultiFab.");
  return -1.0;
}

// static member definitions
const double AmrSim::DELTA[NDIMS][NDIMS] = { {1.0 / NMODES, 0.0, 0.0},
                                           {0.0, 1.0 / NMODES, 0.0},
                                           {0.0, 0.0, 1.0 / NMODES} };

const double AmrSim::MODE_MATRIX[NMODES][NMODES] =
{
	{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
	{0, 1, -1, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1},
	{0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1},
	{0, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1},
	{-1./3., 2./3., 2./3., -1./3., -1./3., -1./3., -1./3., 2./3., 2./3., 2./3., 2./3., 2./3., 2./3., 2./3., 2./3.},
	{0, 0, 0, 0, 0, 0, 0, 1./1., 1./1., -1./1., -1./1., -1./1., -1./1., 1./1., 1./1.},
	{0, 0, 0, 0, 0, 0, 0, 1./1., -1./1., 1./1., -1./1., -1./1., 1./1., -1./1., 1./1.},
	{-1./3., -1./3., -1./3., 2./3., 2./3., -1./3., -1./3., 2./3., 2./3., 2./3., 2./3., 2./3., 2./3., 2./3., 2./3.},
	{0, 0, 0, 0, 0, 0, 0, 1./1., -1./1., -1./1., 1./1., 1./1., -1./1., -1./1., 1./1.},
	{-1./3., -1./3., -1./3., -1./3., -1./3., 2./3., 2./3., 2./3., 2./3., 2./3., 2./3., 2./3., 2./3., 2./3., 2./3.},
	{-2, 1, 1, 1, 1, 1, 1, -2, -2, -2, -2, -2, -2, -2, -2},
	{0, 1, -1, 0, 0, 0, 0, -2, -2, -2, -2, 2, 2, 2, 2},
	{0, 0, 0, 1, -1, 0, 0, -2, -2, 2, 2, -2, -2, 2, 2},
	{0, 0, 0, 0, 0, 1, -1, -2, 2, -2, 2, -2, 2, -2, 2},
	{0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1, -1, 1, 1, -1}
};

const double AmrSim::MODE_MATRIX_INVERSE[NMODES][NMODES] =
{
	{2./9., 0, 0, 0, -1./3., 0, 0, -1./3., 0, -1./3., -2./9., 0, 0, 0, 0},
	{1./9., 1./3., 0, 0, 1./3., 0, 0, -1./6., 0, -1./6., 1./18., 1./6., 0, 0, 0},
	{1./9., -1./3., 0, 0, 1./3., 0, 0, -1./6., 0, -1./6., 1./18., -1./6., 0, 0, 0},
	{1./9., 0, 1./3., 0, -1./6., 0, 0, 1./3., 0, -1./6., 1./18., 0, 1./6., 0, 0},
	{1./9., 0, -1./3., 0, -1./6., 0, 0, 1./3., 0, -1./6., 1./18., 0, -1./6., 0, 0},
	{1./9., 0, 0, 1./3., -1./6., 0, 0, -1./6., 0, 1./3., 1./18., 0, 0, 1./6., 0},
	{1./9., 0, 0, -1./3., -1./6., 0, 0, -1./6., 0, 1./3., 1./18., 0, 0, -1./6., 0},
	{1./72., 1./24., 1./24., 1./24., 1./24., 1./8., 1./8., 1./24., 1./8., 1./24., -1./72., -1./24., -1./24., -1./24., 1./8.},
	{1./72., 1./24., 1./24., -1./24., 1./24., 1./8., -1./8., 1./24., -1./8., 1./24., -1./72., -1./24., -1./24., 1./24., -1./8.},
	{1./72., 1./24., -1./24., 1./24., 1./24., -1./8., 1./8., 1./24., -1./8., 1./24., -1./72., -1./24., 1./24., -1./24., -1./8.},
	{1./72., 1./24., -1./24., -1./24., 1./24., -1./8., -1./8., 1./24., 1./8., 1./24., -1./72., -1./24., 1./24., 1./24., 1./8.},
	{1./72., -1./24., 1./24., 1./24., 1./24., -1./8., -1./8., 1./24., 1./8., 1./24., -1./72., 1./24., -1./24., -1./24., -1./8.},
	{1./72., -1./24., 1./24., -1./24., 1./24., -1./8., 1./8., 1./24., -1./8., 1./24., -1./72., 1./24., -1./24., 1./24., 1./8.},
	{1./72., -1./24., -1./24., 1./24., 1./24., 1./8., -1./8., 1./24., -1./8., 1./24., -1./72., 1./24., 1./24., -1./24., 1./8.},
	{1./72., -1./24., -1./24., -1./24., 1./24., 1./8., 1./8., 1./24., 1./8., 1./24., -1./72., 1./24., 1./24., 1./24., -1./8.}
};
