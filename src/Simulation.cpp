#include <cstring>
#include <iostream>
#include <array>
#include "Simulation.h"

void Simulation::collide() {
  amrex::IntVect pos(0);
  std::array<double, NMODES> mode;
  double usq, TrS;
  double stress[NMODES][NMODES];
  int a, b, m, p;

  for (amrex::MFIter mfi(dist_fn); mfi.isValid(); ++mfi) {
    amrex::FArrayBox& fab_dist_fn = dist_fn[mfi];

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

          density(pos) = mode[0];

          // no forcing is currently present in the model,
          // so we disregard uDOTf for now
          usq = 0.0;
          for (a = 0; a < NDIMS; ++a) {
            velocity(pos, a) = mode[a+1] / density(pos);
            usq += velocity(pos, a) * velocity(pos, a);
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
          TrS -= OMEGA_B * (TrS - density(pos)*usq);
          // Relax the traceless part
          for (a = 0; a < NDIMS; ++a) {
            for (b = 0; b < NDIMS; ++b) {
              stress[a][b] -= OMEGA_S * (stress[a][b] - density(pos)
                                      * ( velocity(pos,a) * velocity(pos,b)
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
              fab_dist_fn(pos, p) = mode[m] * MODE_MATRIX_INVERSE[p][m];
            }
          }

        } // i
      } // j
    } // k
  } // MFIter

  return;
}

void Simulation::propagate() {
  return;
}

Simulation::Simulation(int const nx, int const ny, int const nz,
  double const tau_s, double const tau_b, int (&periodicity)[NDIMS])
: NX(nx), NY(ny), NZ(nz), NUMEL(nx*ny*nz), COORD_SYS(0),
  PERIODICITY{ periodicity[0], periodicity[1], periodicity[2] },
  TAU_S(tau_s), TAU_B(tau_b), OMEGA_S(1.0/(tau_s+0.5)),
  OMEGA_B(1.0/(tau_b+0.5)),
  idx_domain( amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
              amrex::IntVect(AMREX_D_DECL(nx-1, ny-1, nz-1)),
              amrex::IndexType( {AMREX_D_DECL(0, 0, 0)} ) ),
  phys_domain( {AMREX_D_DECL(0.0, 0.0, 0.0)},
               {AMREX_D_DECL( (double) nx-1, (double) ny-1, (double) nz-1 )} ),
  geometry(idx_domain, &phys_domain, COORD_SYS, periodicity),
  ba_domain(idx_domain), dm(ba_domain), density(idx_domain),
  velocity(idx_domain, NDIMS), force(idx_domain, NDIMS),
  dist_fn(ba_domain, dm, NMODES, HALO_DEPTH) {};

void Simulation::setDensity(double const uniform_density) {
  // set density array to one value everywhere
  density.setVal(uniform_density);
  return;
}

void Simulation::setDensity(double const * const rho) {
  // could this be done by getting a position from a boxiter?
  amrex::IntVect pos(0);
  // set density to match provided array
  for (int i = 0; i < NX; ++i) {
    pos.setVal(0, i);
    for (int j = 0; j < NY; ++j) {
      pos.setVal(1, j);
      for (int k = 0; k < NZ; ++k) {
        pos.setVal(2, k);
        density(pos) = rho[i*NY*NZ + j*NZ + k];
      }
    }
  }
  return;
}

void Simulation::setVelocity(double const uniform_velocity) {
  // set velocity to one value everywhere
  velocity.setVal(uniform_velocity);
  return;
}

void Simulation::setVelocity(double const * const u) {
  amrex::IntVect pos(0);
  // set velocity to match provided array
  for (int i = 0; i < NX; ++i) {
    pos.setVal(0, i);
    for (int j = 0; j < NY; ++j) {
      pos.setVal(1, j);
      for (int k = 0; k < NZ; ++k) {
        pos.setVal(2, k);
        for (int n = 0; n < NDIMS; ++n) {
          velocity(pos, n) = u[i*NY*NZ*NDIMS + j*NZ*NDIMS + k*NDIMS + n];
        }
      }
    }
  }
  return;
}

void Simulation::calcEquilibriumDist() {
  double u2[NDIMS], mod_sq, u_cs2[NDIMS], u2_2cs4[NDIMS], uv_cs4, vw_cs4, uw_cs4, mod_sq_2;
  double u[NDIMS], rho, rho_w[NDIMS];
  amrex::IntVect pos(0);

  for (amrex::MFIter mfi(dist_fn); mfi.isValid(); ++mfi) {
    amrex::FArrayBox& fab_dist_fn = dist_fn[mfi];

    // will need to replace NX/NY/NZ with local domain sizes
    for (int k = 0; k < NZ; ++k) {
      pos.setVal(2, k);
      for (int j = 0; j < NY; ++j) {
        pos.setVal(1, j);
        for (int i = 0; i < NX; ++i) {
          pos.setVal(0, i);

          // get density and velocity at this point in space
          rho = density(pos);
          u[0] = velocity(pos, 0);
          u[1] = velocity(pos, 1);
          u[2] = velocity(pos, 2);

          // calculate coefficients
          rho_w[0] = rho * 2.0 / 9.0;
          rho_w[1] = rho / 9.0;
          rho_w[2] = rho / 72.0;

          u2[0] = u[0]*u[0];
          u2[1] = u[1]*u[1];
          u2[2] = u[2]*u[2];

          u_cs2[0] = u[0] / CS2;
          u_cs2[1] = u[1] / CS2;
          u_cs2[2] = u[2] / CS2;

          u2_2cs4[0] = u2[0] / (2.0 * CS2 * CS2);
          u2_2cs4[1] = u2[1] / (2.0 * CS2 * CS2);
          u2_2cs4[2] = u2[2] / (2.0 * CS2 * CS2);

          uv_cs4 = u_cs2[0] * u_cs2[1];
          vw_cs4 = u_cs2[1] * u_cs2[2];
          uw_cs4 = u_cs2[0] * u_cs2[2];

          mod_sq = (u2[0] + u2[1] + u2[2]) / (2.0 * CS2);
          mod_sq_2 = (u2[0] + u2[1] + u2[2]) * (1 - CS2) / (2.0 * CS2 * CS2);

          // set distribution function
          fab_dist_fn(pos, 0) = rho_w[0] * (1.0 - mod_sq);

          fab_dist_fn(pos, 1) = rho_w[1] * (1.0 - mod_sq + u_cs2[0] + u2_2cs4[0]);
          fab_dist_fn(pos, 2) = rho_w[1] * (1.0 - mod_sq - u_cs2[0] + u2_2cs4[0]);
          fab_dist_fn(pos, 3) = rho_w[1] * (1.0 - mod_sq + u_cs2[1] + u2_2cs4[1]);
          fab_dist_fn(pos, 4) = rho_w[1] * (1.0 - mod_sq - u_cs2[1] + u2_2cs4[1]);
          fab_dist_fn(pos, 5) = rho_w[1] * (1.0 - mod_sq + u_cs2[2] + u2_2cs4[2]);
          fab_dist_fn(pos, 6) = rho_w[1] * (1.0 - mod_sq - u_cs2[2] + u2_2cs4[2]);

          fab_dist_fn(pos, 7) = rho_w[2] * (1.0 + u_cs2[0] + u_cs2[1] + u_cs2[2]
                                        + uv_cs4 + vw_cs4 + uw_cs4 + mod_sq_2);
          fab_dist_fn(pos, 8) = rho_w[2] * (1.0 + u_cs2[0] + u_cs2[1] - u_cs2[2]
                                        + uv_cs4 - vw_cs4 - uw_cs4 + mod_sq_2);
          fab_dist_fn(pos, 9) = rho_w[2] * (1.0 + u_cs2[0] - u_cs2[1] + u_cs2[2]
                                        - uv_cs4 - vw_cs4 + uw_cs4 + mod_sq_2);
          fab_dist_fn(pos, 10) = rho_w[2] * (1.0 + u_cs2[0] - u_cs2[1] - u_cs2[2]
                                         - uv_cs4 + vw_cs4 - uw_cs4 + mod_sq_2);
          fab_dist_fn(pos, 11) = rho_w[2] * (1.0 - u_cs2[0] + u_cs2[1] + u_cs2[2]
                                         - uv_cs4 + vw_cs4 - uw_cs4 + mod_sq_2);
          fab_dist_fn(pos, 12) = rho_w[2] * (1.0 - u_cs2[0] + u_cs2[1] - u_cs2[2]
                                         - uv_cs4 - vw_cs4 + uw_cs4 + mod_sq_2);
          fab_dist_fn(pos, 13) = rho_w[2] * (1.0 - u_cs2[0] - u_cs2[1] + u_cs2[2]
                                         + uv_cs4 - vw_cs4 - uw_cs4 + mod_sq_2);
          fab_dist_fn(pos, 14) = rho_w[2] * (1.0 - u_cs2[0] - u_cs2[1] - u_cs2[2]
                                         + uv_cs4 + vw_cs4 + uw_cs4 + mod_sq_2);
        } // i
      } // j
    } // k
  } // MultiFabIter

  return;
}

int Simulation::iterate(int const nsteps) {
  for (int t = 0; t < nsteps; ++t){
    collide();
    propagate();
    time_step++;
  }

  return time_step;
}

void Simulation::printDensity() {
  amrex::IntVect pos(0);

  for (int k = 0; k < NZ; ++k) {
    pos.setVal(2, k);
    for (int j = 0; j < NY; ++j) {
      pos.setVal(1, j);
      for (int i = 0; i < NX; ++i) {
        pos.setVal(0, i);
        printf("%g ", density(pos));
      }
      printf("\n");
    }
    printf("\n");
  }

  return;
}

// static member definitions
const double Simulation::DELTA[NDIMS][NDIMS] = { {1.0 / NMODES, 0.0, 0.0},
                                           {0.0, 1.0 / NMODES, 0.0},
                                           {0.0, 0.0, 1.0 / NMODES} };

const double Simulation::MODE_MATRIX[NMODES][NMODES] =
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

const double Simulation::MODE_MATRIX_INVERSE[NMODES][NMODES] =
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
