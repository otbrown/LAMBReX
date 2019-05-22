#include <iostream>
#include <array>
#include "AMRSim.h"
#include "AMReX_filcc_f.H"
#include "AMReX_FillPatchUtil.H"

void AmrSim::SwapElements(double * const arr, int const offset_a,
int const offset_b) {
  double tmp = arr[offset_a];
  arr[offset_a] = arr[offset_b];
  arr[offset_b] = tmp;
  return;
}

void AmrSim::UpdateBoundaries(int const LEVEL) {
  dist_fn.at(LEVEL).FillBoundary(geom[LEVEL].periodicity());
  return;
}

void AmrSim::Collide(int const LEVEL) {
  amrex::IntVect pos(0);
  amrex::IntVect lo(0);
  amrex::IntVect hi(0);
  std::array<double, NMODES> mode;
  double usq, TrS;
  double stress[NMODES][NMODES];
  int a, b, m, p;

  for (amrex::MFIter mfi(dist_fn.at(LEVEL)); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    amrex::FArrayBox& fab_dist_fn = dist_fn.at(LEVEL)[mfi];
    amrex::FArrayBox& fab_density = density.at(LEVEL)[mfi];
    amrex::FArrayBox& fab_velocity = velocity.at(LEVEL)[mfi];

    lo = box.smallEnd();
    hi = box.bigEnd();

    for (int k = lo[2]; k <= hi[2]; ++k) {
      pos.setVal(2, k);
      for (int j = lo[1]; j <= hi[1]; ++j) {
        pos.setVal(1, j);
        for (int i = lo[0]; i <= hi[0]; ++i) {
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

void AmrSim::Propagate(int const LEVEL) {
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

  for (amrex::MFIter mfi(dist_fn.at(LEVEL)); mfi.isValid(); ++mfi) {
    data = dist_fn.at(LEVEL)[mfi].dataPtr();
    dims = dist_fn.at(LEVEL)[mfi].length();

    // this looping is bad for FAB, but will need to reorder swaps too...
    for (int i = 0; i < dims[0]; ++i) {
      for (int j = 0; j < dims[1]; ++j) {
        for (int k = 0; k < dims[2]; ++k) {

          // [1,0,0]
          if (i < dims[0]-1)
          SwapElements(data, FLindex(i,j,k,1,dims), FLindex(i+1,j,k,2,dims));
          // [0,1,0]
          if (j < dims[1]-1)
          SwapElements(data, FLindex(i,j,k,3,dims), FLindex(i,j+1,k,4,dims));
          // [0,0,1]
          if (k < dims[2]-1)
          SwapElements(data, FLindex(i,j,k,5,dims), FLindex(i,j,k+1,6,dims));

          // [1,1,1]
          if (i < dims[0]-1 && j < dims[1]-1 && k < dims[2]-1)
          SwapElements(data, FLindex(i,j,k,7,dims), FLindex(i+1,j+1,k+1,14,dims));
          // [1,1,-1]
          if (i < dims[0]-1 && j < dims[1]-1 && k > 0)
          SwapElements(data, FLindex(i,j,k,8,dims), FLindex(i+1,j+1,k-1,13,dims));
          // [1,-1,1]
          if (i < dims[0]-1 && j > 0 && k < dims[2]-1)
          SwapElements(data, FLindex(i,j,k,9,dims), FLindex(i+1,j-1,k+1,12,dims));
          // [1,-1,-1]
          if (i < dims[0]-1 && j > 0 && k > 0)
          SwapElements(data, FLindex(i,j,k,10,dims), FLindex(i+1,j-1,k-1,11,dims));

          // reorder
          SwapElements(data, FLindex(i,j,k,1,dims), FLindex(i,j,k,2,dims));
          SwapElements(data, FLindex(i,j,k,3,dims), FLindex(i,j,k,4,dims));
          SwapElements(data, FLindex(i,j,k,5,dims), FLindex(i,j,k,6,dims));

          SwapElements(data, FLindex(i,j,k,7,dims), FLindex(i,j,k,14,dims));
          SwapElements(data, FLindex(i,j,k,8,dims), FLindex(i,j,k,13,dims));
          SwapElements(data, FLindex(i,j,k,9,dims), FLindex(i,j,k,12,dims));
          SwapElements(data, FLindex(i,j,k,10,dims), FLindex(i,j,k,11,dims));
        }
      }
    }

  }

  return;
}

void AmrSim::InitDensity(int const LEVEL) {
  amrex::IntVect dims(NX, NY, NZ);
  amrex::IntVect lo(0);
  amrex::IntVect hi(0);
  amrex::IntVect pos(0);
  amrex::IntVect pos_oob(0);
  int lindex;

  if (!LEVEL) {
    for (amrex::MFIter mfi(density.at(LEVEL)); mfi.isValid(); ++mfi) {
      const amrex::Box& box = mfi.validbox();
      amrex::FArrayBox& rho = density.at(LEVEL)[mfi];

      lo = box.smallEnd();
      hi = box.bigEnd();

      for (int k = lo[2]; k <= hi[2] && k < NZ; ++k) {
        pos.setVal(2, k);
        for (int j = lo[1]; j <= hi[1] && j < NY; ++j) {
          pos.setVal(1, j);
          for (int i = lo[0]; i <= hi[0] && i < NX; ++i) {
            pos.setVal(0, i);
            lindex = CLindex(i, j, k, 0, dims, 1);
            rho(pos) = initial_density.at(lindex);
          }
        }
      }

      // deal with possibly empty boundary cells
      if (NX < hi[0]) {
        pos.setVal(0, NX);
        pos_oob.setVal(0, hi[0]);
        for (int j = lo[1]; j <= hi[1]; ++j) {
          pos.setVal(1, j);
          pos_oob.setVal(1, j);
          for (int k = lo[2]; k <= hi[2]; ++k) {
            pos.setVal(2, k);
            pos_oob.setVal(2, k);
            rho(pos_oob) = rho(pos);
          }
        }
      }

      if (NY < hi[1]) {
        pos.setVal(1, NY);
        pos_oob.setVal(1, hi[1]);
        for (int i = lo[0]; i <= hi[0]; ++i) {
          pos.setVal(0, i);
          pos_oob.setVal(0, i);
          for (int k = lo[2]; k <= hi[2]; ++k) {
            pos.setVal(2, k);
            pos_oob.setVal(2, k);
            rho(pos_oob) = rho(pos);
          }
        }
      }

      if (NZ < hi[2]) {
        pos.setVal(2, NZ);
        pos_oob.setVal(2, hi[2]);
        for (int i = lo[0]; i <= hi[0]; ++i) {
          pos.setVal(0, i);
          pos_oob.setVal(0, i);
          for (int j = lo[1]; j <= hi[1]; ++j) {
            pos.setVal(1, j);
            pos_oob.setVal(1, j);
            rho(pos_oob) = rho(pos);
          }
        }
      }
    }
  } else {
    amrex::Abort("Only level 0 should be initialised from scratch currently.");
  }

  return;
}

void AmrSim::InitVelocity(int const LEVEL) {
  amrex::IntVect dims(NX, NY, NZ);
  amrex::IntVect lo(0);
  amrex::IntVect hi(0);
  amrex::IntVect pos(0);
  amrex::IntVect pos_oob(0);
  int lindex;

  if (!LEVEL) {
    for (amrex::MFIter mfi(velocity.at(LEVEL)); mfi.isValid(); ++mfi) {
      const amrex::Box& box = mfi.validbox();
      amrex::FArrayBox& u = velocity.at(LEVEL)[mfi];

      lo = box.smallEnd();
      hi = box.bigEnd();

      for (int k = lo[2]; k <= hi[2] && k < NZ; ++k) {
        pos.setVal(2, k);
        for (int j = lo[1]; j <= hi[1] && j < NY; ++j) {
          pos.setVal(1, j);
          for (int i = lo[0]; i <= hi[0] && i < NX; ++i) {
            pos.setVal(0, i);
            for (int n = 0; n < NDIMS; ++n) {
              lindex = CLindex(i, j, k, n, dims, NDIMS);
              u(pos, n) = initial_velocity.at(lindex);
            }
          }
        }
      }

      // deal with possibly empty boundary cells
      if (NX < hi[0]) {
        pos.setVal(0, NX);
        pos_oob.setVal(0, hi[0]);
        for (int j = lo[1]; j <= hi[1]; ++j) {
          pos.setVal(1, j);
          pos_oob.setVal(1, j);
          for (int k = lo[2]; k <= hi[2]; ++k) {
            pos.setVal(2, k);
            pos_oob.setVal(2, k);
            for (int n = 0; n < NDIMS; ++n) u(pos_oob, n) = u(pos, n);
          }
        }
      }

      if (NY < hi[1]) {
        pos.setVal(1, NY);
        pos_oob.setVal(1, hi[1]);
        for (int i = lo[0]; i <= hi[0]; ++i) {
          pos.setVal(0, i);
          pos_oob.setVal(0, i);
          for (int k = lo[2]; k <= hi[2]; ++k) {
            pos.setVal(2, k);
            pos_oob.setVal(2, k);
            for (int n = 0; n < NDIMS; ++n) u(pos_oob, n) = u(pos, n);
          }
        }
      }

      if (NZ < hi[2]) {
        pos.setVal(2, NZ);
        pos_oob.setVal(2, hi[2]);
        for (int i = lo[0]; i <= hi[0]; ++i) {
          pos.setVal(0, i);
          pos_oob.setVal(0, i);
          for (int j = lo[1]; j <= hi[1]; ++j) {
            pos.setVal(1, j);
            pos_oob.setVal(1, j);
            for (int n = 0; n < NDIMS; ++n) u(pos_oob, n) = u(pos, n);
          }
        }
      }
    }
  } else {
    amrex::Abort("Only LEVEL 0 should be initialised from scratch currently.");
  }

  return;
}

void AmrSim::ComputeDt(int const LEVEL) {
  // we will likely need something like this (from AmrCoreAdv) to calculate the
  // increase in sim time at different refinement LEVELs
  dt.at(LEVEL) = 1.0;
  return;
}

void AmrSim::IterateLevel(int const LEVEL) {
  // iterates one LEVEL by one time *step*
  Collide(LEVEL);
  UpdateBoundaries(LEVEL);
  Propagate(LEVEL);

  // update time and step count
  sim_time.at(LEVEL) += dt.at(LEVEL);
  time_step.at(LEVEL)++;

  return;
}

void AmrSim::DistFnFillShim(double * data, const int * lo, const int * hi,
  const int * dom_lo, const int * dom_hi, const double * dx,
  const double * grd_lo, const double * time, const int * bc) {
  const int nq = NMODES;
  // if you look at the call stack in the tutorial this looks like it shouldn't
  // work as bc is an int not an array, but operator() is overloaded in
  // amrex::BndryFuncArray

  amrex_fab_filcc(data, lo, hi, &nq, dom_lo, dom_hi, dx, grd_lo, bc);

  return;
}

void AmrSim::DistFnFillPatch(const int level, amrex::MultiFab& dest_mf) {
  // based on AmrCoreAdv::FillPatch from AmrCore tutorial
  // just selects FillPatchSingleLevel if at level 0, or FillPatchTwoLevels, to
  // fill from coarser level otherwise

  if(!level) {
    const amrex::Vector<amrex::MultiFab*> source_mf{&dist_fn[level]};
    const amrex::Vector<double> source_time{sim_time[level]};
    amrex::PhysBCFunct<amrex::BndryFuncArray> physbc(geom[level], f_bndry,
      bfunc);

    amrex::FillPatchSingleLevel(dest_mf, sim_time[level], source_mf,
      source_time, 0, 0, NMODES, geom[level], physbc, 0);
  } else {
    const amrex::Vector<amrex::MultiFab*> coarse_mf{&dist_fn[level-1]};
    const amrex::Vector<amrex::MultiFab*> fine_mf{&dist_fn[level]};
    const amrex::Vector<double> coarse_time{sim_time[level-1]};
    const amrex::Vector<double> fine_time{sim_time[level]};
    amrex::PhysBCFunct<amrex::BndryFuncArray> coarse_physbc(geom[level-1],
      f_bndry, bfunc);
    amrex::PhysBCFunct<amrex::BndryFuncArray> fine_physbc(geom[level],
      f_bndry, bfunc);

    amrex::FillPatchTwoLevels(dest_mf, sim_time[level], coarse_mf, coarse_time,
      fine_mf, fine_time, 0, 0, NMODES, geom[level-1], geom[level],
      coarse_physbc, 0, fine_physbc, 0, refRatio(level-1), mapper, f_bndry, 0);
  }

  return;
}

void AmrSim::DistFnFillFromCoarse(const int level, amrex::MultiFab& fine_mf) {
  // Again based on AmrCore tutorial, this time AmrCoreAdv::FillCoarsePatch
  // fills the distribution function multifab at the specified level from the
  // level above
  // Different from FillPatchTwoLevels which does some sort of rough temporal
  // interpolatation too, apparently...
  if (!level) amrex::Abort("Cannot fill level 0 from coarse.");

  amrex::PhysBCFunct<amrex::BndryFuncArray> coarse_physbc(geom[level-1],
    f_bndry, bfunc);
  amrex::PhysBCFunct<amrex::BndryFuncArray> fine_physbc(geom[level],
    f_bndry, bfunc);

  amrex::InterpFromCoarseLevel(fine_mf, sim_time[level], dist_fn[level-1], 0, 0,
    NMODES, geom[level-1], geom[level], coarse_physbc, 0, fine_physbc, 0,
    refRatio(level-1), mapper, f_bndry, 0);

  return;
}

bool AmrSim::TagCell(const amrex::FArrayBox& f, const amrex::IntVect& pos) {
  // should this cell be tagged for refinement?
  // for now the answer is always no.
  return false;
}

void AmrSim::ErrorEst(int level, amrex::TagBoxArray& tba, double time, int ngrow) {
  // implements control logic for levels which require refinement
  amrex::IntVect pos;
  amrex::IntVect lo;
  amrex::IntVect hi;

  for (amrex::MFIter mfi(dist_fn.at(level)); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    amrex::FArrayBox& fab_dist_fn = dist_fn.at(level)[mfi];
    amrex::TagBox& tagfab = tba[mfi];

    lo = box.smallEnd();
    hi = box.bigEnd();

    for (int k = lo[2]; k <= hi[2]; ++k) {
      pos.setVal(2, k);
      for (int j = lo[1]; j <= hi[1]; ++j) {
        pos.setVal(1, j);
        for (int i = lo[0]; i <= hi[0]; ++i) {
          pos.setVal(0, i);
          if (TagCell(fab_dist_fn, pos)) {
            // tag the cell for refinement
            tagfab(pos) = amrex::TagBox::SET;
          } else tagfab(pos) = amrex::TagBox::CLEAR;
        }
      }
    }
  }

  return;
}

void AmrSim::MakeNewLevelFromScratch(int level, double time, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) {
  // define MultiFabs
  density[level].define(ba, dm, 1, 0);
  velocity[level].define(ba, dm, NDIMS, 0);
  dist_fn[level].define(ba, dm, NMODES, HALO_DEPTH);

  // set up simulation timings
  sim_time.at(level) = time;
  ComputeDt(level);
  time_step.at(level) = 0;

  // copy user-defined initial values to MultiFabs
  InitDensity(level);
  InitVelocity(level);

  // Calculate distribution function from initial density and velocity
  CalcEquilibriumDist(level);

  return;
}

void AmrSim::MakeNewLevelFromCoarse(int level, double time, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) {
  if (!level) amrex::Abort("Cannot construct level 0 from a coarser level.");

  // define MultiFabs
  density[level].define(ba, dm, 1, 0);
  velocity[level].define(ba, dm, NDIMS, 0);
  dist_fn[level].define(ba, dm, NMODES, HALO_DEPTH);

  // set up simulation timings
  sim_time.at(level) = time;
  ComputeDt(level);
  time_step.at(level) = 0;

  // fill distribution function MF with data interpolated from coarse level
  DistFnFillFromCoarse(level, dist_fn.at(level));

  // update boundary conditions
  UpdateBoundaries(level);

  // calculate hydrodynamic variables at this level
  CalcHydroVars(level);

  return;
}

void AmrSim::RemakeLevel(int level, double time, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) {
  amrex::MultiFab new_rho(ba, dm, 1, 0);
  amrex::MultiFab new_u(ba, dm, NDIMS, 0);
  amrex::MultiFab new_f(ba, dm, NMODES, HALO_DEPTH);

  // fill new distribution function with (possibly interpolated data) from old
  // one
  DistFnFillPatch(level, new_f);

  // swap new MFs in, swapping ensures old ones are appropriately destructed at
  // the return of this function
  std::swap(new_f, dist_fn[level]);
  std::swap(new_rho, density[level]);
  std::swap(new_u, velocity[level]);

  // update boundaries of distribution function
  UpdateBoundaries(level);

  // calculate hydrodynamic variables on this level using new dist_fn
  CalcHydroVars(level);

  return;
}

void AmrSim::ClearLevel(int level) {
  density.at(level).clear();
  velocity.at(level).clear();
  dist_fn.at(level).clear();

  return;
}

AmrSim::AmrSim(double const tau_s, double const tau_b)
  : NX(geom[0].Domain().length(0)), NY(geom[0].Domain().length(1)),
    NZ(geom[0].Domain().length(2)), NUMEL(NX*NY*NZ), COORD_SYS(0),
    PERIODICITY{ geom[0].period(0), geom[0].period(1), geom[0].period(2) },
    TAU_S(tau_s), TAU_B(tau_b), OMEGA_S(1.0/(tau_s+0.5)),
    OMEGA_B(1.0/(tau_b+0.5)), bfunc(DistFnFillShim), sim_time(1, 0.0), dt(1, 0.0),
    time_step(1, 0)
{
  std::cout << "NX: " << NX << " NY: " << NY << " NZ: " << NZ << std::endl;
  // resize vectors
  int num_levels = max_level + 1;
  density.resize(num_levels);
  velocity.resize(num_levels);
  dist_fn.resize(num_levels);
  f_bndry.resize(NMODES);

  for (int i = 0; i < NDIMS; ++i) {
    if (PERIODICITY[i]) {
      for (int comp = 0; comp < NMODES; ++comp) {
        f_bndry[comp].setLo(i, amrex::BCType::int_dir);
        f_bndry[comp].setHi(i, amrex::BCType::int_dir);
      }
    } else {
      amrex::Abort("Currently only periodic boundary conditions allowed.");
    }
  }
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
  initial_velocity.assign(3*NUMEL, u_init);
  return;
}

void AmrSim::SetInitialVelocity(const std::vector<double> u_init) {
  initial_velocity = u_init;
  return;
}

double AmrSim::GetDensity(const int i, const int j, const int k,
const int LEVEL) const {
  amrex::IntVect pos(i,j,k);
  for (amrex::MFIter mfi(density.at(LEVEL)); mfi.isValid(); ++mfi) {
    if (density.at(LEVEL)[mfi].box().contains(pos))
      return density.at(LEVEL)[mfi](pos);
  }
  amrex::Abort("Invalid index to density MultiFab.");
  return -1.0;
}

double AmrSim::GetVelocity(const int i, const int j, const int k,
const int n, const int LEVEL) const {
  amrex::IntVect pos(i,j,k);
  for (amrex::MFIter mfi(velocity.at(LEVEL)); mfi.isValid(); ++mfi) {
    if (velocity.at(LEVEL)[mfi].box().contains(pos))
      return velocity.at(LEVEL)[mfi](pos, n);
  }
  amrex::Abort("Invalid index to density MultiFab.");
  return -9999.0;
}

void AmrSim::CalcEquilibriumDist(int const LEVEL) {
  double u2[NDIMS], mod_sq, u_cs2[NDIMS], u2_2cs4[NDIMS], uv_cs4, vw_cs4, uw_cs4, mod_sq_2;
  double u[NDIMS], rho, rho_w[NDIMS];
  amrex::IntVect pos(0);
  amrex::IntVect lo(0);
  amrex::IntVect hi(0);

  for (amrex::MFIter mfi(dist_fn.at(LEVEL)); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    lo = box.smallEnd();
    hi = box.bigEnd();

    amrex::FArrayBox& fab_dist_fn = dist_fn.at(LEVEL)[mfi];
    amrex::FArrayBox& fab_density = density.at(LEVEL)[mfi];
    amrex::FArrayBox& fab_velocity = velocity.at(LEVEL)[mfi];


    for (int k = lo[2]; k <= hi[2]; ++k) {
      pos.setVal(2, k);
      for (int j = lo[1]; j <= hi[1]; ++j) {
        pos.setVal(1, j);
        for (int i = lo[0]; i <= hi[0]; ++i) {
          pos.setVal(0, i);

          // get density and velocity at this point in space
          rho = fab_density(pos);
          u[0] = fab_velocity(pos, 0);
          u[1] = fab_velocity(pos, 1);
          u[2] = fab_velocity(pos, 2);

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

  UpdateBoundaries(LEVEL);

  return;
}

void AmrSim::CalcHydroVars(int const LEVEL) {
  amrex::IntVect pos(0);
  amrex::IntVect lo(0);
  amrex::IntVect hi(0);
  std::array<double, NMODES> mode;

  for (amrex::MFIter mfi(dist_fn.at(LEVEL)); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    amrex::FArrayBox& fab_dist_fn = dist_fn.at(LEVEL)[mfi];
    amrex::FArrayBox& fab_density = density.at(LEVEL)[mfi];
    amrex::FArrayBox& fab_velocity = velocity.at(LEVEL)[mfi];

    lo = box.smallEnd();
    hi = box.bigEnd();

    for (int k = lo[2]; k <= hi[2]; ++k) {
      pos.setVal(2, k);
      for (int j = lo[1]; j <= hi[1]; ++j) {
        pos.setVal(1, j);
        for (int i = lo[0]; i <= hi[0]; ++i) {
          pos.setVal(0, i);

          // it seems like only the first 4 elements of mode are used, even with
          // a force present (which there is not here...)
          for (int m = 0; m < NMODES; ++m) {
            mode[m] = 0.0;
            for (int p = 0; p < NMODES; ++p) {
              mode[m] += fab_dist_fn(pos, p) * MODE_MATRIX[m][p];
            }
          }

          fab_density(pos) = mode[0];
          for (int a = 0; a < NDIMS; ++a) {
            fab_velocity(pos, a) = mode[a+1] / mode[0];
          }

        } // i
      } // j
    } // k
  } // MFIter

  return;
}

void AmrSim::Iterate(int const nsteps) {
  // once there are multiple levels this will need to handle synchronisation
  // between them
  for (int t = 0; t < nsteps; ++t) {
    IterateLevel(0);
  }
  return;
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
