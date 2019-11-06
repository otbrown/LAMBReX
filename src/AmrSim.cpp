#include <iostream>
#include <array>
#include "AMRSim.h"
#include "AMReX_filcc_f.H"
#include "AMReX_FillPatchUtil.H"
#include "AMReX_Interpolater.H"

#include "amr_help.h"

void AmrSim::SwapElements(double * const arr, int const offset_a,
int const offset_b) {
  double tmp = arr[offset_a];
  arr[offset_a] = arr[offset_b];
  arr[offset_b] = tmp;
  return;
}

void AmrSim::UpdateBoundaries(int const LEVEL) {
  // TODO: check time
  levels.at(LEVEL).next.get<DistFn>().FillBoundary(geom[LEVEL].periodicity());
  return;
}

void AmrSim::Collide(const amrex::MultiFab& f_old, amrex::MultiFab& f_new,
  const double OMEGA_S, const double OMEGA_B) {

  for_point_in(f_old, [&](auto accessor) {
    auto mode = std::array<double, DistFn::NV>{};

    for (int m = 0; m < DistFn::NV; ++m) {
      for (int p = 0; p < DistFn::NV; ++p) {
        mode[m] += accessor(f_old, p) * MODE_MATRIX[m][p];
      }
    }

    const auto& density = mode[0];
    // no forcing is currently present in the model,
    // so we disregard uDOTf for now
    double velocity[NDIMS];
    double usq = 0.0;
    for (int a = 0; a < NDIMS; ++a) {
      velocity[a] = mode[a+1] / density;
      usq += velocity[a] * velocity[a];
    }

    double stress[NDIMS][NDIMS] = {
      {mode[4], mode[5], mode[6]},
      {mode[5], mode[7], mode[8]},
      {mode[6], mode[8], mode[9]}
    };

    // Form the trace
    double TrS = 0.0;
    for (int a = 0; a < NDIMS; ++a) {
	     TrS += stress[a][a];
    }

    // Form the traceless part
    for (int a = 0; a < NDIMS; ++a) {
      stress[a][a] -= (TrS / NDIMS);
    }

    // Relax the trace
    TrS -= OMEGA_B * (TrS - density*usq);

    // Relax the traceless part
    for (int a = 0; a < NDIMS; ++a) {
      for (int b = 0; b < NDIMS; ++b) {
        stress[a][b] -= OMEGA_S * (stress[a][b] - density
				                * ( velocity[a]
					              * velocity[b]
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
    for (int p = 0; p < NMODES; ++p) {
      double fp = 0;
      for (int m = 0; m < NMODES; ++m) {
	       fp += mode[m] * MODE_MATRIX_INVERSE[p][m];
      }
      accessor(f_new, p) = fp;
    }
    });

  return;
}

void AmrSim::Stream(int const LEVEL) {
  auto& f_nxt = levels[LEVEL].next.get<DistFn>();
  auto f_prop = field_traits<DistFn>::MakeLevelData(f_nxt.boxArray(),
    f_nxt.DistributionMap());

  for_point_in(f_nxt, [&f_nxt, &f_prop](const auto& dest) {
    DistFn::PropagatePoint(dest, f_nxt, f_prop);
  });

  // swap propagated distribution function in to next
  std::swap(f_nxt, f_prop);

  return;
}

void AmrSim::CollideLevel(const int LEVEL) {
  const double OMEGA_S = 1.0 / (tau_s.at(LEVEL)+0.5);
  const double OMEGA_B = 1.0 / (tau_b.at(LEVEL)+0.5);
  const auto& f_old = levels.at(LEVEL).now.get<DistFn>();
  auto& f_new = levels.at(LEVEL).next.get<DistFn>();

  Collide(f_old, f_new, OMEGA_S, OMEGA_B);

  f_new.FillBoundary(geom[LEVEL].periodicity());

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
    auto& rho_mf = levels.at(LEVEL).now.get<Density>();
    for (amrex::MFIter mfi(rho_mf); mfi.isValid(); ++mfi) {
      const amrex::Box& box = mfi.validbox();
      amrex::FArrayBox& rho = rho_mf[mfi];

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
  int r = 1;

  // dt_0 = 1, m_0 = 1
  // dt_n = (1/(r_n)) dt_(n-1), m_n = (1/r_n) m_(n-1)
  // this clearly requires dt to be defined already at LEVEL-1, but this should
  // always be the case, and there's no way around it unless we assume it's
  // constant. Which to be fair, it is.
  // tau_0 is initialised in AmrSim constructor
  // tau_n = r_n(tau_(n-1) - 0.5) + 0.5
  auto& fine_time = levels[LEVEL].time;
  if (LEVEL) {
    const auto& coarse_time = levels[LEVEL-1].time;
    r = refRatio(LEVEL-1)[0];
    fine_time.delta = coarse_time.delta / r;
    mass.at(LEVEL) = mass.at(LEVEL-1) / r;
    tau_s.at(LEVEL) = r * (tau_s.at(LEVEL-1) - 0.5) + 0.5;
    tau_b.at(LEVEL) = r * (tau_b.at(LEVEL-1) - 0.5) + 0.5;
  }
  else {
    fine_time.delta = 1.0;
    mass.at(LEVEL) = 1.0;
  }

  return;
}

void AmrSim::IterateLevel(int const LEVEL) {
  // iterates one LEVEL by one time *step*
  CollideAndStream(LEVEL);
  // update time and step count
  auto& time = levels.at(LEVEL).time;
  time.current += time.delta;
  ++time.step;

  return;
}

void AmrSim::SubCycle(int const BASE_LEVEL, int const NUM_STEPS) {
  if (BASE_LEVEL == finest_level)
    for (int iter = 0; iter < NUM_STEPS; ++iter) IterateLevel(BASE_LEVEL);
  else {
    IterateLevel(BASE_LEVEL);
    SubCycle(BASE_LEVEL+1, refRatio(BASE_LEVEL)[0]);
  }

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
    auto& lvl = levels[level];
    const amrex::Vector<amrex::MultiFab*> source_mf{&lvl.now.get<DistFn>()};
    const amrex::Vector<double> source_time{lvl.time.current};
    amrex::PhysBCFunct<amrex::BndryFuncArray> physbc(geom[level], f_bndry,
      bfunc);

    amrex::FillPatchSingleLevel(dest_mf, lvl.time.current, source_mf,
      source_time, 0, 0, NMODES, geom[level], physbc, 0);
  } else {
    auto& coarse = levels[level-1];
    auto& fine = levels[level];
    const amrex::Vector<amrex::MultiFab*> coarse_mf{&coarse.now.get<DistFn>()};
    const amrex::Vector<amrex::MultiFab*> fine_mf{&fine.now.get<DistFn>()};
    const amrex::Vector<double> coarse_time{coarse.time.current};
    const amrex::Vector<double> fine_time{fine.time.current};
    amrex::PhysBCFunct<amrex::BndryFuncArray> coarse_physbc(geom[level-1],
      f_bndry, bfunc);
    amrex::PhysBCFunct<amrex::BndryFuncArray> fine_physbc(geom[level],
      f_bndry, bfunc);

    amrex::FillPatchTwoLevels(dest_mf, fine.time.current, coarse_mf, coarse_time,
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
  // interpolatation too
  if (!level) amrex::Abort("Cannot fill level 0 from coarse.");

  amrex::PhysBCFunct<amrex::BndryFuncArray> coarse_physbc(geom[level-1],
    f_bndry, bfunc);
  amrex::PhysBCFunct<amrex::BndryFuncArray> fine_physbc(geom[level],
    f_bndry, bfunc);

  amrex::InterpFromCoarseLevel(fine_mf, levels[level].time.current, levels[level-1].now.get<DistFn>(), 0, 0,
    NMODES, geom[level-1], geom[level], coarse_physbc, 0, fine_physbc, 0,
    refRatio(level-1), mapper, f_bndry, 0);

  return;
}

void AmrSim::RohdeCycle(int const COARSE_LEVEL) {
  // Following algorithm for advancing LB simulation of a fine and coarse grid
  // as per section 2.1 of Rohde, Kandhai, Derksen, and van den Akker (2006).
  // This method operates on pairs of levels, and includes subcycling for
  // the refined level.
  // The method has been modified to account for the fact that our grid is not
  // locally refined, but instead multiple grids exist

  // calculate relaxation frequency at coarse and fine levels
  const double OMEGA_S_C = 1.0 / (tau_s.at(COARSE_LEVEL)+0.5);
  const double OMEGA_B_C = 1.0 / (tau_b.at(COARSE_LEVEL)+0.5);
  const double OMEGA_S_F = 1.0 / (tau_s.at(COARSE_LEVEL+1)+0.5);
  const double OMEGA_B_F = 1.0 / (tau_b.at(COARSE_LEVEL+1)+0.5);

  // build post collision multifabs
  const int REF_RATIO = refRatio(COARSE_LEVEL)[0];
  const auto& f_C = levels[COARSE_LEVEL].now.get<DistFn>();
  const auto& f_F = levels[COARSE_LEVEL+1].now.get<DistFn>();

  amrex::MultiFab f_C_pc(f_C.boxArray(), f_C.DistributionMap(), NMODES, 1);
  amrex::MultiFab f_F_pc(f_F.boxArray(), f_F.DistributionMap(), NMODES,
    2*REF_RATIO);

  // initialise post collision multifabs
  f_C_pc.setVal(0.0);
  f_F_pc.setVal(0.0);

  // step 1: Collide on coarse level
  Collide(f_C, f_C_pc, OMEGA_S_C, OMEGA_B_C);

  // branch here for subcycling -- if fine level is *finest* level just do
  // collide and stream, otherwise call RohdeCycle on fine level
  // in order to maintain a constant speed of sound the fine level should be
  // iterated as many times as the refinement ratio
  if (COARSE_LEVEL + 1 == finest_level) {
    for (int iter = 0; iter < REF_RATIO; ++iter) {
      // step 2: Collide on fine level
      Collide(f_F, f_F_pc, OMEGA_S_F, OMEGA_B_F);
      // step 3: Fill fine boundary from coarse
      Explode(COARSE_LEVEL);
      // step 4: Stream on fine level
      Stream(COARSE_LEVEL+1);
    }
    // we should make sure actual level data is kept up-to-date for finest level
    // note that this assumes the data in levels[FINE_LEVEL].next is meaningful
    CompleteTimeStep(COARSE_LEVEL+1);
  } else {
    for (int iter = 0; iter < REF_RATIO; ++iter) {
      RohdeCycle(COARSE_LEVEL+1);
    }
  }

  // step 5: Stream on coarse level
  // Note that we actually have to make sure we don't pull from cells which are covered by the fine grid
  // probably best to zero during collision
  StreamInterior(f_C_pc, levels[COARSE_LEVEL].next.get<DistFn>());

  // step 6: Sum from fine level in to ghost level
  SumFromFine(COARSE_LEVEL);

  // mark completion of time step for coarse level (includes now <-> next)
  CompleteTimeStep(COARSE_LEVEL);

  return;
}

void AmrSim::Explode(int const COARSE_LEVEL) {
  // Should take the coarse value and copy it to all fine cells covered by it.
  // They should have the *same* value as the coarse cell, it should *not* be evenly
  // divided amongst the fine cells.
  return;
}

void AmrSim::SumFromFine(int const COARSE_LEVEL) {
  return;
}

void AmrSim::StreamInterior(const amrex::MultiFab& f_SRC,
amrex::MultiFab& f_dest) {
  for_point_in(f_SRC, [&f_SRC, &f_dest](const auto& dest_acc) {
    DistFn::PropagatePoint(dest_acc, f_SRC, f_dest);
  });
  return;
}

void AmrSim::CompleteTimeStep(int const LEVEL) {
  auto& lvl = levels[LEVEL];
  lvl.time.current += lvl.time.delta;
  ++lvl.time.step;
  lvl.UpdateNow();
  return;
}

bool AmrSim::TagCell(int const level, const amrex::IntVect& pos) {
  // should this cell be tagged for refinement?
  if (static_tags.at(level).contains(pos)) return true;
  else return false;
}

void AmrSim::ErrorEst(int level, amrex::TagBoxArray& tba, double time, int ngrow) {
  // implements control logic for levels which require refinement
  amrex::IntVect pos;
  amrex::IntVect lo;
  amrex::IntVect hi;

  const auto& f = levels[level].now.get<DistFn>();
  for (amrex::MFIter mfi(f); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    amrex::TagBox& tagfab = tba[mfi];

    lo = box.smallEnd();
    hi = box.bigEnd();

    for (int k = lo[2]; k <= hi[2]; ++k) {
      pos.setVal(2, k);
      for (int j = lo[1]; j <= hi[1]; ++j) {
        pos.setVal(1, j);
        for (int i = lo[0]; i <= hi[0]; ++i) {
          pos.setVal(0, i);
          if (TagCell(level, pos)) {
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
  auto& lvl = levels[level];
  // define MultiFabs
  velocity[level].define(ba, dm, NDIMS, 0);
  lvl.Define(ba, dm);

  // set up simulation timings

  lvl.time.current = time;
  ComputeDt(level);
  lvl.time.step = 0;

  // Only try and fill the level if this is the coarsest level
  if (!level) {
    // copy user-defined initial values to MultiFabs
    InitDensity(level);
    InitVelocity(level);
    // Calculate distribution function from initial density and velocity
    CalcEquilibriumDist(level);
  }

  return;
}

void AmrSim::MakeNewLevelFromCoarse(int level, double time, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) {
  if (!level) amrex::Abort("Cannot construct level 0 from a coarser level.");

  // define MultiFabs
  velocity[level].define(ba, dm, NDIMS, 0);
  auto& lvl = levels[level];
  lvl.Define(ba, dm);
  // set up simulation timings
  lvl.time.current = time;
  ComputeDt(level);
  lvl.time.step = 0;

  // fill distribution function MF with data interpolated from coarse level
  DistFnFillFromCoarse(level, lvl.now.get<DistFn>());

  // update boundary conditions
  UpdateBoundaries(level);

  // calculate hydrodynamic variables at this level
  CalcHydroVars(level);

  return;
}

void AmrSim::RemakeLevel(int level, double time, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) {
  //amrex::MultiFab new_rho(ba, dm, 1, 0);
  amrex::MultiFab new_u(ba, dm, NDIMS, 0);
  auto new_f = field_traits<DistFn>::MakeLevelData(ba, dm);
  auto new_rho = field_traits<Density>::MakeLevelData(ba, dm);
  // fill new distribution function with (possibly interpolated data) from old
  // one
  DistFnFillPatch(level, new_f);

  // swap new MFs in, swapping ensures old ones are appropriately destructed at
  // the return of this function
  //std::swap(new_f, dist_fn[level]);
  auto& state = levels[level].now;
  std::swap(new_f, state.get<DistFn>());
  std::swap(new_rho, state.get<Density>());
  std::swap(new_u, velocity[level]);

  // update sim time
  levels[level].time.current = time;

  // update boundaries of distribution function
  UpdateBoundaries(level);

  // calculate hydrodynamic variables on this level using new dist_fn
  CalcHydroVars(level);

  return;
}

void AmrSim::ClearLevel(int level) {
  velocity.at(level).clear();
  levels.at(level).Clear();

  return;
}

AmrSim::AmrSim(double const tau_s_0, double const tau_b_0)
  : NX(geom[0].Domain().length(0)), NY(geom[0].Domain().length(1)),
    NZ(geom[0].Domain().length(2)), NUMEL(NX*NY*NZ), COORD_SYS(0),
    PERIODICITY{ geom[0].period(0), geom[0].period(1), geom[0].period(2) },
    bfunc(DistFnFillShim), levels(max_level + 1)
{
  std::cout << "NX: " << NX << " NY: " << NY << " NZ: " << NZ << std::endl;
  // resize vectors
  int num_levels = max_level + 1;
  velocity.resize(num_levels);

  f_bndry.resize(NMODES);
  tau_s.resize(num_levels);
  tau_b.resize(num_levels);
  mass.resize(num_levels);
  static_tags.resize(num_levels);

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

  // initialise relaxation time vectors at level 0
  tau_s.at(0) = tau_s_0;
  tau_b.at(0) = tau_b_0;
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

// TODO: have this return an optional
double AmrSim::GetDensity(const int i, const int j, const int k,
const int LEVEL) const {
  amrex::IntVect pos(i,j,k);
  const auto& rho_mf = levels[LEVEL].now.get<Density>();
  for (amrex::MFIter mfi(rho_mf); mfi.isValid(); ++mfi) {
    if (rho_mf[mfi].box().contains(pos))
      return rho_mf[mfi](pos);
  }
  return NL_DENSITY;
}

double AmrSim::GetVelocity(const int i, const int j, const int k,
const int n, const int LEVEL) const {
  amrex::IntVect pos(i,j,k);
  for (amrex::MFIter mfi(velocity.at(LEVEL)); mfi.isValid(); ++mfi) {
    if (velocity.at(LEVEL)[mfi].box().contains(pos))
      return velocity.at(LEVEL)[mfi](pos, n);
  }
  return NL_VELOCITY;
}

void AmrSim::CalcEquilibriumDist(int const LEVEL) {
  double u2[NDIMS], mod_sq, u_cs2[NDIMS], u2_2cs4[NDIMS], uv_cs4, vw_cs4, uw_cs4, mod_sq_2;
  double u[NDIMS], rho_w[NDIMS];
  amrex::IntVect pos(0);
  amrex::IntVect lo(0);
  amrex::IntVect hi(0);
  auto& state = levels[LEVEL].now;
  auto& f = state.get<DistFn>();
  const auto& rho = state.get<Density>();

  for (amrex::MFIter mfi(f); mfi.isValid(); ++mfi) {
    const auto& box = mfi.validbox();
    const auto& lo = box.smallEnd();
    const auto& hi = box.bigEnd();

    amrex::FArrayBox& fab_dist_fn = f[mfi];
    auto& fab_density = rho[mfi];
    auto& fab_velocity = velocity.at(LEVEL)[mfi];


    for (int k = lo[2]; k <= hi[2]; ++k) {
      pos.setVal(2, k);
      for (int j = lo[1]; j <= hi[1]; ++j) {
        pos.setVal(1, j);
        for (int i = lo[0]; i <= hi[0]; ++i) {
          pos.setVal(0, i);

          // get density and velocity at this point in space
          double rho = fab_density(pos);
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
  auto& state = levels[LEVEL].now;
  const auto& f = state.get<DistFn>();
  auto& rho = state.get<Density>();

  for (amrex::MFIter mfi(f); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    auto& fab_dist_fn = f[mfi];
    auto& fab_density = rho[mfi];
    amrex::FArrayBox& fab_velocity = velocity.at(LEVEL)[mfi];

    const auto& lo = box.smallEnd();
    const auto& hi = box.bigEnd();

    for (int k = lo[2]; k <= hi[2]; ++k) {
      for (int j = lo[1]; j <= hi[1]; ++j) {
        for (int i = lo[0]; i <= hi[0]; ++i) {
	  const amrex::IntVect pos = {i, j, k};

	  const auto mode = [&](){
	    auto mode = std::array<double, DistFn::NV>{};
	    for (int m = 0; m < DistFn::NV; ++m) {
	      mode[m] = 0.0;
	      for (int p = 0; p < DistFn::NV; ++p) {
		mode[m] += fab_dist_fn(pos, p) * MODE_MATRIX[m][p];
	      }
	    }
	    return mode;
	  }();

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
  if (!finest_level) {
    for (int t = 0; t < nsteps; ++t) {
      IterateLevel(0);
    }
  } else {
    for (int t = 0; t < nsteps; ++t) {
      //SubCycle(0, 1);  // old method
      RohdeCycle(0);
    }
  }

  return;
}

void AmrSim::SetStaticRefinement(int const level, const std::array<int, NDIMS>&
  lo_corner, const std::array<int, NDIMS>& hi_corner) {
  const amrex::IntVect LO(lo_corner.data());
  const amrex::IntVect HI(hi_corner.data());
  const amrex::Box tag_area(LO, HI);

  static_tags.at(level).define(tag_area);

  regrid(level, GetTime(level));

  return;
}

void AmrSim::UnsetStaticRefinement(int const level) {
  static_tags.at(level).clear();
  regrid(level, GetTime(level));
  return;
}

std::pair<std::array<int,NDIMS>, std::array<int,NDIMS>>
AmrSim::GetExtent(int const LEVEL) const {
  const auto& f = levels[LEVEL].now.get<DistFn>();
  const amrex::Box& domain = f.boxArray().minimalBox();
  amrex::IntVect lo = domain.smallEnd();
  amrex::IntVect hi = domain.bigEnd();

  std::array<int,NDIMS> lo_arr{lo[0], lo[1], lo[2]};
  std::array<int,NDIMS> hi_arr{hi[0], hi[1], hi[2]};

  return std::make_pair(lo_arr, hi_arr);
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
