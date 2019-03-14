#include <iostream>
#include <array>
#include "AMRSim.h"

void AmrSim::ErrorEst(int level, amrex::TagBoxArray& tba, double time, int ngrow) { return; }

void AmrSim::MakeNewLevelFromScratch(int level, double time, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) { return ; }

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
OMEGA_B(1.0/(tau_b+0.5)),
idx_domain( amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
            amrex::IntVect(AMREX_D_DECL(nx-1, ny-1, nz-1)),
            amrex::IndexType( {AMREX_D_DECL(0, 0, 0)} ) ),
phys_domain( domain ), ba_domain(idx_domain), dm(ba_domain),
density(ba_domain, dm, 1, 0), velocity(ba_domain, dm, NDIMS, 0),
dist_fn(ba_domain, dm, NMODES, HALO_DEPTH)
{
  printf("Initialising AmrSim\n");
  // rewrite dummy domain here
  verbose = 1;
  SetMaxGridSize(amrex::IntVect(AMREX_D_DECL(nx-1,ny-1,nz-1)));
  SetBlockingFactor(amrex::IntVect(AMREX_D_DECL(nx-1,ny-1,nz-1)));
  SetDistributionMap(0, dm);
  SetBoxArray(0, ba_domain);
  printGridSummary(std::cout, 0, 0);
  std::cout << Geom(0).isAllPeriodic() << std::endl;
};
