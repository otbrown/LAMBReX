#include <cstring>
#include "Lamb.h"
#include "AMReX_IntVect.H"
#include "AMReX_IndexType.H"
#include "AMReX_Box.H"
#include "AMReX_FArrayBox.H"

void Lamb::buildGeometry() {
  const amrex::IntVect LO_CORNER( AMREX_D_DECL(0, 0, 0) );
  const amrex::IntVect HI_CORNER( AMREX_D_DECL(NX-1, NY-1, NZ-1) );
  const amrex::IndexType IDX_TYPE( {AMREX_D_DECL(0, 0, 0)} );

  idx_domain.setSmall(LO_CORNER);
  idx_domain.setBig(HI_CORNER);
  idx_domain.setType(IDX_TYPE);

  phys_domain.setLo( {AMREX_D_DECL(0.0, 0.0, 0.0)} );
  phys_domain.setHi( {AMREX_D_DECL((double) NX-1, (double) NY-1, (double) NZ-1)} );

  geometry.define(idx_domain, &phys_domain, COORD_SYS, PERIODICITY);

  return;
}

Lamb::Lamb(int nx, int ny, int nz, double tau_s, double tau_b, int * periodicity)
: NX(nx), NY(ny), NZ(nz), TAU_S(tau_s), TAU_B(tau_b) {
  for (int dim = 0; dim < NDIMS; ++dim) {
    PERIODICITY[dim] = periodicity[dim];
  }

  buildGeometry();

  density = new amrex::FArrayBox(idx_domain);
  velocity = new amrex::FArrayBox(idx_domain, NDIMS);
  force = new amrex::FArrayBox(idx_domain, NDIMS);

  dist_fn = new amrex::FArrayBox(idx_domain, NUM_VELOCITIES);

  return;
};

Lamb::~Lamb() {
  delete density;
  delete velocity;
  delete force;
  delete dist_fn;

  return;
}

void Lamb::setDensity(double uniform_density) {
  // set density array to one value everywhere
  density->setVal(uniform_density);
  return;
}

void Lamb::setDensity(double * rho) {
  amrex::IntVect pos(0);
  // set density to match provided array
  for (int i = 0; i < NX; ++i) {
    pos.setVal(0, i);
    for (int j = 0; j < NY; ++j) {
      pos.setVal(1, j);
      for (int k = 0; k < NZ; ++k) {
        pos.setVal(2, k);
        (*density)(pos) = rho[i*NY*NZ + j*NZ + k];
      }
    }
  }
  return;
}

void Lamb::setVelocity(double uniform_velocity) {
  // set velocity to one value everywhere
  velocity->setVal(uniform_velocity);
  return;
}

void Lamb::setVelocity(double * u) {
  amrex::IntVect pos(0);
  // set velocity to match provided array
  for (int i = 0; i < NX; ++i) {
    pos.setVal(0, i);
    for (int j = 0; j < NY; ++j) {
      pos.setVal(1, j);
      for (int k = 0; k < NZ; ++k) {
        pos.setVal(2, k);
        for (int n = 0; n < NDIMS; ++n) {
          (*velocity)(pos, n) = u[i*NY*NZ*NDIMS + j*NZ*NDIMS + k*NDIMS + n];
        }
      }
    }
  }
  return;
}
