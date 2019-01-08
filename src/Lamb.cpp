#include <cstring>
#include "Lamb.h"
#include "AMReX_IntVect.H"
#include "AMReX_IndexType.H"
#include "AMReX_Box.H"
#include "AMReX_FArrayBox.H"

void Lamb::buildGeometry() {
  const amrex::IntVect lo_corner( AMREX_D_DECL(0, 0, 0) );
  const amrex::IntVect hi_corner( AMREX_D_DECL(NX-1, NY-1, NZ-1) );
  const amrex::IndexType idx_type( {AMREX_D_DECL(0, 0, 0)} );

  idx_domain.setSmall(lo_corner);
  idx_domain.setBig(hi_corner);
  idx_domain.setType(idx_type);

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

  density = new double[NUMEL];
  velocity = new double[NUMEL * NDIMS];
  force = new double[NUMEL * NDIMS];

  buildGeometry();

  dist_fn = new amrex::FArrayBox(idx_domain, 1);

  return;
};

Lamb::~Lamb() {
  delete[] density;
  delete[] velocity;
  delete[] force;
  delete dist_fn;

  return;
}

void Lamb::setDensity(double const_density) {
  // set density array to one value everywhere
  for (int i = 0; i < NUMEL; ++i) density[i] = const_density;
  return;
}

void Lamb::setDensity(double * rho) {
  // set density to match provided array
  size_t count = NUMEL * sizeof(double);
  memcpy(density, rho, count);
  return;
}

void Lamb::setVelocity(double const_velocity) {
  // set velocity to one value everywhere
  for (int i = 0; i < NUMEL * NDIMS; ++i) velocity[i] = const_velocity;
  return;
}

void Lamb::setVelocity(double * u) {
  // set velocity to match provided array
  size_t count = NUMEL * NDIMS * sizeof(double);
  memcpy(velocity, u, count);
  return;
}
