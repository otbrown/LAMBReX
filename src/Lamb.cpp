#include "Lamb.h"
#include "AMReX_IntVect.H"
#include "AMReX_IndexType.H"
#include "AMReX_Box.H"

void Lamb::buildGeo() {
  const amrex::IntVect lo_corner( AMREX_D_DECL(0, 0, 0) );
  const amrex::IntVect hi_corner( AMREX_D_DECL(_NX-1, _NY-1, _NZ-1) );
  const amrex::IndexType idx_type( {AMREX_D_DECL(0, 0, 0)} );

  idx_domain.setSmall(lo_corner);
  idx_domain.setBig(hi_corner);
  idx_domain.setType(idx_type);

  phys_domain.setLo( {AMREX_D_DECL(0.0, 0.0, 0.0)} );
  phys_domain.setHi( {AMREX_D_DECL((double) _NX-1, (double) _NY-1, (double) _NZ-1)} );

  geometry.define(idx_domain, &phys_domain, COORD_SYS, _periodicity);

  return;
}

Lamb::Lamb(int nx, int ny, int nz, double tau_s, double tau_b, int * periodicity)
: _NX(nx), _NY(ny), _NZ(nz), _TAU_S(tau_s), _TAU_B(tau_b) {
  for (int dim = 0; dim < _NDIMS; ++dim) {
    _periodicity[dim] = periodicity[dim];
  }

  density = new double[_NUMEL];
  velocity = new double[_NUMEL];
  force = new double[_NUMEL];

  buildGeo();

  return;
};

Lamb::~Lamb() {
  delete[] density;
  delete[] velocity;
  delete[] force;

  return;
}
