#include <cstring>
#include <iostream>
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

void Lamb::calcEquilibriumDist() {
  double u2[3], mod_sq, u_cs2[3], u2_2cs4[3], uv_cs4, vw_cs4, uw_cs4, mod_sq_2;
  double u[3], rho, rho_w[3];
  amrex::IntVect pos(0);

  for (int i = 0; i < NX; ++i) {
    pos.setVal(0, i);
    for (int j = 0; j < NY; ++j) {
      pos.setVal(1, j);
      for (int k = 0; k < NZ; ++k) {
        pos.setVal(2, k);

        // get density and velocity at this point in space
        rho = (*density)(pos);
        u[0] = (*velocity)(pos, 0);
        u[1] = (*velocity)(pos, 1);
        u[2] = (*velocity)(pos, 2);

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
        (*dist_fn)(pos, 0) = rho_w[0] * (1.0 - mod_sq);

        (*dist_fn)(pos, 1) = rho_w[1] * (1.0 - mod_sq + u_cs2[0] + u2_2cs4[0]);
        (*dist_fn)(pos, 2) = rho_w[1] * (1.0 - mod_sq - u_cs2[0] + u2_2cs4[0]);
        (*dist_fn)(pos, 3) = rho_w[1] * (1.0 - mod_sq + u_cs2[1] + u2_2cs4[1]);
        (*dist_fn)(pos, 4) = rho_w[1] * (1.0 - mod_sq - u_cs2[1] + u2_2cs4[1]);
        (*dist_fn)(pos, 5) = rho_w[1] * (1.0 - mod_sq + u_cs2[2] + u2_2cs4[2]);
        (*dist_fn)(pos, 6) = rho_w[1] * (1.0 - mod_sq - u_cs2[2] + u2_2cs4[2]);

        (*dist_fn)(pos, 7) = rho_w[2] * (1.0 + u_cs2[0] + u_cs2[1] + u_cs2[2]
                                      + uv_cs4 + vw_cs4 + uw_cs4 + mod_sq_2);
        (*dist_fn)(pos, 8) = rho_w[2] * (1.0 + u_cs2[0] + u_cs2[1] - u_cs2[2]
                                      + uv_cs4 - vw_cs4 - uw_cs4 + mod_sq_2);
        (*dist_fn)(pos, 9) = rho_w[2] * (1.0 + u_cs2[0] - u_cs2[1] + u_cs2[2]
                                      - uv_cs4 - vw_cs4 + uw_cs4 + mod_sq_2);
        (*dist_fn)(pos, 10) = rho_w[2] * (1.0 + u_cs2[0] - u_cs2[1] - u_cs2[2]
                                       - uv_cs4 + vw_cs4 - uw_cs4 + mod_sq_2);
        (*dist_fn)(pos, 11) = rho_w[2] * (1.0 - u_cs2[0] + u_cs2[1] + u_cs2[2]
                                       - uv_cs4 + vw_cs4 - uw_cs4 + mod_sq_2);
        (*dist_fn)(pos, 12) = rho_w[2] * (1.0 - u_cs2[0] + u_cs2[1] - u_cs2[2]
                                       - uv_cs4 - vw_cs4 + uw_cs4 + mod_sq_2);
        (*dist_fn)(pos, 13) = rho_w[2] * (1.0 - u_cs2[0] - u_cs2[1] + u_cs2[2]
                                       + uv_cs4 - vw_cs4 - uw_cs4 + mod_sq_2);
        (*dist_fn)(pos, 14) = rho_w[2] * (1.0 - u_cs2[0] - u_cs2[1] - u_cs2[2]
                                       + uv_cs4 + vw_cs4 + uw_cs4 + mod_sq_2);
      }
    }
  }

  return;
}

void Lamb::printDensity() {
  amrex::IntVect pos(0);

  for (int k = 0; k < NZ; ++k) {
    pos.setVal(2, k);
    for (int j = 0; j < NY; ++j) {
      pos.setVal(1, j);
      for (int i = 0; i < NX; ++i) {
        pos.setVal(0, i);
        printf("%g ", (*density)(pos));
      }
      printf("\n");
    }
    printf("\n");
  }

  return;
}
