/*----------------------------------------------------------------------------*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"
#include "bft_printf.h"

#include "fvm_writer.h"

#include "cs_base.h"
#include "cs_field.h"
#include "cs_gui_util.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_math.h"
#include "cs_mesh.h"
#include "cs_mesh_location.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_multigrid.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_prototypes.h"
#include "cs_sles.h"
#include "cs_sles_it.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_selector.h"

#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "prototypes.h"

/*----------------------------------------------------------------------------*/

cs_real_t phi(int f_id,
              cs_real_t xx,
              cs_real_t yy,
              cs_real_t tcur)
{
  cs_var_cal_opt_t var_cal_opt_u, var_cal_opt_sca;
  int keysca = cs_field_key_id("scalar_id");
  int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  cs_real_t omega = 2., pi = acos(-1), coef = 1.e-3, ro1 = 5., ro2 = 1.;
  cs_real_t mu = 1.e-3, kk= 2., uf = 0.5, vf = 0.5;
  cs_field_t *sca1 = NULL;

  for (int fid = 0; fid < cs_field_n_fields(); fid++) {
    cs_field_t *f = cs_field_by_id(fid);
    if (cs_field_get_key_int(f, keysca) == 1) {
      sca1 = f;
      break;
    }
  }

  cs_field_get_key_struct(CS_F_(u), key_cal_opt_id, &var_cal_opt_u);

  if (sca1 != NULL)
    cs_field_get_key_struct(sca1, key_cal_opt_id, &var_cal_opt_sca);

  cs_real_t xh, yh, u, v, p, sca, rho, q_z, q_u, q_v, s1, s2, s3, s4, s5;
  cs_real_t phi_, s6, s7, s8, s9, s10, s11, s12, s13, s14, dzdx, dzdy;
  cs_real_t d2rhodxdy, d2zdxdy, d2ud2x, d2vd2x, d2vd2y, d2ud2y,d2rhodzdx,d2rhodzdy;
  cs_real_t d2udxdy,d2vdxdy, d2rhod2x, d2rhod2y;
  cs_real_t dudx,dudy,dudt,dvdx,dvdy,dvdt,dpdx,dpdy, tvx, tvy;
  cs_real_t dro, sro, d2zd2x,d2zd2y,dzdt,drhodz,drhodx,drhody,drhodt;
  int iconvu, iconvs, idiffu, idiffs;

  iconvu = var_cal_opt_u.iconv;
  idiffu = var_cal_opt_u.idiff;

  if (sca1 != NULL) {
    /* For the coupled cases */
    iconvs = var_cal_opt_sca.iconv;
    idiffs = var_cal_opt_sca.idiff;
  } else {
    /* For the uncoupled cases (Stokes or Euler) */
    iconvs = var_cal_opt_u.iconv;
    idiffs = var_cal_opt_u.idiff;
  }

  /*-----------------------------*
   *     Case u_f = v_f = 0. *
   *-----------------------------*/

  dro = ro2-ro1;
  sro = ro1+ro2;

  xh = xx - uf*tcur;
  yh = yy - vf*tcur;

  s1 = cos(pi*kk*xh);
  s2 = sin(pi*kk*xh);
  s3 = cos(pi*kk*yh);
  s4 = sin(pi*kk*yh);
  s5 = cos(pi*omega*tcur);
  s6 = sin(pi*omega*tcur);
  s7 = ro1 / ro2;
  s8 = s2*s4*s5;
  s9 = s1*s4*s6;
  s10 = s2*s3*s6;
  s11 = 1+s7+(1-s7)*s8;
  s12 = 1-s7;
  s13 = s1*s4*s5/s11;
  s14 = s2*s3*s5/s11;

  /* Variables */
  phi_ = 0.;

  /* Scalar */
  sca = (1+s8)/s11;

  /* Density */
  rho = 1./(sca/ro2 + (1-sca)/ro1);

  /* Velocity */
  u = uf - omega*s9/(4.*kk)*dro/rho;

  v = vf - omega*s10/(4.*kk)*dro/rho;

  /* Pressure */
  p = rho*u*v/2.;

  /* Derivatives */

  /* Z derivatives */
  dzdx = pi*kk*s13*(1-s12*sca);
  dzdy = pi*kk*s14*(1-s12*sca);
  d2zd2x = pi*kk*(-s13*dzdx*s12+(1-s12*sca)*(pi*kk*(-s8/s11 - pow(s13,2)*s12)));
  d2zd2y = pi*kk*(-s14*dzdy*s12+(1-s12*sca)*(pi*kk*(-s8/s11 - pow(s14,2)*s12)));
  dzdt = pi*omega*s2*s4*s6/s11*(s12*sca-1);
  d2zdxdy = pi*kk*(-s12*s13*dzdy+(1-s12*sca)*(pi*kk*s1*s3*s5*s11-pi*kk*s1*
            s4*s5*(1-s7)*s2*s3*s5)/pow(s11,2));


  /* Rho derivatives */
  drhodz = -(1/ro2 - 1/ro1)*pow(rho,2);
  drhodx = drhodz * dzdx;
  drhody = drhodz * dzdy;
  drhodt = drhodz * dzdt;
  d2rhodzdx = -2*(1/ro2 - 1/ro1)*drhodx*rho;
  d2rhodzdy = -2*(1/ro2 - 1/ro1)*drhody*rho;
  d2rhod2x = drhodz * d2zd2x + dzdx*d2rhodzdx;
  d2rhod2y = drhodz * d2zd2y + dzdx*d2rhodzdy;
  d2rhodxdy = drhodz * d2zdxdy + dzdy * d2rhodzdx;

  /* U derivatives */
  dudx = -dro*omega/4/kk*(-pi*kk*s2*s4*s6*rho-s9*drhodx)/pow(rho,2);
  dudy = -dro*omega/4/kk*(pi*kk*s1*s3*s6*rho-s9*drhody)/pow(rho,2);
  dudt = -dro*omega/4/kk*(pi*omega*s1*s4*s5*rho-s1*s4*s6*drhodt)/pow(rho,2);

  d2ud2x = -dro*omega/4/kk/pow(rho,4)*(pow(rho,2)*s9*(-pow(pi*kk,2)*rho - d2rhod2x)
         - (-pi*kk*s2*s4*s6*rho-s9*drhodx)*2*rho*drhodx);
  d2ud2y = -dro*omega/4/kk/pow(rho,4)*(pow(rho,2)*s9*(-pow(pi*kk,2)*rho - d2rhod2y)
         - (pi*kk*s1*s3*s6*rho-s9*drhody)*2*rho*drhody);

  d2udxdy = -dro*omega/4/kk/pow(rho,4)*s6*((-pow(pi*kk,2)*s2*s3*rho+(pi*kk)*s1*
           s3*drhodx+pi*kk*s2*s4*drhody-s1*s4*d2rhodxdy)*pow(rho,2) - (pi*kk*
           s1*s3*rho - s1*s4*drhody)*2*rho*drhodx);

  /* V derivatives */
  dvdx = -dro*omega/4/kk*(pi*kk*s1*s3*s6*rho-s10*drhodx)/pow(rho,2);
  dvdy = -dro*omega/4/kk*(-pi*kk*s2*s4*s6*rho-s10*drhody)/pow(rho,2);
  dvdt = -dro*omega/4/kk*(pi*omega*s2*s3*s5*rho-s2*s3*s6*drhodt)/pow(rho,2);

  d2vd2x = -dro*omega/4/kk/pow(rho,4)*(pow(rho,2)*s10*(-pow(pi*kk,2)*rho - d2rhod2x)
         - (pi*kk*s1*s3*s6*rho-s10*drhodx)*2*rho*drhodx);
  d2vd2y = -dro*omega/4/kk/pow(rho,4)*(pow(rho,2)*s10*(pow(pi*kk,2)*rho - d2rhod2y)
         - (-pi*kk*s2*s4*s6*rho-s10*drhody)*2*rho*drhody);

  d2vdxdy = -dro*omega/4/kk/pow(rho,4)*s6*((-pow(pi*kk,2)*s1*s4*rho+(pi*kk)*s1*
            s3*drhody+pi*kk*s2*s4*drhodx-s2*s3*d2rhodxdy)*pow(rho,2) - (pi*kk*
            s1*s3*rho - s2*s3*drhodx)*2*rho*drhody);

  /* Pressure gradient */
  dpdx = 0.5*(rho*u*dvdx+rho*v*dudx+u*v*drhodx);
  dpdy = 0.5*(rho*u*dvdy+rho*v*dudy+u*v*drhody);

  /* Viscous terms */
  tvx = mu*((d2ud2x + d2ud2y) + 1./3. *(d2ud2x + d2vdxdy));

  tvy = mu*((d2vd2x + d2vd2y) + 1./3. *(d2vd2x + d2udxdy));

  /*-----------------------------------------------------------------------------*
   * Scalar source term                                                          *
   * q_z = d(rhoZ)/dt + d(rho u Z)/dx + d(rho v Z)/dy - D(d2(Z)/dx2 + d2(Z)/dy2) *
   * q_z = rho(dZ/dt + (u-uf)dZ/dx + (v-vf)dZ/dy) - D(d2(Z)/dx2 + d2(Z)/dy2)     *
   *-----------------------------------------------------------------------------*/

  q_z = rho*(dzdt + (iconvs*u-uf)*dzdx + (iconvs*v-vf)*dzdy) - idiffs*coef*(d2zd2x + d2zd2y);

  /*-----------------------------------------------------------------------------------*
   * X Velocity source term                                                            *
   * q_u = d(rho u)/dt + d(rho u u)/dx + d(rho u v)/dy + dp/dx - (d(t11)/dx+d(t12)/dy) *
   * q_u = rho(du/dt + (u-uf)du/dx + (v-vf)du/dy) + dp/dx - tvx                        *
   *-----------------------------------------------------------------------------------*/

  q_u = rho*(dudt + (iconvu*u-uf)*dudx + (iconvu*v-vf)*dudy) + dpdx - idiffu*tvx;

  /*-----------------------------------------------------------------------------------*
   * Y Velocity source term                                                            *
   * q_v = d(rho v)/dt + d(rho u v)/dx + d(rho v v)/dy + dp/dy - (d(t12)/dx+d(t22)/dy) *
   * q_v = rho(dv/dt + (u-uf)dv/dx + (v-vf)dv/dy) + dp/dy - tvy                        *
   *-----------------------------------------------------------------------------------*/

  q_v = rho*(dvdt + (iconvu*u-uf)*dvdx + (iconvu*v-vf)*dvdy) + dpdy - idiffu*tvy;

  /*-----------------------------------------------------------------------------------*/

  if (f_id == CS_F_(p)->id) {

    phi_ = p;

  } else if (f_id == CS_F_(u)->id) {

    phi_ = u;

  } else if (f_id == CS_F_(u)->id + 100) {

    phi_ = v;

  } else if (f_id == CS_F_(u)->id + 1000) {

    phi_ = 0.;

  } else if (sca1 != NULL && f_id == sca1->id) {

    phi_ = sca;

  } else if (f_id == -1) {

    phi_ = rho;

  } else if (f_id == -(CS_F_(u)->id+10)) {

    phi_ = q_u;

  } else if (f_id == -(CS_F_(u)->id + 100)) {

    phi_ = q_v;

  } else if (f_id == -(CS_F_(u)->id + 1000)) {

    phi_ = 0.;

  } else if (sca1 != NULL && f_id == -sca1->id) {

    phi_ = q_z;

  } else if (f_id == -50) {

    phi_ = q_z + idiffu*coef*(d2zd2x + d2zd2y);

  } else {
    bft_printf("Unexpected variable : %d in function call\n", f_id);
    cs_exit(EXIT_FAILURE);
  }

  return phi_;
}

/*----------------------------------------------------------------------------*/
