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
              cs_real_t yy)
{
  cs_var_cal_opt_t var_cal_opt_u, var_cal_opt_sca;
  const cs_fluid_properties_t *fp = cs_glob_fluid_properties;
  const cs_physical_constants_t *pc = cs_glob_physical_constants;
  int keysca = cs_field_key_id("scalar_id");
  int key_cal_opt_id = cs_field_key_id("var_cal_opt");
  int kvisls0 = cs_field_key_id("scalar_diffusivity_ref");
  cs_field_t *sca1 = NULL;

  for (int fid = 0; fid < cs_field_n_fields(); fid++) {
    cs_field_t *f = cs_field_by_id(fid);
    if (cs_field_get_key_int(f, keysca) == 1) {
      sca1 = f;
      break;
    }
  }

  cs_field_get_key_struct(CS_F_(u), key_cal_opt_id, &var_cal_opt_u);
  cs_field_get_key_struct(sca1, key_cal_opt_id, &var_cal_opt_sca);

  cs_real_t visls0 = cs_field_get_key_double(sca1, kvisls0);

  cs_real_t omega = 0.1, pi = acos(-1);
  cs_real_t mu, mv, uu, vv, pp, sca, ro1, ro2, rho, vmass, beta, phi_;
  cs_real_t dudx, dudy, dvdx, dvdy, dpdx, dpdy, dmudx, dmudy, dmvdx, dmvdy;
  cs_real_t d2udx2, d2udy2, d2udxdy, d2vdx2, d2vdy2, d2vdxdy;
  cs_real_t d2mudx2, d2mudy2, d2mudxdy, d2mvdx2, d2mvdy2, d2mvdxdy;
  cs_real_t dscadx, dscady, d2scadx2, d2scady2, d2scadxdy;
  int iconvu = var_cal_opt_u.iconv, iconvs = var_cal_opt_sca.iconv;
  int idiffu = var_cal_opt_u.idiff, idiffs = var_cal_opt_sca.idiff;

  /* Variables */
  phi_ = 0.;

  pp = sin(pi*(xx+0.5))*sin(pi*(yy+0.5));

  mu = (1.-cos(2.*pi*xx))*sin(2.*pi*yy);
  mv = sin(2.*pi*xx)*(cos(2.*pi*yy)-1.);

  sca = 0.5*(sin(2.*pi*xx)*sin(2.*pi*yy)+1.);

  ro1 = 1.;
  ro2 = omega*ro1;

  vmass  = sca/ro1 + (1.-sca)/ro2;
  rho = 1./vmass;
  beta= 1./ro1 - 1./ro2;

  uu = mu/rho;
  vv = mv/rho;

  /* Derivatives */

  dmudx    =  2.*pi*sin(2.*pi*xx)*sin(2.*pi*yy);
  dmudy    =  2.*pi*(1.-cos(2.*pi*xx))*cos(2.*pi*yy);
  d2mudx2  =  4.*pow(pi,2)*cos(2.*pi*xx)*sin(2.*pi*yy);
  d2mudy2  = -4.*pow(pi,2)*(1.-cos(2.*pi*xx))*sin(2.*pi*yy);
  d2mudxdy =  4.*pow(pi,2)*sin(2.*pi*xx)*cos(2.*pi*yy);

  dmvdx    =  2.*pi*cos(2.*pi*xx)*(cos(2.*pi*yy)-1.);
  dmvdy    = -2.*pi*sin(2.*pi*xx)*sin(2.*pi*yy);
  d2mvdx2  = -4.*pow(pi,2)*sin(2.*pi*xx)*(cos(2.*pi*yy)-1.);
  d2mvdy2  = -4.*pow(pi,2)*sin(2.*pi*xx)*cos(2.*pi*yy);
  d2mvdxdy = -4.*pow(pi,2)*cos(2.*pi*xx)*sin(2.*pi*yy);

  dscadx    =  pi*cos(2.*pi*xx)*sin(2.*pi*yy);
  dscady    =  pi*sin(2.*pi*xx)*cos(2.*pi*yy);
  d2scadx2  = -2.*pow(pi,2)*sin(2.*pi*xx)*sin(2.*pi*yy);
  d2scady2  = -2.*pow(pi,2)*sin(2.*pi*xx)*sin(2.*pi*yy);
  d2scadxdy =  2.*pow(pi,2)*cos(2.*pi*xx)*cos(2.*pi*yy);

  dpdx = pi*cos(pi*(xx+0.5))*sin(pi*(yy+0.5));
  dpdy = pi*sin(pi*(xx+0.5))*cos(pi*(yy+0.5));

  dudx = vmass*dmudx + mu*beta*dscadx;
  dudy = vmass*dmudy + mu*beta*dscady;
  dvdx = vmass*dmvdx + mv*beta*dscadx;
  dvdy = vmass*dmvdy + mv*beta*dscady;

  d2udx2 = vmass*d2mudx2  + beta*(2.*dmudx*dscadx + mu*d2scadx2);
  d2udy2 = vmass*d2mudy2  + beta*(2.*dmudy*dscady + mu*d2scady2);
  d2udxdy= vmass*d2mudxdy + beta*(dmudx*dscady + dmudy*dscadx + mu*d2scadxdy);
  d2vdx2 = vmass*d2mvdx2  + beta*(2.*dmvdx*dscadx + mv*d2scadx2);
  d2vdy2 = vmass*d2mvdy2  + beta*(2.*dmvdy*dscady + mv*d2scady2);
  d2vdxdy= vmass*d2mvdxdy + beta*(dmvdx*dscady + dmvdy*dscadx + mv*d2scadxdy);

  if (f_id == CS_F_(p)->id) {

    phi_ = pp;

  } else if (f_id == CS_F_(u)->id) {

    phi_ = uu;

  } else if (f_id == CS_F_(u)->id + 100) {

    phi_ = vv;

  } else if (f_id == CS_F_(u)->id + 1000) {

    phi_ = 0.;

  } else if (sca1 != NULL && f_id == sca1->id) {

    phi_ = sca;

  } else if (f_id == -1) {

    phi_ = rho;

  /* source terms (mu=rho*u, mv=rho*v, vm=1/rho) */
  } else if (f_id == -(CS_F_(u)->id+10)) {

    phi_ = iconvu * (mu*dudx + mv*dudy)
         - idiffu * fp->viscl0 * (d2udx2 + d2udy2 + 1./3.*(d2udx2 + d2vdxdy))
         + dpdx - (rho-fp->ro0)*pc->gravity[0];

  } else if (f_id == -(CS_F_(u)->id + 100)) {

    phi_ = iconvu * (mu*dvdx + mv*dvdy)
         - idiffu * fp->viscl0 * (d2vdx2 + d2vdy2 + 1./3.*(d2vdy2 + d2udxdy))
         + dpdy - (rho-fp->ro0)*pc->gravity[1];

  } else if (f_id == -(CS_F_(u)->id + 1000)) {

    phi_ = 0.;

  } else if (sca1 != NULL && f_id == -sca1->id) {

    phi_ = iconvs * (mu*dscadx + mv*dscady)- idiffs
         * visls0 * (d2scadx2 + d2scady2);

  } else {
    bft_printf("Unexpected variable : %d in function call\n", f_id);
    cs_exit(EXIT_FAILURE);
  }

  return phi_;
}

/*----------------------------------------------------------------------------*/
