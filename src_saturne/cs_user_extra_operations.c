/*============================================================================
 * This function is called at the end of each time step, and has a very
 *  general purpose
 *  (i.e. anything that does not have another dedicated user subroutine)
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_halo.h"
#include "cs_halo_perio.h"
#include "cs_log.h"
#include "cs_parall.h"
#include "cs_parameters.h"
#include "cs_physical_constants.h"
#include "cs_prototypes.h"
#include "cs_rotation.h"
#include "cs_time_moment.h"
#include "cs_time_step.h"
#include "cs_turbomachinery.h"
#include "cs_selector.h"

#include "cs_post.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_extra_operations.c
 *
 * \brief This function is called at the end of each time step, and has a very
 * general purpose (i.e. anything that does not have another dedicated
 * user subroutine)
 *
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is called at the end of each time step.
 *
 * It has a very general purpose, although it is recommended to handle
 * mainly postprocessing or data-extraction type operations.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations(void)
{
  FILE *file = NULL;
  int n_fields = cs_field_n_fields();
  const cs_lnum_t n_cells     = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const cs_real_t *volume = cs_glob_mesh_quantities->cell_vol;
  const cs_real_3_t *cell_cen = (const cs_real_3_t *)cs_glob_mesh_quantities->cell_cen;
  cs_time_step_t *ts = cs_get_glob_time_step();
  cs_real_t Pth, Uth, Vth, Tth, xx, yy, zz;
  cs_real_t *varres;
  cs_field_t *sca1 = NULL;
  cs_solving_info_t sinfo;
  int ntlist = 1;
  int *ids;

  /*===============================================================================
   * 1. Convergence test
   *===============================================================================*/

  cs_real_t residu = 0.;

  /* Allocate a temporary array for cells or interior/boundary faces selection */
  BFT_MALLOC(varres, n_fields, cs_real_t);
  BFT_MALLOC(ids, n_fields, int);

  bool cved = true;
  int ivar = 0, nvar;

  for (int f_id = 0 ; f_id < n_fields ; f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (f->type & CS_FIELD_VARIABLE) {
      if (f->id == CS_F_(u)->id) {
        cs_field_get_key_struct(f, cs_field_key_id("solving_info"), &sinfo);
        cved = cved && (sinfo.l2residual  <= 1e-04);
        varres[ivar] = sinfo.l2residual;
        for (int isou = 0 ; isou < 3 ; isou++)
          residu += sinfo.l2residual;
        ids[ivar] = f_id;
        ivar++;
      } else {
        cs_field_get_key_struct(f, cs_field_key_id("solving_info"), &sinfo);
        cved = cved && (sinfo.l2residual  <= 1e-04);
        varres[ivar] = sinfo.l2residual;
        residu += sinfo.l2residual;
        ids[ivar] = f_id;
        ivar++;
      }
    }
  }

  nvar = ivar;
  BFT_REALLOC(varres, nvar, cs_real_t);
  BFT_REALLOC(ids, nvar, int);

  /* Warning: L2 time residual is computed after extra_operations */
  cved = cved && (ts->nt_cur > 1);

  if (cved) {
    ts->nt_max = ts->nt_cur;
    bft_printf("Converged at %d\n",ts->nt_max);
  }

  if (cs_glob_rank_id <= 0 && (ts->nt_cur-ts->nt_prev) == 1) {
    file = fopen("residuals.dat","w");
    fprintf(file,"#    Time step          Time     Residuals");
    for (ivar = 0 ; ivar < nvar ; ivar++) {
      int f_id = ids[ivar];
      if (ivar < nvar -1)
        fprintf(file,"          %s",cs_field_by_id(f_id)->name);
      else
        fprintf(file,"          %s\n",cs_field_by_id(f_id)->name);
    }
    fclose(file);
  }

  if (cs_glob_rank_id <= 0 && (ts->nt_cur == 1 || ts->nt_cur % ntlist == 0 || ts->nt_cur == ts->nt_max)) {
    file = fopen("residuals.dat","a");
    fprintf(file,"%14d %13.5e %13.5e", ts->nt_cur, ts->t_cur, residu);
    for (ivar = 0 ; ivar < nvar ; ivar++) {
      if (ivar < nvar -1)
        fprintf(file,"    %14.5e",varres[ivar]);
      else
        fprintf(file,"    %14.5e\n",varres[ivar]);
    }
    fclose(file);
  }

  BFT_FREE(varres);
  BFT_FREE(ids);

  /*===============================================================================
   * 2. Compute error
   *===============================================================================*/

  if (ts->nt_cur == ts->nt_max) {

    cs_real_t Per, Uer, Ver, Ter, Pcs, Ucs, Vcs, Tcs;

    /* Error */
    Per = 0.;
    Ter = 0.;
    Ver = 0.;
    Uer = 0.;

    cs_real_3_t *vel = (cs_real_3_t *)CS_F_(u)->val;

    int keysca = cs_field_key_id("scalar_id");

    for (int f_id = 0; f_id < cs_field_n_fields(); f_id++) {
      cs_field_t *f = cs_field_by_id(f_id);
      if (cs_field_get_key_int(f, keysca) == 1) {
        sca1 = f;
        break;
      }
    }

    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
      xx = cell_cen[iel][0];
      yy = cell_cen[iel][1];

      Pth = phi(CS_F_(p)->id, xx, yy);
      Pcs = CS_F_(p)->val[iel];

      Per += pow(Pcs-Pth,2)*volume[iel];

      Uth = phi(CS_F_(u)->id, xx, yy);
      Ucs = vel[iel][0];

      Uer += pow(Ucs-Uth,2)*volume[iel];

      Vth = phi(CS_F_(u)->id + 100, xx, yy);
      Vcs = vel[iel][1];

      Ver += pow(Vcs-Vth,2)*volume[iel];

      if (sca1 != NULL) {
        Tth = phi(sca1->id, xx, yy);
        Tcs = sca1->val[iel];
        Ter += pow(Tcs-Tth,2)*volume[iel];
      }
    }

    if (cs_glob_rank_id >= 0) {
      cs_parall_sum(1, CS_DOUBLE, &Per);
      cs_parall_sum(1, CS_DOUBLE, &Uer);
      cs_parall_sum(1, CS_DOUBLE, &Ver);
      cs_parall_sum(1, CS_DOUBLE, &Ter);
    }

    cs_real_t voltot = cs_glob_mesh_quantities->tot_vol;

    Per = sqrt(Per/voltot);
    Uer = sqrt(Uer/voltot);
    Ver = sqrt(Ver/voltot);
    if (sca1 != NULL)
      Ter = sqrt(Ter/voltot);

    cs_gnum_t gncel = n_cells;
    if(cs_glob_rank_id >= 0)
      cs_parall_counter(&gncel,1);

    /* Writing */
    if (cs_glob_rank_id <= 0) {
      file = fopen("L2error.dat","w");
      if (sca1 != NULL) {
        fprintf(file,"#       Ncells  P Error norm  U Error norm  V Error norm  T Error norm\n");
        fprintf(file,"%14lu %14.5e %14.5e %14.5e %14.5e\n",
                     gncel, Per, Uer, Ver, Ter);
      } else {
        fprintf(file,"#       Ncells  P Error norm  U Error norm  V Error norm\n");
        fprintf(file,"%14lu %14.5e %14.5e %14.5e\n",
                     gncel, Per, Uer, Ver);
      }
      fclose(file);
    }

    /* Writing CS profile */
    if (cs_glob_rank_id <= 0) {
      file = fopen("profile.dat","w");
      if (sca1 != NULL)
        fprintf(file,"# x U V P T\n");
      else
        fprintf(file,"# x U V P\n");
    }

    cs_lnum_t npoint = ceil(sqrt((double)gncel));
    cs_lnum_t iel1 = -999, iel;
    int irang1 = -999, irangv;
    cs_real_t xabs, xu, xv, xp, xt;

    for (cs_lnum_t ii = 0 ; ii < npoint ; ii++) {
      xx = (double)ii/(double)(npoint-1)*1.;
      yy = 0.8;
      zz = 0.;

      CS_PROCF(findpt, FINDPT)(&n_cells_ext, &n_cells, (cs_real_t *)cell_cen,
                               &xx, &yy, &zz,
                               &iel,    &irangv);

      if (iel != iel1 || irangv != irang1) {
        iel1 = iel;
        irang1 = irangv;

        if (cs_glob_rank_id == irangv) {
          iel--;
          xabs = cell_cen[iel][0];
          xu = vel[iel][0];
          xv = vel[iel][1];
          xp = CS_F_(p)->val[iel];
          if (sca1 != NULL)
            xt = sca1->val[iel];
          else
            xt = 0.;
        } else {
          xabs = 0.;
          xu = 0.;
          xv = 0.;
          xp = 0.;
          xt = 0.;
        }

        /* Broadcast to other ranks in parallel */
        if (cs_glob_rank_id >= 0) {
          cs_parall_bcast(irangv, 1, CS_DOUBLE, &xabs);
          cs_parall_bcast(irangv, 1, CS_DOUBLE, &xu);
          cs_parall_bcast(irangv, 1, CS_DOUBLE, &xv);
          cs_parall_bcast(irangv, 1, CS_DOUBLE, &xp);
          cs_parall_bcast(irangv, 1, CS_DOUBLE, &xt);
        }

        if (cs_glob_rank_id <= 0) {
          if (sca1 != NULL) {
          fprintf(file,"%14.5e %14.5e %14.5e %14.5e %14.5e\n",
                       xabs, xu, xv, xp, xt);
          } else {
          fprintf(file,"%14.5e %14.5e %14.5e %14.5e\n",
                       xabs, xu, xv, xp);
          }
        }
      }
    }

    if (cs_glob_rank_id <= 0)
      fclose(file);
  }

  /*===============================================================================
   * 3. Save analytical solution
   *===============================================================================*/

  if (cs_glob_rank_id <= 0 && ts->nt_cur == ts->nt_max) {
    file = fopen("theory.dat","w");
    fprintf(file,"# x Uth Vth Pth Tth\n");

    /* Analytical solution */
    cs_lnum_t npoint = 100;

    for (cs_lnum_t ii = 0 ; ii < npoint ; ii++) {
      xx = (double)ii/(double)(npoint-1)*1.;
      yy = 0.8;
      Pth = phi(CS_F_(p)->id, xx, yy);
      Uth = phi(CS_F_(u)->id, xx, yy);
      Vth = phi(CS_F_(u)->id + 100, xx, yy);

      if (sca1 != NULL) {
        Tth = phi(sca1->id, xx, yy);
        fprintf(file,"%14.5e %14.5e %14.5e %14.5e %14.5e\n",
                     xx, Uth, Vth, Pth, Tth);
      } else {
        fprintf(file,"%14.5e %14.5e %14.5e %14.5e\n",
                     xx, Uth, Vth, Pth);
      }
    }
    fclose(file);
  }
}

END_C_DECLS
