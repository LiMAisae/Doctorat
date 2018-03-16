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
#include "cs_blas.h"
#include "cs_field.h"
#include "cs_field_pointer.h"
#include "cs_field_operator.h"
#include "cs_math.h"
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
  FILE *file = NULL, *f1 = NULL, *f2 = NULL;
  const cs_lnum_t n_cells     = cs_glob_mesh->n_cells;
  const cs_lnum_t n_cells_ext = cs_glob_mesh->n_cells_with_ghosts;
  const cs_real_t *volume = cs_glob_mesh_quantities->cell_vol;
  const cs_real_3_t *cell_cen = (const cs_real_3_t *)cs_glob_mesh_quantities->cell_cen;
  cs_time_step_t *ts = cs_get_glob_time_step();
  cs_real_t *Per, *Pthvol, *Ter, *Tthvol, xx, yy, zz, Per1, Per2, Perror;
  cs_real_t Ter1, Ter2, Terror, Uerror, Uthx, Uthy, Uerx, Uery;
  cs_real_t Pth, Uth, Vth, Tth;
  cs_real_t *Uer, *Uthvol, *Ver, *Vthvol;
  cs_field_t *sca1 = NULL;
  cs_real_t uf = 0.5;

  /*===============================================================================
   * 2. Compute error
   *===============================================================================*/

  BFT_MALLOC(Per, n_cells_ext, cs_real_t);
  BFT_MALLOC(Pthvol, n_cells_ext, cs_real_t);
  BFT_MALLOC(Ter, n_cells_ext, cs_real_t);
  BFT_MALLOC(Tthvol, n_cells_ext, cs_real_t);
  BFT_MALLOC(Uer, n_cells_ext, cs_real_t);
  BFT_MALLOC(Uthvol, n_cells_ext, cs_real_t);
  BFT_MALLOC(Ver, n_cells_ext, cs_real_t);
  BFT_MALLOC(Vthvol, n_cells_ext, cs_real_t);

  int keysca = cs_field_key_id("scalar_id");

  for (int f_id = 0; f_id < cs_field_n_fields(); f_id++) {
    cs_field_t *f = cs_field_by_id(f_id);
    if (cs_field_get_key_int(f, keysca) == 1) {
      sca1 = f;
      break;
    }
  }

  cs_real_3_t *vel = (cs_real_3_t *)CS_F_(u)->val;

  if (ts->t_cur >= 0.4096) {
    ts->nt_max = ts->nt_cur;
  }

  if (ts->nt_cur == ts->nt_max) {
    for (cs_lnum_t iel = 0 ; iel < n_cells ; iel++) {
      xx = cell_cen[iel][0];
      yy = cell_cen[iel][1];

      if (fabs(uf - 0.5) < cs_math_epzero) {
        Pth = phi(CS_F_(p)->id, xx, yy, ts->t_cur) - 0.375;
      } else {
        Pth = phi(CS_F_(p)->id, xx, yy, ts->t_cur);
      }

      Per[iel] = (CS_F_(p)->val[iel]-Pth)*volume[iel];
      Pthvol[iel] = phi(CS_F_(p)->id, xx, yy, ts->t_cur)*volume[iel];

      if (sca1 != NULL) {
        Tth = phi(sca1->id, xx, yy, ts->t_cur);
        Ter[iel] = (sca1->val[iel]-Tth)*volume[iel];
        Tthvol[iel] = phi(sca1->id, xx, yy, ts->t_cur)*volume[iel];
      }

      Uer[iel] = (vel[iel][0]-phi(CS_F_(u)->id, xx, yy, ts->t_cur))*volume[iel];
      Ver[iel] = (vel[iel][1]-phi(CS_F_(u)->id+100, xx, yy, ts->t_cur))*volume[iel];

      Uthvol[iel] = phi(CS_F_(u)->id, xx, yy, ts->t_cur)*volume[iel];
      Vthvol[iel] = phi(CS_F_(u)->id+100, xx, yy, ts->t_cur)*volume[iel];
    }

    Per1 = sqrt(cs_gres(n_cells, volume, Per, Per));
    Per2 = sqrt(cs_gres(n_cells, volume, Pthvol, Pthvol));

    Perror = Per1 / Per2;

    if (sca1 != NULL) {
      Ter1 = sqrt(cs_gres(n_cells, volume, Ter, Ter));
      Ter2 = sqrt(cs_gres(n_cells, volume, Tthvol, Tthvol));
      Terror = Ter1 / Ter2;
    }

    Uerx = cs_gres(n_cells, volume, Uer, Uer);
    Uery = cs_gres(n_cells, volume, Ver, Ver);
    Uthx = cs_gres(n_cells, volume, Uthvol, Uthvol);
    Uthy = cs_gres(n_cells, volume, Vthvol, Vthvol);

    Uerror = sqrt((Uerx + Uery)/(Uthx + Uthy));

    cs_gnum_t gncel = n_cells;
    if(cs_glob_rank_id >= 0)
      cs_parall_counter(&gncel,1);

    /* Writing */
    if (cs_glob_rank_id <= 0) {
      file = fopen("L2error.dat","w");
      if (sca1 != NULL) {
        fprintf(file,"#       Ncells  P Error norm  U Error norm  T Error norm\n");
        fprintf(file,"%14lu %14.5e %14.5e %14.5e\n", gncel, Perror, Uerror, Terror);
      } else {
        fprintf(file,"#       Ncells  P Error norm  U Error norm\n");
        fprintf(file,"%14lu %14.5e %14.5e\n", gncel, Perror, Uerror);
      }
      fclose(file);
    }

    /* Writing CS profile */
    if (cs_glob_rank_id <= 0) {
      f1 = fopen("profile.dat","w");
      f2 = fopen("solth_discrete.dat","w");
      if (sca1 != NULL) {
        fprintf(f1,"# x U P T\n");
        fprintf(f2,"# x U P T\n");
      } else {
        fprintf(f1,"# x U P\n");
        fprintf(f2,"# x U P\n");
      }
    }

    cs_lnum_t npoint = ceil(sqrt((double)gncel));
    cs_lnum_t iel1 = -999, iel;
    int irang1 = -999, irangv;
    cs_real_t xabs, xu, xv, xp, xt;

    for (cs_lnum_t ii = 0 ; ii < npoint ; ii++) {
      xx = (double)ii/(double)(npoint-1)*2.;
      yy = 1.8;
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

          Pth = phi(CS_F_(p)->id, xabs, 1.8, ts->t_cur);
          Uth = phi(CS_F_(u)->id, xabs, 1.8, ts->t_cur);
          Vth = phi(CS_F_(u)->id+100, xabs, 1.8, ts->t_cur);

          if (sca1 != NULL)
            Tth = phi(sca1->id, xabs, 1.8, ts->t_cur);
          else
            Tth = 0.;

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
          cs_parall_bcast(irangv, 1, CS_DOUBLE, &Uth);
          cs_parall_bcast(irangv, 1, CS_DOUBLE, &Vth);
          cs_parall_bcast(irangv, 1, CS_DOUBLE, &Pth);
          cs_parall_bcast(irangv, 1, CS_DOUBLE, &Tth);
        }

        if (cs_glob_rank_id <= 0) {
          if (sca1 != NULL) {
                  fprintf(f1,"%14.5e %14.5e %14.5e %14.5e %14.5e\n",
                             xabs, xu, xv, xp, xt);
                  fprintf(f2,"%14.5e %14.5e %14.5e %14.5e %14.5e\n",
                             xabs, Uth, Vth, Pth, Tth);
          } else {
                  fprintf(f1,"%14.5e %14.5e %14.5e %14.5e\n",
                             xabs, xu, xv, xp);
                  fprintf(f2,"%14.5e %14.5e %14.5e %14.5e\n",
                             xabs, Uth, Vth, Pth);
          }
        }
      }
    }

    if (cs_glob_rank_id <= 0) {
      fclose(f1);
      fclose(f2);
    }
  }

  /*===============================================================================
   * 3. Save analytical solution
   *===============================================================================*/

  bool soltheo = true;

  if (cs_glob_rank_id <= 0 && ts->nt_cur == ts->nt_max && soltheo) {
    file = fopen("theory.dat","w");
    if (sca1 != NULL)
      fprintf(file,"# x Uth Vth Pth Tth\n");
    else
      fprintf(file,"# x Uth Vth Pth\n");

    /* Analytical solution */
    cs_lnum_t npoint = 400;

    for (cs_lnum_t ii = 0 ; ii < npoint ; ii++) {
      xx = (double)ii/(double)(npoint-1)*2.;
      yy = 1.8;
      Pth = phi(CS_F_(p)->id, xx, yy, ts->t_cur);
      Uth = phi(CS_F_(u)->id, xx, yy, ts->t_cur);
      Vth = phi(CS_F_(u)->id + 100, xx, yy, ts->t_cur);

      if (sca1 != NULL) {
        Tth = phi(sca1->id, xx, yy, ts->t_cur);
        fprintf(file,"%14.5e %14.5e %14.5e %14.5e %14.5e\n",
                     xx, Uth, Vth, Pth, Tth);
      } else {
        fprintf(file,"%14.5e %14.5e %14.5e %14.5e\n",
                     xx, Uth, Vth, Pth);
      }
    }
    fclose(file);
  }

  BFT_FREE(Per);
  BFT_FREE(Ter);
  BFT_FREE(Uer);
  BFT_FREE(Ver);
  BFT_FREE(Pthvol);
  BFT_FREE(Uthvol);
  BFT_FREE(Vthvol);
  BFT_FREE(Tthvol);
}

END_C_DECLS
