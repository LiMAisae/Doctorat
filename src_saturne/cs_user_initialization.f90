!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

!===============================================================================
! Purpose:
! -------

!> \file cs_user_initialization.f90
!>
!> \brief Initialize variables
!>
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine cs_user_f_initialization &
 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use atincl
use ctincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use mesh
use field
use turbomachinery
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          iflmas, iflmab
integer          iflmb0
integer          init, inc, itypfl
integer          nswrgp, imligp, iwarnp

double precision dt(ncelet)

! Local variables

integer          iel
double precision xx, yy

integer, allocatable, dimension(:) :: lstelt

double precision, dimension(:,:), pointer :: cvar_vel
double precision, dimension(:), pointer :: cvar_pr, cvar_scal

double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: brom, crom
double precision, dimension(:,:), pointer :: coefau, cofafu
double precision, dimension(:,:,:), pointer :: coefbu, cofbfu
double precision  epsrgp, climgp, extrap

type(var_cal_opt) :: vcopt

!===============================================================================
! Interfaces
!===============================================================================

interface

  function phi(f_id, x, y) result(phi_) &
    bind(C, name='phi')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(kind=c_int), value :: f_id
    real(kind=c_double), value :: x, y
    real(kind=c_double) :: phi_
  end function phi

end interface

!===============================================================================
! Initialization
!===============================================================================

allocate(lstelt(ncel)) ! temporary array for cells selection

call field_get_val_s(ivarfl(ipr), cvar_pr)
call field_get_val_v(ivarfl(iu), cvar_vel)
call field_get_val_s(ivarfl(isca(1)), cvar_scal)

do iel = 1, ncel

  xx = xyzcen(1,iel)
  yy = xyzcen(2,iel)

  cvar_pr(iel) = phi(ivarfl(ipr),xx,yy)
  cvar_vel(1,iel) = phi(ivarfl(iu) ,xx,yy)
  cvar_vel(2,iel) = phi(ivarfl(iu)+100 ,xx,yy)
  cvar_vel(3,iel) = phi(ivarfl(iu)+1000 ,xx,yy)

  if (nscal.eq.1)  cvar_scal(iel) = phi(ivarfl(isca(1)),xx,yy)

enddo

!===============================================================================
! BS. mass fluxes imposed
!===============================================================================

  call field_get_key_int(ivarfl(iu), kimasf, iflmas)
  call field_get_key_int(ivarfl(iu), kbmasf, iflmab)

  ! Pointers to the mass fluxes
  call field_get_val_s(iflmas, imasfl)
  call field_get_val_s(iflmab, bmasfl)

  call field_get_val_s(icrom, crom)
  call field_get_val_s(ibrom, brom)

  call field_get_coefa_v(ivarfl(iu), coefau)
  call field_get_coefb_v(ivarfl(iu), coefbu)
  call field_get_coefaf_v(ivarfl(iu), cofafu)
  call field_get_coefbf_v(ivarfl(iu), cofbfu)

  itypfl = 1
  init   = 1
  inc    = 1
  iflmb0 = 1

  call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt)

  nswrgp = vcopt%nswrgr
  imligp = vcopt%imligr
  iwarnp = vcopt%iwarni
  epsrgp = vcopt%epsrgr
  climgp = vcopt%climgr
  extrap = vcopt%extrag

  call inimav                                                     &
  !==========
 ( ivarfl(iu)      , itypfl ,                                     &
   iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
   iwarnp ,                                                       &
   epsrgp , climgp ,                                              &
   crom, brom,                                                    &
   cvar_vel    ,                                                  &
   coefau , coefbu ,                                              &
   imasfl , bmasfl )

deallocate(lstelt) ! temporary array for cells selection

return
end subroutine cs_user_f_initialization
