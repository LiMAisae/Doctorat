!-------------------------------------------------------------------------------

!                      Code_Saturne version 4.0-alpha
!                      --------------------------
! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     mbrom         indicator of filling of romb array
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine usphyv &
 ( nvar   , nscal  ,                                              &
   mbrom  ,                                                       &
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
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use field
use mesh
use user_module

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          mbrom

double precision dt(ncelet)

! Local variables

integer          ivart, iel, ifac, izero

double precision, dimension(:), pointer :: coefap, coefbp
double precision, dimension(:), pointer :: bfpro_rom, cpro_rom, cproa_rom,cproaa_rom
double precision, dimension(:), pointer :: cpro_drhodt,cpro_ts_scalt
double precision, dimension(:), pointer :: cvar_scalt

double precision xx, yy, ro1, ro2, xrtp, rho, drhods, rhoaa, rhoa

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

if (irovar.eq.1) then

  izero = 0

  ! Position of variables, coefficients
  ! -----------------------------------

  ! --- Number of the helium variable
  !       (and of its boundary conditions)

  ivart = isca(1)

  call field_get_val_s(ivarfl(ivart), cvar_scalt)

  ! --- Position of boundary conditions for variable 'ivart'

  call field_get_coefa_s(ivarfl(ivart), coefap)
  call field_get_coefb_s(ivarfl(ivart), coefbp)

  ! --- Rank of density
  !     in 'propce', physical properties at element centers:       'ipcrom'
  !     in 'propfb', physical properties at boundary face centers: 'ipbrom'

  call field_get_val_s(icrom, cpro_rom)
  call field_get_val_prev_s(icrom, cproa_rom)
  call field_get_val_s(ibrom, bfpro_rom)

  if (idilat.ge.4) then
    call field_get_val_s(icroaa, cproaa_rom)
  endif

  ! Density at cell centers
  !------------------------

  ro1 = 1.d0
  ro2 = omega*ro1

  if (nscal.eq.1) then
    do iel = 1, ncel
      xrtp = cvar_scalt(iel)
      cpro_rom(iel) = 1.d0/(xrtp/ro1 + (1.d0-xrtp)/ro2)
    enddo
  else
    do iel = 1, ncel
      xx = xyzcen(1,iel)
      yy = xyzcen(2,iel)
      cpro_rom(iel) = phi(-1,xx,yy)
    enddo
  endif

  ! reference density
  ro0 = ro1

  if (ntcabs.eq.1.and.idilat.ge.4) then
    do iel = 1, ncel
      cproaa_rom(iel) = cpro_rom(iel)
      cproa_rom(iel) = cpro_rom(iel)
    enddo
  endif

  ! Density at boundary faces
  !---------------------------

  mbrom = 1

  ! Boundary condition coefficients
  call field_get_coefa_s(ivarfl(ivart), coefap)
  call field_get_coefb_s(ivarfl(ivart), coefbp)

  if (nscal.eq.1) then
    do ifac = 1, nfabor
      ! ifabor(ifac) is the cell adjacent to the boundary face
      iel = ifabor(ifac)
      xrtp = coefap(ifac) + cvar_scalt(iel)*coefbp(ifac)
      bfpro_rom(ifac) = 1.d0/(xrtp/ro1 + (1.d0-xrtp)/ro2)
    enddo
  else
    do ifac = 1, nfabor
      xx = cdgfbo(1,ifac)
      yy = cdgfbo(2,ifac)
      bfpro_rom(ifac) = phi(-1,xx,yy)
    enddo
  endif

  !Cas idilat = 4 ou 5 : calcul de Drho/Dt

  if (idilat.ge.4) then
    call field_get_val_s(iustdy(itsrho), cpro_drhodt)
    if (nscal.eq.1) then
      call field_get_val_s(iustdy(1), cpro_ts_scalt)
      do iel = 1, ncel
        rho    = cpro_rom(iel)
        rhoa   = cproa_rom(iel)
        rhoaa  = cproaa_rom(iel)
        drhods = (1.d0/ro2 - 1.d0/ro1)*rho*rho
        ! Drho/Dt = drho/ds * Dh/Dt
        !         = drho/ds * 1/rho^n-2 [div(rho D grad h)]
        cpro_drhodt(iel) = drhods / rhoaa * cpro_ts_scalt(iel)
        cpro_ts_scalt(iel) = 0.d0
      enddo
    else
      do iel = 1, ncel
        xx = xyzcen(1,iel)
        yy = xyzcen(2,iel)
        rho    = cpro_rom(iel)
        rhoaa  = cproaa_rom(iel)
        drhods = 1.d0/ro2 - 1.d0/ro1
        cpro_drhodt(iel) = drhods / rhoaa * phi(-15,xx,yy)*volume(iel)
      enddo
    endif
  endif

endif

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine usphyv
