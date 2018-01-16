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
!
!-------------------------------------------------------------------------------

subroutine cs_f_user_boundary_conditions &
 ( nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   dt     ,                                                       &
   rcodcl )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ihmpre
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use atincl
use atsoil
use ctincl
use cs_fuel_incl
use mesh
use field
use turbomachinery

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          icodcl(nfabor,nvarcl)
integer          itrifb(nfabor), itypfb(nfabor)
integer          izfppp(nfabor)

double precision dt(ncelet)
double precision rcodcl(nfabor,nvarcl,3)

! Local variables

integer          nlelt, ilelt, ifac
double precision xx, yy

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

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

allocate(lstelt(nfabor))  ! temporary array for boundary faces selection

!===============================================================================
! Assign boundary conditions to boundary faces here

! For each subset:
! - use selection criteria to filter boundary faces of a given subset
! - loop on faces from a subset
!   - set the boundary condition for each face
!===============================================================================

call getfbr('Symmetric_walls',nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  itypfb(ifac) = isymet

enddo


call getfbr('not Symmetric_walls',nlelt, lstelt)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  itypfb(ifac) = iindef

  xx = cdgfbo(1,ifac)
  yy = cdgfbo(2,ifac)

  icodcl(ifac,ipr) = 3
  rcodcl(ifac,ipr,3) = 0.d0

  icodcl(ifac,iu) = 1
  icodcl(ifac,iv) = 1
  icodcl(ifac,iw) = 1
  if (nscal.eq.1) icodcl(ifac,isca(1)) = 1

  rcodcl(ifac,iu ,1) = phi (ivarfl(iu) ,xx,yy)
  rcodcl(ifac,iv ,1) = phi (ivarfl(iu)+100 ,xx,yy)
  rcodcl(ifac,iw ,1) = phi (ivarfl(iu)+1000 ,xx,yy)
  if (nscal.eq.1) rcodcl(ifac,isca(1) ,1) = phi (ivarfl(isca(1)) ,xx,yy)

enddo

deallocate(lstelt)  ! temporary array for boundary faces selection

return
end subroutine cs_f_user_boundary_conditions
