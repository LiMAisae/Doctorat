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
!> \param[in]     ipart         number of the post-processing mesh (< 0 or > 0)
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     nvlsta        number of Lagrangian statistical variables
!> \param[in]     ncelps        number of cells in post-processing mesh
!> \param[in]     nfacps        number of interior faces in post-process. mesh
!> \param[in]     nfbrps        number of boundary faces in post-process. mesh
!> \param[in]     itypps        global presence flag (0 or 1) for cells (1),
!>                              interior faces (2), or boundary faces (3) in
!>                              post-processing mesh
!> \param[in]     lstcel        list of cells in post-processing mesh
!> \param[in]     lstfac        list of interior faces in post-processing mesh
!> \param[in]     lstfbr        list of boundary faces in post-processing mesh
!_______________________________________________________________________________

subroutine usvpst &
 ( ipart  ,                                                       &
   nvar   , nscal  , nvlsta ,                                     &
   ncelps , nfacps , nfbrps ,                                     &
   itypps ,                                                       &
   lstcel , lstfac , lstfbr )

!===============================================================================

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use pointe
use entsor
use optcal
use numvar
use parall
use period
use mesh
use field
use post
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          ipart
integer          nvar,   nscal , nvlsta
integer          ncelps, nfacps, nfbrps

integer          itypps(3)
integer          lstcel(ncelps), lstfac(nfacps), lstfbr(nfbrps)


! Local variables

integer          iel, iloc
integer          idimt, ii
logical          ientla, ivarpr
double precision xx, yy
double precision rvoid(1)

double precision, dimension(:), allocatable :: scel
double precision, dimension(:,:), allocatable :: vcel

character*32     namevr

integer          intpst
data             intpst /0/
save             intpst

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
! Increment call counter once per time step (possibly used in some tests)
!===============================================================================

if (ipart .eq. -1) then
  intpst = intpst + 1
endif

if (ipart .eq. -1) then

  allocate(vcel(3,ncelps),scel(ncelps))

  ! Pressure
  ! ********

  ! Initialize variable name
  do ii = 1, 32
    namevr(ii:ii) = ' '
  enddo

  ! Variable name
  namevr = 'p_0'

  do iloc = 1, ncelps
    iel = lstcel(iloc)

    xx = xyzcen(1,iel)
    yy = xyzcen(2,iel)

    scel(iloc) = phi(ivarfl(ipr),xx,yy)
  enddo

  idimt = 1        ! 1: scalar, 3: vector, 6/9: symm/non-symm tensor
  ientla = .true.  ! dimension 1 here, so no effect
  ivarpr = .false. ! defined on the work array, not on the parent

  ! Output values; as we have no face values, we can pass a
  ! trivial array rvoid instead of trafac and trafbr.
  call post_write_var(ipart, namevr, idimt, ientla, ivarpr,  &
                      ntcabs, ttcabs, scel, rvoid, rvoid)

  ! Velocity
  ! ********

  ! Initialize variable name
  do ii = 1, 32
    namevr(ii:ii) = ' '
  enddo

  ! Variable name
  namevr = 'U_0'

  do iloc = 1, ncelps
    iel = lstcel(iloc)

    xx = xyzcen(1,iel)
    yy = xyzcen(2,iel)

    vcel(1,iloc) = phi(ivarfl(iu),xx,yy)
    vcel(2,iloc) = phi(ivarfl(iu)+100,xx,yy)
    vcel(3,iloc) = phi(ivarfl(iu)+1000,xx,yy)
  enddo

  idimt = 3        ! 1: scalar, 3: vector, 6/9: symm/non-symm tensor
  ientla = .true.  ! Values are interlaced
  ivarpr = .false. ! defined on the work array, not on the parent

  ! Output values; as we have no face values, we can pass a
  ! trivial array rvoid instead of trafac and trafbr.
  call post_write_var(ipart, namevr, idimt, ientla, ivarpr,  &
                      ntcabs, ttcabs, vcel, rvoid, rvoid)


  ! Scalar theta
  ! ************

  ! Initialize variable name
  do ii = 1, 32
    namevr(ii:ii) = ' '
  enddo

  ! Variable name
  namevr = 'Theta_0'

  do iloc = 1, ncelps
    iel = lstcel(iloc)

    xx = xyzcen(1,iel)
    yy = xyzcen(2,iel)

    scel(iloc) = phi(ivarfl(isca(1)),xx,yy)
  enddo

  idimt = 1        ! 1: scalar, 3: vector, 6/9: symm/non-symm tensor
  ientla = .true.  ! dimension 1 here, so no effect
  ivarpr = .false. ! defined on the work array, not on the parent

  ! Output values; as we have no face values, we can pass a
  ! trivial array rvoid instead of trafac and trafbr.
  call post_write_var(ipart, namevr, idimt, ientla, ivarpr,  &
                      ntcabs, ttcabs, scel, rvoid, rvoid)
  deallocate(vcel,scel)

endif ! end of test on post-processing mesh number

return

end subroutine usvpst
