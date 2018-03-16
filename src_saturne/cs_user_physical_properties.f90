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
integer          ith, iscal, ii, ifcvsl
double precision vara, varb, varc, varam, varbm, varcm, vardm
double precision                   varal, varbl, varcl, vardl
double precision                   varac, varbc

double precision, dimension(:), pointer :: coefap, coefbp
double precision, dimension(:), pointer :: bfpro_rom, cpro_rom, cproa_rom,cproaa_rom
double precision, dimension(:), pointer :: cpro_viscl, cpro_vscalt, cpro_cp
double precision, dimension(:), pointer :: cvar_scalt, cpro_drhodt,cpro_ts_scalt

double precision xx, yy, xrtp, rho, drhods, rhoaa, rhoa

!===============================================================================
! Interfaces
!===============================================================================

interface

  function phi(f_id, x, y, t) result(phi_) &
    bind(C, name='phi')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(kind=c_int), value :: f_id
    real(kind=c_double), value :: x, y, t
    real(kind=c_double) :: phi_
  end function phi

end interface

!===============================================================================

if (irovar.eq.1) then

  izero = 0
  ivart = isca(1)

  if (nscal.eq.1) call field_get_val_s(ivarfl(ivart), cvar_scalt)
  if (nscal.eq.1) call field_get_coefa_s(ivarfl(ivart), coefap)
  if (nscal.eq.1) call field_get_coefb_s(ivarfl(ivart), coefbp)
  call field_get_val_s(icrom, cpro_rom)
  call field_get_val_prev_s(icrom, cproa_rom)
  call field_get_val_s(ibrom, bfpro_rom)

  if (idilat.ge.4) then
    call field_get_val_s(icroaa, cproaa_rom)
  endif

  if (nscal.eq.1) then
    do iel = 1, ncel
      xrtp = cvar_scalt(iel)
      cpro_rom(iel) = 1.d0/(xrtp/ro2 + (1.d0-xrtp)/ro1)
    enddo
  else
    do iel = 1, ncel
      xx = xyzcen(1,iel)
      yy = xyzcen(2,iel)
      cpro_rom(iel) = phi(-1,xx,yy,ttcabs)
    enddo
  endif
  if (ntcabs.eq.1.and.idilat.ge.2) then
    do iel = 1, ncel
      cproa_rom(iel) = cpro_rom(iel)!phi(-1,xx,yy,0.d0)
      if (idilat.ge.4) cproaa_rom(iel) = cpro_rom(iel)
    enddo
  endif


  mbrom = 1
  print*,"coefap_sca, coefbp_sca", coefap(50), coefbp(50)
  if (nscal.eq.1) then
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      xrtp = coefap(ifac) + cvar_scalt(iel)*coefbp(ifac)
      bfpro_rom(ifac) = 1.d0/(xrtp/ro2 + (1.d0-xrtp)/ro1)
    enddo
  else
    do ifac = 1, nfabor
      xx = cdgfbo(1,ifac)
      yy = cdgfbo(2,ifac)
      bfpro_rom(ifac) = phi(-1,xx,yy,ttcabs)
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
        drhods = (1.d0/ro1 - 1.d0/ro2)*rho*rho
        cpro_drhodt(iel) = drhods / rhoaa * cpro_ts_scalt(iel)
        cpro_ts_scalt(iel) = 0.d0
      enddo
    else
      do iel = 1, ncel
        xx = xyzcen(1,iel)
        yy = xyzcen(2,iel)
        rho    = cpro_rom(iel)
        rhoaa  = cproaa_rom(iel)
        drhods = (1.d0/ro1 - 1.d0/ro2)*rho*rho
        !phi(-10) = q_z + div(K grad(Z))
        cpro_drhodt(iel) = drhods / rhoaa * phi(-50,xx,yy,ttcabs)*volume(iel)
      enddo
    endif
  endif

endif

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    DONNEES DE CALCUL INCOHERENTES                          ',/,&
'@                                                            ',/,&
'@      usipph indique que la chaleur specifique est uniforme ',/,&
'@        ICP = ',I10   ,' alors que                          ',/,&
'@      usphyv impose une chaleur specifique variable.        ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier usipph ou usphyv.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    DONNEES DE CALCUL INCOHERENTES                          ',/,&
'@                                                            ',/,&
'@    Pour le scalaire ', i10                                  ,/,&
'@      la diffusivite est uniforme alors que                 ',/,&
'@      usphyv impose une diffusivite variable.               ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier usipsu ou usphyv.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    APPEL A csexit DANS LE SOUS PROGRAMME usphyv            ',/,&
'@                                                            ',/,&
'@    La variable dont dependent les proprietes physiques ne  ',/,&
'@      semble pas etre une variable de calcul.               ',/,&
'@    En effet, on cherche a utiliser la temperature alors que',/,&
'@      ISCALT = ',I10                                         ,/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Verifier le codage de usphyv (et le test lors de la     ',/,&
'@      definition de IVART).                                 ',/,&
'@    Verifier la definition des variables de calcul dans     ',/,&
'@      usipsu. Si un scalaire doit jouer le role de la       ',/,&
'@      temperature, verifier que ISCALT a ete renseigne.     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:  stop when computing physical quantities',/,       &
'@    =======',/,                                                 &
'@    Inconsistent calculation data',/,                           &
'@',/,                                                            &
'@      usipph specifies that the specific heat is uniform',/,    &
'@        icp = ',i10   ,' while',/,                              &
'@      usphyv prescribes a variable specific heat.',/,           &
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Modify usipph or usphyv.',/,                                &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 1010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:  stop when computing physical quantities',/,       &
'@    =======',/,                                                 &
'@    Inconsistent calculation data',/,                           &
'@',/,                                                            &
'@    For scalar', i10,/,                                         &
'@      the diffusivity is uniform while',/,                      &
'@      usphyv prescribes a variable diffusivity.',/,             &
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Modify usipsu or usphyv.',/,                                &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 9010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:  stop when computing physical quantities',/,       &
'@    =======',/,                                                 &
'@',/,                                                            &
'@    The variable on which physical properties depend does',/,   &
'@      seem to be a calculation variable.',/,                    &
'@    Indeed, we are trying to use the temperature while',/,      &
'@      iscalt = ',i10                                  ,/,&
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Check the programming in usphyv (and the test when',/,      &
'@      defining ivart).',/,                                      &
'@    Check the definition of calculation variables in',/,        &
'@      usipsu. If a scalar should represent the,',/,             &
'@      temperature, check that iscalt has been defined',/,       &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

#endif

!----
! End
!----

return
end subroutine usphyv
