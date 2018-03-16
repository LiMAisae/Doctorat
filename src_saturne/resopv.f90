!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2017 EDF S.A.
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
! Function:
! ---------

!> \file resopv.f90
!>
!> \brief This subroutine performs the pressure correction step of the Navier
!> Stokes equations for incompressible or slightly compressible flows for
!> the coupled velocity components solver.
!>
!> This function solves the following Poisson equation on the pressure:
!> \f[
!>     D \left( \Delta t, \delta p \right) =
!> \divs \left( \rho \vect{\widetilde{u}}\right)
!>     - \Gamma^n
!>     + \dfrac{\rho^n - \rho^{n-1}}{\Delta t}
!> \f]
!> The mass flux is then updated as follows:
!> \f[
!>  \dot{m}^{n+1}_\ij = \dot{m}^{n}_\ij
!>                    - \Delta t \grad_\fij \delta p \cdot \vect{S}_\ij
!> \f]
!>
!> Remarks:
!> - an iterative process is used to solve the Poisson equation.
!> - if the coefficient arak is set to 1, the the Rhie & Chow filter is
!>   activated.
!>
!> Please refer to the
!> <a href="../../theory.pdf#resopv"><b>resopv</b></a>
!> section of the theory guide for more informations.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     nfbpcd        number of faces with condensation source term
!> \param[in]     ncmast        number of cells with condensation source terms
!> \param[in]     icetsm        index of cells with mass source term
!> \param[in]     ifbpcd        index of faces with condensation source term
!> \param[in]     ltmast        index of cells with condensation source terms
!> \param[in]     isostd        indicator of standard outlet and index
!>                               of the reference outlet face
!> \param[in]     dt            time step (per cell)
!> \param[in]     vel           velocity
!> \param[in]     coefav        boundary condition array for the variable
!>                               (explicit part)
!> \param[in]     coefbv        boundary condition array for the variable
!>                               (implicit part)
!> \param[in]     coefa_dp      boundary conditions for the pressure increment
!> \param[in]     coefb_dp      boundary conditions for the pressure increment
!> \param[in]     smacel        variable value associated to the mass source
!>                               term (for ivar=ipr, smacel is the mass flux
!>                               \f$ \Gamma^n \f$)
!> \param[in]     spcond        variable value associated to the condensation
!>                               source term (for ivar=ipr, spcond is the flow rate
!>                               \f$ \Gamma_{s,cond}^n \f$)
!> \param[in]     svcond        variable value associated to the condensation
!>                              source term (for ivar=ipr, svcond is the flow rate
!>                              \f$ \Gamma_{v, cond}^n \f$)
!> \param[in]     frcxt         external forces making hydrostatic pressure
!> \param[in]     dfrcxt        variation of the external forces
!>                              composing the hydrostatic pressure
!> \param[in]     tpucou        non scalar time step in case of
!>                               velocity pressure coupling
!> \param[in]     trav          right hand side for the normalizing
!>                               the residual
!> \param[in]     viscf         visc*surface/dist aux faces internes
!> \param[in]     viscb         visc*surface/dist aux faces de bord
!> \param[in]     dpvar         tableau de travail pour increment
!> \param[in]     tslagr        coupling term for the Lagrangian module
!> \param[in]     trava         tableau de travail pour couplage
!_______________________________________________________________________________

subroutine resopv &
 ( nvar   , ncesmp , nfbpcd , ncmast ,                            &
   icetsm , ifbpcd , ltmast , isostd ,                            &
   dt     , vel    ,                                              &
   coefav , coefbv , coefa_dp        , coefb_dp ,                 &
   smacel , spcond , svcond ,                                     &
   frcxt  , dfrcxt , tpucou , trav   ,                            &
   viscf  , viscb  ,                                              &
   dpvar  , tslagr ,                                              &
   trava  )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use numvar
use entsor
use cstphy
use cstnum
use optcal
use pointe, only: itypfb, b_head_loss, gamcav, dgdpca
use albase
use parall
use period
use ppincl, only: icondv
use lagran
use cplsat
use mesh
use field
use field_operator
use cavitation
use vof
use cs_f_interfaces
use cs_c_bindings
use cs_tagms, only:s_metal

!===============================================================================

implicit none

! Arguments

integer          nvar
integer          ncesmp, nfbpcd, ncmast

integer          icetsm(ncesmp), ifbpcd(nfbpcd)
integer          ltmast(ncelet)
integer          isostd(nfabor+1)

double precision, dimension (1:ncelet), target :: dt
double precision smacel(ncesmp,nvar), spcond(nfbpcd,nvar)
double precision svcond(ncelet,nvar)
double precision frcxt(3,ncelet), dfrcxt(3,ncelet)
double precision, dimension (1:6,1:ncelet), target :: tpucou
double precision trav(3,ncelet)
double precision viscf(nfac), viscb(nfabor)
double precision dpvar(ncelet)
double precision tslagr(ncelet,*)
double precision trava(ndim,ncelet)
double precision coefav(3  ,nfabor)
double precision coefbv(3,3,nfabor)
double precision vel   (3  ,ncelet)
double precision coefa_dp(nfabor)
double precision coefb_dp(nfabor)

! Local variables

character(len=80) :: chaine
integer          lchain
integer          iccocg, inc   , iprev, init  , isym
integer          ii, jj, iel   , ifac  , ifac0 , iel0
integer          nswmpr
integer          isweep, niterf
integer          iflmb0, ifcsor
integer          nswrgp, imligp, iwarnp
integer          iflmas, iflmab
integer          idiffp, iconvp, ndircp
integer          iinvpe, indhyd
integer          itypfl
integer          isou  , ibsize, iesize
integer          imucpp, idftnp, iswdyp
integer          iescap, ircflp, ischcp, isstpp, ivar, f_id0
integer          nswrsp
integer          imvisp
integer          iflid, iflwgr, f_dim, imasac

integer          icvflb
integer          ivoid(1)

double precision residu, phydr0
double precision ardtsr, arsr  , unsara, thetap
double precision dtsrom, unsvom, romro0
double precision epsrgp, climgp, extrap, epsilp
double precision drom  , dronm1, relaxp
double precision hint, qimp, qimpv(3), epsrsp, blencp
double precision ressol, rnorm2
double precision nadxkm1, nadxk, paxm1ax, paxm1rk, paxkrk, alph, beta
double precision visci(3,3), fikis, viscis, distfi
double precision cfl, kpdc, rho, pimp, bpmasf

type(solving_info) sinfo
type(var_cal_opt) :: vcopt_p, vcopt_u

double precision rvoid(1)

double precision, allocatable, dimension(:) :: dam, xam
double precision, allocatable, dimension(:) :: res, divu, presa, divu_prev
double precision, dimension(:,:), allocatable :: gradp
double precision, allocatable, dimension(:) :: coefaf_dp, coefbf_dp
double precision, allocatable, dimension(:) :: coefap, coefbp, coefa_dp2
double precision, allocatable, dimension(:) :: coefa_rho, coefb_rho
double precision, allocatable, dimension(:) :: cofafp, cofbfp, coefaf_dp2
double precision, allocatable, dimension(:) :: rhs, rovsdt
double precision, allocatable, dimension(:) :: hydro_pres
double precision, allocatable, dimension(:) :: velflx, velflb, ddpvar
double precision, allocatable, dimension(:,:) :: coefar, cofafr
double precision, allocatable, dimension(:,:,:) :: coefbr, cofbfr
double precision, allocatable, dimension(:) :: adxk, adxkm1, dpvarm1, rhs0
double precision, allocatable, dimension(:,:) :: weighf
double precision, allocatable, dimension(:) :: weighb
double precision, allocatable, dimension(:,:) :: frchy, dfrchy
double precision, dimension(:), pointer :: coefa_p, coefb_p
double precision, dimension(:), pointer :: coefaf_p, coefbf_p
double precision, allocatable, dimension(:) :: iflux, bflux
double precision, allocatable, dimension(:) :: xunsro
double precision, allocatable, dimension(:), target :: xdtsro
double precision, allocatable, dimension(:,:), target :: tpusro
double precision, dimension(:), pointer :: viscap
double precision, dimension(:,:), pointer :: vitenp
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: i_mass_flux_prev, b_mass_flux_prev

double precision, dimension(:), pointer :: brom, crom, croma, crom_prev2
double precision, dimension(:), pointer :: cvar_pr, cvara_pr
double precision, dimension(:,:), pointer :: cpro_wgrec_v
double precision, dimension(:), pointer :: cpro_wgrec_s
double precision, dimension(:), pointer :: c_estim_der
double precision, dimension(:), pointer :: cpro_tsrho
double precision, allocatable, dimension(:) :: surfbm

!===============================================================================

!===============================================================================
! 1. Initialisations
!===============================================================================

! Initializations to avoid compiler warnings
rnorm2 = 0.d0
niterf = 0

call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt_u)
call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt_p)

! Allocate temporary arrays
allocate(dam(ncelet), xam(nfac))
allocate(res(ncelet), presa(ncelet), divu(ncelet),divu_prev(ncelet))
allocate(rhs(ncelet), rovsdt(ncelet))
allocate(iflux(nfac), bflux(ndimfb))
iswdyp = vcopt_p%iswdyn
if (iswdyp.ge.1) allocate(adxk(ncelet), adxkm1(ncelet),   &
                          dpvarm1(ncelet), rhs0(ncelet))
if (icalhy.eq.1) allocate(frchy(ndim,ncelet),             &
                          dfrchy(ndim,ncelet), hydro_pres(ncelet))

! Diffusive flux Boundary conditions for delta P
allocate(coefaf_dp(ndimfb), coefbf_dp(ndimfb))

! Associate pointers to pressure diffusion coefficient
viscap => dt(:)
if (vcopt_p%idften.eq.6)  vitenp => tpucou(:,:)

! Index of the field
iflid = ivarfl(ipr)
call field_get_key_struct_solving_info(iflid, sinfo)


! --- Writing
call field_get_name(ivarfl(ipr), chaine)
lchain = 16

f_id0 = -1

! --- Boundary conditions

call field_get_coefa_s(ivarfl(ipr), coefa_p)
call field_get_coefb_s(ivarfl(ipr), coefb_p)
call field_get_coefaf_s(ivarfl(ipr), coefaf_p)
call field_get_coefbf_s(ivarfl(ipr), coefbf_p)

! --- Physical quantities
call field_get_val_s(icrom, crom)
if (icalhy.eq.1.or.idilat.gt.1) then
  call field_get_val_prev_s(icrom, croma)
endif
call field_get_val_s(ibrom, brom)

call field_get_val_s(ivarfl(ipr), cvar_pr)
call field_get_val_prev_s(ivarfl(ipr), cvara_pr)

call field_get_key_int(ivarfl(ipr), kimasf, iflmas)
call field_get_key_int(ivarfl(ipr), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

! --- Solving options
isym  = 1
if(vcopt_p%iconv.gt.0 ) then
  isym  = 2
endif

! Matrix block size
ibsize = 1
iesize = 1
! Initialization dedicated to VOF algo.


! Calculation of dt/rho


!===============================================================================
! 3. Compute a approximated pressure increment if needed
!    that is when there is buoyancy terms (gravity and variable density)
!    with a free outlet.
!===============================================================================

! Standard initialization
do ifac = 1, nfac
  iflux(ifac) = 0.d0
enddo

do ifac = 1, nfabor
  coefa_dp(ifac) = 0.d0
  coefaf_dp(ifac) = 0.d0
  coefb_dp(ifac) = coefb_p(ifac)
  coefbf_dp(ifac) = coefbf_p(ifac)
  bflux(ifac) = 0.d0
enddo
! Compute a pseudo hydrostatic pressure increment stored
! in hydro_pres(.) with Homogeneous Neumann BCs everywhere

! Compute the BCs for the pressure increment
! (first we set the BCs of a standard pressure increment,
!  that are (A = 0, B_dp = B_p) for the gradient BCs
!  Then the A_dp is set thank to the pre-computed hydrostatic pressure
!  so that the pressure increment will be 0 on the reference outlet face.


!===============================================================================
! 4. Building of the linear system to solve
!===============================================================================

! ---> Implicit term

do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo
! Implicit part of the cavitation source

! Strengthen the diagonal for Low Mach Algorithm

! ---> Face diffusivity
if (vcopt_p%idiff.ge.1) then

  ! Scalar diffusivity
 

    if (ivofmt.ge.0) then
      imvisp = 1  ! VOF algorithm: continuity of the flux across internal faces
    else
      imvisp = imvisf
    endif

    call viscfa &
    !==========
   ( imvisp ,            &
     viscap ,            &
     viscf  , viscb  )

    if (vcopt_p%iwgrec.eq.1) then
      ! Weighting for gradient
      do iel = 1, ncel
        cpro_wgrec_s(iel) = viscap(iel)
      enddo
      call synsca(cpro_wgrec_s)
    endif

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

endif

iconvp = vcopt_p%iconv
idiffp = vcopt_p%idiff
ndircp = ndircl(ipr)

thetap = 1.d0
imucpp = 0

call matrix &
!==========
 ( iconvp , idiffp , ndircp , isym   ,                            &
   thetap , imucpp ,                                              &
   coefb_dp , coefbf_dp     , rovsdt ,                            &
   imasfl , bmasfl , viscf  , viscb  ,                            &
   rvoid  , dam    , xam    )

! Free memory
deallocate(iflux, bflux)
!===============================================================================
! 5. Mass flux initialization
!===============================================================================

! --- Flux de masse predit et premiere composante Rhie et Chow

! Allocate a work array for the gradient calculation
allocate(gradp(3,ncelet))

iccocg = 1
iprev  = 1
inc    = 1
 !print*, "cvar_pr,cvara_pr",cvar_pr(50),cvara_pr(50),"In resopv pressure start"
if (ivofmt.lt.0) then
  call field_gradient_potential(ivarfl(ipr), iprev, imrgra, inc,    &
                                iccocg, iphydr,                     &
                                frcxt, gradp)

endif

do iel = 1, ncelet
  do isou = 1, 3
    trav(isou,iel) = gradp(isou,iel)
  enddo
enddo


! Standard algorithm


do iel = 1, ncel
  ardtsr  = arak*(dt(iel)/crom(iel))
  do isou = 1, 3
    trav(isou,iel) = vel(isou,iel) + ardtsr*trav(isou,iel)
  enddo
enddo

! ---> Traitement du parallelisme et de la periodicite

if (irangp.ge.0.or.iperio.eq.1) then
  call synvin(trav)
endif

init   = 1
inc    = 1
! BCs will be taken into account after in idilat>=4
iflmb0 = 1
if (iale.eq.1.or.imobil.eq.1) iflmb0 = 0
nswrgp = vcopt_u%nswrgr
imligp = vcopt_u%imligr
iwarnp = vcopt_p%iwarni
epsrgp = vcopt_u%epsrgr
climgp = vcopt_u%climgr
itypfl = 1


call inimav &
!==========
 ( f_id0  , itypfl ,                                              &
   iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
   iwarnp ,                                                       &
   epsrgp , climgp ,                                              &
   crom, brom   ,                                                 &
   trav   ,                                                       &
   coefav , coefbv ,                                              &
   imasfl , bmasfl )


init   = 0
inc    = 1
iccocg = 1

!----------------------
! Rhie and Choow filter
!----------------------

if (arak.gt.0.d0) then

! --- Prise en compte de Arak : la viscosite face est multipliee
!       Le pas de temps aussi. On retablit plus bas.
  do ifac = 1, nfac
    viscf(ifac) = arak*viscf(ifac)
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = arak*viscb(ifac)
  enddo

  ! On annule la viscosite facette pour les faces couplees pour ne pas modifier
  ! le flux de masse au bord dans le cas d'un dirichlet de pression: la correction
  ! de pression et le filtre sont annules.
  if (nbrcpl.gt.0) then
    do ifac = 1, nfabor
      if (itypfb(ifac).eq.icscpd) then
        viscb(ifac) = 0.d0
      endif
    enddo
  endif

  ! Scalar diffusivity
  !-------------------
    do iel = 1, ncel
      viscap(iel) = arak*viscap(iel)
    enddo

    nswrgp = vcopt_p%nswrgr
    imligp = vcopt_p%imligr
    iwarnp = vcopt_p%iwarni
    epsrgp = vcopt_p%epsrgr
    climgp = vcopt_p%climgr
    extrap = vcopt_p%extrag

    call itrmas &
 ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr ,     &
   0      , iwarnp ,                                                           &
   epsrgp , climgp , extrap ,                                                  &
   frcxt  ,                                                                    &
   cvara_pr          ,                                                         &
   coefa_p , coefb_p , coefaf_p , coefbf_p ,                                   &
   viscf  , viscb  ,                                                           &
   viscap ,                                                                    &
   imasfl , bmasfl )

    ! --- Correction du pas de temps
    unsara = 1.d0/arak
    do iel = 1, ncel
      viscap(iel) = viscap(iel)*unsara
    enddo


  ! --- Correction de la viscosite aux faces
  do ifac = 1, nfac
    viscf(ifac) = viscf(ifac)*unsara
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = viscb(ifac)*unsara
  enddo
  
endif
!===============================================================================
! 6. Solving (Loop over the non-orthogonalities)
!===============================================================================

! --- Number of sweeps
nswmpr = vcopt_p%nswrsm

! --- Variables are set to 0
!       cvar_pr    is the increment of the pressure
!       dpvar      is the increment of the increment between sweeps
!       divu       is the initial divergence of the predicted mass flux
do iel = 1, ncel
  cvar_pr(iel) = 0.d0
  dpvar(iel) = 0.d0
  presa(iel) = 0.d0
enddo

relaxp = vcopt_p%relaxv


! --- Initial divergence
init = 1

call divmas(init, imasfl , bmasfl , divu)

! --- Source term associated to the mass aggregation
if (istmpf.eq.2) then
  call field_get_key_int(ivarfl(iu), kimasf, iflmas)
  call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
  call field_get_val_prev_s(iflmas, i_mass_flux_prev)
  call field_get_val_prev_s(iflmab, b_mass_flux_prev)
  init = 1
  call divmas(init, i_mass_flux_prev, b_mass_flux_prev, divu_prev)
  print*,"imasflu_prev",i_mass_flux_prev(50)
endif

if (idilat.eq.2.or.idilat.eq.3) then
  if (ischtp .eq.2) then
    do iel = 1, ncel
      drom = crom(iel) - croma(iel) !First order dicretization at time
      divu(iel) = divu(iel) + drom*cell_f_vol(iel)/dt(iel)*3.d0/2.d0 + 1.d0/2.d0*divu_prev(iel)
    enddo
  else
    do iel = 1, ncel
      drom = crom(iel) - croma(iel) !First order dicretization at time
      divu(iel)= divu(iel) + drom*cell_f_vol(iel)/dt(iel)
    enddo
  endif
  print*,'crom, croma', crom(50), croma(50),"source term mass aggregation in resopv"
endif
 
! --- Initial right hand side
do iel = 1, ncel
  rhs(iel) = - divu(iel) - rovsdt(iel)*cvar_pr(iel)
enddo


! --- Right hand side residual
residu = sqrt(cs_gdot(ncel,rhs,rhs))

sinfo%rnsmbr = residu

! Pressure derive for the log
if (rnormp.lt.epzero) then
  sinfo%dervar = - sinfo%rnsmbr
else
  sinfo%dervar = sinfo%rnsmbr/rnormp
endif

isweep = 1

! Writing
if (vcopt_p%iwarni.ge.2) then
  write(nfecra,1400)chaine(1:16),isweep,residu, relaxp
endif

! Reconstruction loop (beginning)
!--------------------------------
do while (isweep.le.nswmpr.and.residu.gt.vcopt_p%epsrsm*rnormp)

  ! Solving on the increment dpvar
  !-------------------------------
  if (iswdyp.eq.0) then
    do iel = 1, ncel
      dpvar(iel) = 0.d0
    enddo
  endif

  iwarnp = vcopt_p%iwarni
  epsilp = vcopt_p%epsilo

  ! The pressure is a scalar => no problem for the periodicity of rotation
  ! (iinvpe=1)
  iinvpe = 1

  ! Solver reisudal
  ressol = residu

  call sles_solve_native(ivarfl(ipr), '',                            &
                         isym, ibsize, iesize, dam, xam, iinvpe,     &
                         epsilp, rnormp, niterf, ressol, rhs, dpvar)
  

  ! Update the increment of pressure
  !---------------------------------
  if (iswdyp.eq.0) then
    if (idtvar.ge.0.and.isweep.le.nswmpr.and.residu.gt.vcopt_p%epsrsm*rnormp) then
      do iel = 1, ncel
        presa(iel) = cvar_pr(iel)
        cvar_pr(iel) = presa(iel) + vcopt_p%relaxv*dpvar(iel)
      enddo
    ! If it is the last sweep, update with the total increment
    else
      do iel = 1, ncel
        presa(iel) = cvar_pr(iel)
        cvar_pr(iel) = presa(iel) + dpvar(iel)
      enddo
    endif
 
  endif
  ! --- Update the right hand side and update the residual
  !      rhs^{k+1} = - div(rho u^n) - D(dt, delta delta p^{k+1})
  !-------------------------------------------------------------

  iccocg = 1
  init = 1
  inc  = 0
  if (iphydr.eq.1.or.iifren.eq.1) inc = 1
  nswrgp = vcopt_p%nswrgr
  imligp = vcopt_p%imligr
  iwarnp = vcopt_p%iwarni
  epsrgp = vcopt_p%epsrgr
  climgp = vcopt_p%climgr
  extrap = vcopt_p%extrag

  if (vcopt_p%idften.eq.1) then

    call itrgrp &
    !==========
 ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr ,     &
   iwarnp ,                                                                    &
   epsrgp , climgp , extrap ,                                                  &
   dfrcxt ,                                                                    &
   cvar_pr   ,                                                                 &
   coefa_dp  , coefb_dp  ,                                                     &
   coefaf_dp , coefbf_dp ,                                                     &
   viscf  , viscb  ,                                                           &
   viscap ,                                                                    &
   rhs    )


  endif

  do iel = 1, ncel
    rhs(iel) = - divu(iel) - rhs(iel) - rovsdt(iel)*cvar_pr(iel)
  enddo

  ! --- Convergence test
  residu = sqrt(cs_gdot(ncel,rhs,rhs))

  ! Writing
  sinfo%nbivar = sinfo%nbivar + niterf

  ! Writing
  if (vcopt_p%iwarni.ge.2) then
    write(nfecra,1400) chaine(1:16), isweep, residu, relaxp
    write(nfecra,1500) chaine(1:16), isweep, residu, rnormp, niterf
  endif

  isweep = isweep + 1

enddo
! --- Reconstruction loop (end)

! For logging
if (abs(rnormp).gt.0.d0) then
  sinfo%resvar = residu/rnormp
else
  sinfo%resvar = 0.d0
endif


! Writing
if(vcopt_p%iwarni.ge.1) then
  if (residu.le.vcopt_p%epsrsm*rnormp) then
    write(nfecra,1500) chaine(1:16), isweep-1, residu, rnormp, niterf

  ! Writing: non-convergence
  else if(isweep.gt.nswmpr) then
    write(nfecra,1600) chaine(1:16),nswmpr
  endif
endif

! Save convergence info
call field_set_key_struct_solving_info(ivarfl(ipr), sinfo)

! --- Compute the indicator, taken the volume into account (L2 norm)
!     or not
if(iescal(iesder).gt.0) then
  call field_get_val_s(iestim(iesder), c_estim_der)
  do iel = 1, ncel
    c_estim_der(iel) = abs(rhs(iel))/volume(iel)
  enddo
  if(iescal(iesder).eq.2) then
    do iel = 1, ncel
      c_estim_der(iel) = c_estim_der(iel)*sqrt(volume(iel))
    enddo
  endif
endif
! Update the mass flux
!---------------------

! On annule la viscosite facette pour les faces couplees pour ne pas modifier
! le flux de masse au bord dans le cas d'un dirichlet de pression: la correction
! de pression et le filtre sont annules.
if (nbrcpl.ge.0) then
  do ifac = 1, nfabor
    if (itypfb(ifac).eq.icscpd) then
      viscb(ifac) = 0.d0
    endif
  enddo
endif

iccocg = 1
init = 0
inc  = 0
! In case of hydrostatic pressure, inc is set to 1 to take explicit
! boundary conditions on the pressure (coefa)
if (iphydr.eq.1.or.iifren.eq.1) inc = 1
nswrgp = vcopt_p%nswrgr
imligp = vcopt_p%imligr
iwarnp = vcopt_p%iwarni
epsrgp = vcopt_p%epsrgr
climgp = vcopt_p%climgr
extrap = vcopt_p%extrag

if (vcopt_p%idften.eq.1) then

  call itrmas &
 ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr ,     &
   0      , iwarnp ,                                                           &
   epsrgp , climgp , extrap ,                                                  &
   dfrcxt ,                                                                    &
   presa  ,                                                                    &
   coefa_dp  , coefb_dp  ,                                                     &
   coefaf_dp , coefbf_dp ,                                                     &
   viscf  , viscb  ,                                                           &
   viscap ,                                                                    &
   imasfl , bmasfl )

  ! The last increment is not reconstructed to fullfill exactly the continuity
  ! equation (see theory guide). The value of dfrcxt has no importance.
  iccocg = 0
  nswrgp = 0
  inc = 0

  call itrmas &
 ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr ,     &
   0      , iwarnp ,                                                           &
   epsrgp , climgp , extrap ,                                                  &
   dfrcxt ,                                                                    &
   dpvar  ,                                                                    &
   coefa_dp  , coefb_dp  ,                                                     &
   coefaf_dp , coefbf_dp ,                                                     &
   viscf  , viscb  ,                                                           &
   viscap ,                                                                    &
   imasfl , bmasfl )
endif

!===============================================================================
! 8. Suppression of the mesh hierarchy
!===============================================================================

call sles_free_native(ivarfl(ipr), '')

!===============================================================================
! 9. Weakly compressible algorithm: semi analytic scheme
!    2nd step solving a convection diffusion equation
!===============================================================================


!===============================================================================
! 10. Update the pressure field
!===============================================================================

if (idtvar.lt.0) then
  do iel = 1, ncel
    cvar_pr(iel) = cvara_pr(iel) + vcopt_p%relaxv*cvar_pr(iel)
  enddo
else
  do iel = 1, ncel
    cvar_pr(iel) = cvara_pr(iel) + cvar_pr(iel)
  enddo
endif

! Transformation of volumic mass fluxes into massic mass fluxes

! Free memory
deallocate(dam, xam)
deallocate(res, divu, presa,divu_prev)
deallocate(gradp)
deallocate(coefaf_dp, coefbf_dp)
deallocate(rhs, rovsdt)
if (allocated(weighf)) deallocate(weighf, weighb)
if (iswdyp.ge.1) deallocate(adxk, adxkm1, dpvarm1, rhs0)
if (icalhy.eq.1) deallocate(frchy, dfrchy, hydro_pres)
if (ivofmt.ge.0.or.idilat.eq.4) then
  if (allocated(xdtsro)) deallocate(xdtsro)
  if (allocated(xunsro)) deallocate(xunsro)
  if (allocated(tpusro)) deallocate(tpusro)
endif

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1200 format ( &
 1X,A16,' Sweep: ',I5,' Dynamic relaxation: alpha = ',E12.5,' beta = ',E12.5,/,&
'    < dI^k  ; R^k > = ',E12.5,' ||dI^k  ||^2 = ',E12.5                     ,/,&
'    < dI^k-1; R^k > = ',E12.5,' ||dI^k-1||^2 = ',E12.5                     ,/,&
'   < dI^k-1; dI^k > = ',E12.5)
 1300 format(1X,A16,' : RESIDU DE NORMALISATION =', E14.6)
 1400 format(1X,A16,' : SWEEP = ',I5,' NORME SECOND MEMBRE = ',E14.6,  &
             ', RELAXP = ',E14.6)
 1500 format ( &
 1X,A16,' : Current reconstruction sweep = ',I5                     ,/,&
'           sweep residual = ',E12.5,', norm = ',E12.5              ,/,&
'           number of sweeps for solver = ',I5)
 1600 format( &
'@'                                                                 ,/,&
'@ @@ ATTENTION : ', A16,' ETAPE DE PRESSION'                       ,/,&
'@    ========='                                                    ,/,&
'@  Nombre d''iterations maximal ',I10   ,' atteint'                ,/,&
'@' )

#else

 1200 format ( &
 1X,A16,' Sweep: ',I5,' Dynamic relaxation: alpha = ',E12.5,' beta = ',E12.5,/,&
'    < dI^k  ; R^k > = ',E12.5,' ||dI^k  ||^2 = ',E12.5                     ,/,&
'    < dI^k-1; R^k > = ',E12.5,' ||dI^k-1||^2 = ',E12.5                     ,/,&
'   < dI^k-1; dI^k > = ',E12.5)
 1300 format(1X,A16,' : NORMED RESIDUALS = ', E14.6)
 1400 format(1X,A16,' : SWEEP = ',I5,' RIGHT HAND SIDE NORM = ',E14.6, &
             ', RELAXP = ',E14.6)
 1500 format ( &
 1X,A16,' : Current reconstruction sweep = ',I5                     ,/,&
'           sweep residual = ',E12.5,', norm = ',E12.5              ,/,&
'           number of sweeps for solver = ',I5)
 1600 format( &
'@'                                                                 ,/,&
'@ @@ WARNING: ', A16,' PRESSURE STEP'                              ,/,&
'@    ========'                                                     ,/,&
'@  Maximum number of iterations ',I10   ,' reached'                ,/,&
'@'                                                              )

#endif

!----
! End
!----

return

end subroutine resopv
