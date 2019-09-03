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

!> \file predvv.f90
!>
!> \brief This subroutine performs the velocity prediction step of the Navier
!> Stokes equations for incompressible or slightly compressible flows for
!> the coupled velocity components solver.
!>
!> - at the first call, the predicted velocities are computed and also
!>   an estimator on the predicted velocity is computed.
!>
!> - at the second call, a global estimator on Navier Stokes is computed.
!>   This second call is done after the correction step (\ref resopv).
!>
!> Please refer to the
!> <a href="../../theory.pdf#predvv"><b>predvv</b></b></a> section
!> of the theory guide for more informations.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     iappel        call number (1 or 2)
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     iterns        index of the iteration on Navier-Stokes
!> \param[in]     ncepdp        number of cells with head loss
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     nfbpcd        number of faces with condensation source terms
!> \param[in]     ncmast        number of cells with condensation source terms
!> \param[in]     icepdc        index of cells with head loss
!> \param[in]     icetsm        index of cells with mass source term
!> \param[in]     ifbpcd        index of faces with condensation source terms
!> \param[in]     ltmast        index of cells with condensation source terms
!> \param[in]     itypsm        type of mass source term for the variables
!> \param[in]     dt            time step (per cell)
!> \param[in]     vel           velocity
!> \param[in]     vela          velocity at the previous time step
!> \param[in]     flumas        internal mass flux (depending on iappel)
!> \param[in]     flumab        boundary mass flux (depending on iappel)
!> \param[in]     tslagr        coupling term for the Lagrangian module
!> \param[in]     coefav        boundary condition array for the variable
!>                               (explicit part)
!> \param[in]     coefbv        boundary condition array for the variable
!>                               (implicit part)
!> \param[in]     cofafv        boundary condition array for the diffusion
!>                               of the variable (explicit part)
!> \param[in]     cofbfv        boundary condition array for the diffusion
!>                               of the variable (implicit part)
!> \param[in]     ckupdc        work array for the head loss
!> \param[in]     smacel        variable value associated to the mass source
!>                               term (for ivar=ipr, smacel is the mass flux
!>                               \f$ \Gamma^n \f$)
!> \param[in]     spcond        variable value associated to the condensation
!>                              source term (for ivar=ipr, spcond is the flow rate
!>                              \f$ \Gamma_{s, cond}^n \f$)
!> \param[in]     svcond        variable value associated to the condensation
!>                              source term (for ivar=ipr, svcond is the flow rate
!>                              \f$ \Gamma_{v, cond}^n \f$)
!> \param[in]     frcxt         external forces making hydrostatic pressure
!> \param[in]     trava         working array for the velocity-pressure coupling
!> \param[in]     ximpa         same
!> \param[in]     uvwk          same (stores the velocity at the previous iteration)
!> \param[in]     dfrcxt        variation of the external forces
!                               making the hydrostatic pressure
!> \param[in]     grdphd        hydrostatic pressure gradient to handle the
!>                              imbalance between the pressure gradient and
!>                              gravity source term
!> \param[in]     tpucou        non scalar time step in case of
!>                              velocity pressure coupling
!> \param[in]     trav          right hand side for the normalizing
!>                              the residual
!> \param[in]     viscf         visc*surface/dist aux faces internes
!> \param[in]     viscb         visc*surface/dist aux faces de bord
!> \param[in]     viscfi        same as viscf for increments
!> \param[in]     viscbi        same as viscb for increments
!> \param[in]     secvif        secondary viscosity at interior faces
!> \param[in]     secvib        secondary viscosity at boundary faces
!> \param[in]     w1            working array
!> \param[in]     w7            working array
!> \param[in]     w8            working array
!> \param[in]     w9            working array
!_______________________________________________________________________________

subroutine predvv &
 ( iappel ,                                                       &
   nvar   , nscal  , iterns ,                                     &
   ncepdp , ncesmp , nfbpcd , ncmast ,                            &
   icepdc , icetsm , ifbpcd , ltmast ,                            &
   itypsm ,                                                       &
   dt     , vel    , vela   ,                                     &
   flumas , flumab ,                                              &
   tslagr , coefav , coefbv , cofafv , cofbfv ,                   &
   ckupdc , smacel , spcond , svcond , frcxt  , grdphd ,          &
   trava  , ximpa  , uvwk   , dfrcxt , tpucou , trav   ,          &
   viscf  , viscb  , viscfi , viscbi , secvif , secvib ,          &
   w1     , w7     , w8     , w9     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use cstphy
use cstnum
use optcal
use parall
use period
use lagran
use ppppar
use ppthch
use ppincl
use cplsat
use ihmpre, only: iihmpr
use mesh
use rotation
use turbomachinery
use cs_f_interfaces
use cs_c_bindings
use cfpoin
use field
use field_operator
use pointe, only: gamcav
use cavitation
use vof
use cs_tagms, only:s_metal
use atincl, only: kopint

!===============================================================================

implicit none

! Arguments

integer          iappel
integer          nvar   , nscal  , iterns
integer          ncepdp , ncesmp , nfbpcd , ncmast

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          ifbpcd(nfbpcd)
integer          ltmast(ncelet)

double precision dt(ncelet)
double precision flumas(nfac), flumab(nfabor)
double precision tslagr(ncelet,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision spcond(nfbpcd,nvar), svcond(ncelet,nvar)
double precision frcxt(3,ncelet), dfrcxt(3,ncelet)
double precision grdphd(3, ncelet)
double precision trava(ndim,ncelet)
double precision ximpa(ndim,ndim,ncelet),uvwk(ndim,ncelet)
double precision tpucou(6, ncelet)
double precision trav(3,ncelet)
double precision viscf(*), viscb(nfabor)
double precision viscfi(*), viscbi(nfabor)
double precision secvif(nfac), secvib(nfabor)
double precision w1(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision coefav(3  ,nfabor)
double precision cofafv(3  ,nfabor)
double precision coefbv(3,3,nfabor)
double precision cofbfv(3,3,nfabor)

double precision vel   (3  ,ncelet)
double precision vela  (3  ,ncelet)

! Local variables

integer          f_id  , iel   , ielpdc, ifac  , isou  , itypfl, n_fans
integer          iccocg, inc   , iprev , init  , ii    , jj
integer          nswrgp, imligp, iwarnp
integer          iswdyp, idftnp
integer          iconvp, idiffp, ndircp, nswrsp
integer          ircflp, ischcp, isstpp, iescap
integer          iflmb0, nswrp
integer          idtva0, icvflb, f_oi_id
integer          jsou  , ivisep, imasac
integer          ivoid(1)

double precision rnorm , vitnor
double precision romvom, drom  , rom
double precision epsrgp, climgp, extrap, relaxp, blencp, epsilp
double precision epsrsp
double precision vit1  , vit2  , vit3, xkb, pip, pfac, pfac1
double precision cpdc11, cpdc22, cpdc33, cpdc12, cpdc13, cpdc23
double precision d2s3  , thetap, thetp1, thets , dtsrom
double precision diipbx, diipby, diipbz
double precision ccorio
double precision dvol

double precision rvoid(1)

! Working arrays
double precision, allocatable, dimension(:,:) :: eswork
double precision, allocatable, dimension(:,:) :: grad
double precision, dimension(:,:), allocatable :: smbr
double precision, dimension(:,:,:), allocatable :: fimp
double precision, dimension(:,:), allocatable :: gavinj
double precision, dimension(:,:), allocatable :: tsexp
double precision, dimension(:,:,:), allocatable :: tsimp
double precision, allocatable, dimension(:,:) :: viscce
double precision, dimension(:,:), allocatable :: vect
double precision, dimension(:), allocatable :: xinvro
double precision, dimension(:), pointer :: brom, crom, croma, pcrom
double precision, dimension(:), pointer :: coefa_k, coefb_k
double precision, dimension(:), pointer :: coefa_p, coefb_p
double precision, dimension(:,:), allocatable :: rij
double precision, dimension(:), pointer :: coef1, coef2, coef3, coef4, coef5, coef6
double precision, dimension(:,:), pointer :: coefap
double precision, dimension(:,:,:), pointer :: coefbp
double precision, dimension(:,:), allocatable :: coefat
double precision, dimension(:,:,:), allocatable :: coefbt
double precision, dimension(:,:), allocatable :: tflmas, tflmab
double precision, dimension(:,:), allocatable :: divt
double precision, dimension(:), allocatable ::xnormp
double precision, dimension(:,:), pointer :: forbr, c_st_vel
double precision, dimension(:), pointer :: cvara_pr, cvara_k
double precision, dimension(:), pointer :: cvara_r11, cvara_r22, cvara_r33
double precision, dimension(:), pointer :: cvara_r12, cvara_r23, cvara_r13
double precision, dimension(:,:), pointer :: cvara_rij
double precision, dimension(:), pointer :: viscl, visct, c_estim
double precision, allocatable, dimension(:) :: surfbm
double precision, dimension(:,:), pointer :: lapla
double precision, dimension(:), pointer :: cpro_tsrho

type(var_cal_opt) :: vcopt_p, vcopt_u, vcopt

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate temporary arrays
allocate(smbr(3,ncelet))
allocate(fimp(3,3,ncelet))
allocate(tsexp(3,ncelet))
allocate(tsimp(3,3,ncelet))
call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt_u)
call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt_p)
if (vcopt_u%idften.eq.6) allocate(viscce(6,ncelet))

! Allocate a temporary array for the prediction-stage error estimator
if (iescal(iespre).gt.0) then
  allocate(eswork(3,ncelet))
endif

! Reperage de rho au bord
call field_get_val_s(ibrom, brom)
! Reperage de rho courant (ie en cas d'extrapolation rho^n+1/2)
call field_get_val_s(icrom, crom)
! Reperage de rho^n en cas d'extrapolation
if (iroext.gt.0.or.idilat.gt.1) then
  call field_get_val_prev_s(icrom, croma)
endif


if (iappel.eq.2) then
  if (iforbr.ge.0 .and. iterns.eq.1 .or. ivofmt.ge.0) then
    call field_get_val_s(ivarfl(ipr), cvara_pr)
  endif
  if(iforbr.ge.0 .and. iterns.eq.1                                          &
     .and. (itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60) .and. igrhok.eq.1) then
    call field_get_val_s(ivarfl(ik), cvara_k)
  endif
  if (itytur.eq.3.and.iterns.eq.1) then
    if (irijco.eq.1) then
      call field_get_val_v(ivarfl(irij), cvara_rij)
    else
      call field_get_val_s(ivarfl(ir11), cvara_r11)
      call field_get_val_s(ivarfl(ir22), cvara_r22)
      call field_get_val_s(ivarfl(ir33), cvara_r33)
      call field_get_val_s(ivarfl(ir12), cvara_r12)
      call field_get_val_s(ivarfl(ir23), cvara_r23)
      call field_get_val_s(ivarfl(ir13), cvara_r13)
    endif
  endif
else
  if (iforbr.ge.0 .and. iterns.eq.1 .or. ivofmt.ge.0) then
    call field_get_val_prev_s(ivarfl(ipr), cvara_pr)
  endif
  if(iforbr.ge.0 .and. iterns.eq.1                                          &
     .and. (itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60) .and. igrhok.eq.1) then
      call field_get_val_prev_s(ivarfl(ik), cvara_k)
  endif
  if (itytur.eq.3.and.iterns.eq.1) then
    if (irijco.eq.1) then
      call field_get_val_v(ivarfl(irij), cvara_rij)
    else
      call field_get_val_s(ivarfl(ir11), cvara_r11)
      call field_get_val_s(ivarfl(ir22), cvara_r22)
      call field_get_val_s(ivarfl(ir33), cvara_r33)
      call field_get_val_s(ivarfl(ir12), cvara_r12)
      call field_get_val_s(ivarfl(ir23), cvara_r23)
      call field_get_val_s(ivarfl(ir13), cvara_r13)
    endif
  endif
endif

if (iforbr.ge.0 .and. iterns.eq.1) call field_get_val_v(iforbr, forbr)

! Theta relatif aux termes sources explicites
thets  = thetsn
if (isno2t.gt.0) then
  call field_get_key_int(ivarfl(iu), kstprv, f_id)
  call field_get_val_v(f_id, c_st_vel)
else
  c_st_vel => null()
endif

! Coefficient of the "Coriolis-type" term
if (icorio.eq.1) then
  ! Relative velocity formulation
  ccorio = 2.d0
else if (iturbo.eq.1) then
  ! Mixed relative/absolute velocity formulation
  ccorio = 1.d0
else
  ccorio = 0.d0
endif

!===============================================================================
! 2. Potential forces (pressure gradient and gravity)
!===============================================================================

!-------------------------------------------------------------------------------
! ---> Pressure gradient

! Allocate a work array for the gradient calculation
allocate(grad(3,ncelet))

iccocg = 1
inc    = 1

! For compressible flows, the new Pressure field is required
if (ippmod(icompf).ge.0) then
  iprev = 0
! For incompressible flows, keep the pressure at time n
! in case of PISO algorithm
else
  iprev = 1
endif

  call field_gradient_potential(ivarfl(ipr), iprev, imrgra, inc,    &
                                iccocg, iphydr,                     &
                                frcxt, grad)



!-------------------------------------------------------------------------------
! ---> RESIDU DE NORMALISATION POUR RESOLP
!     Test d'un residu de normalisation de l'etape de pression
!       plus comprehensible = div(rho u* + dt gradP^(n))-Gamma
!       i.e. second membre du systeme en pression hormis la partie
!       pression (sinon a convergence, on tend vers 0)
!       Represente les termes que la pression doit equilibrer
!     On calcule ici div(rho dt/rho gradP^(n)) et on complete a la fin
!       avec  div(rho u*)
!     Pour grad P^(n) on suppose que des CL de Neumann homogenes
!       s'appliquent partout : on peut donc utiliser les CL de la
!       vitesse pour u*+dt/rho gradP^(n). Comme on calcule en deux fois,
!       on utilise les CL de vitesse homogenes pour dt/rho gradP^(n)
!       ci-dessous et les CL de vitesse completes pour u* a la fin.

if (iappel.eq.1.and.irnpnw.eq.1) then

!     Calcul de dt/rho*grad P
  do iel = 1, ncel
    dtsrom = dt(iel)/crom(iel)
    trav(1,iel) = grad(1,iel)*dtsrom
    trav(2,iel) = grad(2,iel)*dtsrom
    trav(3,iel) = grad(3,iel)*dtsrom
  enddo


!     Calcul de rho dt/rho*grad P.n aux faces
!       Pour gagner du temps, on ne reconstruit pas.
  itypfl = 1
  ! VOF algorithm: the pressure step corresponds to the
  ! correction of the volumetric flux, not the mass flux
  if (ivofmt.ge.0)  itypfl = 0
  init   = 1
  inc    = 0
  iflmb0 = 1
  nswrp  = 1
  iwarnp = vcopt_p%iwarni
  imligp = vcopt_u%imligr
  epsrgp = vcopt_u%epsrgr
  climgp = vcopt_u%climgr

  call inimav                                                     &
 ( ivarfl(iu)      , itypfl ,                                     &
   iflmb0 , init   , inc    , imrgra , nswrp  , imligp ,          &
   iwarnp ,                                                       &
   epsrgp , climgp ,                                              &
   crom   , brom   ,                                              &
   trav   ,                                                       &
   coefav , coefbv ,                                              &
   viscf  , viscb  )

  ! Compute div(rho dt/rho*grad P)
  allocate(xnormp(ncelet))

  init = 1
  call divmas(init,viscf,viscb,xnormp)


  ! Dilatable mass conservative algorithm
  if (idilat.eq.2) then
    do iel = 1, ncel
      drom = crom(iel) - croma(iel)
      xnormp(iel) = xnormp(iel) + drom*cell_f_vol(iel)/dt(iel)
    enddo

  endif

!     On conserve XNORMP, on complete avec u* a la fin et
!       on le transfere a resopv

endif


!     Au premier appel, TRAV est construit directement ici.
!     Au second  appel (estimateurs), TRAV contient deja
!       l'increment temporel.
!     On pourrait fusionner en initialisant TRAV a zero
!       avant le premier appel, mais ca fait des operations en plus.

!     Remarques :
!       rho g sera a l'ordre 2 s'il est extrapole.
!       si on itere sur navsto, ca ne sert a rien de recalculer rho g a
!         chaque fois (ie on pourrait le passer dans trava) mais ce n'est
!         pas cher.
if (iappel.eq.1) then
  do iel = 1, ncel
    drom = (crom(iel)-ro0)
    trav(1,iel) = (drom*gx - grad(1,iel) ) * cell_f_vol(iel)
    trav(2,iel) = (drom*gy - grad(2,iel) ) * cell_f_vol(iel)
    trav(3,iel) = (drom*gz - grad(3,iel) ) * cell_f_vol(iel)
  enddo

else if(iappel.eq.2) then

  
  do iel = 1, ncel
    drom = (crom(iel)-ro0)
    trav(1,iel) = trav(1,iel) + (drom*gx - grad(1,iel))*cell_f_vol(iel)
    trav(2,iel) = trav(2,iel) + (drom*gy - grad(2,iel))*cell_f_vol(iel)
    trav(3,iel) = trav(3,iel) + (drom*gz - grad(3,iel))*cell_f_vol(iel)
  enddo
 

endif

! Free memory
deallocate(grad)


!   Pour IAPPEL = 1 (ie appel standard sans les estimateurs)
!     TRAV rassemble les termes sources  qui seront recalcules
!       a toutes les iterations sur navsto
!     Si on n'itere pas sur navsto et qu'on n'extrapole pas les
!       termes sources, TRAV contient tous les termes sources
!       jusqu'au basculement dans SMBR
!     A ce niveau, TRAV contient -grad P et rho g
!       P est suppose pris a n+1/2
!       rho est eventuellement interpole a n+1/2


!-------------------------------------------------------------------------------
! ---> INITIALISATION DU TABLEAU TRAVA et terme source AU PREMIER PASSAGE
!     (A LA PREMIERE ITER SUR NAVSTO)

!     TRAVA rassemble les termes sources qu'il suffit de calculer
!       a la premiere iteration sur navsto quand il y a plusieurs iter.
!     Quand il n'y a qu'une iter, on cumule directement dans TRAV
!       ce qui serait autrement alle dans TRAVA
!     Les termes sources explicites serviront
!       pour le pas de temps suivant en cas d'extrapolation (plusieurs
!       iter sur navsto ou pas)

!     A la premiere iter sur navsto
if (iterns.eq.1) then

  ! Si on   extrapole     les T.S. : -theta*valeur precedente

  if (isno2t.gt.0) then
    ! S'il n'y a qu'une    iter : TRAV  incremente
    if (nterup.eq.1) then
      do iel = 1, ncel
        do ii = 1, ndim
          trav (ii,iel) = trav (ii,iel) - thets*c_st_vel(ii,iel)
        enddo
      enddo
      ! S'il   y a plusieurs iter : TRAVA initialise
    else
      do iel = 1, ncel
        do ii = 1, ndim
          trava(ii,iel) = - thets*c_st_vel(ii,iel)
        enddo
      enddo
    endif
    ! Et on initialise le terme source pour le remplir ensuite
    do iel = 1, ncel
      do ii = 1, ndim
        c_st_vel(ii,iel) = 0.d0
      enddo
    enddo

  ! Si on n'extrapole pas les T.S.
  else
    ! S'il   y a plusieurs iter : TRAVA initialise
    !  sinon TRAVA n'existe pas
    if(nterup.gt.1) then
      do ii = 1, ndim
        do iel = 1, ncel
          trava(ii,iel)  = 0.d0
        enddo
      enddo
    endif
  endif

endif

!-------------------------------------------------------------------------------
! Initialization of the implicit terms

if (iappel.eq.1) then

  ! Low Mach compressible Algos
  if (idilat.gt.1.or.ippmod(icompf).ge.0) then
    call field_get_val_prev_s(icrom, pcrom)

  ! Standard algo
  else

    call field_get_val_s(icrom, pcrom)
  endif

  do iel = 1, ncel
    do isou = 1, 3
      fimp(isou,isou,iel) = vcopt_u%istat*pcrom(iel)/dt(iel)*cell_f_vol(iel)
      do jsou = 1, 3
        if(jsou.ne.isou) fimp(isou,jsou,iel) = 0.d0
      enddo
    enddo
  enddo

!     Le remplissage de FIMP est toujours indispensable,
!       meme si on peut se contenter de n'importe quoi pour IAPPEL=2.
else
  do iel = 1, ncel
    do isou = 1, 3
      do jsou = 1, 3
        fimp(isou,jsou,iel) = 0.d0
      enddo
    enddo
  enddo
endif

!-------------------------------------------------------------------------------
! ---> 2/3 RHO * GRADIENT DE K SI k-epsilon ou k-omega
!      NB : ON NE PREND PAS LE GRADIENT DE (RHO K), MAIS
!           CA COMPLIQUERAIT LA GESTION DES CL ...
!     On peut se demander si l'extrapolation en temps sert a
!       quelquechose

!     Ce terme explicite est calcule une seule fois,
!       a la premiere iter sur navsto : il est stocke dans un champ si on
!       doit l'extrapoler en temps ; il va dans TRAVA si on n'extrapole
!       pas mais qu'on itere sur navsto. Il va dans TRAV si on
!       n'extrapole pas et qu'on n'itere pas sur navsto.


!-------------------------------------------------------------------------------
! ---> Transpose of velocity gradient in the diffusion term

!     These terms are taken into account in bilscv.
!     We only compute here the secondary viscosity.

if (ivisse.eq.1) then

  call visecv(secvif, secvib)

endif




!-------------------------------------------------------------------------------
! ---> Face diffusivity for the velocity
if (vcopt_u%idiff.ge. 1) then

  call field_get_val_s(iviscl, viscl)
  call field_get_val_s(ivisct, visct)

  do iel = 1, ncel
    w1(iel) = viscl(iel) + vcopt_u%idifft*visct(iel)
  enddo


  ! Scalar diffusivity (Default)


    call viscfa &
   ( imvisf ,                                                       &
     w1     ,                                                       &
     viscf  , viscb  )



! --- If no diffusion, viscosity is set to 0.
else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

  if(itytur.eq.3.and.irijnu.eq.1) then
    do ifac = 1, nfac
      viscfi(ifac) = 0.d0
    enddo
    do ifac = 1, nfabor
      viscbi(ifac) = 0.d0
    enddo
  endif

endif

!-------------------------------------------------------------------------------
! ---> Take external forces partially equilibrated with the pressure gradient
!      into account (only for the first call, the second one is dedicated
!      to error estimators)



!===============================================================================
! 3. Solving of the 3x3xNcel coupled system
!===============================================================================


! ---> AU PREMIER APPEL,
!      MISE A ZERO DE L'ESTIMATEUR POUR LA VITESSE PREDITE
!      S'IL DOIT ETRE CALCULE

if (iappel.eq.1) then
  if (iestim(iespre).ge.0) then
    call field_get_val_s(iestim(iespre), c_estim)
    do iel = 1, ncel
      c_estim(iel) =  0.d0
    enddo
  endif
endif

! ---> AU DEUXIEME APPEL,
!      MISE A ZERO DE L'ESTIMATEUR TOTAL POUR NAVIER-STOKES
!      (SI ON FAIT UN DEUXIEME APPEL, ALORS IL DOIT ETRE CALCULE)

if (iappel.eq.2) then
  call field_get_val_s(iestim(iestot), c_estim)
  do iel = 1, ncel
    c_estim(iel) =  0.d0
  enddo
endif

!-------------------------------------------------------------------------------
! ---> User source terms

do iel = 1, ncel
  do isou = 1, 3
    tsexp(isou,iel) = 0.d0
    do jsou = 1, 3
      tsimp(isou,jsou,iel) = 0.d0
    enddo
  enddo
enddo

! The computation of esplicit and implicit source terms is performed
! at the first iteration only.
if (iterns.eq.1) then

  if (iihmpr.eq.1) then
    call uitsnv (vel, tsexp, tsimp)
  endif

  call ustsnv &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iu   ,                                                         &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel , tsexp  , tsimp  )

  if (vcopt_u%ibdtso.gt.1.and.ntcabs.gt.ntinit &
      .and.(idtvar.eq.0.or.idtvar.eq.1)) then
    ! TODO: remove test on ntcabs and implemente a "proper" condition for
    ! initialization.
    f_id = ivarfl(iu)
    call cs_backward_differentiation_in_time(f_id, tsexp, tsimp)
  endif
  ! Skip first time step after restart if previous values have not been read.
  if (vcopt_u%ibdtso.lt.0) then
    vcopt_u%ibdtso = iabs(vcopt_u%ibdtso)
    call field_set_key_struct_var_cal_opt(ivarfl(iu), vcopt_u)
  endif

  ! Nudging towards optimal interpolation for velocity
  if (ippmod(iatmos).ge.0) then
    call field_get_key_int(ivarfl(iu), kopint, f_oi_id)
    if (f_oi_id.ge.0) then
      call cs_at_data_assim_source_term(ivarfl(iu), tsexp, tsimp)
    endif
  endif


endif

! if PISO sweeps are expected, implicit user sources terms are stored in ximpa
if (iterns.eq.1.and.nterup.gt.1) then
  do iel = 1, ncel
    do isou = 1, 3
      do jsou = 1, 3
        ximpa(isou,jsou,iel) = tsimp(isou,jsou,iel)
      enddo
    enddo
  enddo
endif

! ---> Explicit contribution due to implicit terms

if (iterns.eq.1) then
  if (nterup.gt.1) then
    do iel = 1, ncel
      do isou = 1, 3
        do jsou = 1, 3
          trava(isou,iel) = trava(isou,iel)                                  &
                          + tsimp(isou,jsou,iel)*vela(jsou,iel)
        enddo
      enddo
    enddo
  else
    do iel = 1, ncel
      do isou = 1, 3
        do jsou = 1, 3
          trav(isou,iel) = trav(isou,iel)                                    &
                         + tsimp(isou,jsou,iel)*vela(jsou,iel)
        enddo
      enddo
    enddo
  endif
endif

! At the first PISO iteration, explicit source terms are added
if (iterns.eq.1.and.(iphydr.ne.1.or.igpust.ne.1)) then
  ! If source terms are time-extrapolated, they are stored in fields
  if (isno2t.gt.0) then
    do iel = 1, ncel
      do isou = 1, 3
        c_st_vel(isou,iel) = c_st_vel(isou,iel) + tsexp(isou,iel)
      enddo
    enddo

  else
    ! If no PISO sweep
    if (nterup.eq.1) then
      do iel = 1, ncel
        do isou = 1, 3
          trav(isou,iel) = trav(isou,iel) + tsexp(isou,iel)
        enddo
      enddo
    ! If PISO sweeps
    else
      do iel = 1, ncel
        do isou = 1, 3
          trava(isou,iel) = trava(isou,iel) + tsexp(isou,iel)
        enddo
      enddo
    endif
  endif
endif

! ---> Implicit terms
if (iappel.eq.1) then
  ! If source terms are time-extrapolated
  if (isno2t.gt.0) then
    thetap = vcopt_u%thetav
    if (iterns.gt.1) then
      do iel = 1, ncel
        do isou = 1, 3
          do jsou = 1, 3
            fimp(isou,jsou,iel) = fimp(isou,jsou,iel)                      &
                                - ximpa(isou,jsou,iel)*thetap
          enddo
        enddo
      enddo
    else
      do iel = 1, ncel
        do isou = 1, 3
          do jsou = 1, 3
            fimp(isou,jsou,iel) = fimp(isou,jsou,iel)                      &
                                - tsimp(isou,jsou,iel)*thetap
          enddo
        enddo
      enddo
    endif
  else
    if (iterns.gt.1) then
      do iel = 1, ncel
        do isou = 1, 3
          do jsou = 1, 3
            fimp(isou,jsou,iel) = fimp(isou,jsou,iel)                      &
                                + max(-ximpa(isou,jsou,iel),zero)
          enddo
        enddo
      enddo
    else
      do iel = 1, ncel
        do isou = 1, 3
          do jsou = 1, 3
            fimp(isou,jsou,iel) = fimp(isou,jsou,iel)                      &
                                + max(-tsimp(isou,jsou,iel),zero)
          enddo
        enddo
      enddo
    endif
  endif
endif

!-------------------------------------------------------------------------------
! --->  Mass source terms


! ---> Right Hand Side initialization

! If source terms are extrapolated in time
if (isno2t.gt.0) then
  thetp1 = 1.d0 + thets
  ! If no PISO iteration: trav
  if (nterup.eq.1) then
    do iel = 1, ncel
      do isou = 1, 3
        smbr(isou,iel) = trav(isou,iel) + thetp1*c_st_vel(isou,iel)
      enddo
    enddo

  else
    do iel = 1, ncel
      do isou = 1, 3
        smbr(isou,iel) = trav(isou,iel) + trava(isou,iel)       &
                       + thetp1*c_st_vel(isou,iel)
      enddo
    enddo
  endif

! No time extrapolation
else
  ! No PISO iteration
  if (nterup.eq.1) then
    do iel = 1, ncel
      do isou = 1, 3
        smbr(isou,iel) = trav(isou,iel)
      enddo
    enddo
  ! PISO iterations
  else
    do iel = 1, ncel
      do isou = 1, 3
        smbr(isou,iel) = trav(isou,iel) + trava(isou,iel)
      enddo
    enddo
  endif
endif


! Solver parameters
iconvp = vcopt_u%iconv
idiffp = vcopt_u%idiff
ndircp = ndircl(iu)
nswrsp = vcopt_u%nswrsm
nswrgp = vcopt_u%nswrgr
imligp = vcopt_u%imligr
ircflp = vcopt_u%ircflu
ischcp = vcopt_u%ischcv
isstpp = vcopt_u%isstpc
idftnp = vcopt_u%idften
iswdyp = vcopt_u%iswdyn
iwarnp = vcopt_u%iwarni
blencp = vcopt_u%blencv
epsilp = vcopt_u%epsilo
epsrsp = vcopt_u%epsrsm
epsrgp = vcopt_u%epsrgr
climgp = vcopt_u%climgr
extrap = vcopt_u%extrag
relaxp = vcopt_u%relaxv
thetap = vcopt_u%thetav

if (ippmod(icompf).ge.0) then
  ! impose boundary convective flux at some faces (face indicator icvfli)
  icvflb = 1
else
  ! all boundary convective flux with upwind
  icvflb = 0
endif

if (iappel.eq.1) then

  iescap = iescal(iespre)

  if (iterns.eq.1) then

    ! Warning: in case of convergence estimators, eswork give the estimator
    ! of the predicted velocity
    call coditv &
 ( idtvar , ivarfl(iu)      , iconvp , idiffp , ndircp ,          &
   imrgra , nswrsp , nswrgp , imligp , ircflp , ivisse ,          &
   ischcp , isstpp , iescap , idftnp , iswdyp ,                   &
   iwarnp ,                                                       &
   blencp , epsilp , epsrsp , epsrgp , climgp ,                   &
   relaxp , thetap ,                                              &
   vela   , vela   ,                                              &
   coefav , coefbv , cofafv , cofbfv ,                            &
   flumas , flumab ,                                              &
   viscfi , viscbi , viscf  , viscb  , secvif , secvib ,          &
   icvflb , icvfli ,                                              &
   fimp   ,                                                       &
   smbr   ,                                                       &
   vel    ,                                                       &
   eswork )

  else if(iterns.gt.1) then

    call coditv &
 ( idtvar , ivarfl(iu)      , iconvp , idiffp , ndircp ,          &
   imrgra , nswrsp , nswrgp , imligp , ircflp , ivisse ,          &
   ischcp , isstpp , iescap , idftnp , iswdyp ,                   &
   iwarnp ,                                                       &
   blencp , epsilp , epsrsp , epsrgp , climgp ,                   &
   relaxp , thetap ,                                              &
   vela   , uvwk   ,                                              &
   coefav , coefbv , cofafv , cofbfv ,                            &
   flumas , flumab ,                                              &
   viscfi , viscbi , viscf  , viscb  , secvif , secvib ,          &
   icvflb , icvfli ,                                              &
   fimp   ,                                                       &
   smbr   ,                                                       &
   vel    ,                                                       &
   eswork )

  endif

  ! Velocity-pression coupling: compute the vector T, stored in tpucou,
  !  coditv is called, only one sweep is done, and tpucou is initialized
  !  by 0. so that the advection/diffusion added by bilscv is 0.
  !  nswrsp = -1 indicated that only one sweep is required and inc=0
  !  for boundary contitions on the weight matrix.
  if (ipucou.eq.1) then

    ! Allocate temporary arrays for the velocity-pressure resolution
    allocate(vect(3,ncelet))

    nswrsp = -1
    do iel = 1, ncel
      do isou = 1, 3
        smbr(isou,iel) = cell_f_vol(iel)
      enddo
    enddo
    do iel = 1, ncelet
      do isou = 1, 3
        vect(isou,iel) = 0.d0
      enddo
    enddo
    iescap = 0

    ! We do not take into account transpose of grad
    ivisep = 0

    call coditv &
 ( idtvar , ivarfl(iu)      , iconvp , idiffp , ndircp ,          &
   imrgra , nswrsp , nswrgp , imligp , ircflp , ivisep ,          &
   ischcp , isstpp , iescap , idftnp , iswdyp ,                   &
   iwarnp ,                                                       &
   blencp , epsilp , epsrsp , epsrgp , climgp ,                   &
   relaxp , thetap ,                                              &
   vect   , vect   ,                                              &
   coefav , coefbv , cofafv , cofbfv ,                            &
   flumas , flumab ,                                              &
   viscfi , viscbi , viscf  , viscb  , secvif , secvib ,          &
   icvflb , ivoid  ,                                              &
   fimp   ,                                                       &
   smbr   ,                                                       &
   vect   ,                                                       &
   rvoid  )

    do iel = 1, ncelet
      rom = crom(iel)
      do isou = 1, 3
        tpucou(isou,iel) = rom*vect(isou,iel)
      enddo
      do isou = 4, 6
        tpucou(isou,iel) = 0.d0
      enddo
    enddo

    ! Free memory
    deallocate(vect)

  endif

  ! ---> The estimator on the predicted velocity is summed up over the components
  if (iestim(iespre).ge.0) then
    call field_get_val_s(iestim(iespre), c_estim)
    do iel = 1, ncel
      do isou = 1, 3
        c_estim(iel) =  c_estim(iel) + eswork(isou,iel)
      enddo
    enddo
  endif


! ---> End of the construction of the total estimator:
!       RHS resiudal of (U^{n+1}, P^{n+1}) + rho*volume*(U^{n+1} - U^n)/dt
else if (iappel.eq.2) then

  inc = 1
  ! Pas de relaxation en stationnaire
  idtva0 = 0
  imasac = 0

  call bilscv &
 ( idtva0 , ivarfl(iu)      , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , ivisse ,                            &
   iwarnp , idftnp , imasac ,                                              &
   blencp , epsrgp , climgp , relaxp , thetap ,                            &
   vel    , vel    ,                                                       &
   coefav , coefbv , cofafv , cofbfv ,                                     &
   flumas , flumab , viscf  , viscb  , secvif , secvib ,                   &
   icvflb , icvfli ,                                                       &
   smbr   )

  call field_get_val_s(iestim(iestot), c_estim)
  do iel = 1, ncel
    do isou = 1, 3
      c_estim(iel) = c_estim(iel) + (smbr(isou,iel)/volume(iel))**2
    enddo
  enddo
endif

!===============================================================================
! 4. Finalize the norm of the pressure step (see resopv)
!===============================================================================

if (iappel.eq.1.and.irnpnw.eq.1) then

  ! Compute div(rho u*)


  ! To save time, no space reconstruction
  itypfl = 1
  ! VOF algorithm: the pressure step corresponds to the
  ! correction of the volumetric flux, not the mass flux
  if (ivofmt.ge.0)  itypfl = 0
  init   = 1
  inc    = 1
  iflmb0 = 1
  nswrp  = 1

  iwarnp = vcopt_p%iwarni
  imligp = vcopt_u%imligr
  epsrgp = vcopt_u%epsrgr
  climgp = vcopt_u%climgr

  call inimav &
 ( ivarfl(iu)      , itypfl ,                                     &
   iflmb0 , init   , inc    , imrgra , nswrp  , imligp ,          &
   iwarnp ,                                                       &
   epsrgp , climgp ,                                              &
   crom   , brom   ,                                              &
   vel    ,                                                       &
   coefav , coefbv ,                                              &
   viscf  , viscb  )

  init = 0
  call divmas(init,viscf,viscb,xnormp)

  ! Compute the norm rnormp used in resopv
  rnormp = sqrt(cs_gdot(ncel,xnormp,xnormp))

  ! Free memory
  deallocate(xnormp)
endif

! ---> Finilaze estimators + Printings

if (iappel.eq.1) then

  ! ---> Estimator on the predicted velocity:
  !      square root (norm) or square root of the sum times the volume (L2 norm)
  if (iestim(iespre).ge.0) then
    call field_get_val_s(iestim(iespre), c_estim)
    if (iescal(iespre).eq.1) then
      do iel = 1, ncel
        c_estim(iel) = sqrt(c_estim(iel))
      enddo
    else if (iescal(iespre).eq.2) then
      do iel = 1, ncel
        c_estim(iel) = sqrt(c_estim(iel)*volume(iel))
      enddo
    endif
  endif

  ! ---> Norm printings
  if (vcopt_u%iwarni.ge.2) then
    rnorm = -1.d0
    do iel = 1, ncel
      vitnor = sqrt(vel(1,iel)**2+vel(2,iel)**2+vel(3,iel)**2)
      rnorm = max(rnorm,vitnor)
    enddo

    if (irangp.ge.0) call parmax (rnorm)

    write(nfecra,1100) rnorm

    do iel = 1, ncel
      vitnor = sqrt(vel(1,iel)**2+vel(2,iel)**2+vel(3,iel)**2)
      rnorm = min(rnorm,vitnor)
    enddo

    if (irangp.ge.0) call parmin (rnorm)

    write(nfecra,1200) rnorm

  endif

! ---> Estimator on the whole Navier-Stokes:
!      square root (norm) or square root of the sum times the volume (L2 norm)
else if (iappel.eq.2) then

  call field_get_val_s(iestim(iestot), c_estim)
  if (iescal(iestot).eq.1) then
    do iel = 1, ncel
      c_estim(iel) = sqrt(c_estim(iel))
    enddo
  else if (iescal(iestot).eq.2) then
    do iel = 1, ncel
      c_estim(iel) = sqrt(c_estim(iel)*volume(iel))
    enddo
  endif

endif

! Free memory
!------------
deallocate(smbr)
deallocate(fimp)
deallocate(tsexp)
deallocate(tsimp)
if (allocated(viscce)) deallocate(viscce)
if (allocated(divt)) deallocate(divt)

!--------
! Formats
!--------
#if defined(_CS_LANG_FR)

 1100 format(/,                                                   &
 1X,'Vitesse maximale apres prediction ',E12.4)

 1200 format(/,                                                   &
 1X,'Vitesse minimale apres prediction ',E12.4)

#else

 1100 format(/,                                                   &
 1X,'Maximum velocity after prediction ',E12.4)

 1200 format(/,                                                   &
 1X,'Minimum velocity after prediction ',E12.4)

#endif

!----
! End
!----

return

end subroutine
