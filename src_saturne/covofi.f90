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

!> \file covofi.f90
!>
!> \brief This subroutine performs the solving the convection/diffusion
!> equation (with eventually source terms and/or drift) for a scalar quantity
!> over a time step.
!>
!> Please refer to the
!> <a href="../../theory.pdf#covofi"><b>covofi</b></a> section of the
!> theory guide for more informations.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     nfbpcd        number of faces with condensation source terms
!> \param[in]     ncmast        number of cells with condensation source terms
!> \param[in]     iscal         scalar number
!> \param[in]     itspdv        indicator to compute production/dissipation
!>                              terms for a variance:
!>                               - 0: no
!>                               - 1: yes
!> \param[in]     icepdc        index of cells with head loss
!> \param[in]     icetsm        index of cells with mass source term
!> \param[in]     ifbpcd        index of faces with condensation source terms
!> \param[in]     ltmast        index of cells with condensation source terms
!> \param[in]     itypsm        type of mass source term for the variables
!> \param[in]     itypcd        type of surface condensation source term
!> \param[in]     itypst        type of volume  condensation source term
!> \param[in]     dt            time step (per cell)
!> \param[in]     tslagr        coupling term for the Lagrangian module
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
!> \param[in]     flxmst        variable value associated to heat transfer flux
!>                              associated to the metal mass condensation
!> \param[in]     viscf         visc*surface/dist at internal faces
!> \param[in]     viscb         visc*surface/dist at boundary faces
!_______________________________________________________________________________

subroutine covofi &
 ( nvar   , nscal  , ncepdp , ncesmp , nfbpcd , ncmast ,          &
   iscal  , itspdv ,                                              &
   icepdc , icetsm , ifbpcd , ltmast ,                            &
   itypsm , itypcd , itypst ,                                     &
   dt     , tslagr ,                                              &
   ckupdc , smacel , spcond , svcond , flxmst ,                   &
   viscf  , viscb, iter_sca,icvrge)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use ppppar
use ppthch
use coincl
use cpincl
use cs_fuel_incl
use ppincl
use ppcpfu
use lagran
use radiat
use field
use field_operator
use ihmpre, only: iihmpr
use mesh
use parall
use period
use cs_f_interfaces
use atchem
use darcy_module
use cs_c_bindings
use pointe, only: itypfb, pmapper_double_r1
use atincl, only: kopint

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp , nfbpcd ,  ncmast
integer          iscal  , itspdv, iter_sca, icvrge

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          ifbpcd(nfbpcd), itypcd(nfbpcd,nvar)
integer          ltmast(ncelet), itypst(ncelet,nvar)

double precision dt(ncelet)
double precision tslagr(ncelet,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision spcond(nfbpcd,nvar)
double precision svcond(ncelet,nvar), flxmst(ncelet)
double precision viscf(nfac), viscb(nfabor)


! Local variables

logical          lprev
character(len=80) :: chaine, fname
integer          ivar
integer          ii, ifac , iel, isou
integer          iprev , inc   , iccocg, iiun, ibcl
integer          ivarsc
integer          iiscav
integer          ifcvsl, iflmas, iflmab, f_oi_id
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          imucpp, idftnp, iswdyp
integer          iflid , f_id, st_prv_id, st_id,  keydri, iscdri
integer          icvflb, f_dim, iflwgr
integer          icla
integer          icrom_scal, isorb, keysrb, igwfpr, keypre

integer          ivoid(1)

double precision epsrgp, climgp, extrap, relaxp, blencp, epsilp
double precision epsrsp
double precision rhovst, xk    , xe    , sclnor
double precision thetv , thets , thetap, thetp1
double precision smbexp, dvar, cprovol, prod, expkdt
double precision temp, idifftp, roskpl, kplskm, ctot_tmp
double precision turb_schmidt
double precision vark(ncelet)
double precision xvarmu, xvarmu0, xdvar, xvartmp
double precision rvoid(1)

double precision, allocatable, dimension(:) :: w1, smbrs, rovsdt
double precision, allocatable, dimension(:,:) :: viscce
double precision, allocatable, dimension(:,:) :: weighf
double precision, allocatable, dimension(:) :: weighb
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: coefa_p, coefb_p
double precision, allocatable, dimension(:) :: dpvar
double precision, allocatable, dimension(:) :: xcpp
double precision, allocatable, dimension(:) :: srcmas
double precision, allocatable, dimension(:) :: srccond
double precision, allocatable, dimension(:) :: srcmst

double precision, dimension(:,:), pointer :: xut, visten
double precision, dimension(:,:), pointer :: cpro_wgrec_v
double precision, dimension(:), pointer :: cpro_wgrec_s
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: crom, croma, pcrom
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:), pointer :: cvara_k, cvara_ep, cvara_omg
double precision, dimension(:), pointer :: cvara_r11, cvara_r22, cvara_r33
double precision, dimension(:,:), pointer :: cvara_rij
double precision, dimension(:), pointer :: visct, cpro_cp, cproa_scal_st
double precision, dimension(:), pointer :: cpro_scal_st
double precision, dimension(:), pointer :: cpro_viscls, cpro_visct
double precision, dimension(:), pointer :: cpro_tsscal
double precision, dimension(:), pointer :: cpro_x2icla
! Darcy arrays
double precision, allocatable, dimension(:) :: diverg, divflu
double precision, dimension(:), pointer :: cpro_delay, cpro_sat
double precision, dimension(:), pointer :: cproa_delay, cproa_sat
double precision, dimension(:), pointer :: cvar_var, cvara_var, cvara_varsca
double precision, dimension(:), pointer :: cpro_rosoil, cpro_sorb, cpro_precip
double precision, dimension(:), pointer :: cpro_kplus, cpro_kminus, cpro_mxsol
! Radiat arrays
double precision, dimension(:), pointer :: cpro_tsre1, cpro_tsre, cpro_tsri1
character(len=80) :: f_name

type(var_cal_opt) :: vcopt, vcopt_varsc
type(gwf_soilwater_partition) :: sorption_scal

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! --- Variable number
ivar   = isca(iscal)

call field_get_val_s(ivarfl(ivar), cvar_var)
call field_get_val_prev_s(ivarfl(ivar), cvara_var)

! Index of the field
iflid = ivarfl(ivar)

! Key id for drift scalar
call field_get_key_id("drift_scalar_model", keydri)

! Id of the mass flux
call field_get_key_int(iflid, kimasf, iflmas) ! interior mass flux
! Pointer to the internal mass flux
call field_get_val_s(iflmas, imasfl)

! Id of the mass flux
call field_get_key_int(iflid, kbmasf, iflmab) ! boundary mass flux
! Pointer to the Boundary mass flux
call field_get_val_s(iflmab, bmasfl)

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)


! Allocate temporary arrays
allocate(w1(ncelet))
allocate(dpvar(ncelet))
allocate(smbrs(ncelet), rovsdt(ncelet))

if (ippmod(idarcy).eq.1) then
  allocate(diverg(ncelet))
endif

! Initialize variables to avoid compiler warnings

xe = 0.d0
xk = 0.d0

! --- Numero du scalaire eventuel associe dans le cas fluctuation
!         et numero de variable de calcul
iiscav = iscavr(iscal)
if (iiscav.gt.0.and.iiscav.le.nscal) then
  ivarsc = isca(iiscav)
else
  ivarsc = 0
endif

! --- Numero des grandeurs physiques
call field_get_key_int (ivarfl(isca(iscal)), kromsl, icrom_scal)

if (icrom_scal.eq.-1) then
  icrom_scal = icrom
endif

call field_get_val_s(icrom_scal, crom)
call field_have_previous(icrom_scal, lprev)
if (lprev) then
  call field_get_val_prev_s(icrom_scal, croma)
endif
call field_get_val_s(ivisct, visct)

call field_get_key_int (ivarfl(isca(iscal)), kivisl, ifcvsl)
if (ifcvsl.ge.0) then
  call field_get_val_s(ifcvsl, cpro_viscls)
endif

if (idilat.ge.4) then
  call field_get_val_s(iustdy(iscal), cpro_tsscal)
endif

! --- Numero de propriété du terme source si extrapolation
call field_get_key_int(iflid, kstprv, st_prv_id)
if (st_prv_id .ge.0) then
  call field_get_val_s(st_prv_id, cproa_scal_st)
else
  cproa_scal_st => null()
endif

! S pour Source, V pour Variable
thets  = thetss(iscal)
thetv  = vcopt%thetav
print*,"thetav in covofi", thetv
call field_get_name(ivarfl(ivar), chaine)

if(vcopt%iwarni.ge.1) then
  write(nfecra,1000) chaine(1:16)
endif

! When solving the Temperature, we solve:
!  rho*cp*Vol*dT/dt + ...

imucpp = 0
if (iscavr(iscal).gt.0) then
  if (abs(iscacp(iscavr(iscal))).eq.1) then
    imucpp = 1
  endif
else
  if (abs(iscacp(iscal)).eq.1) then
    imucpp = 1
  endif
endif

allocate(xcpp(ncelet))

if (imucpp.eq.0) then
  do iel = 1, ncel
    xcpp(iel) = 1.d0
  enddo
elseif (imucpp.eq.1) then
  if (icp.ge.0) then
    call field_get_val_s(icp, cpro_cp)
    do iel = 1, ncel
      xcpp(iel) = cpro_cp(iel)
    enddo
  else
    do iel = 1, ncel
      xcpp(iel) = cp0
    enddo
  endif
endif

! Handle parallelism and periodicity
if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(xcpp)
endif

! Retrieve turbulent Schmidt value for current scalar
call field_get_key_double(ivarfl(isca(iscal)), ksigmas, turb_schmidt)

!===============================================================================
! 2. Source terms
!===============================================================================

! --> Initialization

do iel = 1, ncel
  rovsdt(iel) = 0.d0
  smbrs(iel) = 0.d0
enddo
if (iihmpr.eq.1) then
  if (iscal.ne.iscalt) then
    call uitssc &
    ( ippmod(idarcy), iflid  , cvar_var , smbrs  , rovsdt )
  else
    call uitsth &
    ( iflid  , cvar_var , smbrs  , rovsdt )
  endif
endif

call ustssc &
!==========
( nvar   , nscal  , ncepdp , ncesmp ,                            &
  iscal  ,                                                       &
  icepdc , icetsm , itypsm ,                                     &
  dt     ,                                                       &
  ckupdc , smacel , smbrs  , rovsdt )

! Store the source terms for convective limiter
call field_get_key_int(iflid, kst, st_id)
if (st_id .ge.0) then
  call field_get_val_s(st_id, cpro_scal_st)

  do iel = 1, ncel
    !Fill the scalar source term field
    cpro_scal_st(iel) = smbrs(iel)
  end do
  ! Handle parallelism and periodicity
  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(cpro_scal_st)
  endif
end if

if (vcopt%ibdtso.gt.ntinit.and.ntcabs.gt.1 &
    .and.(idtvar.eq.0.or.idtvar.eq.1)) then
  ! TODO: remove test on ntcabs and implemente a "proper" condition for
  ! initialization.
  f_id = ivarfl(ivar)
  call cs_backward_differentiation_in_time(f_id, smbrs, rovsdt)
endif
! Skip first time step after restart if previous values have not been read.
if (vcopt%ibdtso.lt.0) vcopt%ibdtso = iabs(vcopt%ibdtso)

! Set ibdtso value
call field_set_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

! Nudging towards optimal interpolation for current scalar
if (ippmod(iatmos).ge.0) then
  call field_get_key_int(ivarfl(ivar), kopint, f_oi_id)
  if (f_oi_id.ge.0) then
    call cs_at_data_assim_source_term(ivarfl(ivar), smbrs, rovsdt)
  endif
endif

! Si on extrapole les TS :
!   SMBRS recoit -theta TS du pas de temps precedent
!     (on aurait pu le faire avant ustssc, mais avec le risque que
!      l'utilisateur l'ecrase)
!   SMBRS recoit la partie du terme source qui depend de la variable
!   A l'ordre 2, on suppose que le ROVSDT fourni par l'utilisateur est <0
!     on implicite le terme (donc ROVSDT*RTPA va dans SMBRS)
!   En std, on adapte le traitement au signe de ROVSDT, mais ROVSDT*RTPA va
!     quand meme dans SMBRS (pas d'autre choix)
!print*,"cproa_scal_st before written", cproa_scal_st(50)
if (st_prv_id .ge. 0) then
  do iel = 1, ncel
    ! Stockage temporaire pour economiser un tableau
    smbexp = cproa_scal_st(iel)
    ! Terme source utilisateur explicite
    cproa_scal_st(iel) = smbrs(iel)
    ! Terme source du pas de temps precedent et
    ! On suppose -ROVSDT > 0 : on implicite
    !    le terme source utilisateur (le reste)
    smbrs(iel) = rovsdt(iel)*cvara_var(iel) - thets*smbexp
    ! Diagonale
    rovsdt(iel) = - thetv*rovsdt(iel)
  enddo
! Si on n'extrapole pas les TS :
else
  do iel = 1, ncel
    ! Terme source utilisateur
    smbrs(iel) = smbrs(iel) + rovsdt(iel)*cvara_var(iel)
    ! Diagonale
    rovsdt(iel) = max(-rovsdt(iel),zero)
  enddo
endif
!print*,"cproa_scal_st after written", cproa_scal_st(50)
! Add thermodynamic pressure variation for the low-Mach algorithm:
! NB: iscalt is the Enthalpy
if ((idilat.eq.3.or.ipthrm.eq.1).and.iscal.eq.iscalt) then
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel) + (pther - pthera)/dt(iel)*cell_f_vol(iel)
  enddo
endif


if (st_prv_id .ge. 0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel) + thetp1 * cproa_scal_st(iel)
  enddo
endif

! Low Mach compressible algos (conservative in time).
! Same algo. for Volume of Fluid method
if (idilat.gt.1 .or. ivofmt.ge.0) then
  call field_get_val_prev_s(icrom_scal, pcrom)
! Standard algo
else
  call field_get_val_s(icrom_scal, pcrom)
endif

! "VITESSE" DE DIFFUSION FACETTE

! On prend le MAX(mu_t,0) car en LES dynamique mu_t peut etre negatif
! (clipping sur (mu + mu_t)). On aurait pu prendre
! MAX(K + K_t,0) mais cela autoriserait des K_t negatif, ce qui est
! considere ici comme non physique.
if (vcopt%idiff.ge.1) then
  ! Scalar diffusivity
  if (vcopt%idften.eq.1) then

    idifftp = vcopt%idifft
    if (ityturt(iscal).eq.3) then
      idifftp = 0
    endif

    if (ifcvsl.lt.0) then
      do iel = 1, ncel
        w1(iel) = visls0(iscal)                                     &
           + idifftp*xcpp(iel)*max(visct(iel),zero)/turb_schmidt
      enddo
    else
      do iel = 1, ncel
        w1(iel) = cpro_viscls(iel)                                &
           + idifftp*xcpp(iel)*max(visct(iel),zero)/turb_schmidt
      enddo
    endif

    if (vcopt%iwgrec.eq.1) then
      ! Weighting for gradient
      do iel = 1, ncel
        cpro_wgrec_s(iel) = w1(iel)
      enddo
      call synsca(cpro_wgrec_s)
      if (irangp.ge.0.or.iperio.eq.1) then
         call synsca(cpro_wgrec_s)
      endif
    endif

    call viscfa &
    !==========
   ( imvisf ,                      &
     w1     ,                      &
     viscf  , viscb  )
  endif

  ! AFM model or DFM models: add div(Cp*rho*T'u') to smbrs
  ! Compute T'u' for GGDH

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

endif


! Not Darcy

if (ippmod(idarcy).eq.-1) then
  ! --> Unsteady term and mass aggregation term
  do iel = 1, ncel
    rovsdt(iel) = rovsdt(iel)                                                 &
                + vcopt%istat*xcpp(iel)*pcrom(iel)*cell_f_vol(iel)/dt(iel)
  enddo
endif

call field_get_key_int(iflid, keydri, iscdri)


!===============================================================================
! 3. Solving
!===============================================================================

iconvp = vcopt%iconv
idiffp = vcopt%idiff
idftnp = vcopt%idften
ndircp = ndircl(ivar)
nswrsp = vcopt%nswrsm
nswrgp = vcopt%nswrgr
imligp = vcopt%imligr
ircflp = vcopt%ircflu
ischcp = vcopt%ischcv
isstpp = vcopt%isstpc
iescap = 0
iswdyp = vcopt%iswdyn
iwarnp = vcopt%iwarni
blencp = vcopt%blencv
epsilp = vcopt%epsilo
epsrsp = vcopt%epsrsm
epsrgp = vcopt%epsrgr
climgp = vcopt%climgr
extrap = vcopt%extrag
relaxp = vcopt%relaxv
! all boundary convective flux with upwind
icvflb = 0

call field_get_coefa_s(ivarfl(ivar), coefap)
call field_get_coefb_s(ivarfl(ivar), coefbp)
call field_get_coefaf_s(ivarfl(ivar), cofafp)
call field_get_coefbf_s(ivarfl(ivar), cofbfp)
if(iter_sca .gt. 1) then

  do iel = 1, ncel
    vark(iel) = cvar_var(iel)
  enddo
endif

print*,"cvar_var,cvara_var, vark",cvar_var(50),cvara_var(50), vark(50),"before codits"
if (iter_sca .eq. 1) then
  call codits &
  !==========
   ( idtvar , ivarfl(ivar)    , iconvp , idiffp , ndircp ,          &
     imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
     ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
     iwarnp ,                                                       &
     blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
     relaxp , thetv  ,                                              &
     cvara_var       , cvara_var       ,                            &
     coefap , coefbp , cofafp , cofbfp ,                            &
     imasfl , bmasfl ,                                              &
     viscf  , viscb  , viscf  , viscb  , viscce ,                   &
     weighf , weighb ,                                              &
     icvflb , ivoid  ,                                              &
     rovsdt , smbrs  , cvar_var        , dpvar  ,                   &
     xcpp   , rvoid  )
else
     call codits &
  !==========
   ( idtvar , ivarfl(ivar)    , iconvp , idiffp , ndircp ,          &
     imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
     ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
     iwarnp ,                                                       &
     blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
     relaxp , thetv  ,                                              &
     cvara_var       , vark       ,                            &
     coefap , coefbp , cofafp , cofbfp ,                            &
     imasfl , bmasfl ,                                              &
     viscf  , viscb  , viscf  , viscb  , viscce ,                   &
     weighf , weighb ,                                              &
     icvflb , ivoid  ,                                              &
     rovsdt , smbrs  , cvar_var        , dpvar  ,                   &
     xcpp   , rvoid  )
endif
print*,"cvar_var,cvara_var,vark",cvar_var(50),cvara_var(50),vark(50),"after codits"


 xvarmu0 = 0.d0
 do iel=1, ncel
   xvarmu0 = xvarmu0 + (cvara_var(iel)**2)*cell_f_vol(iel)
 enddo


if (nterup.gt.1) then

  icvrge = 1

  xvartmp = 0.d0
  do iel = 1,ncel
    xdvar = cvar_var(iel) - vark(iel)
    xvartmp = xvartmp +(xdvar**2) * cell_f_vol(iel)
  enddo
  xvarmu = sqrt(xvartmp)
  if (xvarmu.ge.epsup*xvarmu0) icvrge = 0

endif


!===============================================================================
! 4. Writing and clipping
!===============================================================================

call clpsca(iscal)




! BILAN EXPLICITE (VOIR CODITS : ON ENLEVE L'INCREMENT)
! Ceci devrait etre valable avec le theta schema sur les Termes source

! Free memory
deallocate(w1)
deallocate(smbrs, rovsdt)
if (allocated(viscce)) deallocate(viscce)
if (allocated(weighf)) deallocate(weighf, weighb)
deallocate(dpvar)
deallocate(xcpp)
if (allocated(diverg)) deallocate(diverg)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'   ** RESOLUTION POUR LA VARIABLE ',A16                       ,/,&
'      ---------------------------                            ',/)
 1200 format(1X,A16,' : BILAN EXPLICITE = ',E14.5)
 9000 format( &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR DANS COVOFI                          ',/,&
'@    =========                                               ',/,&
'@    IVARSC DOIT ETRE UN ENTIER POSITIF STRICTEMENT          ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(/,                                                   &
'   ** SOLVING VARIABLE ',A16                                  ,/,&
'      ----------------'                                       ,/)
 1200 format(1X,A16,' : EXPLICIT BALANCE = ',E14.5)
 9000 format( &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ERROR IN COVOFI'                                ,/,&
'@    ========'                                                ,/,&
'@    IVARSC MUST BE A STRICTLY POSITIVE INTEGER'              ,/,&
'@    ITS VALUE IS ',I10                                       ,/,&
'@'                                                            ,/,&
'@  The calculation will not be run.'                          ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!----
! End
!----

return

end subroutine
