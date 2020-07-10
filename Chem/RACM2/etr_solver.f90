!!-----------------------------------------------------------------------
!!     Copyright (C) 2020 Laboratoire Chimie Environnement - CEREA (ENPC)
!!     The H2I is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------
module globals
  implicit none
  double precision, parameter:: TINYM=1.d-20
  double precision, parameter:: DTAEROMIN=1.d-5
  double precision, parameter:: EPSER=1.d-4
end module globals


subroutine Etr_solver(DLconc_tmpL, DLconc_tmpS, DLconc_L, DLconc_S, ns,& 
     current_sub_time, dt_transport, kBOX_L, kBOX_S, kAER, Fbuilding, Oconc, EpaintL, EpaintS, Eroom)
  !------------------------------------------------------------------------
  !
  !     -- DESCRIPTION
  !     This subroutine solves a system of Ordinary Differential Equations
  !     provided by the GDE for aerosols with the Explicit Trapezoidal Rule
  !     algorithm (ETR).
  !
  !------------------------------------------------------------------------
  !
  !     -- INPUT VARIABLES
  !
  !     ns : number of RACM2 species
  !     EpainL: paint emission rates in the light (ug/m^3/s)
  !     EpainS: paint emission rates in the shadow (ug/m^3/s)
  !     Eroom:  room emission rates (ug/m^3/s)
  !     kBOX : air exchange rate between the boxes (/s)
  !
  !     -- INPUT/OUTPUT VARIABLES
  !
  !     DLconc_L: concentration in the light box of second order evaluation (ug/m^3)
  !     DLconc_S: concentration in the shade box of second order evaluation (ug/m^3)
  !
  !     -- OUTPUT VARIABLES
  !
  !     DLconc_tmpL: concentration in the light box of first order evaluation (ug/m^3)
  !     DLconc_tmpS: concentration in the shade box of first order evaluation (ug/m^3)
  !
  !------------------------------------------------------------------------
  use globals
  implicit none

  integer:: j, ns
  double precision:: dCdt_tmpL(ns), dCdt_tmpS(ns) !1st order derivation
  double precision:: dCdt_L(ns), dCdt_S(ns) !2nd order derivation
  double precision:: DLconc_tmpL(ns), DLconc_tmpS(ns)!1st order concentration
  double precision:: DLconc_L(ns), DLconc_S(ns)!2nd order number concentration
  double precision:: current_sub_time,dt_transport,dtetr,tmp
  double precision:: EpaintL(ns), EpaintS(ns), Eroom(ns), Oconc(ns)
  double precision:: kBOX_L, kBOX_S, kAER, Fbuilding

  dcdt_L = 0
  dcdt_S = 0
  dcdt_tmpL = 0
  dcdt_tmpS = 0

  !     First step
  call fgde(DLconc_L, DLconc_S, dcdt_tmpL, dcdt_tmpS, kBOX_L, kBOX_S, kAER, Fbuilding, Oconc, EpaintL, EpaintS, Eroom, ns)
  do j=1,ns

     ! Light
     if(DLconc_L(j)+dcdt_tmpL(j)*dt_transport.GE.TINYM) then
        DLconc_tmpL(j)=DLconc_L(j)+dt_transport*dcdt_tmpL(j)
     else 
        DLconc_tmpL(j) = 0.0
     endif

     ! Shade
     if(DLconc_S(j)+dcdt_tmpS(j)*dt_transport.GE.TINYM) then
        DLconc_tmpS(j)=DLconc_S(j)+dt_transport*dcdt_tmpS(j)
     else 
        DLconc_tmpS(j) = 0.0
     endif
    
  enddo

  !     Second step
  call fgde(DLconc_tmpL, DLconc_tmpS, dcdt_L, dcdt_S, kBOX_L, kBOX_S, kAER, Fbuilding, Oconc, EpaintL, EpaintS, Eroom, ns)
  dtetr = dt_transport*5.0D-01
  current_sub_time = current_sub_time + dt_transport

  do j=1,ns
     ! Light
     tmp = dtetr * (dcdt_tmpL(j) + dcdt_L(j))
     if(DLconc_L(j)+tmp .GE. TINYM) then
        DLconc_L(j) = DLconc_L(j)+tmp
     else
        DLconc_L(j) = 0.0
     endif

     ! Shade
     tmp = dtetr * (dcdt_tmpS(j) + dcdt_S(j))
     if(DLconc_S(j)+tmp .GE. TINYM) then
        DLconc_S(j) = DLconc_S(j)+tmp
     else
        DLconc_S(j) = 0.0
     endif
  enddo

end subroutine Etr_solver

subroutine initstep(DLconc_L, DLconc_S, kBOX_L, kBOX_S, EpaintL, EpaintS, Eroom, ns,&
     initial_time_splitting, time_splitting, t_total)
  !------------------------------------------------------------------------
  !
  !     -- DESCRIPTION
  !     This subroutine performs the initialization of the timestep
  !     for the integration of the GDE. The criterion is related to
  !     the timescales of the aerosol processes.
  !
  !------------------------------------------------------------------------
  !
  !     -- INPUT VARIABLES
  !
  !     ns : number of RACM2 species
  !     t_total: total time limitation in current loop(s)
  !     DLconc_L: concentration in the light box of second order evaluation (ug/m^3)
  !     DLconc_S: concentration in the shade box of second order evaluation (ug/m^3)
  !
  !     -- OUTPUT VARIABLES
  !
  !     time_splitting: splitting time step
  !
  !------------------------------------------------------------------------
  use globals
  implicit none

  integer j, ns
  double precision:: tmp, tscale, time_splitting, initial_time_splitting
  double precision:: t_total 
  double precision:: DLconc_L(ns), DLconc_S(ns)
  double precision:: dcdt_tmpL(ns), dcdt_tmpS(ns)
  double precision:: EpaintL(ns), EpaintS(ns), Eroom(ns), Oconc(ns)
  double precision:: kBOX_L, kBOX_S, kAER, Fbuilding

  call fgde(DLconc_L, DLconc_S, dcdt_tmpL, dcdt_tmpS, kBOX_L, kBOX_S, kAER, Fbuilding, Oconc, EpaintL, EpaintS, Eroom, ns)

  time_splitting = t_total - initial_time_splitting
  do j=1,ns
     tmp = DLconc_L(j)*dcdt_tmpL(j)
     if (tmp.ne.0.D0) then
        tscale=DLconc_L(j)/DABS(dcdt_tmpL(j))
        time_splitting=DMIN1(time_splitting,tscale)
     endif
     tmp = DLconc_S(j)*dcdt_tmpS(j)
     if (tmp.ne.0.D0) then
        tscale=DLconc_S(j)/DABS(dcdt_tmpS(j))
        time_splitting=DMIN1(time_splitting,tscale)
     endif
  end do

  time_splitting=DMIN1(time_splitting,t_total-initial_time_splitting)
  time_splitting=DMAX1(time_splitting,DTAEROMIN)

end subroutine initstep

subroutine adaptime(DLconc_tmpL, DLconc_tmpS, DLconc_L, DLconc_S, ns,&
     T_dt,current_sub_time,final_sub_time)
  !------------------------------------------------------------------------
  !
  !     -- DESCRIPTION
  !     This subroutine computes N_aerosol time step for time
  !     integration, based on the difference between the
  !     first and second order evaluations
  !
  !------------------------------------------------------------------------
  !
  !     -- INPUT VARIABLES
  !
  !     ns : number of RACM2 species
  !     DLconc_tmpL: concentration in the light box of first order evaluation (ug/m^3)
  !     DLconc_tmpS: concentration in the shade box of first order evaluation (ug/m^3)
  !     DLconc_L: concentration in the light box of second order evaluation (ug/m^3)
  !     DLconc_S: concentration in the shade box of second order evaluation (ug/m^3)
  !
  !     -- OUTPUT VARIABLES
  !
  !     T_dt: dynamic time step for N_aerosol integration
  !
  !------------------------------------------------------------------------
  use globals
  implicit none

  integer j, ns
  double precision:: DLconc_tmpL(ns), DLconc_tmpS(ns)!1st order concentration
  double precision:: DLconc_L(ns), DLconc_S(ns)!2nd order number concentration
  double precision :: tmp,n2err
  double precision :: T_dt
  double precision :: current_sub_time,final_sub_time
  !     ******zero init
  n2err=0.d0
  !     ******local error estimation
  do j=1,ns
     if(DLconc_L(j).gt.TINYM) then
        tmp=(DLconc_L(j)-DLconc_tmpL(j))/(DLconc_L(j)+TINYM)
        n2err=n2err+tmp*tmp
     endif

     if(DLconc_S(j).gt.TINYM) then
        tmp=(DLconc_S(j)-DLconc_tmpS(j))/(DLconc_S(j)+TINYM)
        n2err=n2err+tmp*tmp
     endif
  end do

  n2err=DSQRT(n2err)

!!!     ******compute new time step
  ! formula to compute new time step
  if(n2err.NE.0.d0) then
     T_dt=T_dt*DSQRT(EPSER/n2err)
  else
     T_dt=final_sub_time-current_sub_time
  endif
  T_dt = DMAX1(DTAEROMIN, T_dt)

end subroutine adaptime

subroutine fgde(DLconcL, DLconcS, dCdt_L, dCdt_S, kBOX_L, kBOX_S, kAER, Fbuilding, Oconc, EpaintL, EpaintS, Eroom, ns)
  !------------------------------------------------------------------------
  !
  !     -- DESCRIPTION
  !     This subroutine provides entries for different aerosol dynamic process
  !
  !------------------------------------------------------------------------
  !
  !     -- INPUT VARIABLES
  !
  !     ns : number of RACM2 species
  !     EpainL: paint emission rates in the light (ug/m^3/s)
  !     EpainS: paint emission rates in the shadow (ug/m^3/s)
  !     Eroom:  room emission rates (ug/m^3/s)
  !     kBOX : air exchange rate between the boxes (/s)
  !     DLconc_L: concentration in the light box of second order evaluation (ug/m^3)
  !     DLconc_S: concentration in the shade box of second order evaluation (ug/m^3)
  !
  !     -- OUTPUT VARIABLES
  !
  !     dCdt_L: concentration derivation in the light box (ug/m^3/s)
  !     dCdt_S: concentration derivation in the shade box (ug/m^3/s)
  !
  !------------------------------------------------------------------------
  implicit none
  double precision :: DLconcS(ns), DLconcL(ns), Oconc(ns)
  double precision :: Eroom(ns), EpaintL(ns), EpaintS(ns)
  double precision :: kBOX_L, kBOX_S, kAER, Fbuilding
  double precision :: dCdt_L(ns), dCdt_S(ns)
  integer c, ns

  do c=1,ns

     if ((c.NE.38).AND.(c.NE.39).AND.(c.NE.40).AND.(c.NE.41)) then 
        dCdt_L(c) = kBOX_L * (-DLconcL(c) + DLconcS(c)) + Eroom(c) + EpaintL(c)&
          - kAER * DLconcL(c) + Fbuilding * kAER * Oconc(c)
        dCdt_S(c) = kBOX_S * ( DLconcL(c) - DLconcS(c)) + Eroom(c) + EpaintS(c)&
          - kAER * DLconcS(c) + Fbuilding * kAER * Oconc(c)
     else ! remove sink and source terms for adsorbed species
        dCdt_L(c) = 0.d0
        dCdt_S(c) = 0.d0
     end if

  end do

end subroutine fgde
