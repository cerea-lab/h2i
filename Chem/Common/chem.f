!!-----------------------------------------------------------------------
!!     Copyright (C) 2020 Laboratoire Chimie Environnement - CEREA (ENPC)
!!     The H2I is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------

      SUBROUTINE chem (ns,nr,nrphot,nreactphot,nemis,nemisspecies,
     $     convers_factor,convers_factor_jac,ts,DLattenuation,DLhumid,
     $     DLtemp,DLpress,DLCsourc,DLCphotolysis_rates,delta_t,
     $     DLattenuationf,DLhumidf,DLtempf, DLpressf,DLCsourcf,
     $     DLCphotolysis_ratesf,ncycle,dlon,dlat,DLconc,
     $     option_adaptive_time_step, ATOL, tstep_min,
     $     option_photolysis, option_chemistry, islight, nhtr, k_eff)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION
C     
C     This routine computes one timestep for gas-phase chemistry RACM2.
C     
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C     
C     TS: initial time (GMT, computed from January 1st, [s]).
C     DLATTENUATION: cloud attenuation at initial time.
C     DLHUMID: specific humidity at initial time ([%]).
C     DLTEMP: temperature at initial time ([K]).
C     DLPRESS: pressure at initial time ([Pa]).
C     DLCSOURC: array of chemical volumic emissions at initial time
C     # ([\mu.g/m^3/s]).
C     DLCPHOTOLYSIS_RATES: photochemical kinetic rates
C     # at initial time ([s^{-1}]).
C     DELTA_T: time step ([s]).
C     option_adaptive_time_step: 1 if adaptive time step.
C     ATOL:  relative tolerance for deciding if the time step is kept.
C     tstep_min: minimum time step.
C     The same variables are defined at final time of the timestep.
C     'f' is then put at the end of the name.
C     
C     -- INPUT/OUTPUT VARIABLES
C     
C     DLCONC: array of chemical concentrations ([\mu.g/m^3]).
C     # Before entry, it is given at initial time of the timestep.
C     # On exit, it is computed at final time of the timestep.
C     
C     
C------------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION ts,delta_t
      DOUBLE PRECISION tschem,tfchem

      integer ns,nr,nrphot,nemis,nhtr

      DOUBLE PRECISION DLconc(ns),ZC(ns)
      DOUBLE PRECISION DLtemp,DLtempf
      DOUBLE PRECISION DLattenuation
      DOUBLE PRECISION DLattenuationf
      DOUBLE PRECISION DLhumid,DLhumidf
      DOUBLE PRECISION DLCsourc(Nemis)
      DOUBLE PRECISION DLCsourcf(Nemis)
      DOUBLE PRECISION ZCsourc(ns)
      DOUBLE PRECISION ZCsourcf(ns)
      DOUBLE PRECISION DLRki(Nr),DLRkf(Nr)
      DOUBLE PRECISION DLpress,DLpressf
      DOUBLE PRECISION DLCphotolysis_rates(NRphot)
      DOUBLE PRECISION DLCphotolysis_ratesf(NRphot)
      DOUBLE PRECISION k_eff(nhtr)

      double precision dlon,dlat

      integer ncycle
      double precision convers_factor(ns)
      double precision convers_factor_jac(ns,ns)

      DOUBLE PRECISION Zangzen,Zangzenf,angzen
      DOUBLE PRECISION Zatt,Zattf

      DOUBLE PRECISION muzero,DLmuzero
      EXTERNAL muzero

      double precision pi

      integer nreactphot(nrphot)
      integer nemisspecies(nemis)

      INTEGER Jt,Jsp,i

      DOUBLE PRECISION DLk1(ns), DLk2(ns),ZC_old(ns)
      DOUBLE PRECISION EPSDLK
      PARAMETER (EPSDLK = 1.D-15)
      DOUBLE PRECISION supEdtstep, Edtstep(ns),ATOL
      DOUBLE PRECISION tstep,tstep_new,tstep_min,tfchem_tmp    

      INTEGER option_adaptive_time_step
      INTEGER option_photolysis, option_chemistry
      INTEGER islight

C     Constants.
      pi = 3.14159265358979323846D0

C     Spatial extraction for volumic sources.

      DO Jsp=1,ns
         ZCsourc(jsp)=0.D0
         ZCsourcf(jsp)=0.d0
      ENDDO

      DO Jsp=1,Nemis
         ZCsourc(nemisspecies(jsp)+1)=DLCsourc(Jsp)
         ZCsourcf(nemisspecies(jsp)+1)=
     $        DLCsourcf(Jsp)
      ENDDO

C     Cloud attenuation.

      Zatt = DLattenuation
      Zattf = DLattenuationf

C     Projection.

      ! ==== Before chem.f ===
      DO Jsp=1,ns
         ZC(Jsp) = DLconc(Jsp)
      ENDDO

C     Integration of chemistry (eventually with subcycling).

      DO Jt=1,Ncycle
         tschem=ts+(Jt-1)*delta_t/Ncycle
         tfchem=tschem+delta_t/Ncycle
         tstep = delta_t
         tfchem_tmp = tfchem


         Do while (tschem.LT.tfchem) 
C     If option_photolysis is 1,
C     photolytic reactions are calculated in kinetic.f
            DLmuzero=muzero(tschem,Dlon,Dlat)
            Zangzen=dabs(DACOS(DLmuzero)*180.D0/PI)
            DLmuzero=muzero(tfchem_tmp,Dlon,Dlat)
            Zangzenf=dabs(DACOS(DLmuzero)*180.D0/PI)


            IF (islight.eq.1) then
               CALL Kinetic_racm2_light(nr,DLRKi,DLtemp,DLhumid,DLpress,
     s           Zangzen,Zatt,option_photolysis,nhtr,k_eff)
               CALL Kinetic_racm2_light(nr,DLRKf,DLtempf,DLhumidf,
     s           DLpressf,Zangzenf,Zattf,option_photolysis,nhtr,k_eff)
            ELSE
               CALL Kinetic_racm2_shade(nr,DLRKi,DLtemp,DLhumid,DLpress,
     s           Zangzen,Zatt,option_photolysis,nhtr,k_eff)
               CALL Kinetic_racm2_shade(nr,DLRKf,DLtempf,DLhumidf,
     s           DLpressf,Zangzenf,Zattf,option_photolysis,nhtr,k_eff)
            END IF


C     If option_photolysis is 2,
C     photolytic reactions may be read.
            IF (option_photolysis.eq.2) then
            DO i=1,Nrphot
               DLRKi(Nreactphot(i)+1) = Zatt *
     $              DLCphotolysis_rates(i)
               DLRKf(Nreactphot(i)+1) = Zattf *
     $              DLCphotolysis_ratesf(i)
            ENDDO
            ENDIF

C     Solve gas-phase chemistry for the time step
            CALL roschem(ns, nr, ZC,ZCsourc,ZCsourcf,
     $           convers_factor, convers_factor_jac,tschem
     $           ,tfchem_tmp,DLRki,DLRkf,ZC_old,DLK1,DLK2,
     $           option_chemistry)

C     Integration of chemistry with adaptive time stepping
            IF(option_adaptive_time_step.EQ.1) then
C     Check that the time step was ok
               supEdtstep = 0.D0
               Do Jsp = 1,ns
                  If((DLK1(Jsp).GT.EPSDLK
     &                 .OR.DLK2(Jsp).GT.EPSDLK)
     &                 .AND.ZC(Jsp).GT.EPSDLK) then
                                ! Estimate the relative error
                     Edtstep(Jsp) = 0.5D0 * 
     &                    dabs(DLk1(Jsp) + DLk2(Jsp))
     &                    / ZC(Jsp)
                     If(Edtstep(Jsp).GT.supEdtstep) then 
                        supEdtstep = Edtstep(Jsp)
                     Endif
                  Endif
               Enddo
               supEdtstep = supEdtstep/ATOL
               If(supEdtstep.GT.1.D0
     &              .AND.tstep.GT.tstep_min) then
                                ! The time step is rejected and the computation
                                ! is redone with a smaller time step
                  tstep_new = tstep * 0.9d0 /dsqrt(supEdtstep)
                  tstep_new = DMAX1(tstep_new,tstep_min)
                  tfchem_tmp = tschem + tstep_new
                  If(tfchem_tmp.GT.tfchem) then
                     tfchem_tmp = tfchem
                     tstep_new = tfchem_tmp - tschem
                  Endif
                  tstep = tstep_new
                  Do Jsp=1,ns
                     ZC(Jsp) = ZC_old(Jsp)
                  Enddo
               Else
                                ! The time step is accepted and the time is incremented
                  tschem = tfchem_tmp
                  tfchem_tmp = tfchem_tmp + tstep
                  If(tfchem_tmp.GT.tfchem) then
                     tfchem_tmp = tfchem
                     tstep = tfchem_tmp - tschem
                  Endif
               Endif
            ELSE
               tschem = tfchem
            ENDIF

         Enddo                  !End loop Do while for time stepping

      ENDDO

C     Storage in the 3D array of chemical concentrations.

      angzen = Zangzenf

      ! === After chem.f ===
      DO i=1,ns
         DLconc(i) = ZC(i)
      ENDDO

      END
