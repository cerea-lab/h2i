!!-----------------------------------------------------------------------
!!     Copyright (C) 2020 Laboratoire Chimie Environnement - CEREA (ENPC)
!!     The H2I is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------
      PROGRAM photochemistry

      use globals 

      IMPLICIT NONE

      double precision, external :: GAUSS, VENLIGHT
      double precision, external :: SENLIGHT_ROOM, SENLIGHT_PAINT, SGAS
      double precision, external :: rho_air, nu_air
      double precision, external :: varIrrad

      integer ns, nr, nrphot
C     ns: number of species
C     nr: number of reactions
C     nrphot : number of photolytic reactions

      parameter (ns = 117, nr = 389, nrphot = 34) 

      integer option_chemistry 
C     Type of kinetic chemical mechanism 

      ! ------ Photolytic variables
      integer photLbox, photSbox
      integer, parameter :: photolysis=1, no_phot=0 
C     photLbox, photSbox: photolysis activation
C     # photolysis = photolytic reactions activated
C     # no_phot = photolytic reactions inactivated
      integer, parameter :: photoLight=1, photoShade=2 
C     Sets of photolysis constants in the sunlit or shaded box

      integer ns2, nr2, nrphot2
      integer nreactphot(nrphot)

      ! ------ Time variables
      DOUBLE PRECISION ts, tf, delta_t, dt, current_t
C     TS: initial time (GMT, computed from January 1st, [s])
C     TF: final time (GMT, computed from January 1st, [s])
C     DELTA_T, DT: time step of the main loop [s]
      INTEGER, PARAMETER ::  ncycle=1
C     NCYCLE: number of sub-cycles 

      integer option_adaptive_time_step
      double precision ATOL, tstep_min
C     tstep_min: minimum time step allowed to solve chemical reactivity [s] 

      ! ------ ETR solver variables
      DOUBLE PRECISION current_sub_time, current_sub_time_save
      DOUBLE PRECISION final_sub_time, dt_transport

      ! ------ Environment variables
      DOUBLE PRECISION DLtemp, DLhumid, DLpress
      DOUBLE PRECISION DLtempf, DLhumidf, DLpressf
C     DLHUMID, DLHUMIDF: specific humidity [kg/kg]
C     DLTEMP, DLTEMPF: temperature [K]
C     DLPRESS, DLPRESSF: pressure [Pa]
      DOUBLE PRECISION relative_humid, Psat, water_frac
C     PSAT: water vapor pressure [Pa]
C     RELATIVE_HUMID: relative humidity (RH)
C     WATER_FRAC: water molar fraction
      DOUBLE PRECISION YlH2O, SumM
C     YLH2O: Number of water molecules computed from the massic fraction (DLHUMID)
C     SUMM: total number of molecules
      DOUBLE PRECISION Tsurf
C     TSURF: temperature of the room surfaces [K]

      DOUBLE PRECISION kAER, Fbuilding
C     KAER: Air exchange rate beween In and Out [/s]
C     FBUILDING: Building filtration factor [-]

      ! ------ Air physical parameters
      double precision nu, rho
C     NU: air kinematic viscosity [m2/s]
C     RHO: air volumetric mass [kg/m3]
      double precision u_friction, dU, u_inf
C     UINF: airstream velocity in the room [m/s]
C     U_FRICTION: friction velocity [m/s]
C     dU: derivative of the mean air speed parallel to the surface [/s]

      double precision dlon, dlat
C     DLON, DLAT: latitude and longitude
      integer ndays
C     NDAYS: number of elapsed days since January 1st
      
      ! ------ Space variables
      DOUBLE PRECISION v_room, v_light, v_shd
      DOUBLE PRECISION s_room, sroom_light, sroom_shd
C     V_ROOM, S_ROOM: total volume [m3]/ surface [m2] of the room
C     V_LIGHT, SROOM_LIGHT: volume [m3]/ surface [m2] of the sunlit part 
C     V_SHD, SROOM_SHD: volume [m3]/ surface [m2] of the shaded part 
      DOUBLE PRECISION s_paint, spaint_light, spaint_shd
C     S_PAINT: room surface covered by paint [m2]
C     SPAINT_LIGHT/SHD: s_paint fraction in the direct/indirect light [m2]
      DOUBLE PRECISION s_rest, srest_light, srest_shd
C     S_REST: room surface not covered by paint [m2]
C     SREST_LIGHT/SHD: s_rest fraction in the direct/indirect light [m2]
      DOUBLE PRECISION s_borders_light
C     S_BORDERS_LIGHT: surface allowing gas transfer between boxes [m2]
      DOUBLE PRECISION kBOX_L, kBOX_S
C     K_BOX_L/S: exchange rate between the boxes [/s]

      ! ------ Chemical species
      character*6 species_name(ns)
      character*6, allocatable :: input_names(:)
C     SPECIES_NAME: names of the RACM2 compounds
C     INPUT_NAMES: names of the input compounds
      integer nic
      integer, allocatable :: ieq(:)
C     NIC: number of input compounds
C     IEQ: index correspondances between RACM2 list (species-list-RACM2.dat) 
C          and input compounds (init.dat)
      integer sim_init
C     SIM_INIT: Number of species with initialized concentration 

      ! ------ Species concentrations
      DOUBLE PRECISION DLconcL(ns) ! L = in the light
      DOUBLE PRECISION DLconcS(ns) ! S = in the shadow
C     DLCONC: arrays of modelled gas concentrations [\mu.g/m^3] 
C     # Before entry, it is given at initial time of the timestep
C     # On exit, it is computed at final time of the timestep
      DOUBLE PRECISION DLconc_tmpL(ns), DLconc_tmpS(ns)

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Omeas
C     OMEAS: boolean array of outdoor measurement availability: 
C            -1 = outdoor record available; constant value otherwise
      DOUBLE PRECISION Oconc(ns)
C     OCONC: outdoor concentration records [\mu.g/m^3]

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Cin_meas
C     Cin_meas: VOCs indoor concentration records [\mu.g/m^3]

      ! ------ Thermodynamical variables
      double precision navogadro
      parameter (navogadro = 6.02213d23)
C     navogadro: Avogadro's number

      double precision R
      parameter (R = 8.314d0)
C     R: perfect gas constant

      DOUBLE PRECISION molecular_weight(ns)
C     MOLECULAR_WEIGHT: molar mass [g/mol] 
      DOUBLE PRECISION convers_factor(ns)
      double precision convers_factor_jac(ns, ns)
      double precision Dc
C     DC: gas-phase diffusivity [m2/s]
      double precision omega
C     OMEGA: quadratic mean molecular or Thermal velocity [m/s]
      double precision vtrd
C     VTRD: transport-limited deposition velocity [m/s]
      double precision, allocatable :: Tc_species(:)
      double precision, allocatable :: Vc_species(:)
C     TC_SPECIES: species critical temperatures [K]
C     VC_SPECIES: species critical molar volumes [cm3/mol]
      double precision, allocatable :: sig(:)
      double precision, allocatable :: eps(:)
C     Lennard-Jones parameters : 
C     SIG: species collision diameters [Angstrom]
C     EPS: species Lennard-Jones potential parameters [K]

      ! ------ Room and paint emissions
      DOUBLE PRECISION Eroom(ns), EpaintS(ns), EpaintL(ns)
C     EROOM: room emission rates [\mu.g/s]
C     EPAINTL: paint surface emission rates in the direct light [\mu.g/m2/s]
C     EPAINTS: paint surface emission rates in the shadow [\mu.g/m2/s]
      INTEGER ptype
C     PTYPE: type of paint (0 = no TiO2; 3 = 3.5% TiO2)

      ! ------ Surface reactivity 
      double precision, allocatable :: k_react(:)
C     K_REACT: types of surface reactivity
C     # -1 = no deposition, or surface reactivity with reaction rate 
C            fixed by user
C     # 0  = depostion fully controlled by transport (inifinite uptake)
C     # 3  = deposition with parameterized uptake

      ! case k_react = 3
      double precision, allocatable :: uptake_ref(:)
C     UPTAKE_REF: uptake value at T=296K, I=20W/m2, RH=40%, [NO2]=40ppb
      double precision C_gas_L, C_gas_S
C     C_GAS_L/S: gas-phase species concentrations in box L/S
      double precision uptake_paint_L, uptake_paint_S
      double precision uptake_rest_L, uptake_rest_S
C     UPTAKE_PAINT: paramaterized paint uptake values 
C     UPTAKE_REST: parameterized uptake values of rest of room

      double precision k_vtrd_L, k_vtrd_S 
C     K_VTRD_L/S: transport-limited deposition rate in box L/S [/s]
      double precision k_collision_L, k_collision_S
C     K_COLLISION_L/S: surface adhesion rate in box L/S [/s]
      double precision, allocatable :: k_eff_L(:)
      double precision, allocatable :: k_eff_S(:)
C     K_EFF_L/S: heterogeneous kinetic constants in box L/S

      ! ------ Heterogeneous reactions
      integer nhtr
      integer, allocatable :: ieq_xn(:)
      character*20, allocatable :: xn_names(:)
C     NHTR: number of reactions with surfaces (deposition + reactivity)
C     IEQ_XN: index correspondance between RACM2 list (species-list-RACM2.dat)
C        and species undergoing heterogeneous reaction (critical_constants.dat)
C     XN_NAMES: names of heterogeneous reactions
      double precision kNO2ad, kHONOdiss
C     kNO2ad: NO2ad -> 0.5 HONOad + 0.5 HNO3ad
C     kHONOdiss: HONOad -> 0.5 NOad + 0.5 NO2ad
      double precision kNOconvert, kHONOconvert
C     kHONOconvert: HONOad + HNO3ad -> 2 NOad
C     kNOconvert: NOad + HNO3ad -> NO2ad + HONOad
      double precision kHONOevap, kNOevap, kNO2evap
C     kHONOevap, kNOevap, kNO2evap: molecular evaporation constants [/s/mlc]

      ! ------ Outdoor variables (unused) 
      integer, parameter ::  nemis = 0
      integer nemisspecies(nemis)
      double precision DLCsourc(nemis), DLCsourcf(nemis)
      double precision DLattenuation, DLattenuationf
      double precision DLCphotolysis_rates(nrphot)
      double precision DLCphotolysis_ratesf(nrphot)

      ! ------ Experimental records
      character*50 VOCs_record
      character*50 VOCs_time

      integer, parameter :: nCOVs=21
      integer :: nmeas
C     nCOVs: number of VOC species
C     nmeas: number of VOC measurements during simulation time

      ! ------ Filenames and folder paths
      character*80 path_config
C     PATH_CONFIG: path to the configuration file
      character*60 fileTRH
      character*80 path_T_RH
C     FILETRH: name of the file containing indoor temperature and RH profiles
C     PATH_T_RH: path to fileTRH
      character*60 fileInit
      character*80 path_simInit
C     FILEINIT: name of the file containing initial concentrations assessed 
C               from simulations
C     PATH_SIMINIT: path to fileInit
      character*80 oconc_path
C     OCONC_PATH: path to the outdoor concentration measurements
      character*80 hetrxn_path
C     HETRXN_PATH: path to the file containing types of reactivity with surfaces
C                  (k_react) and uptake values
      character*20 output_dir
C     OUTPUT_DIR: directory of the simulation outputs
      character*120 output_path_L
      character*120 output_path_S
      character*120 output_path_param
C     OUTPUT_PATH: paths to the output files

      ! ------ Run parameters
      character*20 nmexp
      character*50 test_name
C     NMEXP: name of the experiment
C     TEST_NAME: name of the simulation run

      ! ------ Reading variables
      character*20 read_name
      character*6 ic_name
      character*3 variable
      integer i, j, t
      integer ind, nt, Jsp
      integer t_read, t_old, nread
      double precision tmp, tmp2, tmp3
      double precision read_gamma_rest, Vc, Tc
      double precision C_read,C_old
      double precision valCin(nCOVs)
      double precision HONOad_L, HNO3ad_L, NO2ad_L, NOad_L
      double precision HONOad_S, HNO3ad_S, NO2ad_S, NOad_S
      logical search_t   
      logical file_exists
      integer errCode


 
      ! ===================
      ! === Input data
      ! ===================

      ! Read input data
      OPEN(12, file='CONTROL/which_exp.dat', status='old')
      READ(12,*)
      READ(12,*) nmexp
      READ(12,*)
      READ(12,*) variable, dt
      READ(12,*)
      READ(12,*) output_dir
      READ(12,*)
      READ(12,*) test_name 
      CLOSE(12)

      path_config = 'CONTROL/'//trim(nmexp)//'/config.dat'
      OPEN(12, file=path_config, status='old')
      READ(12,*)
      READ(12,*) variable, dlat
      READ(12,*) variable, dlon
      READ(12,*) variable, ndays
      READ(12,*)
      READ(12,*) variable, ts
      READ(12,*)
      READ(12,*) variable, tf
      READ(12,*)
      READ(12,*) fileTRH
      READ(12,*)
      READ(12,*) ptype
      READ(12,*)
      READ(12,*) variable, kAER
      READ(12,*) variable, v_room
      READ(12,*) variable, s_room
      READ(12,*) variable, s_paint
      READ(12,*)
      READ(12,*) fileInit
      CLOSE(12)

      path_config = 'CONTROL/'//trim(nmexp)//'/param.dat'
      OPEN(12, file=path_config, status='old')
      READ(12,*) 
      READ(12,*) variable, Fbuilding  
      READ(12,*) variable, u_inf
      READ(12,*) 
      READ(12,*) variable, kNO2ad 
      READ(12,*) variable, kNOconvert 
      READ(12,*) variable, kHONOconvert 
      READ(12,*) variable, kHONOevap 
      READ(12,*) variable, kNOevap 
      READ(12,*) variable, kNO2evap 
      READ(12,*) variable, kHONOdiss 
      READ(12,*) 
      READ(12,*) variable, HONOad_L
      READ(12,*) variable, HONOad_S
      READ(12,*) variable, HNO3ad_L
      READ(12,*) variable, HNO3ad_S
      READ(12,*) variable, NO2ad_L
      READ(12,*) variable, NO2ad_S
      READ(12,*) variable, NOad_L
      READ(12,*) variable, NOad_S
      CLOSE(12)

      kAER = kAER/3600.0d0

      write(*,*) "Latitude:", dlat, "Longitude:", dlon
      write(*,*) "Simulation duration (in s):", tf-ts
      write(*,*) "Time step (in s):", dt

      option_adaptive_time_step = 1
      ATOL = 0.001d0
      tstep_min = 1.d0 

      option_chemistry = 2 ! RACM2 

      call dimensions_racm2(ns2, nr2, nrphot2)

      do i = 1, nrphot
         nreactphot(i) = i - 1 ! no of photolysis ranges from 0 to nrphot 
      end do

      photLbox = photolysis
      photSbox = photolysis

      OPEN(12,file='CONTROL/species-list-racm2.dat',status='old')
      Do Jsp=1,ns
         READ(12,*) species_name(Jsp), molecular_weight(Jsp)
         ! The vectors given to etr_solver are initialized to zero
         DLconcL(Jsp) = 0.d0
         DLconcS(Jsp) = 0.d0
         DLconc_tmpL(Jsp) = 0.d0
         DLconc_tmpS(Jsp) = 0.d0
         Eroom(Jsp) = 0.d0
         EpaintL(Jsp) = 0.d0
         EpaintS(Jsp) = 0.d0
         Oconc(Jsp) = 0.d0
      Enddo
      CLOSE(12)

      ! Compute conversion factor
      Do Jsp=1,ns
         convers_factor(jsp) = navogadro * 1.d-12 / 
     $        molecular_weight(jsp)
      end do
      
      do i = 1, ns
         do j = 1, ns
            convers_factor_jac(i, j) = molecular_weight(i)
     $           / molecular_weight(j)
         end do
      end do

      ! Read meteo data.
      DLattenuation = 1.0d0
      DLattenuationf = DLattenuation

      DLpress = 101325.0d0
      DLpressf = DLpress

      write(*,*) "Attenuation:", DLattenuation
      write(*,*) "Pressure (kPa):", DLpress
      write(*,*) "File for temp(K) and relative humidity:", fileTRH

      do i = 1, nrphot
         DLCphotolysis_rates(i) = 0.d0
         DLCphotolysis_ratesf(i) = 0.d0
      end do


      ! Compute Lennard-Jones parameters
      ind = 0
      OPEN(12,file='CONTROL/critical_constants.dat',status='old')
      read(12,*)
      read(12,*) nhtr
      ALLOCATE(ieq_xn(nhtr))
      ALLOCATE(Tc_species(nhtr))
      ALLOCATE(Vc_species(nhtr))
      ALLOCATE(eps(nhtr))
      ALLOCATE(sig(nhtr))
      READ(12,*)
      do i = 1,nhtr
         read(12,*) ic_name, Tc, Vc
         Tc_species(i) = Tc
         Vc_species(i) = Vc
         do j = 1,ns
            if (species_name(j) .eq. ic_name) then
               ind = 1
               ieq_xn(i) = j 
            end if
         end do
         if (ind .eq. 0) then
            write(*,*) "Error"
            stop
         end if
         ! Force constants of the Lennard-Jones potential
         ! Ref : Flynn & Thodos 1962
         eps(i) = 2.36d0 * Tc / (Vc**0.25d0)    
         sig(i) = 0.618d0 * Vc**(1.d0/3.d0) * Tc**(1.d0/18.d0)
      end do
      CLOSE(12)

      ! Surface reactivity
      ALLOCATE(k_react(nhtr))
      ALLOCATE(k_eff_L(nhtr))
      ALLOCATE(k_eff_S(nhtr))
      ALLOCATE(uptake_ref(nhtr))
      ALLOCATE(xn_names(nhtr))
      hetrxn_path = 'CONTROL/'//trim(nmexp)//'/surface-react.dat'
      OPEN(12,file=hetrxn_path,status='old')
      read(12,*)
      do i = 1,nhtr
         read(12,*) read_name, tmp, tmp2
         xn_names(i) = read_name
         k_react(i) = tmp
         uptake_ref(i) = tmp2
      end do
      CLOSE(12)
 

      ! Initial indoor concentrations assessed from simulations
      path_simInit = 'CONTROL/'//trim(nmexp)//'/'//trim(fileInit)
     &  //'_L.dat'
      OPEN(12,file=path_simInit, status='old')
      READ(12,*) 
      READ(12,*) sim_init
      READ(12,*) 
      do i = 1, sim_init
         read(12,*) ic_name, tmp
         do j = 1,ns
            if (species_name(j) .eq. ic_name) then
               DLconcL(j) = tmp + 0.0001d0*tmp
            end if
         end do
      end do
      CLOSE(12)

      path_simInit = 'CONTROL/'//trim(nmexp)//'/'//trim(fileInit)
     &  //'_S.dat'
      OPEN(12,file=path_simInit, status='old')
      READ(12,*) 
      READ(12,*) sim_init
      READ(12,*) 
      do i = 1, sim_init
         read(12,*) ic_name, tmp
         do j = 1,ns
            if (species_name(j) .eq. ic_name) then
               DLconcS(j) = tmp - 0.0001d0*tmp
            end if
         end do
      end do
      CLOSE(12)


      ! Initial indoor concentrations measured in the room
      OPEN(12,file='CONTROL/'//trim(nmexp)//'/init.dat',status='old')
      READ(12,*) 
      READ(12,*) nic
      ALLOCATE(input_names(nic))
      ALLOCATE(ieq(nic))
      READ(12,*) 
c$$$      write(*,*) "number of ic species:", nic
      do i = 1, nic
         read(12,*) ic_name, tmp
         input_names(i) = ic_name
         do j = 1,ns
            if (species_name(j) .eq. ic_name) then
               DLconcL(j) = tmp + 0.0001d0*tmp
               DLconcS(j) = tmp - 0.0001d0*tmp
               ieq(i) = j 
            end if
         end do
      end do
      CLOSE(12)

      ! Outdoor concentrations availability
      ALLOCATE(Omeas(nic), STAT=errCode)
      OPEN(12,file='CONTROL/'//trim(nmexp)//'/outdoor.dat',status='old')
      READ(12,*) 
      do i = 1, nic
         read(12,*) ic_name, tmp
         Omeas(i) = tmp
      end do
      CLOSE(12)

      ! Number of time steps calculation
      if (dt < tstep_min) then
         write(*,*) "Warning: time step is set to ",tstep_min
         dt = tstep_min
      end if

      if (tf .lt. dt) then
         write(*,*) "tf should be greater than dt."
         stop
      end if

      nt = int((tf-ts) / dt)

      ! VOCs experimental records
      VOCs_record = 'CONTROL/'//trim(nmexp)//'/records/VOCs_conc.txt'
      VOCs_time = 'CONTROL/'//trim(nmexp)//'/records/VOCs_time.txt'
      nmeas = ceiling((tf-ts)/900.d0) ! 900s=15min

      ALLOCATE(Cin_meas(nCOVs,nmeas), STAT=errCode)
      OPEN(15,file=VOCs_record,status='old')
      do Jsp=1,nmeas
         read(15, *) valCin
         do i = 1,nCOVs ! COVs only
            Cin_meas(i,Jsp) = valCin(i)
         end do
      end do
      CLOSE(15)


      ! Creates output files 
      output_path_L = 'RESULTS/'//trim(output_dir)//'/'//trim(nmexp)//
     &   '_light_'//trim(test_name)//'.dat'
      inquire(file=output_path_L, exist=file_exists)
      OPEN(12, file=output_path_L)
      if (file_exists) then
            close(12, status='delete')
      else
            close(12)
      end if

      output_path_S = 'RESULTS/'//trim(output_dir)//'/'//trim(nmexp)//
     &   '_shade_'//trim(test_name)//'.dat'
      inquire(file=output_path_S, exist=file_exists)
      OPEN(12, file=output_path_S)
      if (file_exists) then
            close(12, status='delete')
      else
            close(12)
      end if

      output_path_param = 'RESULTS/'//trim(output_dir)//'/'//trim(nmexp)
     &   //'_param_'//trim(test_name)//'.dat'
      inquire(file=output_path_param, exist=file_exists)
      OPEN(12, file=output_path_param)
      if (file_exists) then
            close(12, status='delete')
      else
            close(12)
      end if


      ! Saves parameters
      OPEN(13, file=output_path_param)
      write(13,*) '--- Experiment : ', nmexp
      write(13,*)
      write(13,*) '---- Test name : ', test_name
      write(13,*) 'Fbuilding    = ', Fbuilding
      write(13,*) 'u_air        = ', u_inf
      write(13,*) 'kNO2ad       = ', kNO2ad 
      write(13,*) 'kNOconvert   = ', kNOconvert
      write(13,*) 'kHONOconvert = ', kHONOconvert
      write(13,*) 'kHONOevap    = ', kHONOevap
      write(13,*) 'kNOevap      = ', kNOevap
      write(13,*) 'kNO2evap     = ', kNO2evap
      write(13,*) 'kHONOdiss    = ', kHONOdiss
      do i=1,nhtr
         if (xn_names(i).EQ.'O3') then
            write(13,*) 'gamma_O3     = ', uptake_ref(i)
         end if
         if (xn_names(i).EQ.'NO2') then
            write(13,*) 'gamma_NO2    = ', uptake_ref(i)
         end if
         if (xn_names(i).EQ.'HONO') then
            write(13,*) 'gamma_HONO   = ', uptake_ref(i)
         end if
         if (xn_names(i).EQ.'NO') then
            write(13,*) 'gamma_NO     = ', uptake_ref(i)
         end if
      end do
      write(13,*) 
      write(13,*) 'Adsorbed species initial concentrations'

      do j=1,ns
         if (species_name(j) .eq. 'HONOad') then
            DLconcL(j) = HONOad_L
            write(13,*) 'initHONO_L   = ', DLconcL(j)
            DLconcS(j) = HONOad_S
            write(13,*) 'initHONO_S   = ', DLconcS(j) 
         end if
         if (species_name(j) .eq. 'HNO3ad') then
            DLconcL(j) = HNO3ad_L 
            write(13,*) 'initHNO3_L   = ', DLconcL(j) 
            DLconcS(j) = HNO3ad_S 
            write(13,*) 'initHNO3_S   = ', DLconcS(j) 
         end if
         if (species_name(j) .eq. 'NO2ad') then
            DLconcL(j) = NO2ad_L 
            write(13,*) 'initNO2_L    = ', DLconcL(j) 
            DLconcS(j) = NO2ad_S 
            write(13,*) 'initNO2_S    = ', DLconcS(j) 
         end if
         if (species_name(j) .eq. 'NOad') then
            DLconcL(j) = NOad_L 
            write(13,*) 'initNO_L     = ', DLconcL(j) 
            DLconcS(j) = NOad_S 
            write(13,*) 'initNO_S     = ', DLconcS(j) 
         end if
      end do
      CLOSE(13) 



      ! ===================
      ! === Main Loop
      ! ===================
      do t = 1, nt
        
         delta_t = dt 
         current_t = ts + delta_t * (t - 1)

C        # Reads temperature and relative humidity
         path_T_RH = 'CONTROL/'//trim(nmexp)//'/'//trim(fileTRH)
         OPEN(12,file=path_T_RH,status='old')
         search_t = .TRUE.
         do while ( search_t )
            READ(12,*) tmp, DLtemp, relative_humid
            if (AINT(tmp/10.0d0)*10 .EQ. current_t) then
               search_t = .FALSE.
            end if
         end do
         CLOSE(12)
         DLtempf = DLtemp

C        # Computes humidity 
         Psat = 611.2d0 * exp( 17.67d0 * ( DLtemp-273.15d0) / 
     &        (DLtemp - 29.65d0) )
         water_frac = relative_humid/100.d0 * Psat /DLpress 
         DLhumid = 0.62197d0 * water_frac / ( 1.0d0 + water_frac *
     &        (0.62197d0 - 1.0d0) )
         DLhumidf = DLhumid

         SumM = DLpress * 7.243D16 / DLtemp
         YlH2O = 29.d0*SumM*DLhumid/(18.d0+11.d0*DLhumid)


C        # Computes friction velocity
         rho = rho_air(DLtemp, DLpress)
         nu = nu_air(DLtemp, DLpress)
         dU = (0.074d0/rho/nu) * (rho*u_inf*u_inf/2.d0) * (u_inf*
     &          v_room**(1.d0/3.d0)/nu)**(-1.d0/5.d0)
         u_friction = SQRT(nu*dU)


C        # Box's surface and volume
         v_light = VENLIGHT(current_t)
         v_shd = v_room - v_light

         spaint_light = SENLIGHT_PAINT(current_t)
         spaint_shd = s_paint - spaint_light 

         sroom_light = SENLIGHT_ROOM(current_t)
         sroom_shd = s_room - sroom_light

         s_rest = s_room - s_paint
         srest_light = sroom_light - spaint_light
         srest_shd = sroom_shd - spaint_shd

         s_borders_light = SGAS(current_t)
         kBOX_L = u_inf * s_borders_light / v_light 
         kBOX_S = u_inf * s_borders_light / v_shd


C        # Prints variables
         write(*,*) '   t = ',(current_t-ndays*3600*24)/3600
         write(*,*) 'T=', DLtemp, ' RH=', relative_humid
         write(*,*) 'Light :  V=', v_light
         write(*,*) '   spaint=', spaint_light, 'sroom=', sroom_light
         write(*,*) 'Shadow :  V=', v_shd
         write(*,*) '   spaint=', spaint_shd, 'sroom=', sroom_shd
         write(*,*) "Gas exchanges : S=", s_borders_light
         write(*,*) "Friction velocity (cm/s) : ", u_friction*100

         Tsurf = (DLtemp - 273.15d0) * 1.2d0 + 273.15d0
         write(*,*) 'Tsurf=', Tsurf 
         write(*,*)


C        # Loop on input compounds
         do i=1,nic

C          ! Injections
           if (nmexp.EQ.'Exp1-Oct27')then
              if (((current_t-delta_t/2.d0).LT.25968600) .AND.
     &             ((current_t+delta_t/2.d0).GT.25968600)) then
                 if (input_names(i).EQ.'NO2') then
                    write(*,*) '15H30 INJECTION NO2'
                    DLconcS(ieq(i)) = DLconcS(ieq(i)) + 120.d0
                 end if
                 if (input_names(i).EQ.'NO') then
                    write(*,*) '15H30 INJECTION NO'
                    DLconcS(ieq(i)) = DLconcS(ieq(i)) + 2.d0
                 end if
              end if
           end if

C          ! Exchanges with outdoor
           if ((Omeas(i).GT.0) .OR. (Omeas(i).EQ.0)) then
              Oconc(ieq(i)) = Omeas(i)
           else
              oconc_path ='CONTROL/'//trim(nmexp)//'/outdoor_conc/'//
     $           trim(input_names(i))//'.txt'
              OPEN(12,file=oconc_path,status='old')
              READ(12,*) t_read, C_read
              if (current_t .EQ. ts) then
                 Oconc(ieq(i)) = C_read
              else
                 do while (t_read .LT. current_t )
                    t_old = t_read
                    C_old = C_read
                    READ(12,*) t_read, C_read
                 end do
                 Oconc(ieq(i))=C_old+(current_t-t_old)/(t_read-t_old)
     $              *(C_read-C_old)
              end if
              CLOSE(12)
           end if


C          ! VOCs assigned to experimental concentrations
           if (i.LT.nCOVs+1) then
              OPEN(12,file=VOCs_time,status='old')
              READ(12,*) t_read
              nread = 1
              if (current_t .EQ. ts) then
                 DLconcL(ieq(i)) = Cin_meas(i,nread)
                 DLconcS(ieq(i)) = Cin_meas(i,nread)
              else
                 do while (t_read .LT. current_t)
                    t_old = t_read
                    nread = nread + 1
                    READ(12,*) t_read
                 end do
                 DLconcL(ieq(i)) = Cin_meas(i,nread-1) + (current_t
     &              - t_old) / (t_read - t_old) * (Cin_meas(i,nread)
     &              - Cin_meas(i,nread-1))
                 DLconcS(ieq(i)) = DLconcL(ieq(i))
              end if
              CLOSE(12)
           end if
         end do


C        # Heterogeneous reactivity
         do i=1,nhtr
           if (k_react(i) .LT. 0) then
C             ! No deposition 
              k_eff_L(i) = 0
              k_eff_S(i) = 0

C             ! Tunable kinetic rates
              if (xn_names(i).EQ.'NO2ad') then
                 k_eff_L(i) = kNO2ad 
                 k_eff_S(i) = k_eff_L(i)
              end if
              if (xn_names(i).EQ.'NOconvert') then
                 k_eff_L(i) = kNOconvert 
                 k_eff_S(i) = k_eff_L(i)
              end if
              if (xn_names(i).EQ.'HONOconvert') then
                 k_eff_L(i) = kHONOconvert  
                 k_eff_S(i) = k_eff_L(i)
              end if
              if (xn_names(i).EQ.'HONOevap') then
                 k_eff_L(i) = kHONOevap*YlH2O 
                 k_eff_S(i) = kHONOevap*YlH2O 
              end if
              if (xn_names(i).EQ.'NOevap') then
                 k_eff_L(i) = kNOevap*YlH2O 
                 k_eff_S(i) = kNOevap*YlH2O 
              end if
              if (xn_names(i).EQ.'NO2evap') then
                 k_eff_L(i) = kNO2evap*YlH2O 
                 k_eff_S(i) = kNO2evap*YlH2O 
              end if
              if (xn_names(i).EQ.'HONOdiss') then
                 k_eff_L(i) = kHONOdiss 
                 k_eff_S(i) = k_eff_L(i) 
              end if

           else
C             ! Computes k_transport 
              call chen_othmer_diffusivity(Dc, DLtemp, DLpress,
     &            Tc_species(i), Vc_species(i),
     &            molecular_weight(ieq_xn(i)))

              call compute_deposition_velocity(u_friction,nu,Dc,4,vtrd)

              k_vtrd_L = vtrd * sroom_light / v_light  
              k_vtrd_S = vtrd * sroom_shd / v_shd 

              k_eff_L(i) = k_vtrd_L
              k_eff_S(i) = k_vtrd_S


              if (k_react(i) .GT. 0) then
                 call compute_quadratic_mean_velocity(DLtemp, 
     &              molecular_weight(ieq_xn(i)), omega)

                 if (k_react(i) .EQ. 3) then 
C                ! Computes k_collision 

                    ! Conversion µg/m3 -> ppb
                    C_gas_L = DLconcL(ieq_xn(i)) 
     &                 / (molecular_weight(ieq_xn(i))
     &                 * 1.d6) * R * DLtemp / DLpress * 1.0d9
                    C_gas_S = DLconcS(ieq_xn(i)) 
     &                 / (molecular_weight(ieq_xn(i))
     &                 * 1.d6) * R * DLtemp / DLpress * 1.0d9

                    ! Uptake of the paint boards
                    call param_gamma_NO2(uptake_ref(i), 
     &                  Tsurf, relative_humid, C_gas_L, 
     &                  ptype, uptake_paint_L)
                    call param_gamma_NO2(uptake_ref(i), DLtemp, 
     &                  relative_humid, C_gas_S, ptype, 
     &                  uptake_paint_S) 

                    ! Uptake of the rest of the room
                    call param_gamma_NO2(uptake_ref(i), Tsurf, 
     &                  relative_humid, C_gas_L, 0,
     &                  uptake_rest_L)
                    call param_gamma_NO2(uptake_ref(i), DLtemp,
     &                  relative_humid, C_gas_S, 0, 
     &                  uptake_rest_S) 


                    ! k_collision calculation
                    k_collision_L = ( (uptake_paint_L*spaint_light) +
     &                 (uptake_rest_L*srest_light) ) * omega /4.d0
     &                 /v_light 
                    k_collision_S = ( (uptake_paint_S*spaint_shd) +
     &                 (uptake_rest_S*srest_shd) ) * omega/4.d0/v_shd 

                 end if

                 if ((k_collision_L.EQ.0).OR.(k_collision_S.EQ.0)) then
                    write(*,*) "WARNING: k_collision = 0 for ", 
     &                 xn_names(i)
                 end if

                 k_eff_L(i) = 1.d0 /(1.d0/k_vtrd_L + 1.d0/k_collision_L)
                 k_eff_S(i) = 1.d0 /(1.d0/k_vtrd_S + 1.d0/k_collision_S)

              end if
           end if

         end do


C        # Small v_light
         if (v_light .LT. 1.0e-6) then
            do i=1,ns
               DLconcL(i) = 0 
            end do
         end if


C        // ------ SOLVER------ //
         current_sub_time = current_t
         final_sub_time = current_t + delta_t
         if (current_sub_time .EQ. ts) then
            dt_transport = DTAEROMIN
            call initstep(DLconcS, DLconcL, kBOX_L, kBOX_S, EpaintL, 
     $           EpaintS, Eroom, ns, ts, dt_transport, delta_t)
         end if
         ! ******Dynamic time loop
         write(*,*) 'Dynamic loop...'
         do while ( current_sub_time .lt. final_sub_time )

            current_sub_time_save = current_sub_time

            ! Computes concentration derivatives with sink and source terms
            call Etr_solver(DLconc_tmpL, DLconc_tmpS, DLconcL, DLconcS,
     $          ns, current_sub_time, dt_transport, kBOX_L, kBOX_S, 
     $          kAER, Fbuilding, Oconc, EpaintL, EpaintS, Eroom) 

            if (current_sub_time.le.final_sub_time) then
               call adaptime(DLconc_tmpL, DLconc_tmpS, DLconcL, DLconcS,
     $              ns, dt_transport, current_sub_time,final_sub_time)
            endif
 
C           ! Computes gas-phase chemistry in the sunlit box
            call chem (ns,nr,nrphot,nreactphot,nemis,nemisspecies,
     $       convers_factor,convers_factor_jac,current_sub_time_save,
     $       DLattenuation,DLhumid,
     $       DLtemp,DLpress,DLCsourc,DLCphotolysis_rates,dt_transport,
     $       DLattenuationf,DLhumidf,DLtempf,DLpressf,DLCsourcf,
     $       DLCphotolysis_ratesf,ncycle,dlon,dlat,DLconcL,
     $       option_adaptive_time_step, ATOL, tstep_min,
     $       photLbox, option_chemistry, photoLight, nhtr, k_eff_L)
       
C           ! Computes gas-phase chemistry in the shaded box
            call chem (ns,nr,nrphot,nreactphot,nemis,nemisspecies,
     $       convers_factor,convers_factor_jac,current_sub_time_save,
     $       DLattenuation,DLhumid,
     $       DLtemp,DLpress,DLCsourc,DLCphotolysis_rates,dt_transport,
     $       DLattenuationf,DLhumidf,DLtempf,DLpressf,DLCsourcf,
     $       DLCphotolysis_ratesf,ncycle,dlon,dlat,DLconcS,
     $       option_adaptive_time_step, ATOL, tstep_min,
     $       photSbox, option_chemistry, photoShade, nhtr, k_eff_S)

         end do


C        # WRITES OUTPUT FILES
         OPEN(12, file=output_path_L, position='append')
         write(12,*) "====== Simulation time :", current_t
         Do Jsp=1,ns
            WRITE(12,*) species_name(jsp), DLconcL(Jsp)
         Enddo
         write(12,*) ""
         CLOSE(12)

         OPEN(12, file=output_path_S, position='append')
         write(12,*) "====== Simulation time :", current_t
         Do Jsp=1,ns
            WRITE(12,*) species_name(jsp), DLconcS(Jsp)
         Enddo
         write(12,*) ""
         CLOSE(12)


      end do


     

      END
      

      DOUBLE PRECISION FUNCTION GAUSS(t,A,B,C)
      IMPLICIT NONE
      DOUBLE PRECISION A, B, C, pi, t
      pi = 4.d0*ATAN(1.d0)
      GAUSS = A/(B*SQRT(2.d0*pi))*EXP(-(t-C)*(t-C)/(2.d0*B*B)) 
      RETURN
      END

      DOUBLE PRECISION FUNCTION VENLIGHT(tu)
C     tu: GMT time [s] starting from January 1st
      IMPLICIT NONE
      DOUBLE PRECISION tu, j, h, A, B, C
      DOUBLE PRECISION, EXTERNAL :: GAUSS
      j = tu/86400.d0 ! time in days
      h = (j - AINT(j)) * 24.d0 ! time in hours
      A = 36.50515d0 
      B = 2.154006d0
      C = 11.19869d0
      VENLIGHT = GAUSS(h,A,B,C)
      RETURN
      END

      DOUBLE PRECISION FUNCTION SGAS(tu)
C     tu: GMT time [s] starting from January 1st
      IMPLICIT NONE
      DOUBLE PRECISION tu, j, h, A, B, C
      DOUBLE PRECISION, EXTERNAL :: GAUSS
      j = tu/86400.d0 ! time in days
      h = (j - AINT(j)) * 24.d0 ! time in hours
      A = 120.03678173d0
      B = 2.4683186d0
      C = 11.15402536d0
      SGAS = GAUSS(h,A,B,C)
      RETURN
      END

      DOUBLE PRECISION FUNCTION SENLIGHT_PAINT(tu)
C     tu: GMT time [s] starting from January 1st
      IMPLICIT NONE
      DOUBLE PRECISION tu, j, h, A1, B1, C1, G1, A2, B2, C2, G2
      DOUBLE PRECISION, EXTERNAL :: GAUSS
      j = tu/86400.d0 ! time in days
      h = (j - AINT(j)) * 24.d0 ! time in hours
      A1 = 9.133593d0
      B1 = 1.141670d0
      C1 = 9.530217d0
      G1 =  GAUSS(h,A1,B1,C1)
      A2 = 7.472820d0
      B2 = 1.523130d0
      C2 = 13.39981d0
      G2 = GAUSS(h,A2,B2,C2)
      SENLIGHT_PAINT = G1 + G2 - 2.103735d0*G1*G2
      RETURN
      END

      DOUBLE PRECISION FUNCTION SENLIGHT_ROOM(tu)
C     tu: GMT time [s] starting from January 1st
      IMPLICIT NONE
      DOUBLE PRECISION tu, j, h, A, B, C
      DOUBLE PRECISION, EXTERNAL :: GAUSS
      j = tu/86400.d0 ! time in days
      h = (j - AINT(j)) * 24.d0 ! time in hours
      A = 36.95750842d0
      B = 2.19501359d0
      C = 11.55451472d0
      SENLIGHT_ROOM = GAUSS(h,A,B,C)
      RETURN
      END

      DOUBLE PRECISION FUNCTION rho_air(temperature, pressure)
      IMPLICIT NONE
      DOUBLE PRECISION temperature, pressure, R, weight_air
      R = 8.314d0
      weight_air = 2.897d1
      rho_air = pressure / (R * temperature) * weight_air / 1000.d0
      RETURN
      END

      DOUBLE PRECISION FUNCTION nu_air(temperature, pressure)
C     Semi-empirical Sutherland relationship
      IMPLICIT NONE
      DOUBLE PRECISION temperature, pressure, eta_air, c, B
      DOUBLE PRECISION, EXTERNAL :: rho_air
      B = 1.46228d-6
      c = 113.d0
      eta_air = B * dsqrt(temperature) / (1.d0 + c / temperature)
      nu_air = eta_air / rho_air(temperature, pressure)
      RETURN
      END

