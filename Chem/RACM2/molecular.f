!!-----------------------------------------------------------------------
!!     Copyright (C) 2020 Laboratoire Chimie Environnement - CEREA (ENPC)
!!     The H2I is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------
C-----------------------------------------------------------------------
c     File: molecular.f (modified from Aerosol.f)
C-----------------------------------------------------------------------


c     Function: chen_othmer_diffusivity
c
c     Computes gas-phase diffusivity in air
c     using Chen & Othmer (1962) formula
c
c     Parameters:
c     temperature - Temperature [K]
c     pressure - Atmospheric pressure [Pa]
c     Tc - Critical temperature [K] 
c     Vc - Critical molar volume [cm3/mol] 
c     M - Molar weight [g/mol]
c
c     Returns:
c     diffusivity - Gas-phase diffusivity [m2/s].
      subroutine chen_othmer_diffusivity(diffusivity, temperature, 
     &   pressure, Tc_species, Vc_species, weight_species)

      double precision temperature
      double precision pressure
      double precision weight_air, weight_species
      double precision Vc_air, Vc_species
      double precision Tc_air, Tc_species
      parameter (Tc_air=132.63d0, Vc_air=92.35, weight_air=2.897d1)

      double precision diffusivity

      diffusivity = 4.3d-5 * (temperature/100.d0)**1.81d0
     $    * dsqrt((weight_air+weight_species)/weight_air/weight_species)
     $    / (pressure/101325.d0 * (Tc_air*Tc_species/10000.d0)**0.1406d0
     $    * ((Vc_air/100.d0)**0.4d0+(Vc_species/100.d0)**0.4d0)**2.d0)

      end



c     Function: compute_quadratic_mean_velocity
c
c     Computes condensation transfer rate for a given species.
c
c     Parameters:
c     temperature - Temperature [K]
c     weight - Molecular weight [g/mol]
c
c     Returns:
c     velocity - Quadratic mean molecular velocity [m/s]
      subroutine compute_quadratic_mean_velocity(temperature, weight,
     $     velocity)

      double precision temperature
      double precision weight
      double precision velocity

      velocity = dsqrt(2.11714271498563d4 * temperature / weight)

      end



c     Function: compute_deposition_velocity
c
c     Computes the transport-limited deposition velocity
c     Ref : Lai & Nazaroff 2000
c
c     Parameters:
c     u_friction - friction velocity on the surface [m/s]
c     D - diffusion coefficient [m2/s]
c
c     Returns:
c     vtrd - deposition velocity [m/s]
      subroutine compute_deposition_velocity(u_friction, nu, Dc, nchunk,
     &    vtrd)

      double precision nu
      double precision B0, B1, B2, B3
      double precision x1_a, x1_b, x2_a, x2_b, x3_a, x3_b
      double precision nt1_a, nt1_b, nt2_a, nt2_b, nt3_a, nt3_b
      double precision I1, I2, I3
      double precision u_friction, Dc, inv_vd, vtrd
      integer i, nchunk

      B0 = 0.d0
      B1 = 4.3d0
      B2 = 12.5d0
      B3 = 30.d0
      I1 = 0.d0
      I2 = 0.d0
      I3 = 0.d0

      do i=1,nchunk
         ! Computes interval boundaries
         x1_a = B0 + (i-1.d0)*(B1-B0)/nchunk
         x1_b = B0 + i*(B1-B0)/nchunk 
         x2_a = B1 + (i-1.d0)*(B2-B1)/nchunk
         x2_b = B1 + i*(B2-B1)/nchunk 
         x3_a = B2 + (i-1.d0)*(B3-B2)/nchunk
         x3_b = B2 + i*(B3-B2)/nchunk 

         ! Computes function values at boundaries
         nt1_a = nu * 7.669d-4 * x1_a**3.d0     
         nt1_b = nu * 7.669d-4 * x1_b**3.d0    
         nt2_a = nu * 1.0d-3 * x2_a**2.8214d0   
         nt2_b = nu * 1.0d-3 * x2_b**2.8214d0  
         nt3_a = nu * 1.07d-2 * x3_a**1.8895d0  
         nt3_b = nu * 1.07d-2 * x3_b**1.8895d0 
       
         ! Computes invd integrals  
         I1 = I1 + 0.5d0*(1.d0/(nt1_a+Dc)+ 1.d0/(nt1_b+Dc))*(x1_b-x1_a)
         I2 = I2 + 0.5d0*(1.d0/(nt2_a+Dc)+ 1.d0/(nt2_b+Dc))*(x2_b-x2_a)
         I3 = I3 + 0.5d0*(1.d0/(nt3_a+Dc)+ 1.d0/(nt3_b+Dc))*(x3_b-x3_a)

      end do

      ! Returns vtrd
      inv_vd = nu * (I1 + I2 + I3)
      vtrd = u_friction / inv_vd
      
      end
