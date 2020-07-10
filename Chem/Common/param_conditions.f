!!-----------------------------------------------------------------------
!!     Copyright (C) 2020 Laboratoire Chimie Environnement - CEREA (ENPC)
!!     The H2I is distributed under the GNU General Public License v3
!!-----------------------------------------------------------------------
c     File: param_cond.f C-----------------------------------------------------------------------


c     Function: gamma_NO2_NO2
c     Parameterization as a function of NO2 concentration
      subroutine gamma_NO2_NO2(gamma_ref, N_ppb, T, gamma_no2)
      double precision gamma_ref, gamma_norm, gamma_no2
      double precision M_NO2, Patm, R, T, a, b, c, d, N_ppb

      if ((N_ppb .LT. 0) .OR. (N_ppb .GT. 180)) then
         write(*,*) 'WARNING: [NO2] out of bounds in gamma_NO2_NO2'
         write(*,*) N_ppb
      end if

      a = 118.06327844053042d0
      b = 20.405614604432444d0
      d = 0.6054352683691557d0
      gamma_norm = a * EXP(-(N_ppb - (b*LOG(5/a)))/b) + d
      gamma_no2 = gamma_ref * gamma_norm
      end


c     Function: gamma_NO2_HUMID
c     Parameterization as a function of relative humidity
      subroutine gamma_NO2_HUMID(gamma_ref, H, ptype, gamma_humid)
      double precision H, gamma_ref, gamma_norm, gamma_humid
      integer ptype

      if ((H .LT. 0) .OR. (T .GT. 70)) then
         write(*,*) 'WARNING: RH out of bounds in gamma_NO2_HUMID'
      end if

      gamma_norm = -2.3108196675682176d-4*H*H + 1.4983184605993079d-2*H
     $       + 7.056900415326914d-1
      gamma_humid = gamma_ref * gamma_norm
      end


c     Function: gamma_NO2_TEMP
c     Parameterization as a function of surface temperature
      subroutine gamma_NO2_TEMP(gamma_ref, Tsurf, ptype, gamma_temp)
      double precision Tsurf, gamma_ref, gamma_norm, gamma_temp
      integer ptype
      if ((Tsurf .LT. 296) .OR. (Tsurf .GT. 313)) then
         write(*,*) 'WARNING: temp out of bounds in gamma_NO2_TEMP (',
     &     Tsurf, ' K )'
      end if
      
      if (ptype .EQ. 0) then
         gamma_norm = 1 
      else if (ptype .EQ. 3) then
         gamma_norm = 0.06250165d0*Tsurf - 17.61740619d0
      else
         write(*,*) 'WARNING: unknown type of paint in gamma_NO2_TEMP'
      end if
      gamma_temp = gamma_ref * gamma_norm
      end


c     Function: param_gamma_NO2
c     
c     Parameters:
c     gamma_ref - reference uptake measured at T=296K, RH=40%,
c                                              [NO2]=40ppb, I=20W/m2
c     T - room temperature [K]
c     H - relative humidity [%]
c     N - NO2 concentration [µg/m3]
c
c     Returns: gamma_par
      subroutine param_gamma_NO2(gamma_ref, T, H, N, ptype, gamma_par)
      double precision gamma_ref, gamma_par, T, H, N
      double precision gamma_temp, gamma_humid
      integer ptype
      call gamma_NO2_HUMID(gamma_ref, H, ptype, gamma_humid)
      call gamma_NO2_TEMP(gamma_humid, T, ptype, gamma_temp)
      call gamma_NO2_NO2(gamma_temp, N, T, gamma_par)
      end







