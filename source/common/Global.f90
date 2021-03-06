!-------------------------------------------------------------------------------
!   Copyright (C) 2017  Whitman T. Dailey
!   
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License version 3 as 
!   published by the Free Software Foundation.
!   
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!   
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
Module Global
    
    Use Kinds, Only: dp
    Implicit None
    Public
    
! Mathematical Constants
    Real(dp), Parameter :: pi = 3.141592653589793238462643383279502884_dp  !36 digits carried in case needed for quad precision
    Real(dp), Parameter :: TwoPi = 2._dp * pi
    Real(dp), Parameter :: FourPi = 4._dp * pi
    Real(dp), Parameter :: HalfPi = 0.5_dp * pi
    Real(dp), Parameter :: inv_TwoPi = 0.5_dp / pi
    Real(dp), Parameter :: inv_FourPi = 0.25_dp / Pi
    Real(dp), Parameter :: SqrtPi = Sqrt(pi)
    Real(dp), Parameter :: X_Hat(1:3) = (/ 1._dp , 0._dp , 0._dp /)
    Real(dp), Parameter :: Y_Hat(1:3) = (/ 0._dp , 1._dp , 0._dp /)
    Real(dp), Parameter :: Z_Hat(1:3) = (/ 0._dp , 0._dp , 1._dp /)
    Complex(dp), Parameter :: imag_i = ( 0._dp , 1._dp )

! Physical Constants
    Real(dp), Parameter :: m0c2 = 510.998928_dp  ![keV] Rest mass energy of an electron
    Real(dp), Parameter :: speed_of_light = 299792.458_dp ![km/s]
    Real(dp), Parameter :: keV_per_Joule = 6.241509126E15_dp  ![keV/J]
    Real(dp), Parameter :: g0 = 9.80665E-3_dp  ![km/s^2] standard acceleration due to gravity
    Real(dp), Parameter :: k_Boltzmann = 1.38064852E-23_dp  ![J/K] Boltzmann constant
    Real(dp), Parameter :: nA = 6.022140857E23_dp  ![1/mol] Avagadro's number
    Real(dp), Parameter :: h_Planck = 6.62607015E-34_dp  ![J*s] Planck constant
    Real(dp), Parameter :: h_bar_Planck = h_Planck / TwoPi  ![J*s] reduced Planck constant
    
! NEUTRON constants
    Real(dp), Parameter :: neutron_mass = 1.674927471E-27_dp  ![kg]
    Real(dp), Parameter :: neutron_mass_E = neutron_mass * speed_of_light**2 * 1000._dp**2 ![J]  Mass-energy of a neutron
    Real(dp), Parameter :: neutron_speed_conversion = &  !When multiplied by Sqrt(Energy in keV), gives speed of neutron in km/s
                               & Sqrt( (2._dp / keV_per_Joule) / (neutron_mass) ) / 1000._dp  ![km/s per Sqrt(keV)]
    Real(dp), Parameter :: neutron_PE_conversion = &  !When multiplied by change in height, gives change in PE in keV
                               & g0 * neutron_mass * keV_per_Joule  ![kg keV m/s^2 per J]
    Real(dp), Parameter :: mfp_per_barn_per_km_at_seaLevel = 5.07_dp
    ! With 1 barn microscopic cross section, a path length of 1 km in sea level air
    !   is an optical path length of 5.07 (mean free paths)
    !   So this is "mean free paths per barn per km at sea level"
    ! Note: EPL (effective path length) is the path length at sea level
    !   that has the same optical thickness as the actual path
    Real(dp), Parameter :: n_lambda = Log(2._dp) / 611._dp  ![1/s] ln(2) divided by neutron half-life
    Real(dp), Parameter :: n_kill_weight = &  !weight of 1 particular neutron having chance to exist from a 1MT source
                               & 1._dp / (0.35_dp * nA * 1000._dp)

!  Earth and Gravitation
    Real(dp), Parameter :: R_Earth = 6371._dp     ![km] Mean radius of earth
    Real(dp), Parameter :: std_grav_parameter = 398600.4418_dp  ![km^3 / s^2]  standard gravitational parameter of Earth
    Real(dp), Parameter :: Escape_speed = &  !minimum speed to escape earth's gravitational influence from the surface
                               & Sqrt(2._dp * std_grav_parameter / R_Earth)  ![km/s]
    Real(dp), Parameter :: rot_Earth = 7.292115E-5_dp  ![rad/s] Mean rotational speed of the Earth (radians per SIDEREAL second)
    Real(dp), Parameter :: Sid_Day_sec = 86164.09053_dp  !number of seconds in a sidereal day

End Module Global
