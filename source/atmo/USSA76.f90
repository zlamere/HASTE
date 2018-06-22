Module US_Std_Atm_1976
    
    Use Kinds, Only: dp
    Implicit None    
    Private
    Public :: find_base_layer
    Public :: T
    Public :: P
    Public :: rho    
    Public :: Zb

    !US Standard Atmosphere 1976 parameters
    !The following constants are defined here to ensure consistency with 1976 atmosphere model definition.
    !There may be more modern values for these constants, but these values ensure agreement with the 1976 US Std Atmosphere model
    Real(dp), Parameter :: g0 = 9.80665_dp  ![m / s^2]  accelleration due to gravity
    Real(dp), Parameter :: R_Earth = 6356.766_dp  ![km]  Radius of earth (nominal) at 45 deg latitude, used to relate geometric and geopotential heights, from US Std Atmosphere 1976
    Real(dp), Parameter :: R_star = 8.31432_dp  ![J / (mol*K)]  Universal gas constant as defined in US Standard Atmosphere 1976
    Real(dp), Parameter :: M0 = 28.964425912034_dp  ![kg / kmol] Average molecular weight of the ten most abundant species in air, weighted by relative abundance, from US Std Atmosphere 1976
    Real(dp), Parameter :: Hb(0:7) = (/ 0._dp, &  ![km] geopotential heights of layer boundaries, US Standard Atmosphere 1976 table 4 
                                      & 11._dp, & 
                                      & 20._dp, & 
                                      & 32._dp, & 
                                      & 47._dp, & 
                                      & 51._dp, & 
                                      & 71._dp, &
                                      & 86._dp * R_Earth / (86._dp + R_Earth) /)
    Real(dp), Parameter :: Zb(0:11) = (/ 0._dp, &  ![km] geometric heights of layer boundaries, US Standard Atmosphere 1976 
                                       & Hb(1) * R_Earth / (R_Earth - Hb(1)), & 
                                       & Hb(2) * R_Earth / (R_Earth - Hb(2)), & 
                                       & Hb(3) * R_Earth / (R_Earth - Hb(3)), & 
                                       & Hb(4) * R_Earth / (R_Earth - Hb(4)), & 
                                       & Hb(5) * R_Earth / (R_Earth - Hb(5)), & 
                                       & Hb(6) * R_Earth / (R_Earth - Hb(6)), &
                                       & 86._dp, &
                                       & 91._dp, &
                                       & 110._dp, &
                                       & 120._dp, &
                                       & 1000._dp /)
    Real(dp), Parameter :: Lb(0:11) = (/ -6.5_dp, &  ![K/km] temperature lapse rates in each layer, US Standard Atmosphere 1976 table 4 
                                       & 0._dp, & 
                                       & 1._dp, & 
                                       & 2.8_dp, & 
                                       & 0._dp, & 
                                       & -2.8_dp, & 
                                       & -2._dp, &
                                       & 0._dp, &
                                       & 0._dp, &
                                       & 12._dp, &
                                       & 0._dp, &
                                       & 0._dp /)
    Real(dp), Parameter :: Tb(0:11) = (/ 288.15_dp, &  ![K] Computed temperature at layer boundaries 
                                      & 216.65_dp, & 
                                      & 216.65_dp, & 
                                      & 228.65_dp, & 
                                      & 270.65_dp, & 
                                      & 270.65_dp, & 
                                      & 214.65_dp, &
                                      !& 186.9459083101885122_dp, & !<--value is for when M0 correction is NOT accounted for
                                      !& 186.9459083101885122_dp, & !<--value is for when M0 correction is NOT accounted for
                                      & 186.8671666936082608_dp, & !<--value is for when M0 correction is accounted for
                                      & 186.8671666936082608_dp, & !<--value is for when M0 correction is accounted for
                                      & 240._dp, &
                                      & 360._dp, &
                                      & 1000._dp /)
    Real(dp), Parameter :: Pb(0:7) = (/ 101325._dp, &   ![Pa] Computed pressure at layer boundaries
                                      & 22632.0336238972840275_dp, & 
                                      & 5474.87437675730708586_dp, & 
                                      & 868.014988510785148131_dp, & 
                                      & 110.905629143702212828_dp, & 
                                      & 66.9384346263881217465_dp, & 
                                      & 3.95638449983647254755_dp, &
                                      !& 0.37337628269333201966_dp /) !<--value is for when M0 correction is NOT accounted for
                                      & 0.37069900314554960416_dp /)  !<--value is for when M0 correction is accounted for
    Logical, Parameter :: Lb_nonzero(0:11) = (/ .TRUE., & 
                                              & .FALSE., & 
                                              & .TRUE., & 
                                              & .TRUE., & 
                                              & .FALSE., & 
                                              & .TRUE., & 
                                              & .TRUE., &
                                              & .FALSE., &
                                              & .FALSE., &
                                              & .TRUE., &
                                              & .FALSE., &
                                              & .FALSE.  /)  !flags indicating non-zero lapse rate
    Logical, Parameter :: T_linear_by_H(0:11) = (/ .TRUE., & 
                                                 & .FALSE., & 
                                                 & .TRUE., & 
                                                 & .TRUE., & 
                                                 & .FALSE., & 
                                                 & .TRUE., & 
                                                 & .TRUE., &
                                                 & .FALSE., &
                                                 & .FALSE., &
                                                 & .FALSE., &
                                                 & .FALSE., &
                                                 & .FALSE.  /)  !flags indicating linear temperature by geopotential height
    Logical, Parameter :: T_elliptical(0:11) = (/ .FALSE., & 
                                                & .FALSE., & 
                                                & .FALSE., & 
                                                & .FALSE., & 
                                                & .FALSE., & 
                                                & .FALSE., & 
                                                & .FALSE., &
                                                & .FALSE., &
                                                & .TRUE., &
                                                & .FALSE., &
                                                & .FALSE., &
                                                & .FALSE.  /)  !flags indicating elliptical temperature by geometric height
    Logical, Parameter :: T_exponential(0:11) = (/ .FALSE., & 
                                                 & .FALSE., & 
                                                 & .FALSE., & 
                                                 & .FALSE., & 
                                                 & .FALSE., & 
                                                 & .FALSE., & 
                                                 & .FALSE., &
                                                 & .FALSE., &
                                                 & .FALSE., &
                                                 & .FALSE., &
                                                 & .TRUE., &
                                                 & .FALSE.  /)  !flags indicating exponential temperature by geometric height
    Logical, Parameter :: P_rho_not_by_N(0:11) = (/ .TRUE., & 
                                                  & .TRUE., & 
                                                  & .TRUE., & 
                                                  & .TRUE., & 
                                                  & .TRUE., & 
                                                  & .TRUE., & 
                                                  & .TRUE., &
                                                  & .FALSE., &
                                                  & .FALSE., &
                                                  & .FALSE., &
                                                  & .FALSE., &
                                                  & .FALSE.  /)  !flags indicating Pressure and Density computed by OTHER than number density
    Real(dp), Parameter :: Tb_minus_LbHb(0:7) = Tb(0:7) - Lb(0:7)*Hb(0:7)  !precomputed quantity for 1976 temperature calculations
    Real(dp), Parameter :: L_star = g0 * M0 / R_star  !precomputed quantity for 1976 pressure calculations
    Real(dp), Parameter :: L_star_Lb(0:7) = (/ L_star / Lb(0), &  !precomputed quantity for 1976 pressure calculations
                                             & 0._dp, & 
                                             & L_star / Lb(2), & 
                                             & L_star / Lb(3), & 
                                             & 0._dp, & 
                                             & L_star / Lb(5), & 
                                             & L_star / Lb(6), &
                                             & 0._dp /)
    Real(dp), Parameter :: Pb_Tb_L_star_Lb(0:7) = (/ Pb(0) * Tb(0)**L_star_Lb(0), &  !precomputed quantity for 1976 pressure calculations
                                                   & Pb(1), & 
                                                   & Pb(2) * Tb(2)**L_star_Lb(2), & 
                                                   & Pb(3) * Tb(3)**L_star_Lb(3), & 
                                                   & Pb(4), & 
                                                   & Pb(5) * Tb(5)**L_star_Lb(5), & 
                                                   & Pb(6) * Tb(6)**L_star_Lb(6), &
                                                   & Pb(7) /)
    Real(dp), Parameter :: L_star_Tb(0:7) = -L_star / Tb(0:7)  !precomputed quantity for 1976 pressure calculations
    Real(dp), Parameter :: rho_star = M0 / R_star  !precomputed quantity for 1976 density calculations
    Real(dp), Parameter :: Tc = (Lb(9) * (Zb(9)-Zb(8)) * Tb(9) + Tb(8)**2 - Tb(9)**2) / &
                              & (Lb(9) * (Zb(9)-Zb(8)) + 2._dp * Tb(8) - 2._dp * Tb(9))  !US Standard Atmosphere 1976 equation B-8
    Real(dp), Parameter :: big_A = Tb(8) - Tc  !US Standard Atmosphere 1976 equation B-5
    Real(dp), Parameter :: little_A = (Zb(9)-Zb(8)) * big_A / Sqrt(big_A**2 - (Tb(9)-Tc)**2)  !US Standard Atmosphere 1976 equation B-9
    Real(dp), Parameter :: T_inf = 1000._dp
    Real(dp), Parameter :: lambda = Lb(9) / (T_inf - Tb(10))  !precomputed quantity for 1976 temperature calculations
    Real(dp), Parameter :: R_Z10 = R_Earth + Zb(10)
    Real(dp), Parameter :: Na = 6.022169E26_dp  ![1/kmol] Avagadro's Number
    Real(dp), Parameter :: inv_Na = 1._dp / Na
    Real(dp), Parameter :: N_star = R_star / Na
    Real(dp), Parameter :: K0 = 1.2E2_dp
    Real(dp), Parameter :: Mi(1:3) = (/ 28.0134_dp, & !N2
                                      & 15.9994_dp, & !O1
                                      & 31.9988_dp /) !O2
    Real(dp), Parameter :: ai(2:3) = (/ 6.986E20_dp, & !O1
                                      & 4.863E20_dp /) !O2
    Real(dp), Parameter :: bi = 0.750_dp  !same for O1 and O2
    Real(dp), Parameter :: bigQi(2:3) = (/ -5.809644E-4_dp, & !O1
                                         &  1.366212E-4_dp /) !O2
    Real(dp), Parameter :: bigUi(2:3) = (/ 56.90311_dp, & !O1
                                         & 86._dp /) !O2
    Real(dp), Parameter :: bigWi(2:3) = (/ 2.706240E-5_dp, & !O1
                                         & 8.333333E-5_dp /) !O2
    Real(dp), Parameter :: littleQi = -3.416248E-3_dp  !only defined for O1
    Real(dp), Parameter :: littleUi = 97._dp  !only defined for O1
    Real(dp), Parameter :: littleWi = 5.008765E-4_dp  !only defined for O1
    Real(dp), Parameter :: N7(1:5) = (/ 1.129794E20_dp, &  !N2
                                      & 8.6E16_dp, &  !O1
                                      & 3.030898E19_dp, &  !O2
                                      & 1.351400E18_dp, &  !Ar
                                      & 7.5817E14_dp /)  !He  !US Standard Atmosphere 1976 table 9    
    !Precomuted parameters for Romberg Quadrature routines
    Real(dp), Parameter :: Romb1(1:10) = (/ 4._dp, &
                                          & 4._dp**2, &
                                          & 4._dp**3, &
                                          & 4._dp**4, &
                                          & 4._dp**5, &
                                          & 4._dp**6, &
                                          & 4._dp**7, &
                                          & 4._dp**8, &
                                          & 4._dp**9, &
                                          & 4._dp**10 /)
    Real(dp), Parameter :: Romb2(1:10) = 1._dp / (Romb1 - 1._dp)
    
Contains

Function Find_Base_Layer(Z,iZb) Result(b)
    Use Kinds, Only: dp
    Use Utilities, Only: Bisection_Search
    Implicit None
    Integer :: b
    Real(dp), Intent(In) :: Z
    Integer, Intent(In), Optional :: iZb(1:3)
    
    If (Present(iZb)) Then
        b = (iZb(1)-1) + Bisection_Search(Z,Zb(iZb(1):iZb(2)),iZb(3)) - 1  !subtract 1 to get index for layer below Z
    Else
        b = Bisection_Search(Z,Zb(1:10),10) - 1  !subtract 1 to get index for layer below Z
    End If
End Function Find_Base_Layer

Function T(Z,layer,layer_range)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: T  ![K]
    Real(dp), Intent(In) :: Z ![km]
    Integer, Intent(In), Optional :: layer  !layer in which Z falls
    Integer, Intent(In), Optional :: layer_range(1:3)  !layer range over which to search
    Integer :: b  !index of base for layer
    
    !find atmospheric base layer
    If (Present(layer)) Then
        b = layer - 1
    Else If (Present(layer_range)) Then
        b = Find_Base_Layer(Z,layer_range)
    Else
        b = Find_Base_Layer(Z)
    End If
    If (Lb_nonzero(b)) Then
        If (T_linear_by_H(b)) Then
            !T = Tb(b) + Lb(b) * (Z_to_H(Z) - Hb(b))  !US Standard Atmosphere 1976 equation 23
            T = Tb_minus_LbHb(b) + Lb(b) * Z_to_H(Z)  !US Standard Atmosphere 1976 equation 23
            If (b.EQ.6 .AND. Z.GT.80._dp) T = T * T_M0_correction(Z)  !US Standard Atmosphere 1976 equation 22
        Else
            T = Tb(b) + Lb(b) * (Z - Zb(b))  !US Standard Atmosphere 1976 equation 29
        End If
    Else If (T_exponential(b)) Then
        T = T_inf - (T_inf - Tb(b)) * Exp(-lambda * (Z - Zb(b)) * R_Z10 / (R_Earth + Z))  !US Standard Atmosphere 1976 equation 31
    Else If (T_elliptical(b)) Then
        T = Tc + big_A * Sqrt(1._dp - ((Z - Zb(b)) / little_A)**2)  !US Standard Atmosphere 1976 equation 27
    Else !zero lapse rate
        T = Tb(b)
    End If
End Function T

Function dT_dZ(Z,layer,layer_range)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: dT_dZ  ![K/km]
    Real(dp), Intent(In) :: Z ![km]
    Integer, Intent(In), Optional :: layer
    Integer, Intent(In), Optional :: layer_range(1:3)
    Integer :: b
    
    !find atmospheric base layer
    If (Present(layer)) Then
        b = layer - 1
    Else If (Present(layer_range)) Then
        b = Find_Base_Layer(Z,layer_range)
    Else
        b = Find_Base_Layer(Z)
    End If
    If (Lb_nonzero(b)) Then
        dT_dZ = Lb(b)
    Else If (T_exponential(b)) Then
        dT_dZ = lambda * (T_inf - Tb(b)) * ((R_Earth + Zb(b)) / (R_Earth + Z))**2 * Exp(-lambda * (Z - Zb(b)) * R_Z10 / (R_Earth + Z))  !US Standard Atmosphere 1976 equation 32
    Else If (T_elliptical(b)) Then
        dT_dZ = -big_A * (Z - Zb(b)) / (little_A * Sqrt(1._dp - ((Z - Zb(b)) / little_A)**2))  !US Standard Atmosphere 1976 equation 28
    Else
        dT_dZ = 0._dp
    End If
End Function dT_dZ

Function T_M0_correction(Z) Result(c)
    Use Kinds, Only: dp
    Use Interpolation, Only: Linear_Interp
    Implicit None
    Real(dp) :: c
    Real(dp), Intent(In) :: Z
    Integer :: i
    Real(dp), Parameter :: Zm_corr(0:12) = (/ 80._dp, &
                                            & 80.5_dp, &
                                            & 81._dp, &
                                            & 81.5_dp, &
                                            & 82._dp, &
                                            & 82.5_dp, &
                                            & 83._dp, &
                                            & 83.5_dp, &
                                            & 84._dp, &
                                            & 84.5_dp, &
                                            & 85._dp, &
                                            & 85.5_dp, &
                                            & 86._dp /)
    Real(dp), Parameter :: M0_corr(0:12) = (/ 1._dp, &
                                            & 0.999996_dp, &
                                            & 0.999989_dp, &
                                            & 0.999971_dp, &
                                            & 0.999941_dp, &
                                            & 0.999909_dp, &
                                            & 0.999870_dp, &
                                            & 0.999829_dp, &
                                            & 0.999786_dp, &
                                            & 0.999741_dp, &
                                            & 0.999694_dp, &
                                            & 0.999641_dp, &
                                            & 0.9995788_dp /)
    
    i = Ceiling(2._dp * (Z - 80._dp))
    c = Linear_Interp(Z,Zm_corr(i-1),Zm_corr(i),M0_corr(i-1),M0_corr(i))
End Function T_M0_correction

Function nN2_power(Z,b) Result(x)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: x
    Real(dp), Intent(In) :: Z
    Integer, Intent(In) :: b
    Real(dp) :: M_over_R
    Logical :: Z_below_100
    Real(dp), Parameter :: xb(7:10) = (/ 0._dp, &  !Z = 86km
                                       & 0.8891736933272194_dp, &  !Z = 91km
                                       & 3.9815991755600772_dp, &  !Z = 110km
                                       & 5.0588189719155327_dp /)  !Z = 120km
    Real(dp), Parameter :: xb_100 = 2.4639385720488779_dp  !Z = 100km
    Logical, Parameter :: no_sublayers(7:10) = (/ .TRUE., &
                                                & .FALSE., &
                                                & .TRUE., &
                                                & .TRUE. /)
    Real(dp), Parameter :: rho_star_N2 = Mi(1) / R_star
    
    If (Z .LE. 100._dp) Then
        M_over_R = rho_star
        Z_below_100 = .TRUE.
    Else
        M_over_R = rho_star_N2
        Z_below_100 = .FALSE.
    End If
    If (no_sublayers(b)) Then
        x = xb(b) + M_over_R * Romberg_Quad_nN2(Zb(b),Z,b)
    Else !b=8
        If (Z_below_100) Then
            x = xb(b) + M_over_R * Romberg_Quad_nN2(Zb(b),Z,b)
        Else
            x = xb_100 + M_over_R * Romberg_Quad_nN2(100._dp,Z,b)
        End If
    End If
End Function nN2_power

Function g(Z)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: g
    Real(dp), Intent(In) :: Z
    
    g = g0 * (R_Earth / (R_Earth + Z))**2  !US Standard Atmosphere 1976 equation 17
End Function g

Function nO1_O2_powers(Z,b) Result(x)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: x(1:2)
    Real(dp), Intent(In) :: Z
    Integer, Intent(In) :: b
    Logical :: Z_below_97
    Real(dp), Parameter :: xb(1:2,7:10) = Reshape( (/ 0._dp, &                   !O1, Z = 86km
                                                    & -1.2335160528821659_dp, &  !O1, Z = 91km
                                                    & -1.2350408308007523_dp, &  !O1, Z = 110km
                                                    & -0.7312343085877384_dp, &  !O1, Z = 120km
                                                    &  0._dp, &                  !O2, Z = 86km
                                                    &  0.8987087875276969_dp, &  !O2, Z = 91km
                                                    &  4.5003520516313137_dp, &  !O2, Z = 110km
                                                    &  5.8804674730121536_dp /), &  !O2, Z = 120km
                                                    & (/2,4/) )
    Real(dp), Parameter :: xb_95(1:2) = (/ -1.6326230599173679_dp, &  !O1, Z = 95km
                                         &  1.6401382531252131_dp /)  !O2, Z = 95km
    Real(dp), Parameter :: xb_97(1:2) = (/ -1.6736400073418874_dp, &  !O1, Z = 97km
                                         &  2.0206761115485280_dp /)  !O2, Z = 97km
    Real(dp), Parameter :: xb_100(1:2) = (/ -1.6519509369598046_dp, &  !O1, Z = 100km
                                          &  2.6026364795542958_dp /)  !O2, Z = 100km
    Real(dp), Parameter :: xb_115(1:2) = (/ -0.9801501587157719_dp, &  !O1, Z = 115km
                                          &  5.2767043940711219_dp /)  !O2, Z = 115km
    Logical, Parameter :: no_sublayers(7:10) = (/ .TRUE., &
                                                & .FALSE., &
                                                & .FALSE., &
                                                & .TRUE. /)
    
    If (Z .LE. 97._dp) Then
        Z_below_97 = .TRUE.
    Else
        Z_below_97 = .FALSE.
    End If
    If (no_sublayers(b)) Then
        If (Z_below_97) Then !b=7
            x = Romberg_Quad_nO1_O2(nO1_O2_integrand1,Zb(7),Z,b)
        Else !b=10
            x = xb(:,10) + Romberg_Quad_nO1_O2(nO1_O2_integrand5,Zb(10),Z,b)
        End If
    Else
        If (Z_below_97) Then !b=8
            If (Z .LT. 95._dp) Then
                x = xb(:,8) + Romberg_Quad_nO1_O2(nO1_O2_integrand1,Zb(8),Z,b)
            Else
                x = xb_95 + Romberg_Quad_nO1_O2(nO1_O2_integrand2,95._dp,Z,b)
            End If
        Else !b=8 OR b=9
            If (Z .LT. 100._dp) Then
                x = xb_97 + Romberg_Quad_nO1_O2(nO1_O2_integrand2,97._dp,Z,b)
            Else If (b .EQ. 8) Then
                x = xb_100 + Romberg_Quad_nO1_O2(nO1_O2_integrand3,100._dp,Z,b)
            Else !b=9
                If (Z .LT. 115._dp) Then
                    x = xb(:,9) + Romberg_Quad_nO1_O2(nO1_O2_integrand3,Zb(9),Z,b)
                Else
                    x = xb_115 + Romberg_Quad_nO1_O2(nO1_O2_integrand4,115._dp,Z,b)
                End If
            End If
        End If
    End If
End Function nO1_O2_powers

Function nO1_O2_integrand1(Z,b)  !for 86 to 95 km
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: nO1_O2_integrand1(1:2)
    Real(dp), Intent(In) :: Z
    Integer, Intent(In) :: b
    Real(dp) :: Tz
    Real(dp) :: D(1:2)
    
    Tz = T(Z,b)
    D = ai * (Tz / 273.15_dp)**bi / (N7(1) * Tb(7) * Exp(-nN2_power(Z,b)) / Tz)
    nO1_O2_integrand1 = g(Z) * D * (Mi(2:3) + M0*K0/D) / (R_star * Tz * (D + K0)) + bigQi * (Z - bigUi)**2 * Exp(-bigWi*(Z - bigUi)**3)
    nO1_O2_integrand1(1) = nO1_O2_integrand1(1) + littleQi * (littleUi - Z)**2 * Exp(-littleWi*(littleUi - Z)**3)
End Function nO1_O2_integrand1

Function nO1_O2_integrand2(Z,b)  !for 95 to 97 km
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: nO1_O2_integrand2(1:2)
    Real(dp), Intent(In) :: Z
    Integer, Intent(In) :: b
    Real(dp) :: Tz
    Real(dp) :: D(1:2)
    Real(dp) :: K
    
    Tz = T(Z,b)
    D = ai * (Tz / 273.15_dp)**bi / (N7(1) * Tb(7) * Exp(-nN2_power(Z,b)) / Tz)
    K = K0 * Exp(1._dp - 400._dp / (400._dp - (Z - 95._dp)**2))
    nO1_O2_integrand2 = g(Z) * D * (Mi(2:3) + M0*K/D) / (R_star * Tz * (D + K)) + bigQi * (Z - bigUi)**2 * Exp(-bigWi*(Z - bigUi)**3)
    nO1_O2_integrand2(1) = nO1_O2_integrand2(1) + littleQi * (littleUi - Z)**2 * Exp(-littleWi*(littleUi - Z)**3)
End Function nO1_O2_integrand2

Function nO1_O2_integrand3(Z,b)  !for 97 to 100 km
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: nO1_O2_integrand3(1:2)
    Real(dp), Intent(In) :: Z
    Integer, Intent(In) :: b
    Real(dp) :: Tz
    Real(dp) :: D(1:2)
    Real(dp) :: K

    Tz = T(Z,b)
    D = ai * (Tz / 273.15_dp)**bi / (N7(1) * Tb(7) * Exp(-nN2_power(Z,b)) / Tz)
    K = 1.2E2_dp * Exp(1._dp - 400._dp / (400._dp - (Z - 95._dp)**2))
    nO1_O2_integrand3 = g(Z) * D * (Mi(2:3) + M0*K/D) / (R_star * Tz * (D + K)) + bigQi * (Z - bigUi)**2 * Exp(-bigWi*(Z - bigUi)**3)
End Function nO1_O2_integrand3

Function nO1_O2_integrand4(Z,b)  !for 100 to 115 km
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: nO1_O2_integrand4(1:2)
    Real(dp), Intent(In) :: Z
    Integer, Intent(In) :: b
    Real(dp) :: Tz
    Real(dp) :: D(1:2)
    Real(dp) :: K
    
    Tz = T(Z,b)
    D = ai * (Tz / 273.15_dp)**bi / (N7(1) * Tb(7) * Exp(-nN2_power(Z,b)) / Tz)
    K = K0 * Exp(1._dp - 400._dp / (400._dp - (Z - 95._dp)**2))
    nO1_O2_integrand4 = g(Z) * D * (Mi(2:3) + Mi(1)*K/D) / (R_star * Tz * (D + K)) + bigQi * (Z - bigUi)**2 * Exp(-bigWi*(Z - bigUi)**3)
End Function nO1_O2_integrand4

Function nO1_O2_integrand5(Z,b)  !for 115 to 1000 km
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: nO1_O2_integrand5(1:2)
    Real(dp), Intent(In) :: Z
    Integer, Intent(In) :: b

    nO1_O2_integrand5 = g(Z) * Mi(2:3) / (R_star * T(Z,b)) + bigQi * (Z - bigUi)**2 * Exp(-bigWi*(Z - bigUi)**3)
End Function nO1_O2_integrand5

Subroutine N_density(Z,Tz,b,N,f)
    !returns total number density and fractional composition of atmosphere above 86km geometric altitude
    Use Kinds, Only: dp
    Implicit None
    Real(dp), Intent(In) :: Z
    Real(dp), Intent(In) :: Tz
    Integer, Intent(In) :: b
    Real(dp), Intent(Out) :: N
    Real(dp), Intent(Out), Optional :: f(1:3)
    Real(dp) :: x(1:3)
    
    f = 0._dp
    !N2 power
    x(1) = nN2_power(Z,b)
    !O1 & O2 powers
    x(2:3) = nO1_O2_powers(Z,b)
    If (Present(f)) Then
        !compute number densities of each species
        f = N7(1:3) * Tb(7) * Exp(x) / Tz
        !compute total number density
        N = Sum(f)
        !convert species number densities to fractions
        f = f / N
    Else
        N = Sum(N7(1:3) * Tb(7) * Exp(x)) / Tz
    End If
End Subroutine N_density

Function Romberg_Quad_nN2(a,b,p) Result(q)
    Use Kinds, Only: dp
    Implicit None
    Real(dp):: q    !the result of the integration
    Real(dp), Intent(In) :: a,b    !limits of integration
    Integer, Intent(In) :: p
    Real(dp) :: R(0:10,0:10)  !Romberg table
    Integer :: n,i,j
    Real(dp) :: h,s,ai
    Real(dp), Parameter :: rTol = 1.E-9_dp

    n = 1
    h = b - a
    s = 0.5_dp * (g(a) / T(a,p) + g(b) / T(b,p))
    R(0,0) = h * s
    Do i = 1,10
        !compute trapezoid estimate for next row of table
        n = n * 2
        h = (b - a) / Real(n,dp)
        Do j = 1,n-1,2  !only odd values of j, these are the NEW points at which to evaluate f
            ai = a + Real(j,dp)*h
            s = s + g(ai) / T(ai,p)
        End Do
        R(0,i) = h * s
        !fill out Romberg table row
        Do j = 1,i
            R(j,i) = Romb2(j) * (Romb1(j) * R(j-1,i) - R(j-1,i-1))
        End Do
        !check for convergence
        If ( Abs(R(i-1,i-1) - R(i,i)) .LE. rTol * Abs(R(i,i)) ) Then
            q = R(i,i)  !R(i,i) is the position of the highest precision converged value
            Return  !Normal exit
        End If
    End Do
    !If we get this far, we did not converge
    Call Continue_Romberg_nN2(a,b,p,s,10,R(:,10),2,q)
End Function Romberg_Quad_nN2

Recursive Subroutine Continue_Romberg_nN2(a,b,p,s,d,R0,level,q)  !adds 10 more rows to the previous Romberg_Quad table
    Use Kinds, Only: dp
    Implicit None
    Real(dp), Intent(In) :: a,b    !limits of integration
    Integer, Intent(In) :: p
    Real(dp), Intent(InOut) :: s  !previous sum of ordinates
    Integer, Intent(In) :: d  !length of final row in OLD Romberg Table
    Real(dp), Intent(In) :: R0(0:d)  !final row of OLD romberg table
    Integer, Intent(In) :: level
    Real(dp), Intent(Out) :: q    !the result of the integration, if convergence attained
    Real(dp) :: R(0:d+10,0:10)  !Romberg table extension
    Integer :: n,i,j
    Real(dp) :: h,ai
    Integer :: fours
    Real(dp), Parameter :: rTol = 1.E-9_dp
    
    R(0:d,0) = R0
    Do i = 1,10
        !compute trapezoid estimate for next row of table
        n = 2**(d+i)
        h = (b - a) / Real(n,dp)
        Do j = 1,n-1,2  !only odd values of j, these are the NEW points at which to evaluate f
            ai = a + Real(j,dp)*h
            s = s + g(ai) / T(ai,p)
        End Do
        R(0,i) = h * s
        !fill out Romberg table row
        fours = 1
        Do j = 1,d+i
            fours = fours * 4
            R(j,i) = (Real(fours,dp) * R(j-1,i) - R(j-1,i-1)) / Real(fours - 1,dp)
            !R(j,i) = (((4._dp)**j) * R(j-1,i) - R(j-1,i-1)) / (((4._dp)**j) - 1._dp)
        End Do
        !check for convergence
        If ( Abs(R(i-1,i-1) - R(i,i)) .LE. rTol * Abs(R(i,i)) ) Then
            q = R(d+i,i)
            Return  !Normal exit
        End If
    End Do
    If (level .GT. 10) Then !max allowed recursion depth, interval has been split 100 times...
        Print *,"ERROR:  US_Std_Atm_1976: Continue_Romberg_nN2:  Failed to converge before reaching max recursion depth."
        ERROR STOP
    End If
    !If we get this far, we did not converge, recurse to add 10 more rows
    Call Continue_Romberg_nN2(a,b,p,s,d+10,R(:,10),level+1,q)
End Subroutine Continue_Romberg_nN2

Function Romberg_Quad_nO1_O2(f,a,b,p) Result(q)
    Use Kinds, Only: dp
    Implicit None
    Real(dp):: q(1:2)    !the result of the integration
    Interface
        Function f(x,i)    !the function to be integrated
            Use Kinds,Only: dp
            Implicit None
            Real(dp) :: f(1:2)
            Real(dp), Intent(In) :: x
            Integer, Intent(In) :: i
        End Function f
    End Interface
    Real(dp), Intent(In) :: a,b    !limits of integration
    Integer, Intent(In) :: p
    Real(dp) :: R(1:2,0:10,0:10)  !Romberg table
    Integer :: n,i,j
    Real(dp) :: h,s(1:2)
    Real(dp), Parameter :: rTol = 1.E-6_dp

    n = 1
    h = b - a
    s = 0.5_dp * (f(a,p) + f(b,p))
    R(:,0,0) = h * s
    Do i = 1,10
        !compute trapezoid estimate for next row of table
        n = n * 2
        h = (b - a) / Real(n,dp)
        Do j = 1,n-1,2  !only odd values of j, these are the NEW points at which to evaluate f
            s = s + f(a + Real(j,dp)*h,p)
        End Do
        R(:,0,i) = h * s
        !fill out Romberg table row
        Do j = 1,i
            R(:,j,i) = Romb2(j) * (Romb1(j) * R(:,j-1,i) - R(:,j-1,i-1))
        End Do
        !check for convergence
        If ( All( Abs(R(:,i-1,i-1) - R(:,i,i)) .LE. rTol * Abs(R(:,i,i)) ) ) Then
            q = R(:,i,i)  !R(i,i) is the position of the highest precision converged value
            Return  !Normal exit
        End If
    End Do
    !If we get this far, we did not converge
    Call Continue_Romberg_nO1_O2(f,a,b,p,s,10,R(:,:,10),2,q)
End Function Romberg_Quad_nO1_O2

Recursive Subroutine Continue_Romberg_nO1_O2(f,a,b,p,s,d,R0,level,q)  !adds 10 more rows to the previous Romberg_Quad table
    Use Kinds, Only: dp
    Implicit None
    Interface
        Function f(x,i)    !the function to be integrated
            Use Kinds,Only: dp
            Implicit None
            Real(dp) :: f(1:2)
            Real(dp), Intent(In) :: x
            Integer, Intent(In) :: i
        End Function f
    End Interface
    Real(dp), Intent(In) :: a,b    !limits of integration
    Integer, Intent(In) :: p
    Real(dp), Intent(InOut) :: s(1:2)  !previous sum of ordinates
    Integer, Intent(In) :: d  !length of final row in OLD Romberg Table
    Real(dp), Intent(In) :: R0(1:2,0:d)  !final row of OLD romberg table
    Integer, Intent(In) :: level
    Real(dp), Intent(Out) :: q(1:2)    !the result of the integration, if convergence attained
    Real(dp) :: R(1:2,0:d+10,0:10)  !Romberg table extension
    Integer :: n,i,j
    Real(dp) :: h
    Integer :: fours
    Real(dp), Parameter :: rTol = 1.E-6_dp
    
    R(:,0:d,0) = R0
    Do i = 1,10
        !compute trapezoid estimate for next row of table
        n = 2**(d+i)
        h = (b - a) / Real(n,dp)
        Do j = 1,n-1,2  !only odd values of j, these are the NEW points at which to evaluate f
            s = s + f(a + Real(j,dp)*h,p)
        End Do
        R(:,0,i) = h * s
        !fill out Romberg table row
        fours = 1
        Do j = 1,d+i
            fours = fours * 4
            R(:,j,i) = (Real(fours,dp) * R(:,j-1,i) - R(:,j-1,i-1)) / Real(fours - 1,dp)
            !R(:,j,i) = (((4._dp)**j) * R(:,j-1,i) - R(:,j-1,i-1)) / (((4._dp)**j) - 1._dp)
        End Do
        !check for convergence
        If ( All( Abs(R(:,i-1,i-1) - R(:,i,i)) .LE. rTol * Abs(R(:,i,i)) ) ) Then
            q = R(:,d+i,i)
            Return  !Normal exit
        End If
    End Do
    If (level .GT. 10) Then !max allowed recursion depth, interval has been split 100 times...
        Print *,"ERROR:  US_Std_Atm_1976: Continue_Romberg_nN2:  Failed to converge before reaching max recursion depth."
        ERROR STOP
    End If
    !If we get this far, we did not converge, recurse to add 10 more rows
    Call Continue_Romberg_nO1_O2(f,a,b,p,s,d+10,R(:,:,10),level+1,q)
End Subroutine Continue_Romberg_nO1_O2

Function P(Z,layer,layer_range)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: P  ![Pa]
    Real(dp), Intent(In) :: Z ![km]
    Integer, Intent(In), Optional :: layer
    Integer, Intent(In), Optional :: layer_range(1:3)
    Real(dp) :: T,N
    Integer :: b
    
    !find atmospheric base layer
    If (Present(layer)) Then
        b = layer - 1
    Else If (Present(layer_range)) Then
        b = Find_Base_Layer(Z,layer_range)
    Else
        b = Find_Base_Layer(Z)
    End If
    If (Lb_nonzero(b)) Then
        If (T_linear_by_H(b)) Then
            !T = Tb(b) + Lb(b) * (Z_to_H(Z) - Hb(b))  !US Standard Atmosphere 1976 equation 23
            T = Tb_minus_LbHb(b) + Lb(b) * Z_to_H(Z)  !US Standard Atmosphere 1976 equation 23
            If (b.EQ.6 .AND. Z.GT.80._dp) T = T * T_M0_correction(Z)  !US Standard Atmosphere 1976 equation 22
        Else
            T = Tb(b) + Lb(b) * (Z - Zb(b))  !US Standard Atmosphere 1976 equation 29
        End If
        If (P_rho_not_by_N(b)) Then
            P = Pb_Tb_L_star_Lb(b) * T**(-L_star_Lb(b))  !US Standard Atmosphere 1976 equation 33a
        Else
            Call N_density(Z,T,b,N)
            P = N * T * N_star  !US Standard Atmosphere 1976 equation 33c
        End If
    Else If (T_exponential(b)) Then
        T = T_inf - (T_inf - Tb(b)) * Exp(-lambda * (Z - Zb(b)) * R_Z10 / (R_Earth + Z))  !US Standard Atmosphere 1976 equation 31
        Call N_density(Z,T,b,N)
        P = N * T * N_star  !US Standard Atmosphere 1976 equation 33c
    Else If (T_elliptical(b)) Then
        T = Tc + big_A * Sqrt(1._dp - ((Z - Zb(b)) / little_A)**2)  !US Standard Atmosphere 1976 equation 27
        Call N_density(Z,T,b,N)
        P = N * T * N_star  !US Standard Atmosphere 1976 equation 33c
    Else !zero lapse rate
        If (P_rho_not_by_N(b)) Then
            P = Pb(b) * Exp( L_star_Tb(b) * (Z_to_H(Z) - Hb(b)) )  !US Standard Atmosphere 1976 equation 33b
        Else
            Call N_density(Z,T,b,N)
            P = N * Tb(b) * N_star  !US Standard Atmosphere 1976 equation 33c
        End If
    End If
End Function P

Function rho(Z,layer,layer_range)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: rho  ![g/m^3]
    Real(dp), Intent(In) :: Z ![km]
    Integer, Intent(In), Optional :: layer
    Integer, Intent(In), Optional :: layer_range(1:3)
    Real(dp) :: T,P,N,f(1:3)
    Integer :: b
    
    !find atmospheric base layer
    If (Present(layer)) Then
        b = layer - 1
    Else If (Present(layer_range)) Then
        b = Find_Base_Layer(Z,layer_range)
    Else
        b = Find_Base_Layer(Z)
    End If
    If (Lb_nonzero(b)) Then
        If (T_linear_by_H(b)) Then
            !T = Tb(b) + Lb(b) * (Z_to_H(Z) - Hb(b))  !US Standard Atmosphere 1976 equation 23
            T = Tb_minus_LbHb(b) + Lb(b) * Z_to_H(Z)  !US Standard Atmosphere 1976 equation 23
            If (b.EQ.6 .AND. Z.GT.80._dp) T = T * T_M0_correction(Z)  !US Standard Atmosphere 1976 equation 22
        Else
            T = Tb(b) + Lb(b) * (Z - Zb(b))  !US Standard Atmosphere 1976 equation 29
        End If
        If (P_rho_not_by_N(b)) Then
            P = Pb_Tb_L_star_Lb(b) * T**(-L_star_Lb(b))  !US Standard Atmosphere 1976 equation 33a
            rho = P * rho_star /  T  !US Standard Atmosphere 1976 equation 42-1
        Else
            Call N_density(Z,T,b,N,f)
            rho = Sum(N * f * Mi(1:3)) * inv_Na  !US Standard Atmosphere 1976 equation 42-3
        End If
    Else If (T_exponential(b)) Then
        T = T_inf - (T_inf - Tb(b)) * Exp(-lambda * (Z - Zb(b)) * R_Z10 / (R_Earth + Z))  !US Standard Atmosphere 1976 equation 31
        Call N_density(Z,T,b,N,f)
        rho = Sum(N * f * Mi(1:3)) * inv_Na  !US Standard Atmosphere 1976 equation 42-3
    Else If (T_elliptical(b)) Then
        T = Tc + big_A * Sqrt(1._dp - ((Z - Zb(b)) / little_A)**2)  !US Standard Atmosphere 1976 equation 27
        Call N_density(Z,T,b,N,f)
        rho = Sum(N * f * Mi(1:3)) * inv_Na  !US Standard Atmosphere 1976 equation 42-3
    Else !zero lapse rate
        T = Tb(b)
        If (P_rho_not_by_N(b)) Then
            P = Pb(b) * Exp( L_star_Tb(b) * (Z_to_H(Z) - Hb(b)) )  !US Standard Atmosphere 1976 equation 33b
            rho = P * rho_star /  T  !US Standard Atmosphere 1976 equation 42-1
        Else
            Call N_density(Z,T,b,N,f)
            rho = Sum(N * f * Mi(1:3)) * inv_Na  !US Standard Atmosphere 1976 equation 42-3
        End If
    End If
End Function rho

Elemental Function Z_to_H(Z) Result(H)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: H
    Real(dp), Intent(In) :: Z
    
    H = Z * R_Earth / (Z + R_Earth)
End Function Z_to_H

Elemental Function H_to_Z(H) Result(Z)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: Z
    Real(dp), Intent(In) :: H
    
    Z = H * R_Earth / (R_Earth - H)
End Function H_to_Z

End Module US_Std_Atm_1976
