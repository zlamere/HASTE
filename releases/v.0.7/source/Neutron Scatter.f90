Module Neutron_Scatter
    
    Use Kinds, Only: dp
    Use Cross_Sections, Only: CS_Type
    Implicit None
    
    Private
    Public :: Neutron_Type
    Public :: Scatter_Data_Type
    Public :: Scatter_Model_Type
    Public :: Setup_Scatter_Model
    Public :: Write_Scatter_Model
    Public :: Air_Velocity
    Public :: Apparent_Energy
    Public :: Scattered_Direction
    Public :: Scattered_Angles

    Type Neutron_Type
        Real(dp) :: r(1:3)  ![km] cartesian coordinates (x,y,z) of neutron
        Real(dp) :: big_r  ![km] magnitude of r
        Real(dp) :: Z  ![km]  geometric altitude
        Real(dp) :: zeta  !cosine of the zenith angle
        Real(dp) :: Omega_hat(1:3)  ![1 km/s] unit vector in direction of neutron travel (unit velocity vector)
        Real(dp) :: E0ef  ![keV]energy at emission in emission frame
        Real(dp) :: s0ef  ![km/s] speed in emission frame at emission
        Real(dp) :: E  ![keV] neutron energy
        Real(dp) :: weight  !weight (probability of neutron existence in this state), reduced for things like absorption supression
        Real(dp) :: t  ![s] simulation time since neutron emission
    End Type
    
    Type Scatter_Data_Type
        Real(dp) :: sig_A  ![km^-1]  absorption and total macroscopic (seal-level) cross sections
        Real(dp) :: sig_T
        Integer :: target_index  !index of the target isotope
        Real(dp) :: An  !mass of the target isotope in neutron masses
        Real(dp) :: vAir(1:3)  !intermediate quantity for computing vA, saved as a side effect in case all materials/mechanims are computed
        Real(dp) :: vA(1:3)  ![km/s] velocity vector of the target nucleus
        Integer :: level  !sampled level of scatter
        Real(dp) :: Q  !Q-value for chosen scatter level
        Logical :: da_is_Legendre
        Integer :: n_a  !number of scatter coeffiecents
        Real(dp), Allocatable :: a(:)  !anisotropic scatter coeffs
        Integer :: n_a1,n_a2  !number of tabulated cosine pdf pairs
        Real(dp) :: a_tab_Econv
        Real(dp), Allocatable :: a_tab1(:,:),a_tab2(:,:)  !tabulated scatter cosine pdf values at energies above and below the incident energy
        Real(dp), Allocatable :: iso_cs(:)  !relative scatter cross section of each isotope at last computed energy, saved as a side effect in case all materials/mechanims are computed
        Real(dp), Allocatable :: lev_cs(:)  !relative scatter cross section of each level at last computed energy, saved as a side effect in case all materials/mechanims are computed
        Real(dp) :: Omega_hat0_cm(1:3)
        Real(dp) :: u(1:3),u_speed
        Real(dp) :: s1cm  !speed after the scatter in CM frame is computed as a side effect of the next-event procedure, saving it reduces operations in following scatter simulation
        Real(dp) :: A_hat(1:3),B_hat(1:3),C_hat(1:3)  !cartesian basis in CM frame is computed as a side effect of the next-event procedure, saving it reduces operations in following scatter simulation
    End Type
    
    Type Scatter_Model_Type
        Integer :: n_scatters  !number of scatters to follow for history, scatter is forced to the detector at the final scatter
        Logical :: direct_contribution
        Logical :: estimate_each_scatter
        Logical :: elastic_only
        Logical :: aniso_dist
        Logical :: suppress_absorption
        Logical :: suppress_leakage
        Logical :: all_mat_mech
        Logical :: fast_analog
        Logical :: roulette
        Real(dp) :: roulette_weight
        Real(dp) :: roulette_rate
        Real(dp) :: roulette_mult
        Integer(8) :: n_kills(1:6)  !kills due to (1)Time, (2)Energy, (3)Weight, (4)Leakage, (5)Absorption. (6)Collision limit
        Integer(8) :: next_events(1:3)  !next-events at detector (1)Attempted, (2)Found, (3)Tallied
        Integer(8) :: n_no_tally(1:3)  !next-events at detector NOT tallied due to (1)Time & Energy, (2)Time, (3)Energy
        Type(CS_Type) :: CS
        Type(Scatter_Data_Type) :: scat  !parameters defining the next scatter
        Logical :: Gravity
        Logical :: Neutron_Decay
        !N2H The logical 'Target_Motion' is currently unused except for informational purposes, consider removing
        Logical :: Target_Motion
            Logical :: Doppler_Broaden
            Logical :: Thermal_Motion
            Logical :: Diatomic_Atm
            Logical :: Rotating_Earth
            Logical :: Wind
    Contains
        Procedure, Pass :: Sample_Scatter
        Procedure, Pass :: Set_Scatter_prep
        Procedure, Pass :: Set_Scatter_iso
        Procedure, Pass :: Set_Scatter_lev
    End Type
    
Contains

Function Setup_Scatter_Model(setup_file_name,resources_directory,cs_setup_file,run_file_name) Result(ScatMod)
    Use Cross_Sections, Only: Setup_Cross_Sections
    Use Global, Only: n_kill_weight
    Use FileIO_Utilities, Only: fSHARE
    Implicit None
    Type(Scatter_Model_Type) :: ScatMod
    Character(*), Intent(In) :: setup_file_name
    Character(*), Intent(In) :: resources_directory
    Character(*), Intent(In) :: cs_setup_file
    Character(*), Intent(In) :: run_file_name
    Integer :: n_scatters
    Logical :: direct_contribution,estimate_each_scatter,elastic_only
    Logical :: suppress_absorption,suppress_leakage,all_mat_mech
    Logical :: Gravity,neutron_decay,doppler_broaden
    Logical :: thermal_motion,diatomic_Atm,Rotating_Earth,wind
    Integer :: setup_unit,stat
    Character(10) :: scatter_model  !IsoCM, AnIsoCM
    Character(10) :: dist_to_next_event  !An-exact,An-fast
    Logical :: roulette
    Real(dp) :: roulette_weight
    Integer :: roulette_ratio
    Real(dp) :: E_min,E_max
    
    NameList /NeutronScatterList/ n_scatters,scatter_model,elastic_only,suppress_absorption, &
                                  & suppress_leakage,all_mat_mech,dist_to_next_event, &
                                  & roulette,roulette_weight,roulette_ratio, &
                                  & direct_contribution,estimate_each_scatter,E_min,E_max, &
                                  & Gravity,neutron_decay,doppler_broaden,thermal_motion, &
                                  & Diatomic_Atm,Rotating_Earth,Wind
    
    Open(NEWUNIT = setup_unit , FILE = setup_file_name , STATUS = 'OLD' , ACTION = 'READ' , IOSTAT = stat , SHARE = fSHARE)
    If (stat .NE. 0) Then
        Print *,'ERROR:  Neutron_Scatter: Setup_Scatter_Model:  File open error, '//setup_file_name//', IOSTAT=',stat
        ERROR STOP
    End If
    Read(setup_unit,NML = NeutronScatterList)
    Close(setup_unit)
    If (n_scatters.EQ.-1 .AND. .NOT.estimate_each_scatter) Then
        Print *,'ERROR:  Neutron_Scatter: Setup_Scatter_Model:  estimate_each_scatter must be .TRUE. for n_scatters = -1 (unlimited)'
        ERROR STOP
    End If
    If (n_scatters.EQ.0 .AND. .NOT.direct_contribution) Then
        Print *,'ERROR:  Neutron_Scatter: Setup_Scatter_Model:  direct_contribution must be .TRUE. for n_scatters = 0'
        ERROR STOP
    End If
    Select Case (scatter_model)
        Case('IsoCM')
            ScatMod%aniso_dist = .FALSE.
        Case('AnIsoCM')
            ScatMod%aniso_dist = .TRUE.
        Case Default
            Print *,'ERROR:  Neutron_Scatter: Setup_Scatter_Model:  Undefined scatter model'
            ERROR STOP
    End Select
    ScatMod%n_scatters = n_scatters
    If (direct_contribution) Then
        If (n_scatters.EQ.0 .OR. estimate_each_scatter) Then
            ScatMod%direct_contribution = direct_contribution
        Else
            Print *,'ERROR:  Neutron_Scatter: Setup_Scatter_Model:  direct_contribution cannot be TRUE w/ n_scatters.NE.0 or estimate_each_scatter=FALSE '
            ERROR STOP
        End If
    End If
    ScatMod%estimate_each_scatter = estimate_each_scatter
    ScatMod%elastic_only = elastic_only
    ScatMod%suppress_absorption = suppress_absorption
    ScatMod%suppress_leakage = suppress_leakage
    ScatMod%all_mat_mech = all_mat_mech
    Select Case (dist_to_next_event)
        Case('An-exact')
            ScatMod%fast_analog = .FALSE.
        Case('An-fast')
            !UNDONE Fast analog option not yet implemented
            !ScatMod%fast_analog = .TRUE.
            Print *,'ERROR:  Neutron_Scatter: Setup_Scatter_Model:  Fast Analog Monte-Carlo game not yet implemented'
            ERROR STOP
        Case Default
            Print *,'ERROR:  Neutron_Scatter: Setup_Scatter_Model:  Undefined Monte-Carlo game'
            ERROR STOP
    End Select
    ScatMod%roulette = roulette
    If (roulette) Then
        ScatMod%roulette_weight = roulette_weight * n_kill_weight
        ScatMod%roulette_rate = 1._dp / Real(roulette_ratio,dp)
        ScatMod%roulette_mult = Real(roulette_ratio,dp)
    Else
        ScatMod%roulette_weight = 0._dp
        ScatMod%roulette_rate = 1._dp
        ScatMod%roulette_mult = 1._dp
    End If
    ScatMod%Gravity = Gravity
    ScatMod%Neutron_Decay = Neutron_Decay
    ScatMod%Doppler_Broaden = Doppler_Broaden
    If (Any( (/Thermal_Motion,Rotating_Earth,Wind/) )) Then
        ScatMod%Target_Motion = .TRUE.
    Else
        ScatMod%Target_Motion = .FALSE.
    End If
    ScatMod%Thermal_Motion = Thermal_Motion
    ScatMod%Rotating_Earth = Rotating_Earth
    ScatMod%Wind = Wind
    ScatMod%n_kills = 0
    ScatMod%next_events = 0
    ScatMod%n_no_tally = 0
    ScatMod%CS = Setup_Cross_Sections(resources_directory,cs_setup_file,elastic_only,ScatMod%aniso_dist,E_min,E_max)
    !initialize scatter parameters for sampled scatter, except the lev_cs array (it is not used for the sampled scatter)
    Allocate(ScatMod%scat%a(0:ScatMod%CS%n_a_max))
    ScatMod%scat%a = 0._dp
    Allocate(ScatMod%scat%a_tab1(1:ScatMod%CS%n_a_tab_max,1:2))
    ScatMod%scat%a_tab1 = 0._dp
    Allocate(ScatMod%scat%a_tab2(1:ScatMod%CS%n_a_tab_max,1:2))
    ScatMod%scat%a_tab2 = 0._dp
    Allocate(ScatMod%scat%iso_cs(1:ScatMod%CS%n_iso))
    ScatMod%scat%iso_cs = 0._dp
    If (this_image() .EQ. 1) Then
        Open(NEWUNIT = setup_unit , FILE = run_file_name , STATUS = 'OLD' , ACTION = 'WRITE' , POSITION = 'APPEND' , IOSTAT = stat)
        If (stat .NE. 0) Then
            Print *,'ERROR:  Neutron_Scatter: Setup_Scatter_Model:  File open error, '//run_file_name//', IOSTAT=',stat
            ERROR STOP
        End If
        Write(setup_unit,NML = NeutronScatterList)
        Write(setup_unit,*)
        Close(setup_unit)
    End If
End Function Setup_Scatter_Model

Subroutine Sample_Scatter(ScatMod,n,atm,RNG)
    Use Kinds, Only: dp
    Use Global, Only: mfp_per_barn_per_km_at_seaLevel
    Use Global, Only: Z_hat,X_hat
    Use Atmospheres, Only: Atmosphere_Type
    Use Random_Numbers, Only: RNG_Type
    Use Cross_Sections, Only: sig_composite
    Use Utilities, Only: Bisection_Search
    Use Utilities, Only: Vector_Length
    Use Utilities, Only: Unit_Vector
    Use Utilities, Only: Cross_Product
    Use Interpolation, Only: Linear_Interp
    Use Target_Motion, Only: Random_Thermal_Velocity_AF
    Use Neutron_Utilities, Only: Neutron_Speed
    Use Neutron_Utilities, Only: Neutron_Energy
    Implicit None
    Class(Scatter_Model_Type), Intent(InOut) :: ScatMod
    Type(Neutron_Type), Intent(In):: n
    Type(Atmosphere_Type), Intent(In) :: atm
    Type(RNG_Type), Intent(InOut) :: RNG
    Integer :: E_index,i,index1,index2
    Real(dp) :: r
    Real(dp) :: E_cm
    Real(dp) :: E1,E2
    Real(dp), Allocatable :: a1(:),a2(:),level_cs(:)
    Integer :: iso,lev
    Integer :: target_el
    Real(dp) :: An_prime,Mn
    Real(dp) :: E_apparent
    Real(dp) :: T
    Real(dp) :: v0(1:3),v0cm(1:3)
    
    !Get atmosphere total/absorption and isotope scatter cross sections at the point of scatter
    If (ScatMod%Rotating_Earth .OR. ScatMod%Wind) Then
        ScatMod%scat%vAir = Air_Velocity(ScatMod%Rotating_Earth,ScatMod%Wind,n,atm)
        E_apparent = Apparent_Energy(n%E,n%Omega_hat,ScatMod%scat%vAir)
    Else
        ScatMod%scat%vAir = 0._dp
        E_apparent = n%E
    End If
    If (ScatMod%Doppler_Broaden) Then
        T = atm%T(n%Z)
        Call ScatMod%CS%sig_T_A_broad(E_apparent,T,ScatMod%scat%sig_T,ScatMod%scat%sig_A)
        Do i = 1,ScatMod%CS%n_iso
            ScatMod%scat%iso_cs(i) = ScatMod%CS%sig_S_iso_broad(i,E_apparent,T)
        End Do
    Else
        Call ScatMod%CS%sig_T_A(E_apparent,ScatMod%scat%sig_T,ScatMod%scat%sig_A,iE_get=E_index)
        Do i = 1,ScatMod%CS%n_iso
            ScatMod%scat%iso_cs(i) = ScatMod%CS%sig_S_iso(i,E_apparent,iE_put=E_index)
        End Do
    End If
    ScatMod%scat%sig_A = ScatMod%scat%sig_A * mfp_per_barn_per_km_at_seaLevel  !convert to macroscopic cross section at sea-level
    ScatMod%scat%sig_T = ScatMod%scat%sig_T * mfp_per_barn_per_km_at_seaLevel  !convert to macroscopic cross section at sea-level
    !Choose a target isotope based on relative SCATTER cross sections
    ScatMod%scat%iso_cs = ScatMod%scat%iso_cs * ScatMod%CS%iso_Fractions
    r = RNG%Get_Random() * Sum(ScatMod%scat%iso_cs)
    Do i = 1,ScatMod%CS%n_iso
        If (r .LT. Sum(ScatMod%scat%iso_cs(1:i))) Then
            ScatMod%scat%target_index = i
            ScatMod%scat%An = ScatMod%CS%An(i)
            Exit
        End If
    End Do
    !Add thermal motion
    If (ScatMod%Thermal_Motion) Then
        If (ScatMod%Diatomic_Atm) Then
            If (atm%diatomic(ScatMod%scat%target_index)) Then !Sample isotope of other nucleus in diatomic molecule
                r = RNG%Get_Random()
                target_el = atm%iso_map(ScatMod%scat%target_index)
                Do i = atm%iso_ind(target_el),atm%iso_ind(target_el+1)-1
                    If (r .LT. Sum(atm%iso_frac(atm%iso_ind(target_el):i))) Then
                        An_prime = ScatMod%CS%An(i)
                        Mn = ScatMod%scat%An + An_prime
                        Exit
                    End If
                End Do
                ScatMod%scat%vA = ScatMod%scat%vAir + Random_Thermal_Velocity_AF(atm%T(n%Z),ScatMod%scat%An,An_prime,Mn,RNG)
            Else !monatomic atmospheric constituent
                ScatMod%scat%vA = ScatMod%scat%vAir + Random_Thermal_Velocity_AF(atm%T(n%Z),ScatMod%scat%An,RNG)
            End If
        Else !monotomic approximation of atmosphere
            ScatMod%scat%vA = ScatMod%scat%vAir + Random_Thermal_Velocity_AF(atm%T(n%Z),ScatMod%scat%An,RNG)
        End If
    Else
        ScatMod%scat%vA = ScatMod%scat%vAir
    End If
    !convert energy to the center of mass of the collision frame, storing side effects for later use as well
    v0 = Neutron_Speed(n%E) * n%Omega_hat
    ScatMod%scat%u = (v0 + ScatMod%scat%An * ScatMod%scat%vA) / (ScatMod%scat%An + 1._dp)
    v0cm = v0 - ScatMod%scat%u
    E_cm = Neutron_Energy(v0cm)
    !find E_cm in the unified energy grid
    E_index = Bisection_Search(E_cm,ScatMod%CS%E_uni,ScatMod%CS%n_E_uni)
    !compute a few more predetermined quantities that are known before the scatter
    ScatMod%scat%u_speed = Vector_Length(ScatMod%scat%u)
    ScatMod%scat%Omega_hat0_cm = Unit_Vector(v0cm)
    ScatMod%scat%A_hat = ScatMod%scat%Omega_hat0_cm
    If (Abs(Dot_Product(ScatMod%scat%A_hat,Z_Hat)) .LT. 0.5_dp) Then
        ScatMod%scat%B_hat = Unit_Vector(Cross_Product(Z_Hat,ScatMod%scat%A_hat))
    Else
        ScatMod%scat%B_hat = Unit_Vector(Cross_Product(X_Hat,ScatMod%scat%A_hat))
    End If
    ScatMod%scat%C_hat = Cross_Product(ScatMod%scat%A_hat,ScatMod%scat%B_hat)
    !choose scatter level
    iso = ScatMod%scat%target_index
    If (ScatMod%elastic_only) Then
        ScatMod%scat%level = 0
        ScatMod%scat%Q = 0._dp
    Else  !sample scatter level (including elastic as level 0)
        If (E_index .LT. ScatMod%CS%lev_cs(iso)%thresh(1)) Then !not enough energy for lowest inelastic level, scatter is elastic
            ScatMod%scat%level = 0
            ScatMod%scat%Q = 0._dp
        Else !sample for scatter level
            Allocate(level_cs(0:ScatMod%CS%lev_cs(iso)%n_lev))
            level_cs = 0._dp
            Do i = 0,ScatMod%CS%lev_cs(iso)%n_lev
                If (E_index .LE. ScatMod%CS%lev_cs(iso)%thresh(i)) Exit  !insufficent energy for this or any higher inelastic level
                level_cs(i) = sig_Composite(E_cm,ScatMod%CS%n_E_uni,ScatMod%CS%E_uni,ScatMod%CS%lnE_uni,E_index,1,1,ScatMod%CS%lev_cs(iso)%thresh(i),ScatMod%CS%lev_cs(iso)%sig(i))
            End Do
            r = RNG%Get_Random() * Sum(level_cs)
            Do i = 0,ScatMod%CS%lev_cs(iso)%n_Lev
                If (r .LT. Sum(level_cs(0:i))) Then
                    ScatMod%scat%level = i
                    ScatMod%scat%Q = ScatMod%CS%lev_cs(ScatMod%scat%target_index)%Q(i)
                    Exit
                End If
            End Do
        End If
    End If
    lev = ScatMod%scat%level
    If (ScatMod%aniso_dist) Then
        index1 = ScatMod%CS%lev_cs(iso)%da(lev)%E_map(E_index) - 1
        index2 = ScatMod%CS%lev_cs(iso)%da(lev)%E_map(E_index)
        E1 = ScatMod%CS%E_uni( ScatMod%CS%lev_cs(iso)%da(lev)%E_key( index1 ) )
        E2 = ScatMod%CS%E_uni( ScatMod%CS%lev_cs(iso)%da(lev)%E_key( index2 ) )
        If (ScatMod%CS%lev_cs(iso)%da(lev)%da(index1)%is_Legendre) Then  !angular dist is expressed in legendre coeffs
            ScatMod%scat%da_is_legendre = .TRUE.
            ScatMod%scat%n_a = MaxVal(ScatMod%CS%lev_cs(iso)%da(lev)%da(index1:index2)%n_a)
            Allocate(a1(0:ScatMod%scat%n_a))
            Allocate(a2(0:ScatMod%scat%n_a))
            a1 = 0._dp
            a2 = 0._dp
            a1(0:ScatMod%CS%lev_cs(iso)%da(lev)%da(index1)%n_a) = ScatMod%CS%lev_cs(iso)%da(lev)%da(index1)%a
            a2(0:ScatMod%CS%lev_cs(iso)%da(lev)%da(index2)%n_a) = ScatMod%CS%lev_cs(iso)%da(lev)%da(index2)%a
            ScatMod%scat%a(0:ScatMod%scat%n_a) = Linear_Interp(E_cm,E1,E2,a1,a2)
        Else !angular dist is expressed in tabulated cosine pdf values
            ScatMod%scat%da_is_legendre = .FALSE.
            ScatMod%scat%n_a1 = ScatMod%CS%lev_cs(iso)%da(lev)%da(index1)%n_a
            ScatMod%scat%n_a2 = ScatMod%CS%lev_cs(iso)%da(lev)%da(index2)%n_a
            ScatMod%scat%a_tab1(1:ScatMod%scat%n_a1,:) = ScatMod%CS%lev_cs(iso)%da(lev)%da(index1)%ua
            ScatMod%scat%a_tab2(1:ScatMod%scat%n_a2,:) = ScatMod%CS%lev_cs(iso)%da(lev)%da(index2)%ua
            ScatMod%scat%a_tab_Econv = (E_cm - E1) / (E2 - E1)
        End If
    End If
    !speed of neutron after scatter is independent of direction in CM frame
    ScatMod%scat%s1cm = Neutron_Speed( ( E_cm + ScatMod%scat%An * Neutron_Energy(ScatMod%scat%vA - ScatMod%scat%u) - ScatMod%scat%Q) * &
                                       & ScatMod%scat%An / (ScatMod%scat%An + 1._dp) )
End Subroutine Sample_Scatter

Subroutine Set_Scatter_prep(ScatMod,scat)
    Use Kinds, Only: dp
    Implicit None
    Class(Scatter_Model_Type), Intent(In) :: ScatMod
    Type(Scatter_Data_Type), Intent(Out) :: scat
    
    !Get total and absorption apparent cross sections from sampled scatter
    scat%sig_A = ScatMod%scat%sig_A
    scat%sig_T = ScatMod%scat%sig_T
    !Get air velocity from sampled scatter
    scat%vAir = ScatMod%scat%vAir
    !Set up arrays for angular distributions
    Allocate(scat%a(0:ScatMod%CS%n_a_max))
    scat%a = 0._dp
    Allocate(scat%a_tab1(1:ScatMod%CS%n_a_tab_max,1:2))
    Allocate(scat%a_tab2(1:ScatMod%CS%n_a_tab_max,1:2))
    scat%a_tab1 = 0._dp
    scat%a_tab2 = 0._dp
    !Get precomputed isotope cross sections from sampled scatter
    Allocate(scat%iso_cs(1:ScatMod%CS%n_iso))
    scat%iso_cs = ScatMod%scat%iso_cs / Sum(ScatMod%scat%iso_cs)
End Subroutine Set_Scatter_prep

Subroutine Set_Scatter_iso(ScatMod,n,atm,RNG,scat,iso,n_lev,E_cm,i_E_cm)
    Use Kinds, Only: dp
    Use Global, Only: X_hat,Z_hat
    Use Atmospheres, Only: Atmosphere_Type
    Use Random_Numbers, Only: RNG_Type
    Use Cross_Sections, Only: sig_composite
    Use Utilities, Only: Bisection_Search
    Use Utilities, Only: Vector_Length
    Use Utilities, Only: Unit_Vector
    Use Utilities, Only: Cross_Product
    Use Target_Motion, Only: Random_Thermal_Velocity_AF
    Use Neutron_Utilities, Only: Neutron_Speed
    Use Neutron_Utilities, Only: Neutron_Energy
    Implicit None
    Class(Scatter_Model_Type), Intent(In) :: ScatMod
    Type(Neutron_Type), Intent(In):: n
    Type(Atmosphere_Type), Intent(In) :: atm
    Type(RNG_Type), Intent(InOut) :: RNG
    Type(Scatter_Data_Type), Intent(InOut) :: scat
    Integer, Intent(In) :: iso
    Integer, Intent(Out) :: n_lev
    Real(dp), Intent(Out) :: E_cm
    Integer, Intent(Out) :: i_E_cm
    Integer :: i
    Real(dp) :: r
    Integer :: target_el
    Real(dp) :: An_prime,Mn
    Real(dp) :: v0(1:3),v0cm(1:3)
    
    !Set target isotope
    scat%target_index = iso
    scat%An = ScatMod%CS%An(iso)
    !Add thermal motion
    If (ScatMod%Thermal_Motion) Then
        If (ScatMod%scat%target_index .EQ. iso) Then !the sampled scatter already has values we can use here
            scat%vA = ScatMod%scat%vA
        Else
            If (ScatMod%Diatomic_Atm) Then
                If (atm%diatomic(iso)) Then !Sample isotope of other nucleus in diatomic molecule
                    r = RNG%Get_Random()
                    target_el = atm%iso_map(iso)
                    Do i = atm%iso_ind(target_el),atm%iso_ind(target_el+1)-1
                        If (r .LT. Sum(atm%iso_frac(atm%iso_ind(target_el):i))) Then
                            An_prime = ScatMod%CS%An(i)
                            Mn = scat%An + An_prime
                            Exit
                        End If
                    End Do
                    scat%vA = scat%vAir + Random_Thermal_Velocity_AF(atm%T(n%Z),scat%An,An_prime,Mn,RNG)
                Else !monatomic atmospheric constituent
                    scat%vA = scat%vAir + Random_Thermal_Velocity_AF(atm%T(n%Z),scat%An,RNG)
                End If
            Else !monotomic approximation of atmosphere
                scat%vA = scat%vAir + Random_Thermal_Velocity_AF(atm%T(n%Z),scat%An,RNG)
            End If
        End If
    Else
        scat%vA = scat%vAir
    End If
    !convert energy to the center of mass of the collision frame, storing side effects for later use as well
    v0 = Neutron_Speed(n%E) * n%Omega_hat
    scat%u = (v0 + scat%An * scat%vA) / (scat%An + 1._dp)
    v0cm = v0 - scat%u
    E_cm = Neutron_Energy(v0cm)
    !find E_cm in the unified energy grid
    i_E_cm = Bisection_Search(E_cm,ScatMod%CS%E_uni,ScatMod%CS%n_E_uni)
    !compute a few more predetermined quantities that are known before the scatter
    scat%u_speed = Vector_Length(scat%u)
    scat%Omega_hat0_cm = Unit_Vector(v0cm)
    scat%A_hat = scat%Omega_hat0_cm
    If (Abs(Dot_Product(scat%A_hat,Z_Hat)) .LT. 0.5_dp) Then
        scat%B_hat = Unit_Vector(Cross_Product(Z_Hat,scat%A_hat))
    Else
        scat%B_hat = Unit_Vector(Cross_Product(X_Hat,scat%A_hat))
    End If
    scat%C_hat = Cross_Product(scat%A_hat,scat%B_hat)
    !Get level cross sections
    Allocate(scat%lev_cs(0:ScatMod%CS%lev_cs(iso)%n_lev))
    scat%lev_cs = 0._dp
    If (.NOT. ScatMod%elastic_only) Then
        Do i = 0,ScatMod%CS%lev_cs(iso)%n_lev
            If (i_E_cm .LE. ScatMod%CS%lev_cs(iso)%thresh(i)) Exit  !insufficent energy for this or any higher inelastic level
            n_lev = i
            scat%lev_cs(i) = sig_Composite(E_cm,ScatMod%CS%n_E_uni,ScatMod%CS%E_uni,ScatMod%CS%lnE_uni,i_E_cm,1,1,ScatMod%CS%lev_cs(iso)%thresh(i),ScatMod%CS%lev_cs(iso)%sig(i))
        End Do
        If (n_lev .GT. 0) Then
            scat%lev_cs = scat%lev_cs / Sum(scat%lev_cs)
        Else
            scat%lev_cs = 1._dp
        End If
    End If
End Subroutine Set_Scatter_iso

Subroutine Set_Scatter_lev(ScatMod,scat,lev,E_cm,i_E_cm)
    Use Kinds, Only: dp
    Use Interpolation, Only: Linear_Interp
    Use Neutron_Utilities, Only: Neutron_Speed
    Use Neutron_Utilities, Only: Neutron_Energy
    Implicit None
    Class(Scatter_Model_Type), Intent(In) :: ScatMod
    Type(Scatter_Data_Type), Intent(InOut) :: scat
    Integer, Intent(In) :: lev
    Real(dp), Intent(In) :: E_cm
    Integer, Intent(In) :: i_E_cm
    Integer :: index1,index2
    Real(dp) :: E1,E2
    Real(dp), Allocatable :: a1(:),a2(:)
    
    !Set scatter level
    scat%level = lev
    scat%Q = ScatMod%CS%lev_cs(scat%target_index)%Q(lev)
    If (ScatMod%aniso_dist) Then  !get angular distribution coeffs
        index1 = ScatMod%CS%lev_cs(scat%target_index)%da(lev)%E_map(i_E_cm) - 1
        index2 = ScatMod%CS%lev_cs(scat%target_index)%da(lev)%E_map(i_E_cm)
        E1 = ScatMod%CS%E_uni( ScatMod%CS%lev_cs(scat%target_index)%da(lev)%E_key( index1 ) )
        E2 = ScatMod%CS%E_uni( ScatMod%CS%lev_cs(scat%target_index)%da(lev)%E_key( index2 ) )
        If (ScatMod%CS%lev_cs(scat%target_index)%da(lev)%da(index1)%is_Legendre) Then  !angular dist is expressed in legendre coeffs
            scat%da_is_legendre = .TRUE.
            scat%n_a = MaxVal(ScatMod%CS%lev_cs(scat%target_index)%da(lev)%da(index1:index2)%n_a)
            Allocate(a1(0:scat%n_a))
            Allocate(a2(0:scat%n_a))
            a1 = 0._dp
            a2 = 0._dp
            a1(0:ScatMod%CS%lev_cs(scat%target_index)%da(lev)%da(index1)%n_a) = ScatMod%CS%lev_cs(scat%target_index)%da(lev)%da(index1)%a
            a2(0:ScatMod%CS%lev_cs(scat%target_index)%da(lev)%da(index2)%n_a) = ScatMod%CS%lev_cs(scat%target_index)%da(lev)%da(index2)%a
            scat%a(0:scat%n_a) = Linear_Interp(E_cm,E1,E2,a1,a2)
        Else !angular dist is expressed in tabulated cosine pdf values
            scat%da_is_legendre = .FALSE.
            scat%n_a1 = ScatMod%CS%lev_cs(scat%target_index)%da(lev)%da(index1)%n_a
            scat%n_a2 = ScatMod%CS%lev_cs(scat%target_index)%da(lev)%da(index2)%n_a
            scat%a_tab1(1:scat%n_a1,:) = ScatMod%CS%lev_cs(scat%target_index)%da(lev)%da(index1)%ua
            scat%a_tab2(1:scat%n_a2,:) = ScatMod%CS%lev_cs(scat%target_index)%da(lev)%da(index2)%ua
            scat%a_tab_Econv = (E_cm - E1) / (E2 - E1)
        End If
    End If
    !speed of neutron after scatter is independent of direction in CM frame
    scat%s1cm = Neutron_Speed( ( E_cm + scat%An * Neutron_Energy(scat%vA - scat%u) - scat%Q) * &
                               & scat%An / (scat%An + 1._dp) )
End Subroutine Set_Scatter_lev

Function Air_Velocity(Rotating_Earth,Wind,n,atm) Result(vA)
    Use Kinds, Only: dp
    Use Global, Only: Z_hat
    Use Utilities, Only: Unit_Vector
    Use Utilities, Only: Cross_Product
    Use Atmospheres, Only: Atmosphere_Type
    Use Target_Motion, Only: Atm_Rotation_Velocity_AF
    Implicit None
    Real(dp) :: vA(1:3)
    Logical, Intent(In) :: Rotating_Earth,Wind
    Type(Neutron_Type), Intent(In) :: n
    Type(Atmosphere_Type), Intent(In) :: atm
    Real(dp) :: vA_AF(1:3)
    Real(dp) :: E_hat(1:3),N_hat(1:3),U_hat(1:3)
    
    vA_AF = 0._dp
    !Add rotation of the atmosphere
    If (Rotating_Earth) vA_AF = vA_AF + Atm_Rotation_Velocity_AF(n%big_r,n%r(3)/n%big_r)
    !Add wind
    If (Wind) vA_AF = vA_AF + atm%wind_AF
    !Convert from Air Frame to interial frame
    U_hat = Unit_Vector(n%r)
    E_hat = Unit_Vector(Cross_Product(Z_hat,U_hat))
    N_hat = Cross_Product(U_hat,E_hat)
    vA = E_hat * vA_AF(1) + &
       & N_hat * vA_AF(2) + &
       & U_hat * vA_AF(3)
End Function Air_Velocity

Function Apparent_Energy(E,Omega_hat,vA) Result(E_app)
    Use Kinds, Only: dp
    Use Neutron_Utilities, Only: Neutron_Energy
    Use Neutron_Utilities, Only: Neutron_Speed
    Implicit None
    Real(dp) :: E_app
    Real(dp), Intent(In) :: E,Omega_hat(1:3),vA(1:3)
    
    E_app = Neutron_Energy( Neutron_Speed(E) * Omega_hat - vA )
End Function Apparent_Energy

Function Scattered_Direction(mu,omega,A_hat,B_hat,C_hat) Result(Omega_hat)
    Use Kinds, Only: dp
    Use Global, Only: X_Hat, Z_Hat
    Use Utilities, Only: Unit_Vector
    Use Utilities, Only: Cross_Product
    Implicit None
    Real(dp) :: Omega_hat(1:3)  ! Direction after scatter
    Real(dp), Intent(In) :: mu  ! Deflection angle cosine
    Real(dp), Intent(In) :: omega  ! Orientation angle around Omega_hat0
    Real(dp), Intent(In) :: A_hat(1:3),B_hat(1:3),C_hat(1:3)
    
    Omega_hat = mu * A_hat + Sqrt(1._dp - mu**2) * (Cos(omega) * B_hat + Sin(omega) * C_hat)
End Function Scattered_Direction

Subroutine Scattered_Angles(Omega_hat0,Omega_hat1,mu,omega,B_hat,C_hat)
    Use Kinds, Only: dp
    Implicit None
    Real(dp), Intent(In) :: Omega_hat0(1:3) ! Direction before scatter
    Real(dp), Intent(In) :: Omega_hat1(1:3) ! Direction after scatter
    Real(dp), Intent(Out) :: mu  ! Deflection angle cosine
    Real(dp), Intent(Out) :: omega  ! Orientation angle around Omega_hat0
    Real(dp), Intent(In) :: B_hat(1:3),C_hat(1:3) 
    
    mu = Dot_Product(Omega_hat0,Omega_hat1)
    omega = Atan2( Dot_Product(Omega_hat1,B_hat) , Dot_Product(Omega_hat1,C_hat) )
End Subroutine Scattered_Angles

Subroutine Write_Scatter_Model(s,file_name)
    Use Cross_Sections, Only: Write_Cross_Sections
    Implicit None
    Type(Scatter_Model_Type), Intent(In) :: s
    Character(*), Intent(In) :: file_name
    Integer :: unit,stat
    
    Open(NEWUNIT = unit , FILE = file_name , STATUS = 'UNKNOWN' , ACTION = 'WRITE' , POSITION = 'APPEND' , IOSTAT = stat)
    If (stat .NE. 0) Then
        Print *,'ERROR:  Neutron_Scatter: Write_Scatter_Model:  File open error, '//file_name//', IOSTAT=',stat
        ERROR STOP
    End If
    Write(unit,'(A)') '--------------------------------------------------'
    Write(unit,'(A)') 'SCATTER MODEL INFORMATION'
    Write(unit,'(A)') '--------------------------------------------------'
    Write(unit,'(A)') '  Scatter Model Fidelity:'
    If (s%aniso_dist) Then
        Write(unit,'(A)') '    Scatter Model:       Anisotropic Center-of-Mass'
    Else
        Write(unit,'(A)') '    Scatter Model:       Isotropic Center-of-Mass'
    End If
    If (s%n_scatters .LT. 0) Then
        Write(unit,'(A)') '    Number of scatters:  UNLIMITED'
    Else
        Write(unit,'(A,I4)') '    Number of scatters: ',s%n_scatters
    End If
    If (s%estimate_each_scatter) Then
        Write(unit,'(A)') '                         Contributions estimated at each scatter'
    Else
        Write(unit,'(A)') '                         Contributions estimated at FINAL scatter only'
    End If
    If (s%direct_contribution) Then
        Write(unit,'(A)') '                         Direct Contribution (first event) SAMPLED'
    Else
        Write(unit,'(A)') '                         Direct Contribution (first event) NOT included'
    End If
    If (s%elastic_only) Then
        Write(unit,'(A)') '    Collision Types:     Elastic only'
    Else
        Write(unit,'(A)') '    Collision Types:     Elastic & Inelastic'
    End If
    If (s%fast_analog) Then
        Write(unit,'(A)') '    Collision Distance Estimator:  Fast Approximation to Analog'
    Else
        Write(unit,'(A)') '    Collision Distance Estimator:  Analog'
    End If
    Write(unit,'(A,L)') '    Absorption suppression: ',s%suppress_absorption
    Write(unit,'(A,L)') '    Leakage suppression:    ',s%suppress_leakage
    Write(unit,'(A,L)') '    All Materials & Mech:   ',s%all_mat_mech
    Write(unit,'(A,L)') '    Gravity:                ',s%Gravity
    Write(unit,'(A,L)') '    Neutron Decay:          ',s%Neutron_Decay
    Write(unit,'(A,L)') '    OTF Doppler Broaden:    ',s%Doppler_Broaden
    Write(unit,'(A,L)') '    Target Motion:          ',s%Target_Motion
    Write(unit,'(A,L)') '      Thermal Motion:         ',s%Thermal_Motion
    Write(unit,'(A,L)') '      Diatomic Atm:           ',s%Diatomic_atm
    Write(unit,'(A,L)') '      Rotating Earth:         ',s%Rotating_Earth
    Write(unit,'(A,L)') '      Global Wind:            ',s%wind
    Write(unit,*)
    Write(unit,'(A,I15)') '  Histories killed on:  Time:     ',s%n_kills(1)
    Write(unit,'(A,I15)') '                        Energy:   ',s%n_kills(2)
    Write(unit,'(A,I15)') '                        Weight:   ',s%n_kills(3)
    Write(unit,'(A,I15)') '                        Leaked:   ',s%n_kills(4)
    Write(unit,'(A,I15)') '                        Absorbed: ',s%n_kills(5)
    Write(unit,'(A,I15)') '                        Coll Lim: ',s%n_kills(6)
    Write(unit,'(A,I15)') '  Next-Events to detector:  Attempted: ',s%next_events(1)
    Write(unit,'(A,I15)') '                            Found:     ',s%next_events(2)
    Write(unit,'(A,I15)') '                            Tallied:   ',s%next_events(3)
    Write(unit,'(A,I15)') '  Next-Events NOT tallied:  Time-Energy: ',s%n_no_tally(1)
    Write(unit,'(A,I15)') '                            Time:        ',s%n_no_tally(2)
    Write(unit,'(A,I15)') '                            Energy:      ',s%n_no_tally(3)
    Write(unit,*)
    Write(unit,*)
    Close(unit)
    Call Write_Cross_Sections(s%CS,s%doppler_broaden,file_name)
End Subroutine Write_Scatter_Model

End Module Neutron_Scatter