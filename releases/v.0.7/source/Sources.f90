Module Sources
    
    Use Kinds, Only: dp
    Implicit None
    Private
    Public :: Source_Type
    Public :: Setup_Source
    Public :: Celest_to_XYZ
    Public :: Check_exo_source
    Public :: Write_Source
    
    Type :: Source_Type
        Real(dp) :: r(1:3)     ![km] cartesian (x,y,z) position of source
        Real(dp) :: v(1:3)     ![km/s] velocity of 
        Real(dp) :: z
        Real(dp) :: big_r
        Real(dp) :: speed
        Logical :: has_velocity  !T indicates source has velocity to be added to new neutrons
        Logical :: exoatmospheric  !T indicates source is outside atmosphere, F indicates source IN the atmosphere
        Integer :: dist_index
        Real(dp) :: E_high     ![keV] energy description, E_high=E_low for 'Point', specifies energy range for 'Uniform', specifies cutoff energies for 'Watt235'
        Real(dp) :: E_low
        Real(dp) :: w  !weight of particles as they are created at the source, always equal to 1 for atmospheric source, varies for exoatmospheric sources
        Real(dp) :: A_hat(1:3),B_hat(1:3),C_hat(1:3)  !basis vectors determining distribution of emission angles for direct contributions
    Contains
        Procedure, Pass :: Sample_Energy
    End Type
    
    !Integer designator for source energy distribution choice, can have value equal to one of the following parameters
    Integer, Parameter :: source_dist_Line = 71
    Integer, Parameter :: source_dist_Unif = 73
    Integer, Parameter :: source_dist_Watt235 = 75
    Integer, Parameter :: source_dist_Watt235_uni = 77
    
Contains

Function Setup_Source(setup_file_name,run_file_name,R_top_atm) Result(s)
    Use Kinds, Only: dp
    Use Global, Only: Z_hat,X_hat
    Use Global, Only: R_earth
    Use Utilities, Only: Unit_Vector
    Use Utilities, Only: Vector_Length
    Use Utilities, Only: Cross_Product
    Use FileIO_Utilities, Only: fSHARE
    Implicit None
    Type(Source_Type) :: s
    Character(*), Intent(In) :: setup_file_name,run_file_name
    Real(dp), Intent(In) :: R_top_atm
    Real(dp) :: x_source,y_source,z_source  ![km]  x,y,z coordinates of source
    Real(dp) :: declination_source,right_ascension_source
    Real(dp) :: v_E_source,v_N_source,v_U_source
    Real(dp) :: E_high,E_low
    Character(10) :: source_dist  !Line, Watt235, Uniform
    Character(9) :: position_geometry
    Integer :: setup_unit,stat
    Real(dp) :: E_hat(1:3),N_hat(1:3),U_hat(1:3)
    
    NameList /NeutronSourceList/ position_geometry,x_source,y_source,z_source, &
                                 & declination_source,right_ascension_source, &
                                 & v_E_source,v_N_source,v_U_source, &
                                 & source_dist,E_high,E_low
    
    Open(NEWUNIT = setup_unit , FILE = setup_file_name , STATUS = 'OLD' , ACTION = 'READ' , IOSTAT = stat , SHARE = fSHARE)
    If (stat .NE. 0) Then
        Print *,'ERROR:  Neutron_Source: Initialize_Neutron_Source:  File open error, '//setup_file_name//', IOSTAT=',stat
        ERROR STOP
    End If
    Read(setup_unit,NML = NeutronSourceList)
    Close(setup_unit)
    Select Case(position_geometry)
        Case('Celestial')
            s%r = Celest_to_XYZ(z_source,right_ascension_source,declination_source)
        Case('Cartesian')
            s%r = (/ x_source, y_source, z_source /)
        Case Default
            Print *,'ERROR:  Sources: Setup_Source:  Unknown position geometry.'
            ERROR STOP
    End Select
    s%big_r = Vector_Length(s%r)
    s%z = s%big_r - R_earth
    !check for exoatmospheric source
    Call Check_exo_source(s%big_r,R_top_atm,s%exoatmospheric,s%w)
    !Set source velocity
    s%v = (/ v_E_source, &
           & v_N_source, &
           & v_U_source /)
    If (Any(s%v .NE. 0._dp)) Then
        s%has_velocity = .TRUE.
        U_hat = Unit_Vector(s%r)
        E_hat = Unit_Vector(Cross_Product(Z_hat,U_hat))
        N_hat = Cross_Product(U_hat,E_hat)
        s%v = E_hat * Dot_Product(s%v,E_hat) + &
            & N_hat * Dot_Product(s%v,N_hat) + &
            & U_hat * Dot_Product(s%v,U_hat)
        s%speed = Vector_Length(s%v)
        s%A_hat = Unit_Vector(s%v)
        If (Abs(Dot_Product(s%A_hat,Z_Hat)) .LT. 0.5_dp) Then
            s%B_hat = Unit_Vector(Cross_Product(Z_Hat,s%A_hat))
        Else
            s%B_hat = Unit_Vector(Cross_Product(X_Hat,s%A_hat))
        End If
        s%C_hat = Cross_Product(s%A_hat,s%B_hat)
    Else
        s%has_velocity = .FALSE.
        s%speed = 0._dp
        s%A_hat = 0._dp
        s%B_hat = 0._dp
        s%C_hat = 0._dp
    End If
    Select Case (source_dist)
        Case ('Line')
            s%dist_index = source_dist_Line
            s%E_high = E_high
            s%E_low = E_high
        Case ('Uniform')
            s%dist_index = source_dist_Unif
            s%E_high = E_high
            s%E_low = E_low
        Case ('Watt235')
            s%dist_index = source_dist_Watt235
            s%E_high = E_high
            s%E_low = 0._dp
        Case ('Watt235uni')
            s%dist_index = source_dist_Watt235_uni
            s%E_high = E_high
            s%E_low = E_low
        Case Default
            Print *,'ERROR:  Neutron_Source: Initialize_Neutron_Source: Undefined source energy distribution'
            ERROR STOP
    End Select
    If (this_image() .EQ. 1) Then
        Open(NEWUNIT = setup_unit , FILE = run_file_name , STATUS = 'OLD' , ACTION = 'WRITE' , POSITION = 'APPEND' , IOSTAT = stat)
        If (stat .NE. 0) Then
            Print *,'ERROR:  Sources: Setup_Source:  File open error, '//run_file_name//', IOSTAT=',stat
            ERROR STOP
        End If
        Write(setup_unit,NML = NeutronSourceList)
        Write(setup_unit,*)
        Close(setup_unit)
    End If
End Function Setup_Source

Function Celest_to_XYZ(alt,RA_deg,DEC_deg) Result(r)
    Use Kinds, Only: dp
    Use Global, Only: Pi
    Use Global, Only: R_Earth
    Implicit None
    Real(dp) :: r(1:3)
    Real(dp), Intent(In) :: alt,RA_deg,DEC_deg
    Real(dp) :: RA,DEC
    
    RA = RA_deg * Pi / 180._dp
    DEC = DEC_deg * Pi / 180._dp
    !convert celestial to cartesian soordinates
    r = (/ (R_Earth + alt) * Cos(RA) * Cos(DEC), &
         & -(R_Earth + alt) * Sin(RA) * Cos(DEC), &
         & (R_Earth + alt) * Sin(DEC) /)
End Function Celest_to_XYZ

Subroutine Check_exo_source(Rs,Ra,exo,w)
    Use Kinds, Only: dp
    Use Global, Only: R_Earth
    Use Utilities, Only: Vector_Length
    Implicit None
    Real(dp), Intent(In) :: Rs,Ra
    Logical, Intent(Out) :: exo
    Real(dp), Intent(Out) :: w
    
    If (Rs .LT. Ra) Then
        exo = .FALSE.
        w = 1._dp
    Else
        exo = .TRUE.
        w = 0.5_dp * (1._dp - Sqrt(Rs**2 - Ra**2) / Rs)
    End If
End Subroutine Check_exo_source

Subroutine Sample_Energy(source,RNG,E,w)
    Use Kinds, Only: dp
    Use Random_Numbers, Only: RNG_Type
    Implicit None
    Class(Source_Type), Intent(In) :: source
    Type(RNG_type), Intent(InOut) :: RNG
    Real(dp), Intent(Out) :: E
    Real(dp), Intent(Out) :: w
    
    Select Case (source%dist_index)
        Case (source_dist_Line)
            E = source%E_high
            w = 1._dp
        Case (source_dist_Unif)
            E = Sample_Uniform(RNG,source%E_low,source%E_high)
            w = 1._dp
        Case (source_dist_Watt235)
            E = Sample_Watt235(RNG,source%E_high)
            w = 1._dp
        Case (source_dist_Watt235_uni)
            E = Sample_Uniform(RNG,source%E_low,source%E_high)
            w = 0.000452463_dp * Watt235func(E)            
        Case Default
            Print *,'ERROR:  Neutron_Source: Sample_Neutron_Source_Energy: Undefined source energy distribution'
            ERROR STOP
    End Select
End Subroutine Sample_Energy

Function Sample_Uniform(RNG,E_min,E_max) Result(E)
    Use Kinds, Only: dp
    Use Random_Numbers, Only: RNG_Type
    Implicit None
    Real(dp) :: E  ![keV]
    Type(RNG_Type), Intent(InOut) :: RNG  !random number generator
    Real(dp), Intent(In) :: E_min,E_max  ![keV]  energy range

    E = E_min + RNG%Get_Random() * (E_max - E_min)
End Function Sample_Uniform

Function Sample_Watt235(RNG,E_max) Result(E)
    Use Kinds, Only: dp
    Use Random_Numbers, Only: RNG_Type
    Implicit None
    Real(dp) :: E  ![keV]
    Type(RNG_Type), Intent(InOut) :: RNG  !random number generator
    Real(dp), Intent(In) :: E_max  ![keV]  max energy allowed
    Real(dp) :: w

    !N2H Investigate faster and/or alternate watt sampling schemes
    Do
        E = E_max * RNG%Get_Random()
        w = 1.26463489196562161942_dp * Watt235func(E)
        If (w .GT. RNG%Get_Random()) Exit
    End Do
End Function Sample_Watt235

Function Watt235func(E) Result(w)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: w
    Real(dp), Intent(In) :: E
    
    w = Exp(-0.001036_dp * E) * Sinh(0.047853944456021595545_dp * Sqrt(E))
End Function Watt235func

Subroutine Write_Source(s,file_name)
    Use Global, Only: Pi
    Use Global, Only: R_Earth
    Use Utilities, Only: Vector_Length
    Implicit None
    Type(Source_Type), Intent(In) :: s
    Character(*), Intent(In) :: file_name
    Integer :: unit,stat
    
    Open(NEWUNIT = unit , FILE = file_name , STATUS = 'UNKNOWN' , ACTION = 'WRITE' , POSITION = 'APPEND' , IOSTAT = stat)
    If (stat .NE. 0) Then
        Print *,'ERROR:  Sources: Write_Source:  File open error, '//file_name//', IOSTAT=',stat
        ERROR STOP
    End If
    Write(unit,'(A)') '--------------------------------------------------'
    Write(unit,'(A)') 'SOURCE INFORMATION'
    Write(unit,'(A)') '--------------------------------------------------'
    Write(unit,'(A)') '  Position:'
    Write(unit,'(A,ES23.16E3,A)') '    x = ',s%r(1),' km'
    Write(unit,'(A,ES23.16E3,A)') '    y = ',s%r(2),' km'
    Write(unit,'(A,ES23.16E3,A)') '    z = ',s%r(3),' km'
    Write(unit,'(A,ES23.16E3,A)') '    Right Ascension = ',Acos(s%r(1)/(Vector_Length(s%r)*Cos(Asin(s%r(3)/Vector_Length(s%r))))) / (Pi/180._dp),' deg'
    Write(unit,'(A,ES23.16E3,A)') '    Declination     = ',Asin(s%r(3)/Vector_Length(s%r)) / (Pi/180._dp),' deg'
    Write(unit,'(A,ES23.16E3,A)') '    Geometric Alt   = ',Vector_Length(s%r)-R_Earth,' km'
    If (s%exoatmospheric) Write(unit,'(A,ES23.16E3,A)') '    Exatmospheric source: ',s%w,' starting particle weight'
    Write(unit,'(A)') '  Velocity:'
    Write(unit,'(A,ES23.16E3,A)') '    x = ',s%v(1),' km/s'
    Write(unit,'(A,ES23.16E3,A)') '    y = ',s%v(2),' km/s'
    Write(unit,'(A,ES23.16E3,A)') '    z = ',s%v(3),' km/s'
    Write(unit,*)
    Write(unit,'(A)') '  Energy Distribution:'
    Select Case (s%dist_index)
        Case (source_dist_Line)
            Write(unit,'(A,ES23.16E3,A)') '    Line, ',s%E_high,' keV'
        Case (source_dist_Unif)
            Write(unit,'(A,ES23.16E3,A,ES23.16E3,A)') '    Uniform, ',s%E_low,' keV to ',s%E_high,' keV'
        Case (source_dist_Watt235)
            Write(unit,'(A,ES22.16E3,A)') '    Watt-235, ',s%E_high,' keV max'
        Case (source_dist_Watt235_uni)
            Write(unit,'(A,ES23.16E3,A,ES23.16E3,A)') '    Watt-235-uniform, ',s%E_low,' keV to ',s%E_high,' keV'
        Case Default
            Print *,'ERROR:  Sources: Write_Source: Undefined source energy distribution'
            ERROR STOP
    End Select
    Write(unit,*)
    Write(unit,*)
    Close(unit)
End Subroutine Write_Source

End Module Sources