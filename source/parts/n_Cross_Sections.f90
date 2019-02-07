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
Module n_Cross_Sections

    Use Kinds, Only: dp
    Use Global, Only: k_Boltzmann
    Implicit None
    Private
    Public :: CS_Type
    Public :: sig_Composite
    Public :: Setup_Cross_Sections
    Public :: Write_Cross_Sections
    Public :: Read_CS_file

    Type :: sig_Type
        Integer :: n_sig
        Real(dp), Allocatable :: sig(:) ![barns]  has dimension 1:n_sig
        Real(dp), Allocatable :: lnsig(:) ![barns]  has dimension 1:n_sig
        Integer, Allocatable :: E_map(:)  !has dimension 1:n_E_uni, indexes for each value in unified energy grid to indexes in sig
        Integer, Allocatable :: E_key(:)  !has dimension 1:n_sig, indexes for each value in sig to an energy in unified energy list
        Integer :: n_interp_r  !number of interpolation ranges
        Integer, Allocatable :: interp(:,:)  !has dimension 1:n_interp_r and 1:2, dim 1 is sig index up to which to use the interpolation method specified in dim 2
    End Type

    Type :: da_List_Type
        Integer :: n_a  !number of coeffs/points
        Logical :: is_Legendre
        Logical :: is_tab
        Real(dp), Allocatable :: a(:)  !has dim 0:n_a, list of legendre coefficients or pdf values
        Real(dp), Allocatable :: ua(:,:)  !has dim 1:n_a,1:2, for tabulated cosine pdf, list of cosines
    End Type

    Type :: da_Type
        Integer :: n_da
        Type(da_List_Type), Allocatable :: da(:) ![barns]  has dimension 1:n_da
        Integer, Allocatable :: E_map(:)  !has dimension 1:n_E_uni, indexes for each value in unified energy grid to indexes in da
        Integer, Allocatable :: E_key(:)  !has dimension 1:n_da, indexes for each value in da to an energy in unified energy list
    End Type

    Type :: lev_sig_Type
        Integer :: n_lev
        Real(dp), Allocatable :: Q(:)  ![positive keV]  has dimension 0:n_lev, Q-value for each level
        Integer, Allocatable :: thresh(:)  !has dimension 0:n_lev, index of level threshold energy in unified energy grid
        Type(sig_Type), Allocatable :: sig(:)  !has dimension 0:n_lev, cross sections for each level
        Type(da_Type), Allocatable :: da(:)  !has dimension 0:n_lev, angular cross sections for each level
    End Type

    Type :: abs_sig_Type
        Integer :: n_modes
        Integer, Allocatable :: thresh(:)  !has dimension 1:modes, index of mode threshold energy in unified energy grid
        Type(sig_Type), Allocatable :: sig(:)  !has dim 1:n_modes, cross sections for each absorption mode
    End Type

    Type :: res_sig_Type
        !HACK Resonance CS is included as linearly interpolable listing rather than by actual resonance reconstruction, this violates the small-memory intent of this cross sections implementation
        Real(dp) :: E_range(1:2)  ![keV] Lower and upper limits of the range over which resonant cross sections are included
        Integer :: n_E
        Real(dp), Allocatable :: E(:)  ![keV] has dimension 1:n_E, energy points at which resonance cross section is given
        Real(dp), Allocatable :: sig(:,:)  ![barns] has dimension 1:n_E and 1:2, dim 2 is 1=n-gamma resonance and 2=elastic scatter resonance
    End Type

    Type ::  res_sig_spin_Type
        Integer :: n_r  !number of resonant energies in this spin
        Real(dp), Allocatable :: ErGnGr(:,:)  !has dimension 1:3 and 1:n_r, dim 1 is 1=Er[keV], dim 2 is neutron width, dim 3 is radiation width
    End Type

    Type ::  res_sig_level_Type
        Integer :: n_J  !number of resonant spins in this level
        Real(dp), Allocatable :: gJ(:)  !total angular momentum for each spin group
        Real(dp), Allocatable :: dJ(:)  !second channel spin correction for each spin group
        Type(res_sig_spin_Type), Allocatable :: J(:)  !has dimension 1:n_J, resonance parameters for each spin group
    End Type

    Type :: res_sig_Type_new
        Real(dp) :: k0  !=2.196771E-3_dp*(AWRI/(AWRI+1._dp)), gives k when miltiplied by Sqrt(E)
        Integer :: n_L  !number of levels in which resonances are grouped
        Real(dp), Allocatable :: AP(:)  !has dimension 1:n_L, scattering radius for each level
        Type(res_sig_level_Type), Allocatable :: L(:)   !has dimension 1:n_L, resonance parameters for each leval group
    End Type

    Type :: CS_Type
        Integer :: n_E_uni  !number of energies in the unified energy list
        Real(dp), Allocatable :: E_uni(:) ![keV]  has dimension 1:n_E_uni, list of energies in the unified energy grid
        Real(dp), Allocatable :: lnE_uni(:) ![keV]  has dimension 1:n_E_uni, list of energies in the unified energy grid
        Integer :: n_iso  !number of isotopes in cross sections structure
        Real(dp), Allocatable :: iso_Fractions(:)  !has dimension 1:n_iso, fractional abundance of each isotope in the total atmosphere
        Real(dp), Allocatable :: An(:)  ![neutron masses] has dimension 1:n_iso, mass (in neutron masses) of isotope nucleus
        Real(dp) :: Mn  ![kg] mean mass (in KILOGRAMS) of nuclei in the atmosphere
        Integer :: n_a_max,n_a_tab_max  !max number of coefficents in angular distribution cross sections
        Type(lev_sig_Type), Allocatable :: lev_cs(:)  !has dimension 1:n_iso, cross sections and angular distributions for each inelastic level for each isotope
        Type(abs_sig_Type), Allocatable :: abs_cs(:)  !has dimension 1:n_iso, cross sections for each absorption type for each isotope
        Logical, Allocatable :: has_res_cs(:)  !has dimension 1:n_iso, TRUE indicates resonance cross sections are included in an isotope cross section representation
        Type(res_sig_Type), Allocatable :: res_cs(:)  !has dimension 1:n_iso, resonance cross sections for each isotope
    Contains
        Procedure, Pass :: sig_T  !given energy, returns total microscopic cross section for the atmosphere
        Procedure, Pass :: sig_S  !given energy, returns microscopic scatter cross section for the atmosphere
        Procedure, Pass :: sig_A  !given energy, returns microscopic absorption cross section for the atmosphere
        Procedure, Pass :: sig_T_A  !given energy, returns total and absorption microscopic cross sections for the atmosphere
        Procedure, Pass :: sig_S_iso  !given energy, returns microscopic scatter cross section for a specific isotope
        Procedure, Pass :: sig_T_broad
        Procedure, Pass :: sig_T_A_broad
        Procedure, Pass :: sig_S_iso_broad
    End Type

    Real(dp), Parameter :: k_B = k_Boltzmann / 1000._dp**2  ![J/K]*[1 km^2 / 1000^2 m^2] Boltzmann constant with convenient units locally
    Real(dp), Parameter :: rTol = 1.E-5_dp  !relative tolerance for convergence of broadening integrals, cross section data has about 5 good digits...
    
    Integer, Parameter :: MT_disappearance(1:15) = (/ 102, &  !n,g
                                                    & 103, &  !n,p
                                                    & 104, &  !n,d
                                                    & 105, &  !n,t
                                                    & 106, &  !n,He3
                                                    & 107, &  !n,a
                                                    & 108, &  !n,2a
                                                    & 109, &  !n,3a
                                                    & 111, &  !n,2p
                                                    & 112, &  !n,p+a
                                                    & 113, &  !n,t+2a
                                                    & 114, &  !n,d+2a
                                                    & 115, &  !n,p+d
                                                    & 116, &  !n,p+t
                                                    & 117  /) !n,d+a
    Integer, Parameter :: MT_inelastic(1:40) = (/ (50+i , i=1,40) /)  !n,n' of the i-th level
    Integer, Parameter :: MT_excluded(1:22) = (/  5, &  !interactions not included in any other MT
                                               & 11, &  !n,2n+d
                                               & 16, &  !n,2n
                                               & 17, &  !n,3n
                                               & 22, &  !n,n+a
                                               & 23, &  !n,n+3a
                                               & 24, &  !n,2n+a
                                               & 25, &  !n,3n+a
                                               & 28, &  !n,n+p
                                               & 29, &  !n,n+2a
                                               & 30, &  !n,2n+2a
                                               & 32, &  !n,n+d
                                               & 33, &  !n,n+t
                                               & 34, &  !n,n+He3
                                               & 35, &  !n,n+d+2a
                                               & 36, &  !n,n+t+2a
                                               & 37, &  !n,4n
                                               & 41, &  !n,2n+p
                                               & 42, &  !n,3n+p
                                               & 44, &  !n,n+2p
                                               & 45, &  !n,n+p+a
                                               & 91  /) !n,n' not included in MTs 51-90

Contains

Function Setup_Cross_Sections(resources_directory,cs_setup_file,elastic_only,aniso_dist,E_min,E_max) Result(CS)
    Use Kinds, Only: dp
    Use Sorting, Only: Union_Sort
    Use FileIO_Utilities, Only: max_path_len
    Use FileIO_Utilities, Only: slash
    Use FileIO_Utilities, Only: Output_Message
    Use Global, Only: neutron_mass
    Implicit None
    Type(CS_Type) :: CS
    Character(*), Intent(In) :: resources_directory
    Character(*), Intent(In) :: cs_setup_file
    Logical, Intent(In) :: elastic_only
    Logical, Intent(In) :: aniso_dist
    Real(dp), Intent(In) :: E_min,E_max
    Integer :: n_elements
    Character(:), Allocatable :: file_name_start,cs_file_name
    Integer, Allocatable :: n_isotopes(:)
    Character(4), Allocatable :: isotope_names(:)
    Real(dp), Allocatable :: el_fractions(:),iso_fractions(:)
    Character(2) :: j_char
    Integer :: n_energies
    Integer :: n_p,n_r,n_start
    Integer :: i,j,k
    Real(dp) :: Q_scratch
    Real(dp) :: An_scratch
    Real(dp), Allocatable :: E_uni_scratch(:)
    Real(dp), Allocatable :: E_scratch(:)
    Real(dp), Allocatable :: CS_scratch(:)
    Integer, Allocatable :: Interp_scratch(:,:)
    Type(da_List_type), Allocatable :: Ang_Dist_scratch(:)
    Integer :: setup_unit,stat
    Integer :: n_absorption_modes,n_inelastic_lev
    Logical, Allocatable :: diatomic(:)
    Integer :: ltt
    Logical :: has_resonance
    Real(dp) :: iso_fraction

    NameList /csSetupList1/ n_elements
    NameList /csSetupList2/ el_fractions,n_isotopes
    NameList /csSetupList3/ isotope_names,diatomic

    !read namelists from cross sections setup file
    Allocate(Character(max_path_len) :: file_name_start)
    file_name_start = resources_directory//'cs'//slash//'n_cs'//slash
    Open(NEWUNIT = setup_unit , FILE = file_name_start//cs_setup_file , STATUS = 'OLD' , ACTION = 'READ' , IOSTAT = stat)
    If (stat .NE. 0) Call Output_Message('ERROR:  Cross_Sections: Setup_Cross_Sections:  File open error, '//file_name_start//cs_setup_file//', IOSTAT=',stat,kill=.TRUE.)
    Read(setup_unit,NML = csSetupList1)
    Allocate(el_fractions(1:n_elements))
    Allocate(n_isotopes(1:n_elements))
    Read(setup_unit,NML = csSetupList2)
    el_fractions = el_fractions / Sum(el_fractions) !make sure el_fractions is normalized to 1
    CS%n_iso = Sum(n_isotopes)
    Allocate(isotope_names(1:CS%n_iso))
    Allocate(iso_fractions(1:CS%n_iso))
    Allocate(diatomic(1:CS%n_iso))
    Read(setup_unit,NML = csSetupList3)
    j = 1
    Do i = 1,n_elements
        !make sure iso_fractions are normalized to 1 within each isotope
        iso_fractions(j:Sum(n_isotopes(1:i))) =  iso_fractions(j:Sum(n_isotopes(1:i))) / Sum(iso_fractions(j:Sum(n_isotopes(1:i))))
        !combine element and isotpe fractions for total atmospheric fraction of each isotope
        iso_fractions(j:Sum(n_isotopes(1:i))) = iso_fractions(j:Sum(n_isotopes(1:i))) * el_fractions(i)
        j = Sum(n_isotopes(1:i)) + 1
    End Do
    Close(setup_unit)
    !Initialize the CS structure
    Allocate(CS%iso_fractions(1:CS%n_iso))
    CS%iso_fractions = iso_fractions
    Allocate(CS%An(1:CS%n_iso))
    CS%An = -1._dp
    Allocate(CS%lev_cs(1:CS%n_iso))
    Allocate(CS%abs_cs(1:CS%n_iso))
    Allocate(CS%has_res_cs(1:CS%n_iso))
    CS%has_res_cs = .FALSE.
    Allocate(CS%res_cs(1:CS%n_iso))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  CREATE A UNIFIED ENERGY LIST FOR THE CROSS SECTION DATA
    !count the number of energies (including duplicates) for MF 3 and 4 (interaction cross sections and angular distributions)
    !also count number of disappearance modes (MF=3 and MTs appearing in MT_disappearance)
    !also count number of inelastic levels (MF=3 and MTs 51-90)
    Allocate(Character(max_path_len) :: cs_file_name)
    n_energies = 0
    n_absorption_modes = 0
    n_inelastic_lev = 0
    Do i = 1,CS%n_iso
        !create file name string and open the ENDF tape for this isotope
        cs_file_name = file_name_start//Trim(isotope_names(i))//'.txt'
        Open(NEWUNIT = cs_unit , FILE = cs_file_name , STATUS = 'OLD' , ACTION = 'READ' , IOSTAT = stat)
        If (stat .NE. 0) Call Output_Message('ERROR:  Cross_Sections: Setup_Cross_Sections:  File open error, '//cs_file_name//', IOSTAT=',stat,kill=.TRUE.)
        !count energies in absorption, elastic, and inelastic files (MF 3 and 4)
        DO_SECTIONS: Do
            !check next section type
            Read(cs_unit,'(A73,I1,I3,I5)') trash_c,MF,MT,line_num
            If (MF .EQ. 0) Then  !this indicates a change in MF type, need to look for next section or end of file
                NEXT_MF: Do
                    Read(cs_unit,'(A73,I1,I3,I5)',IOSTAT=stat) trash_c,MF,MT,line_num
                    If (stat .LT. 0) Exit DO_SECTIONS !end of file
                    If (MF .GT. 0) Exit NEXT_MF
                End Do NEXT_MF
            End If
            If (MF.EQ.3 .OR. MF.EQ.4) Then
                If (Any(MT_excluded .EQ. MT)) Then !this interaction type is excluded, advance past it
                    !advance to end of section
                    NEXT_MT: Do
                        Read(cs_unit,'(A73,I1,I3,I5)') trash_c,MF,MT,line_num
                        If (line_num = 99999) Exit NEXT_MT  !this record indicates end of section
                    End Do NEXT_MT
                    !cycle to start next section
                    Cycle DO_SECTIONS
                End If
                !get number of energies in this MF 3 or 4 section
                Select Case (MF)
                    Case (3)
                        If (Any(MT_dissappearance .EQ. MT) n_abs_modes = n_abs_modes + 1
                        If (Any(MT_inelastic .EQ. MT) n_inel_lev = n_inel_lev + 1
                        Read(cs_unit,'(A55,I11)') trash_c, n_p
                    Case (4)
                        !need to read the first line again to get LTT
                        Backspace(cs_unit)
                        Read(cs_unit,'(A33,I11)') trash_c, LTT
                        Read(cs_unit,*)
                        Read(cs_unit,'(I11)') n_p
                        If (LTT .EQ. 3) Then !there is a second range of energies later in the section
                            !advance in the file to the end of the Legendre section
                            Do j = 1,n_p
                                !The first line in each energy contains the energy in eV in the second position and the number of Legendre coefficents in the 5th position
                                Read(cs_unit,'(A44,I11)') trash_c, n_a
                                n_a_lines = (n_a / 6)  !integer divide
                                If (Mod(n_a,6) .GT. 0) n_a_lines = n_a_lines + 1
                                !advance to the next energy
                                Do k = 1,n_a_lines
                                    Read(cs_unit,*)
                                End Do
                            End Do
                            Read(cs_unit,*)
                            !first entry of next line is number of additional energy points
                            Read (cs_unit,'(I11)') n_p_2
                        Else
                            n_p_2 = 0
                        End If
                        np = np + n_p_2
                End Select
                n_energies = n_energies + n_p
            End If
            !advance to end of section
            NEXT_MT: Do
                Read(cs_unit,'(A73,I1,I3,I5)') trash_c,MF,MT,line_num
                If (line_num = 99999) Exit NEXT_MT  !this record indicates end of section
            End Do NEXT_MT
        End Do DO_SECTIONS
        Close(cs_unit)
    End Do
    !Allocate scratch arrays for ordering absorption modes
    Allocate(abs_modes(1:n_abs_modes))
    abs_modes = -1
    Allocate(abs_thresh(1:n_abs_modes))
    abs_thresh = Huge(abs_thresh)
    !Allocate scratch arrays for energy
    Allocate(E_uni_scratch(1:n_energies))
    E_uni_scratch = -1._dp
    !Read in ALL (including duplicates) the energies to be used
    n_start = 0
    Do i = 1,CS%n_iso
        !create file name string and open the ENDF tape for this isotope
        cs_file_name = file_name_start//Trim(isotope_names(i))//'.txt'
        Open(NEWUNIT = cs_unit , FILE = cs_file_name , STATUS = 'OLD' , ACTION = 'READ' , IOSTAT = stat)
        If (stat .NE. 0) Call Output_Message('ERROR:  Cross_Sections: Setup_Cross_Sections:  File open error, '//cs_file_name//', IOSTAT=',stat,kill=.TRUE.)
        !read energies in absorption, elastic, and inelastic files (MF 3 and 4)
        DO_SECTIONS: Do
            !check next section type
            Read(cs_unit,'(A73,I1,I3,I5)') trash_c,MF,MT,line_num
            If (MF .EQ. 0) Then  !this indicates a change in MF type, need to look for next section or end of file
                NEXT_MF: Do
                    Read(cs_unit,'(A73,I1,I3,I5)',IOSTAT=stat) trash_c,MF,MT,line_num
                    If (stat .LT. 0) Exit DO_SECTIONS !end of file
                    If (MF .GT. 0) Exit NEXT_MF
                End Do NEXT_MF
            End If
            If (MF.EQ.3 .OR. MF.EQ.4) Then
                If (Any(MT_excluded .EQ. MT)) Then !this interaction type is excluded, advance past it
                    !advance to end of section
                    NEXT_MT: Do
                        Read(cs_unit,'(A73,I1,I3,I5)') trash_c,MF,MT,line_num
                        If (line_num = 99999) Exit NEXT_MT  !this record indicates end of section
                    End Do NEXT_MT
                    !cycle to start next section
                    Cycle DO_SECTIONS
                End If
                !get number of energies in this MF 3 or 4 section
                Select Case (MF)
                    Case (3)
                        Read(cs_unit,'(A55,I11)') trash_c, n_p
                        Read(cs_unit,*)
                        Do j = 1,n_p,3
                            !Each line in the file has 3 pairs of (eV,barns)
                            If (n_p-j .GT. 1) Then
                                Read(cs_unit,'(5E11.6E1)') E_uni_scratch(n_start+j), trash_r, E_uni_scratch(n_start+j+1), trash_r, E_uni_scratch(n_start+j+2)
                            Else If (n_p-j .EQ. 1) Then
                                Read(cs_unit,'(3E11.6E1)') E_uni_scratch(n_start+j), trash_r, E_uni_scratch(n_start+j+1)
                            Else
                                Read(cs_unit,'(1E11.6E1)') E_uni_scratch(n_start+j)
                            End If
                            If (j .EQ. 1) Then !check if this is an absorption mode
                                If (Any(MT_disappearance .EQ. MT)) Then  !this is an absorption mode
                                    !add it to the list of absorption modes, inserting in order of acending threshold energy
                                    Do k = 1,n_abs_modes
                                        If (E_uni_scratch(n_start+j) .LT. abs_thresh(k)) Then !this k is where this mode goes in the list
                                            If (k .LT. n_abs_modes) Then  !make room by shifting the rest of the list down
                                                abs_thresh(k+1:n_abs_modes) = abs_thresh(k:n_abs_modes-1)
                                                abs_modes(k+1:n_abs_modes) = abs_modes(k:n_abs_modes-1)
                                            End If
                                            !insert
                                            abs_thresh(k) = E_uni_scratch(n_start+j)
                                            abs_modes(k) = MT
                                        End If
                                    End Do
                                End If
                            End If
                        End Do
                    Case (4)
                        !need to read the first line again to get LTT
                        Backspace(cs_unit)
                        Read(cs_unit,'(A33,I11)') trash_c, LTT
                        Read(cs_unit,*)
                        Read(cs_unit,'(I11)') n_p
                        If (LTT .EQ. 3) Then !there is a second range of energies later in the section
                            !advance in the file to the end of the Legendre section
                            Do j = 1,n_p
                                !The first line in each energy contains the energy in eV in the second position and the number of Legendre coefficents in the 5th position
                                Read(cs_unit,'(A44,I11)') trash_c, n_a
                                n_a_lines = (n_a / 6)  !integer divide
                                If (Mod(n_a,6) .GT. 0) n_a_lines = n_a_lines + 1
                                !advance to the next energy
                                Do k = 1,n_a_lines
                                    Read(cs_unit,*)
                                End Do
                            End Do
                            Read(cs_unit,*)
                            !first entry of next line is number of additional energy points
                            Read (cs_unit,'(I11)') n_p_2
                        Else
                            n_p_2 = 0
                        End If
                        If (LTT .EQ. 1) Then  !da is lists of legendre coeffs
                            Do j = 1,n_p
                                !The first line in each energy contains the energy in eV in the second position and the number of Legendre coefficents in the 5th position
                                Read(cs_unit,'(A11,E11.6E1,A22,I11)') trash_c, E_uni_scratch(n_start+j), trash_c, n_a
                                n_a_lines = (n_a / 6)  !integer divide
                                If (Mod(n_a,6) .GT. 0) n_a_lines = n_a_lines + 1
                                !advance to the next energy
                                Do k = 1,n_a_lines
                                    Read(cs_unit,*)
                                End Do
                            End Do
                        Else If (LTT .EQ. 2) Then  !da is tabulated
                            Do j = 1,n_p
                                !The first line in each energy contains the energy in eV in the second position
                                Read(cs_unit,'(A11,E11.6E1)') trash_c, E_uni_scratch(n_start+j)
                                !the next line contains the number of tabulation points in the first position
                                Read(cs_unit,'(I11)') n_a
                                n_a_lines = (n_a / 3)  !integer divide
                                If (Mod(n_a,3) .GT. 0) n_a_lines = n_a_lines + 1
                                !advance to the next energy
                                Do k = 1,n_a_lines
                                    Read(cs_unit,*)
                                End Do
                            End Do
                        Else If (LTT .EQ. 3) Then  !da is tabulated for high energies but legendre for low energies
                            !Read in low energy Legendre points
                            Do j = 1,n_p
                                !The first line in each energy contains the energy in eV in the second position and the number of Legendre coefficents in the 5th position
                                Read(cs_unit,'(A11,E11.6E1,A22,I11)') trash_c, E_uni_scratch(n_start+j), trash_c, n_a
                                n_a_lines = (n_a / 6)  !integer divide
                                If (Mod(n_a,6) .GT. 0) n_a_lines = n_a_lines + 1
                                !advance to the next energy
                                Do k = 1,n_a_lines
                                    Read(cs_unit,*)
                                End Do
                            End Do
                            Read(cs_unit,*)
                            Read(cs_unit,*)
                            !Read in high energy tabulated cosine points
                            Do j = n_p+1,n_p+n_p_2
                                !The first line in each energy contains the energy in eV in the second position
                                Read(cs_unit,'(A11,E11.6E1)') trash_c, E_uni_scratch(n_start+j)
                                !the next line contains the number of tabulation points in the first position
                                Read(cs_unit,'(I11)') n_a
                                n_a_lines = (n_a / 3)  !integer divide
                                If (Mod(n_a,3) .GT. 0) n_a_lines = n_a_lines + 1
                                !advance to the next energy
                                Do k = 1,n_a_lines
                                    Read(cs_unit,*)
                                End Do
                            End Do
                            !update total n_p
                            np = np + n_p_2
                        End If
                End Select
                n_start = n_start + n_p
            End If
            !advance to end of section
            NEXT_MT: Do
                Read(cs_unit,'(A73,I1,I3,I5)') trash_c,MF,MT,line_num
                If (line_num = 99999) Exit NEXT_MT  !this record indicates end of section
            End Do NEXT_MT
        End Do DO_SECTIONS
        Close(cs_unit)    
    End Do
    E_uni_scratch = E_uni_scratch / 1000._dp  !convert to keV
    !E_uni_scratch is now a HUGE list of all the energies in all the files we are going to use, sort and eliminate duplicates
    Call Union_Sort(E_uni_scratch,n_energies,E_min,E_max)
    If (n_energies .GT. Huge(CS%n_E_uni)) Call Output_Message('ERROR:  Cross_Sections: Setup_Cross_Sections:  Length of unified energy grid exceeds available index',kill=.TRUE.)
    !Allocate and fill the unified energy list
    CS%n_E_uni = n_energies
    Allocate(CS%E_uni(1:n_energies))
    CS%E_uni = E_uni_scratch(1:n_energies)
    Deallocate(E_uni_scratch)
    Allocate(CS%lnE_uni(1:n_energies))
    CS%lnE_uni = Log(CS%E_uni)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  READ IN INTERACTION CROSS SECTIONS AND ANGULAR DISTRIBUTIONS
    !Now read each file (yes, again) into its appropriate place in the cross section structure, constructing the integer maps and keys as we go
    CS%n_a_max = 0
    CS%n_a_tab_max = 0
    Do i = 1,CS%n_iso
        !create file name string and open the ENDF tape for this isotope
        cs_file_name = file_name_start//Trim(isotope_names(i))//'.txt'
        Open(NEWUNIT = cs_unit , FILE = cs_file_name , STATUS = 'OLD' , ACTION = 'READ' , IOSTAT = stat)
        If (stat .NE. 0) Call Output_Message('ERROR:  Cross_Sections: Setup_Cross_Sections:  File open error, '//cs_file_name//', IOSTAT=',stat,kill=.TRUE.)
        !read in absorption, resonance, elastic, and inelastic interaction cross sections and angular distributions
        Do j = 1,n_abs_modes
            !Find this interaction in the ENDF tape (MF=3, MT=abs_modes(j))
            Rewind(cs_unit)  !return to the start of the file
            Do
                Read(cs_unit,'(A73,I1,I3,I5)') trash_c,MF,MT,line_num
                If (MF.EQ.3 .AND. MT.EQ.abs_modes(j)) Then
                    Backspace(cs_unit)  !go back one line
                    Exit
                End If
            End Do
            !the next read statement on cs_unit will read the first line of MF=3, MT=abs_modes(j)
            !UNDONE
            Call Read_CS_sect(cs_unit,Q_scratch,An_scratch,E_scratch,CS_scratch,Interp_scratch,n_p,n_r)
            If (Trim_CS_for_E(n_p,E_scratch,CS_scratch,n_r,Interp_scratch,E_min,E_max)) Then
                Call Map_and_Store_CS(CS%n_E_uni,CS%E_uni,n_p,E_scratch,CS_scratch,n_r,Interp_Scratch,CS%abs_cs(i)%sig(j),CS%abs_cs(i)%thresh(j))
            End If
            Deallocate(E_scratch,CS_scratch,Interp_scratch)
        End Do
        If (elastic_only) Then
            CS%lev_cs(i)%n_lev = 0
            Allocate(CS%lev_cs(i)%Q(0:0))
            Allocate(CS%lev_cs(i)%thresh(0:0))
            Allocate(CS%lev_cs(i)%sig(0:0))
            If (aniso_dist) Allocate(CS%lev_cs(i)%da(0:0))
        Else
            CS%lev_cs(i)%n_lev = n_inel_lev
            Allocate(CS%lev_cs(i)%Q(0:n_inel_lev))
            Allocate(CS%lev_cs(i)%thresh(0:n_inel_lev))
            Allocate(CS%lev_cs(i)%sig(0:n_inel_lev))
            If (aniso_dist) Allocate(CS%lev_cs(i)%da(0:n_inel_lev))
        End If
        Rewind(cs_unit)  !return to the start of the file
        Read(cs_unit,*)
        Read(cs_unit,'(A22,I11)') trash_c,LRP
        If (LRP .EQ. 1) Then  !resonance parameters are included in the ENDF tape
            CS%has_res_cs(i) = .TRUE.
            !Find this interaction in the ENDF tape (MF=2, MT=151)
            Rewind(cs_unit)  !return to the start of the file
            Do
                Read(cs_unit,'(A73,I1,I3,I5)') trash_c,MF,MT,line_num
                If (MF.EQ.2 .AND. MT.EQ.151) Then
                    Backspace(cs_unit)  !go back one line
                    Exit
                End If
            End Do
            !the next read statement on cs_unit will read the first line of MF=2, MT=151
            !UNDONE
            Call Read_res_sect(cs_unit,)
            !UNDONE
        Else  !no resonance parameters
            CS%has_res_cs(i) = .FALSE.
        End If
        !Find this interaction in the ENDF tape (MF=3, MT=2)
        Rewind(cs_unit)  !return to the start of the file
        Do
            Read(cs_unit,'(A73,I1,I3,I5)') trash_c,MF,MT,line_num
            If (MF.EQ.3 .AND. MT.EQ.2) Then
                Backspace(cs_unit)  !go back one line
                Exit
            End If
        End Do
        !the next read statement on cs_unit will read the first line of MF=3, MT=2
        !UNDONE
        Call Read_CS_sect(cs_unit,Q_scratch,CS%An(i),E_scratch,CS_scratch,Interp_scratch,n_p,n_r)
        If (Trim_CS_for_E(n_p,E_scratch,CS_scratch,n_r,Interp_scratch,E_min,E_max)) Then
            Call Map_and_Store_CS(CS%n_E_uni,CS%E_uni,n_p,E_scratch,CS_scratch,n_r,Interp_Scratch,CS%abs_cs(i)%sig(j),CS%abs_cs(i)%thresh(j))
        End If
        Deallocate(E_scratch,CS_scratch,Interp_scratch)
        If (aniso_dist) Then  !need elastic ang dist file
            !Find this interaction in the ENDF tape (MF=4, MT=2)
            Rewind(cs_unit)  !return to the start of the file
            Do
                Read(cs_unit,'(A73,I1,I3,I5)') trash_c,MF,MT,line_num
                If (MF.EQ.4 .AND. MT.EQ.2) Then
                    Backspace(cs_unit)  !go back one line
                    Exit
                End If
            End Do
            !the next read statement on cs_unit will read the first line of MF=4, MT=2
            !UNDONE
            Call Read_Ang_Dist_sect(cs_unit,E_scratch,Ang_dist_scratch,n_p,ltt)
            If (Trim_AD_for_E(n_p,E_scratch,Ang_dist_scratch,E_min,E_max)) Then
                Call Map_and_Store_AD(CS%n_E_uni,CS%E_uni,n_p,E_scratch,Ang_dist_scratch,CS%lev_cs(i)%da(0))
            End If
            Deallocate(E_scratch,Ang_dist_scratch)
            If (ltt .EQ. 1) Then
                If (MaxVal(CS%lev_cs(i)%da(0)%da(:)%n_a) .GT. CS%n_a_max) CS%n_a_max = MaxVal(CS%lev_cs(i)%da(0)%da(:)%n_a)
            Else If (ltt .EQ. 2) Then
                If (MaxVal(CS%lev_cs(i)%da(0)%da(:)%n_a) .GT. CS%n_a_tab_max) CS%n_a_tab_max = MaxVal(CS%lev_cs(i)%da(0)%da(:)%n_a)
            Else If (ltt .EQ. 3) Then
                Do k = 1,CS%lev_cs(i)%da(0)%n_da
                    If (CS%lev_cs(i)%da(0)%da(k)%is_legendre) Then
                        If (CS%lev_cs(i)%da(0)%da(k)%n_a .GT. CS%n_a_max) CS%n_a_max = CS%lev_cs(i)%da(0)%da(k)%n_a
                    Else If (CS%lev_cs(i)%da(0)%da(k)%is_tab) Then
                        If (CS%lev_cs(i)%da(0)%da(k)%n_a .GT. CS%n_a_tab_max) CS%n_a_tab_max = CS%lev_cs(i)%da(0)%da(k)%n_a
                    End If
                End Do
            End If
        End If
        If (.NOT. elastic_only) Then  !need inelastic level cross section files
            CS%lev_cs(i)%Q = -1._dp
            CS%lev_cs(i)%Q(0) = 0._dp
            CS%lev_cs(i)%thresh = -1
            CS%lev_cs(i)%thresh(0) = 0
            Do j = 1,n_inel_lev
                !Find this interaction in the ENDF tape (MF=3, MT=50+j)
                Rewind(cs_unit)  !return to the start of the file
                Do
                    Read(cs_unit,'(A73,I1,I3,I5)') trash_c,MF,MT,line_num
                    If (MF.EQ.3 .AND. MT.EQ.50+j) Then
                        Backspace(cs_unit)  !go back one line
                        Exit
                    End If
                End Do
                !the next read statement on cs_unit will read the first line of MF=3, MT=50+j
                !UNDONE
                Call Read_CS_sect(cs_unit,Q_scratch,An_scratch,E_scratch,CS_scratch,Interp_scratch,n_p,n_r)
                If (Trim_CS_for_E(n_p,E_scratch,CS_scratch,n_r,Interp_scratch,E_min,E_max)) Then
                    CS%lev_cs(i)%Q(j) = Q_scratch
                    Call Map_and_Store_CS(CS%n_E_uni,CS%E_uni,n_p,E_scratch,CS_scratch,n_r,Interp_Scratch,CS%lev_cs(i)%sig(j),CS%lev_cs(i)%thresh(j))
                End If
                Deallocate(E_scratch,CS_scratch,Interp_scratch)
                If (aniso_dist) Then  !need inelastic level ang dist files
                    !Find this interaction in the ENDF tape (MF=4, MT=50+j)
                    Rewind(cs_unit)  !return to the start of the file
                    Do
                        Read(cs_unit,'(A73,I1,I3,I5)') trash_c,MF,MT,line_num
                        If (MF.EQ.4 .AND. MT.EQ.50+j) Then
                            Backspace(cs_unit)  !go back one line
                            Exit
                        End If
                    End Do
                    !the next read statement on cs_unit will read the first line of MF=4, MT=50+j
                    !UNDONE
                    Call Read_Ang_Dist_sect(cs_unit,E_scratch,Ang_dist_scratch,n_p,ltt)
                    If (Trim_AD_for_E(n_p,E_scratch,Ang_dist_scratch,E_min,E_max)) Then
                        Call Map_and_Store_AD(CS%n_E_uni,CS%E_uni,n_p,E_scratch,Ang_dist_scratch,CS%lev_cs(i)%da(j),CS%lev_cs(i)%thresh(j))
                    End If
                    Deallocate(E_scratch,Ang_dist_scratch)
                    If (ltt .EQ. 1) Then
                        If (MaxVal(CS%lev_cs(i)%da(j)%da(:)%n_a) .GT. CS%n_a_max) CS%n_a_max = MaxVal(CS%lev_cs(i)%da(j)%da(:)%n_a)
                    Else If (ltt .EQ. 2) Then
                        If (MaxVal(CS%lev_cs(i)%da(j)%da(:)%n_a) .GT. CS%n_a_tab_max) CS%n_a_tab_max = MaxVal(CS%lev_cs(i)%da(j)%da(:)%n_a)
                    Else If (ltt .EQ. 3) Then
                        Do k = 1,CS%lev_cs(i)%da(j)%n_da
                            If (CS%lev_cs(i)%da(j)%da(k)%is_legendre) Then
                                If (CS%lev_cs(i)%da(j)%da(k)%n_a .GT. CS%n_a_max) CS%n_a_max = CS%lev_cs(i)%da(j)%da(k)%n_a
                            Else If (CS%lev_cs(i)%da(j)%da(k)%is_tab) Then
                                If (CS%lev_cs(i)%da(j)%da(k)%n_a .GT. CS%n_a_tab_max) CS%n_a_tab_max = CS%lev_cs(i)%da(j)%da(k)%n_a
                            End If
                        End Do
                    End If
                End If
            End Do
        End If
    End Do
    CS%Mn = neutron_mass * Sum(CS%An) / CS%n_iso
End Function Setup_Cross_Sections

Subroutine Read_CS_file(CS_file_name,Q,An,E_list,CS_list,Int_list,n_p,n_r,just_n)
    Use Kinds, Only: dp
    Use FileIO_Utilities, Only: Output_Message
    Implicit None
    Character(*), Intent(In) :: CS_file_name
    Real(dp), Intent(Out) :: Q
    Real(dp), Intent(Out) :: An
    Real(dp), Allocatable, Intent(Out) :: E_list(:)
    Real(dp), Allocatable, Intent(Out) :: CS_list(:)
    Integer, Allocatable, Intent(Out) :: Int_list(:,:)
    Integer, Intent(Out) :: n_p,n_r
    Logical, Intent(In), Optional :: just_n
    Integer :: i
    Real(dp) :: trash
    Integer :: cs_unit, stat

    Open(NEWUNIT = cs_unit , FILE = CS_file_name , ACTION = 'READ' , IOSTAT = stat)
    If (stat .NE. 0) Call Output_Message('ERROR:  Cross_Sections: Read_CS_file:  File open error, '//CS_file_name//', IOSTAT=',stat,kill=.TRUE.)
    !the second entry of the first line is the mass of the target in neutron masses
    Read(cs_unit,'(2E11.6E1)') trash, An
    !the second entry of the next line is the Q-value for the reaction, the fifth and sixth entries are the number of interpolation ranges and energy levels in the file
    Read(cs_unit,'(4E11.6E1,2I11)') trash, Q, trash, trash, n_r, n_p
    If (n_r .GT. 3) Call Output_Message('ERROR:  Cross_Sections: Read_CS_file:  Number of interpolation ranges greater than 3: '//CS_file_name//', n_r=',n_r,kill=.TRUE.)
    If (Present(just_n)) Then
        If (just_n) Then
            Close(cs_unit)
            Return
        End If
    End If
    !The next line has interpolation methods and ranges, read and fill interpolation list
    Allocate(Int_list(1:n_r,1:2))
    Int_list = 0
    If (n_r .EQ. 1) Then
        Read(cs_unit,'(2I11)') Int_list(1,1), Int_list(1,2)
    Else If (n_r .EQ. 2) Then
        Read(cs_unit,'(4I11)') Int_list(1,1), Int_list(1,2), Int_list(2,1), Int_list(2,2)
    Else !n_r=3
        Read(cs_unit,'(6I11)') Int_list(1,1), Int_list(1,2), Int_list(2,1), Int_list(2,2), Int_list(3,1), Int_list(3,2)
    End If
    If (Any(Int_list(:,2).GT.5) .OR. Any(Int_list(:,2).LT.1)) Call Output_Message('ERROR:  Cross_Sections: Read_CS_file:  Unknown interpolation scheme: '//CS_file_name,kill=.TRUE.)
    !Allocate lists
    Allocate(E_list(1:n_p))
    Allocate(CS_list(1:n_p))
    !Read the rest of the file
    Do i = 1,n_p,3
        !Each line in the file has 3 pairs of (eV,barns)
        If (n_p-i .GT. 1) Then
            Read(cs_unit,'(6E11.6E1)') E_list(i), CS_list(i), E_list(i+1), CS_list(i+1), E_list(i+2), CS_list(i+2)
        Else If (n_p-i .EQ. 1) Then
            Read(cs_unit,'(4E11.6E1)') E_list(i), CS_list(i), E_list(i+1), CS_list(i+1)
        Else
            Read(cs_unit,'(2E11.6E1)') E_list(i), CS_list(i)
        End If
    End Do
    Close(cs_unit)
    !Convert E from eV to keV
    Q = Abs(Q) / 1000._dp  !also store Q as positive value, negativity of Q for inelastic scatter is handled by energy book-keeeping formulations during transport
    E_list = E_list / 1000._dp
End Subroutine Read_CS_file

Subroutine Read_Ang_Dist_file(Coeff_file_name,E_list,da_list,n_p,LTT,just_n)
    Use Kinds, Only: dp
    Use FileIO_Utilities, Only: Output_Message
    Implicit None
    Character(*), Intent(In) :: Coeff_file_name
    Real(dp), Allocatable, Intent(Out) :: E_list(:)
    Type(da_List_type), Allocatable, Intent(Out) :: da_list(:)
    Integer, Intent(Out):: n_p,LTT
    Logical, Intent(In), Optional :: just_n
    Integer :: i,j
    Real(dp) :: trash
    Integer :: n_a_scratch
    Integer :: coeff_unit,stat
    Integer :: LCT
    Integer :: n_p_add
    Logical :: new_line

    Open(NEWUNIT = coeff_unit , FILE = Coeff_file_name , ACTION = 'READ' , IOSTAT = stat)
    If (stat .NE. 0) Call Output_Message('ERROR:  Cross_Sections: Read_Ang_Dist_file:  File open error, '//Coeff_file_name//', IOSTAT=',stat,kill=.TRUE.)
    !check the LTT value on the first line to ensure Legendre Coeffs format
    Read(coeff_unit,'(3E11.6E1,I11)') trash, trash, trash, LTT
    Read(coeff_unit,'(3E11.6E1,I11)') trash, trash, trash, LCT
    !values must be in CM frame, check and throw error if they are not
    If (LCT .NE. 2) Call Output_Message('ERROR:  Cross_Sections: Read_Ang_Dist_file:  Incorrectly formatted file, '//Coeff_file_name//', LCT=',LCT,kill=.TRUE.)
    !Skip the next line
    Read (coeff_unit,*)
    !the next line begins with the number of energy levels in the file
    Read (coeff_unit,'(I11)') n_p
    If (LTT .EQ. 3) Then !need to add more energies from later in the file
        !advance in the file to the end of the Legendre section
        Do i = 1,n_p
            !The first line in each energy contains the energy in eV in the second position and the number of Legendre coefficents in the 5th position
            Read(coeff_unit,'(4E11.6E1,I11)') trash, trash, trash, trash, n_a_scratch
            Do j = 1,n_a_scratch
                If (mod(j,6).EQ.0 .OR. j.EQ.n_a_scratch) Then !this is the last entry on a line, read and advance
                    Read(coeff_unit,'(E11.6E1)', ADVANCE = 'YES') trash
                    new_line = .TRUE.
                Else !read the entry without advancing
                    Read(coeff_unit,'(E11.6E1)', ADVANCE = 'NO') trash
                    new_line = .FALSE.
                End If
            End Do
        End Do
        !advance one or two lines based on where Legendre Coeff reading finished
        If (.NOT. new_line) Read(coeff_unit,*)
        Read(coeff_unit,*)
        !first entry of next line is number of additional energy points
        Read (coeff_unit,'(I11)') n_p_add
    Else
        n_p_add = 0
    End If
    If (Present(just_n)) Then
        If (just_n) Then  !this is the number of energies we need
            Close(coeff_unit)
            n_p = n_p + n_p_add
            Return
        End If
    End If
    !Allocate the data array to the number of provided energy levels
    Allocate(E_list(1:n_p+n_p_add))
    Allocate(da_list(1:n_p+n_p_add))
    If (LTT .EQ. 1) Then  !da is lists of legendre coeffs
        da_list%is_legendre = .TRUE.
        da_list%is_tab = .FALSE.
        Do i = 1,n_p
            !The first line in each energy contains the energy in eV in the second position and the number of Legendre coefficents in the 5th position
            Read(coeff_unit,'(4E11.6E1,I11)') trash, E_list(i), trash, trash, da_list(i)%n_a
            Allocate(da_list(i)%a(0:da_list(i)%n_a))
            da_list(i)%a(0) = 0.5_dp
            Do j = 1,da_list(i)%n_a
                If (mod(j,6).EQ.0 .OR. j.EQ.da_list(i)%n_a) Then !this is the last entry on a line, read and advance
                    Read(coeff_unit,'(E11.6E1)', ADVANCE = 'YES') da_list(i)%a(j)
                Else !read the entry without advancing
                    Read(coeff_unit,'(E11.6E1)', ADVANCE = 'NO') da_list(i)%a(j)
                End If
                !multiply the coeff by (2k+1)/2
                da_list(i)%a(j) = 0.5_dp * da_list(i)%a(j) * Real(2*j + 1,dp)
            End Do
        End Do
    Else If (LTT .EQ. 2) Then  !da is tabulated
        da_list%is_legendre = .FALSE.
        da_list%is_tab = .TRUE.
        Do i = 1,n_p
            !The first line in each energy contains the energy in eV in the second position
            Read(coeff_unit,'(2E11.6E1)') trash, E_list(i)
            !the next line contains the number of tabulation points int he first position
            Read(coeff_unit,'(I11)') da_list(i)%n_a
            Allocate(da_list(i)%a(0:da_list(i)%n_a))
            da_list(i)%a = 0._dp
            Allocate(da_list(i)%ua(1:da_list(i)%n_a,1:2))
            da_list(i)%ua = 0._dp
            Do j = 1,da_list(i)%n_a
                If (mod(j,3).EQ.0 .OR. j.EQ.da_list(i)%n_a) Then !this is the last entry on a line, read and advance
                    Read(coeff_unit,'(2E11.6E1)', ADVANCE = 'YES') da_list(i)%ua(j,1), da_list(i)%ua(j,2)
                Else !read the entry without advancing
                    Read(coeff_unit,'(2E11.6E1)', ADVANCE = 'NO') da_list(i)%ua(j,1), da_list(i)%ua(j,2)
                End If
            End Do
        End Do
    Else If (LTT .EQ. 3) Then  !da is tabulated for high energies but legendre for low energies
        !Read in low energy Legendre points
        Do i = 1,n_p
            da_list(i)%is_legendre = .TRUE.
            da_list(i)%is_tab = .FALSE.
            !The first line in each energy contains the energy in eV in the second position and the number of Legendre coefficents in the 5th position
            Read(coeff_unit,'(4E11.6E1,I11)') trash, E_list(i), trash, trash, da_list(i)%n_a
            Allocate(da_list(i)%a(0:da_list(i)%n_a))
            da_list(i)%a(0) = 0.5_dp
            Do j = 1,da_list(i)%n_a
                If (mod(j,6).EQ.0 .OR. j.EQ.da_list(i)%n_a) Then !this is the last entry on a line, read and advance
                    Read(coeff_unit,'(E11.6E1)', ADVANCE = 'YES') da_list(i)%a(j)
                    new_line = .TRUE.
                Else !read the entry without advancing
                    Read(coeff_unit,'(E11.6E1)', ADVANCE = 'NO') da_list(i)%a(j)
                    new_line = .FALSE.
                End If
                !multiply the coeff by (2k+1)/2
                da_list(i)%a(j) = 0.5_dp * da_list(i)%a(j) * Real(2*j + 1,dp)
            End Do
        End Do
        !advance one or two lines based on where Legendre Coeff reading finished
        If (.NOT. new_line) Read(coeff_unit,*)
        Read(coeff_unit,*)
        !first entry of next line is number of additional energy points
        Read (coeff_unit,'(I11)') n_p_add
        !Read in high energy tabulated cosine points
        Do i = n_p+1,n_p+n_p_add
            da_list(i)%is_legendre = .FALSE.
            da_list(i)%is_tab = .TRUE.
            !The first line in each energy contains the energy in eV in the second position
            Read(coeff_unit,'(2E11.6E1)') trash, E_list(i)
            !the next line contains the number of tabulation points in the first position
            Read(coeff_unit,'(I11)') da_list(i)%n_a
            Allocate(da_list(i)%a(0:da_list(i)%n_a))
            da_list(i)%a = 0._dp
            Allocate(da_list(i)%ua(1:da_list(i)%n_a,1:2))
            da_list(i)%ua = 0._dp
            Do j = 1,da_list(i)%n_a
                If (mod(j,3).EQ.0 .OR. j.EQ.da_list(i)%n_a) Then !this is the last entry on a line, read and advance
                    Read(coeff_unit,'(2E11.6E1)', ADVANCE = 'YES') da_list(i)%ua(j,1), da_list(i)%ua(j,2)
                Else !read the entry without advancing
                    Read(coeff_unit,'(2E11.6E1)', ADVANCE = 'NO') da_list(i)%ua(j,1), da_list(i)%ua(j,2)
                End If
            End Do
        End Do
        !update total n_p
        n_p = n_p + n_p_add
    End If
    Close(coeff_unit)
    !Convert E from eV to keV
    E_list = E_list / 1000._dp
End Subroutine Read_Ang_Dist_file

Subroutine Read_Res_File(Res_file_name,res_List)
    Use Kinds, Only: dp
    Use FileIO_Utilities, Only: Output_Message
    Use Sorting, Only: Union_Sort
    Implicit None
    Character(*), Intent(In) :: Res_file_name
    Type(res_sig_Type_new) :: res_list

    Integer :: i
    Real(dp) :: trash
    Integer :: res_unit,stat
    Integer :: LRF,NRO,NAPS,SPI,AP
    Integer :: nL,nR,nJ
    Real(dp) :: awri
    Real(dp), Allocatable :: res_scratch(:,:),Js(:)

    Open(NEWUNIT = res_unit , FILE = Res_file_name , ACTION = 'READ' , IOSTAT = stat)
    If (stat .NE. 0) Call Output_Message('ERROR:  Cross_Sections: Read_Res_file:  File open error, '//Res_file_name//', IOSTAT=',stat,kill=.TRUE.)
    !check the LRF value on the fourth line to ensure Reich-Moore format
    !check the NRO and NAPS value on the fourth line to ensure energy independent scattering radius and use of SPI instead of channel radius
    Do i = 1,3
        Read(res_unit,*)
    End Do
    Read(res_unit,'(3E11.6E1,3I11)') trash, trash, trash, LRF, NRO, NAPS
    If (LRF .NE. 3) Call Output_Message('ERROR:  Cross_Sections: Read_Res_file:  Incorrectly formatted file, '//Res_file_name//', LRF=',LRF,kill=.TRUE.)
    If (NRO .NE. 0) Call Output_Message('ERROR:  Cross_Sections: Read_Res_file:  Incorrectly formatted file, '//Res_file_name//', NRO=',NRO,kill=.TRUE.)
    If (NAPS .NE. 1) Call Output_Message('ERROR:  Cross_Sections: Read_Res_file:  Incorrectly formatted file, '//Res_file_name//', NAPS=',NAPS,kill=.TRUE.)
    !read AP and n_L from the next line
    Read(res_unit,'(4E11.6E1,I11)') SPI, AP, trash, trash, nL
    !allocate levels
    res_List%n_L = nL
    Allocate(res_List%AP(1:nL))
    res_List%AP = AP  !default value
    Allocate(res_List%L(1:nL))
    Do l = 1,nL
        !read AWRI, APL and number of resonances in first layer from next line
        Read(res_unit,'(5E11.6E1,I11)') AWRI, AP, trash, trash, trash, nR
        If (AP .NE. 0._dp) res_List%AP(l) = AP  !level specific scattering radius
        !Allocate a scratch array
        If (l .GT. 1) Deallocate(res_scratch)
        Allocate(res_scratch(1:6,1:nR))
        res_scratch = -1._dp
        Allocate(Js(1:nR))
        Js = -1._dp
        !read resonance parameters for this level
        Do r = 1,nR
            Read(res_unit,'(6E11.6E1)') res_scratch(1,r),res_scratch(2,r),res_scratch(3,r),res_scratch(4,r),res_scratch(5,r),res_scratch(6,r)
        End Do
        !check for fission resonances
        If (Any(res_scratch(5:6,:).NE.0._dp)) Call Output_Message('ERROR:  Cross_Sections: Read_Res_file:  Dude, a fissionable atmospere is just ridiculous. '//Res_file_name,kill=.TRUE.)
        !count unique spin values
        Js = res_scratch(2,:)
        Call Union_Sort(Js,nJ)
        !Allocate and fill spin groups for this level
        res_list%L(l)%n_J = nJ
        !total angular momentum
        Allocate(res_list%L(l)%gJ(1:nJ))
        res_list%L(l)%gJ = 0.5_dp * (2._dp * Js(1:nJ) + 1._dp) / (2._dp * SPI + 1._dp)
        !second channel spin flags
        Allocate(res_list%L(l)%dJ(1:nJ))
        res_list%L(l)%dJ = .FALSE.
        s_min = Abs(SPI - 0.5_dp)
        s_max = SPI + 0.5_dp
        nS = 1
        Do
            If (s_min + Real(nS,dp)*0.5_dp .GT. s_max) Exit
            nS = nS + 1
        End Do
        If (nS .GT. 1) Then !second channel flags must be set
            !UNDONE
            !UNDONE
            !UNDONE
        End If
        !spin-grouped resonance lists
        Allocate(res_list%L(l)%J(1:nJ))
        Do J = 1,nJ
            !count number of resonance parameters for this spin
            nRj = 0
            Do i = 1,nR
                If (res_scratch(2,i) .EQ. Js(J)) nRj = nRj + 1
            End Do
            res_list%L(l)%J(J)%n_r = nRj
            Allocate(res_list%L(l)%J(J)%ErGnGg(1:3,1:nRj))
            k = 1
            Do i = 1,nR
                If (res_scratch(2,i) .EQ. Js(J)) Then
                    res_list%L(l)%J(J)%ErGnGg(1,k) = res_scratch(1,i) / 1000._dp  !convert to keV
                    res_list%L(l)%J(J)%ErGnGg(2,k) = res_scratch(3,i)
                    res_list%L(l)%J(J)%ErGnGg(3,k) = res_scratch(4,i)
                    res_scratch(4,i) = 0._dp
                    If (k .GE. nRj) Exit
                    k = k + 1
                End If
            End Do
        End Do
    End Do

!UNDONE
!UNDONE
!UNDONE
!UNDONE
    !compute k0
    res_List%k0 = 2.196771E-3_dp*(AWRI/(AWRI+1._dp))
    Close(res_unit)
End Subroutine Read_Res_File

Function Trim_CS_for_E(n_p,E_list,CS_list,n_r,Int_list,E_min,E_max) Result(bingo)
    Use Kinds, Only: dp
    Use Interpolation, Only: Linear_Interp
    Implicit None
    Logical :: bingo
    Integer, Intent(InOut) :: n_p
    Real(dp), Allocatable, Intent(InOut) :: E_list(:)
    Real(dp), Allocatable, Intent(InOut) :: CS_list(:)
    Integer, Intent(InOut) :: n_r
    Integer, Allocatable, Intent(InOut) :: Int_list(:,:)
    Real(dp), Intent(In) :: E_min,E_max
    Real(dp), Allocatable :: E_swap(:)
    Real(dp), Allocatable :: CS_swap(:)
    Integer, Allocatable :: Int_swap(:,:)
    Integer :: i,j

    bingo = .TRUE.
    If (E_list(1).GE.E_min .AND. E_list(n_p).LE.E_max) Return  !no trimming required
    If (E_list(1).GT.E_max .OR. E_list(n_p).LT.E_min) Then !there are no values in range
        bingo = .FALSE.
        n_p = 0
        Return
    End If
    !otherwise, need to find the range of E_list to keep
    i = 1
    Do
        If (E_list(i) .LT. E_min) Then
            i = i + 1
        Else
            Exit
        End If
    End Do
    j = n_p
    Do
        If (E_list(j) .GT. E_max) Then
            j = j - 1
        Else
            Exit
        End If
    End Do
    !compute new number of elements
    n_p = j - i + 1
    !resize E_list
    Allocate(E_swap(1:n_p))
    E_swap = E_list(i:j)
    Deallocate(E_list)
    Allocate(E_list(1:n_p))
    E_list = E_swap
    Deallocate(E_swap)
    !resize CS_list
    Allocate(CS_swap(1:n_p))
    CS_swap = CS_list(i:j)
    Deallocate(CS_list)
    Allocate(CS_list(1:n_p))
    CS_list = CS_swap
    Deallocate(CS_swap)
    !adjust Int_list for new indexes
    Int_list(:,1) = Int_list(:,1) - (i-1)
    If (Any(Int_list(:,1) .LT. 2)) Then !interpolation ranges at the beginning of the list need to be removed
        Do i = 2,n_r
            If (Int_list(i,1) .GE. 2) Then !remove the ranges before this one
                Allocate(Int_swap(1:n_r-(i-1),1:2))
                Int_swap = Int_list(i:n_r,:)
                Deallocate(Int_list)
                n_r = n_r - (i-1)
                Allocate(Int_list(1:n_r,1:2))
                Int_list = Int_swap
                Deallocate(Int_swap)
                Exit
            End If
        End Do
    End If
    If (Any(Int_list(:,1) .GT. n_p)) Then !interpolation ranges at the end of the list need to be removed
        Do i = 1,n_r-1
            If (Int_list(i,1) .GE. n_p) Then !remove the ranges after this one
                Allocate(Int_swap(1:i,1:2))
                Int_swap = Int_list(1:i,:)
                Deallocate(Int_list)
                n_r = i
                Allocate(Int_list(1:n_r,1:2))
                Int_list = Int_swap
                Int_list(n_r,1) = n_p
                Deallocate(Int_swap)
                Exit
            End If
        End Do
    End If
End Function Trim_CS_for_E

Function Trim_AD_for_E(n_p,E_list,AD_list,E_min,E_max) Result(bingo)
    Use Kinds, Only: dp
    Implicit None
    Logical :: bingo
    Integer, Intent(InOut) :: n_p
    Real(dp), Allocatable, Intent(InOut) :: E_list(:)
    Type(da_List_Type), Allocatable, Intent(InOut) :: AD_list(:)
    Real(dp), Intent(In) :: E_min,E_max
    Real(dp), Allocatable :: E_swap(:)
    Type(da_List_Type), Allocatable :: AD_swap(:)
    Integer :: i,j,m,k

    bingo = .TRUE.
    If (E_list(1).GE.E_min .AND. E_list(n_p).LE.E_max) Return  !no trimming required
    If (E_list(1).GT.E_max .OR. E_list(n_p).LT.E_min) Then !there are no values in range
        bingo = .FALSE.
        n_p = 0
        Return
    End If
    !otherwise, need to find the range of E_list to keep
    i = 1
    Do
        If (E_list(i) .LT. E_min) Then
            i = i + 1
        Else
            Exit
        End If
    End Do
    j = n_p
    Do
        If (E_list(j) .GT. E_max) Then
            j = j - 1
        Else
            Exit
        End If
    End Do
    !compute new number of elements
    n_p = j - i + 1
    !resize E_list
    Allocate(E_swap(1:n_p))
    E_swap = E_list(i:j)
    Deallocate(E_list)
    Allocate(E_list(1:n_p))
    E_list = E_swap
    Deallocate(E_swap)
    !resize CS_list
    Allocate(AD_swap(1:n_p))
    m = 1
    Do k = i,j
        AD_swap(m)%n_a = AD_list(k)%n_a
        Allocate(AD_swap(m)%a(0:AD_swap(m)%n_a))
        AD_swap(m)%a = AD_list(k)%a
        m = m + 1
    End Do
    Deallocate(AD_list)
    Allocate(AD_list(1:n_p))
    Do k = 1,n_p
        AD_list(k)%n_a = AD_swap(k)%n_a
        Allocate(AD_list(k)%a(0:AD_list(k)%n_a))
        AD_list(k)%a = AD_swap(k)%a
    End Do
End Function Trim_AD_for_E

Subroutine Map_and_Store_CS(n_E_uni,E_uni,n_p,E_list,CS_list,n_r,Int_list,cs,i_thresh)
    Use Kinds, Only: dp
    Implicit None
    Integer, Intent(In) :: n_E_uni
    Real(dp), Intent(In) :: E_uni(1:n_E_uni)
    Integer, Intent(In) :: n_p
    Real(dp), Intent(In) :: E_list(1:n_p)
    Real(dp), Intent(In) :: CS_list(1:n_p)
    Integer, Intent(In) :: n_r
    Integer, Intent(In) :: Int_list(1:n_r,1:2)
    Type(sig_Type), Intent(Out) :: cs
    Integer, Optional, Intent(Out) :: i_thresh
    Integer :: i,j,t

    !Determine index of threshold energy
    t = 1
    If (Present(i_thresh)) Then  !determine index of the threshold energy, index of map will start at the threshold energy
        If (CS_list(1) .EQ. 0._dp) Then
            Do i = 1,n_E_uni
                If (E_list(1) .GT. E_uni(i)) Cycle
                i_thresh = i
                Exit
            End Do
            t = i_thresh
        Else
            i_thresh = 0
        End If
    End If
    !Create map
    Allocate(cs%E_map(t:n_E_uni))
    cs%E_map = -1
    j = 2
    Do i = t,n_E_uni
        If (E_uni(i) .GT. E_list(j)) j = j + 1
        If (j .GT. n_p) Then
            cs%E_map(i:n_E_uni) = j
            Exit
        End If
        cs%E_map(i) = j
        If (E_uni(i) .EQ. E_list(j)) j = j + 1
    End Do
    !Create key
    Allocate(cs%E_key(1:n_p))
    cs%E_key = -1
    j = 1
    Do i = 1,n_p
        Do j = 1,n_E_uni
            If (E_uni(j) .EQ. E_list(i)) Then
                cs%E_key(i) = j
            End If
        End Do
    End Do
    !Store cross sections
    cs%n_sig = n_p
    Allocate(cs%sig(1:n_p))
    cs%sig = CS_list
    !store interpolation ranges
    cs%n_interp_r = n_r
    Allocate(cs%interp(1:n_r,1:2))
    cs%interp = Int_list
    !Store natural logarithms of cross sections if needed based on interpolation schemes
    If (Any(Int_list(:,2).EQ.4) .OR. Any(Int_list(:,2).EQ.5)) Then
        Allocate(cs%lnsig(1:n_p))
        cs%lnsig = -Huge(cs%lnsig)
        Where (CS_list .GT. 0._dp) cs%lnsig = Log(CS_list)
    End If
End Subroutine Map_and_Store_CS

Subroutine Map_and_Store_AD(n_E_uni,E_uni,n_p,E_list,AD_list,ad,i_thresh)
    Use Kinds, Only: dp
    Implicit None
    Integer, Intent(In) :: n_E_uni
    Real(dp), Intent(In) :: E_uni(1:n_E_uni)
    Integer, Intent(In) :: n_p
    Real(dp), Intent(In) :: E_list(1:n_p)
    Type(da_List_Type), Intent(In) :: AD_list(1:n_p)
    Type(da_Type), Intent(Out) :: ad
    Integer, Optional, Intent(In) :: i_thresh
    Integer :: i,j,t

    !Determine index of threshold energy
    t = 1
    If (Present(i_thresh)) t = i_thresh
    !Create map
    Allocate(ad%E_map(t:n_E_uni))
    ad%E_map = -1
    j = 2
    Do i = t,n_E_uni
        If (E_uni(i) .GT. E_list(j)) j = j + 1
        If (j .GT. n_p) Then
            ad%E_map(i:n_E_uni) = j
            Exit
        End If
        ad%E_map(i) = j
        If (E_uni(i) .EQ. E_list(j)) j = j + 1
    End Do
    !Create key
    Allocate(ad%E_key(1:n_p))
    ad%E_key = -1
    j = 1
    Do i = 1,n_p
        Do j = 1,n_E_uni
            If (E_uni(j) .EQ. E_list(i)) Then
                ad%E_key(i) = j
            End If
        End Do
    End Do
    !Store cross sections
    ad%n_da = n_p
    Allocate(ad%da(1:n_p))
    ad%da = AD_list
End Subroutine Map_and_Store_AD

Subroutine sig_R_iso(r,E,sT,Ss)
    Use Kinds, Only(dp)
    Use Global, Only: TwoPi
    Implicit None
    Type(sig_res_Type_new), Intent(In) :: r
    Real(dp), Intent(In) :: E
    Real(dp), Intent(Out) :: sT
    Real(dp), Intent(Out) :: sS
    Real(dp) :: k,rho
    Integer :: l,J
    Real(dp) :: two_phi
    Complex(dp) :: Rnn,Unn
    Real(dp) :: Dj
    Complex(dp), Parameter :: half_i = CMPLX(0._dp,0.5_dp,KIND=dp)
    Complex(dp), Parameter :: cmplx1 = CMPLX(1._dp,KIND=dp)
    Complex(dp), Parameter :: cmplx2 = CMPLX(2._dp,KIND=dp)
    
    sT = 0._dp
    sS = 0._dp
    k = r%k0 * Sqrt(E)
    Do l = 1,r%n_L  !sum over levels (l)
        rho = k * r%AP(l)
        two_phi = rho
        If (l .GT. 1) Then
            Select Case (l)
                Case (2)
                    two_phi = two_phi - ATAN(rho)
                Case (3)
                    two_phi = two_phi - ATAN(3._dp * rho / (3._dp - rho**2))
                Case (4)
                    two_phi = two_phi - ATAN(rho * (15._dp - rho)**2 / (15._dp - 6._dp * rho**2))
                Case (5)
                    two_phi = two_phi - ATAN(rho * (105._dp - 10._dp * rho**2) / (105._dp - rho**2 * (45._dp + rho**2)))
            End Select
        End If
        two_phi = 2._dp * two_phi
        !TODO See if the following loop can be vectorized (for speed, or at least compactness of notation)
        Do J = 1,r%nJ(l)  !sum over spins (J)
            !compute the R-function
            Rnn = -half_i * Sum( &  !sum over r
                                & CMPLX(r%L(l)%J(j)%ErGnGg(2,:),KIND=dp) / & 
                                & (CMPLX(r%L(l)%J(j)%ErGnGg(1,:) - E,KIND=dp) - half_i*CMPLX(r%L(l)%J(j)%ErGnGg(3,:),KIND=dp)) &
                                & )
            !compute the element of the scattering matrix
            Unn = Exp(CMPLX(0._dp,two_phi,KIND=dp)) * (cmplx2/Rnn - cmplx1)
            !compute second channel contribution
            If (r%L(l)%dj(J)) Then
                Dj = 2._dp * (1._dp - Cos(two_phi))
            Else
                Dj = 0._dp
            End If
            sT = sT + r%L(l)%gj(J) * (1._dp - Real(Unn,dp) + Dj)
            sS = sS + r%L(l)%gj(J) * (Abs(cmplx1 - Unn)**2 + Dj)
        End Do
    End Do
    sT = sT * TwoPi / k**2
    sS = sS * TwoPi / k**2
End Subroutine sig_R_iso

Subroutine sig_T_A(CS,E,sT,sA,iE_get,iE_put)
    Use Kinds, Only: dp
    Implicit None
    Class(CS_Type), Intent(In) :: CS
    Real(dp), Intent(In) :: E
    Real(dp), Intent(Out) :: sT,sA
    Integer, Intent(Out), Optional :: iE_get
    Integer, Intent(In), Optional :: iE_put
    Integer :: E_index
    Real(dp) :: sS

    If (Present(iE_put)) Then
        sA = CS%sig_A(E,iE_put=iE_put)
        sS = CS%sig_S(E,iE_put=iE_put)
    Else
        sA = CS%sig_A(E,iE_get=E_index)
        sS = CS%sig_S(E,iE_put=E_index)
        If (Present(iE_get)) iE_get = E_index
    End If
    sT = sA + sS
End Subroutine sig_T_A

Function sig_T(CS,E,iE_get,iE_put)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: sig_T
    Class(CS_Type), Intent(In) :: CS
    Real(dp), Intent(In) :: E
    Integer, Intent(Out), Optional :: iE_get
    Integer, Intent(In), Optional :: iE_put
    Integer :: E_index
    Real(dp) :: sA,sS

    If (Present(iE_put)) Then
        sA = CS%sig_A(E,iE_put=iE_put)
        sS = CS%sig_S(E,iE_put=iE_put)
    Else
        sA = CS%sig_A(E,iE_get=E_index)
        sS = CS%sig_S(E,iE_put=E_index)
        If (Present(iE_get)) iE_get = E_index
    End If
    sig_T = sA + sS
End Function sig_T

Function sig_S(CS,E,iE_get,iE_put)
    Use Kinds, Only: dp
    Use Utilities, Only: Bisection_Search
    Implicit None
    Real(dp) :: sig_S
    Class(CS_Type), Intent(In) :: CS
    Real(dp), Intent(In) :: E
    Integer, Intent(Out), Optional :: iE_get
    Integer, Intent(In), Optional :: iE_put
    Integer :: E_index
    Integer :: i

    If (Present(iE_put)) Then
        E_index = iE_put
    Else
        E_index = Bisection_Search(E,CS%E_uni,CS%n_E_uni)
    End If
    sig_S = 0._dp
    Do i = 1,CS%n_iso
        sig_S = sig_S + CS%iso_Fractions(i) * sig_Composite(E,CS%n_E_uni,CS%E_uni,CS%lnE_uni,E_index,0,CS%lev_cs(i)%n_lev,CS%lev_cs(i)%thresh,CS%lev_cs(i)%sig)
        If (CS%has_res_cs(i)) Then
            If (E.GT.CS%res_cs(i)%E_range(1) .AND. E.LT.CS%res_cs(i)%E_range(2)) sig_S = sig_S + CS%iso_Fractions(i) * sig_res(E,CS%res_cs(i)%n_E,CS%res_cs(i)%E,CS%res_cs(i)%sig(:,2))
        End If
    End Do
    If (Present(iE_get)) iE_get = E_index
End Function sig_S

Function sig_S_iso(CS,iso,E,iE_get,iE_put)
    Use Kinds, Only: dp
    Use Utilities, Only: Bisection_Search
    Implicit None
    Real(dp) :: sig_S_iso
    Class(CS_Type), Intent(In) :: CS
    Integer, Intent(In) :: iso
    Real(dp), Intent(In) :: E
    Integer, Intent(Out), Optional :: iE_get
    Integer, Intent(In), Optional :: iE_put
    Integer :: E_index

    If (Present(iE_put)) Then
        E_index = iE_put
    Else
        E_index = Bisection_Search(E,CS%E_uni,CS%n_E_uni)
    End If
    sig_S_iso = sig_Composite(E,CS%n_E_uni,CS%E_uni,CS%lnE_uni,E_index,0,CS%lev_cs(iso)%n_lev,CS%lev_cs(iso)%thresh,CS%lev_cs(iso)%sig)
    If (CS%has_res_cs(iso)) Then
        If (E.GT.CS%res_cs(iso)%E_range(1) .AND. E.LT.CS%res_cs(iso)%E_range(2)) sig_S_iso = sig_S_iso + sig_res(E,CS%res_cs(iso)%n_E,CS%res_cs(iso)%E,CS%res_cs(iso)%sig(:,2))
    End If
    If (Present(iE_get)) iE_get = E_index
End Function sig_S_iso

Function sig_A(CS,E,iE_get,iE_put)
    Use Kinds, Only: dp
    Use Utilities, Only: Bisection_Search
    Implicit None
    Real(dp) :: sig_A
    Class(CS_Type), Intent(In) :: CS
    Real(dp), Intent(In) :: E
    Integer, Intent(Out), Optional :: iE_get
    Integer, Intent(In), Optional :: iE_put
    Integer :: E_index
    Integer :: i

    If (Present(iE_put)) Then
        E_index = iE_put
    Else
        E_index = Bisection_Search(E,CS%E_uni,CS%n_E_uni)
    End If
    sig_A = 0._dp
    Do i = 1,CS%n_iso
        sig_A = sig_A + CS%iso_Fractions(i) * sig_Composite(E,CS%n_E_uni,CS%E_uni,CS%lnE_uni,E_index,1,CS%abs_cs(i)%n_modes,CS%abs_cs(i)%thresh,CS%abs_cs(i)%sig)
        If (CS%has_res_cs(i)) Then
            If (E.LT.CS%res_cs(i)%E_range(2)) sig_A = sig_A + CS%iso_Fractions(i) * sig_res(E,CS%res_cs(i)%n_E,CS%res_cs(i)%E,CS%res_cs(i)%sig(:,1))
        End If
    End Do
    If (Present(iE_get)) iE_get = E_index
End Function sig_A

Function sig_Composite(E,n_E,E_list,lnE_list,E_index,n1,n2,t_list,sig_list) Result(sig)
    Use Kinds, Only: dp
    Use Interpolation, Only: Linear_Interp
    Use FileIO_Utilities, Only: Output_Message
    Implicit None
    Real(dp) :: sig
    Real(dp), Intent(In) :: E
    Integer, Intent(In) :: n_E
    Real(dp), Intent(In) :: E_list(1:n_E)
    Real(dp), Intent(In) :: lnE_list(1:n_E)
    Integer, Intent(In) :: E_index
    Integer, Intent(In) :: n1,n2
    Integer, Intent(In) :: t_list(n1:n2)
    Type(sig_Type), Intent(In) :: sig_list(n1:n2)
    Integer :: i,j
    Integer :: index1,index2
    Real(dp) :: E1,E2
    Real(dp) :: sig1,sig2
    Integer :: method
    Integer, Parameter :: Hist_interpolation = 1  !y is constant in x (constant, histogram)
    Integer, Parameter :: LinLin_interpolation = 2  !y is linear in x (linear-linear)
    Integer, Parameter :: LinLog_interpolation = 3  !y is linear in ln(x) (linear-log)
    Integer, Parameter :: LogLin_interpolation = 4  !ln(y) is linear in x (log-linear)
    Integer, Parameter :: LogLog_interpolation = 5  !ln(y) is linear in ln(x) (log-log)

    sig = 0._dp
    Do i = n1,n2
        If (E_index .LE. t_list(i)) Exit
        index1 = sig_list(i)%E_map(E_index) - 1
        index2 = sig_list(i)%E_map(E_index)
        Do j = 1,sig_list(i)%n_interp_r
            If (index2 .LE. sig_list(i)%interp(j,1)) Then
                method = sig_list(i)%interp(j,2)
                Exit
            End If
        End Do
        Select Case (method)
            Case(Hist_interpolation)
                sig = sig + sig_list(i)%sig( index1 )
            Case(LinLin_interpolation)
                E1 = E_list( sig_list(i)%E_key( index1 ) )
                E2 = E_list( sig_list(i)%E_key( index2 ) )
                sig1 = sig_list(i)%sig( index1 )
                sig2 = sig_list(i)%sig( index2 )
                sig = sig + Linear_Interp(E,E1,E2,sig1,sig2)
            Case(LinLog_interpolation)
                E1 = lnE_list( sig_list(i)%E_key( index1 ) )
                E2 = lnE_list( sig_list(i)%E_key( index2 ) )
                sig1 = sig_list(i)%sig( index1 )
                sig2 = sig_list(i)%sig( index2 )
                sig = sig + Linear_Interp(Log(E),E1,E2,sig1,sig2)
            Case(LogLin_interpolation)
                E1 = E_list( sig_list(i)%E_key( index1 ) )
                E2 = E_list( sig_list(i)%E_key( index2 ) )
                sig1 = sig_list(i)%lnsig( index1 )
                sig2 = sig_list(i)%lnsig( index2 )
                sig = sig + Exp(Linear_Interp(E,E1,E2,sig1,sig2))
            Case(LogLog_interpolation)
                E1 = lnE_list( sig_list(i)%E_key( index1 ) )
                E2 = lnE_list( sig_list(i)%E_key( index2 ) )
                sig1 = sig_list(i)%lnsig( index1 )
                sig2 = sig_list(i)%lnsig( index2 )
                sig = sig + Exp(Linear_Interp(Log(E),E1,E2,sig1,sig2))
            Case Default
                Call Output_Message('ERROR:  Cross_Sections: sig_Composite:  Undefined interpolation method',kill=.TRUE.)
        End Select
    End Do
End Function sig_Composite

Function sig_res(E,n_E,E_list,sig_list) Result(sig)
    Use Kinds, Only: dp
    Use Utilities, Only: Bisection_Search
    Use Interpolation, Only: Linear_Interp
    Implicit None
    Real(dp) :: sig
    Real(dp), Intent(In) :: E
    Integer, Intent(In) :: n_E
    Real(dp), Intent(In) :: E_list(1:n_E)
    Real(dp), Intent(In) :: sig_list(1:n_E)
    Integer :: E_i
    Real(dp) :: E1,E2
    Real(dp) :: sig1,sig2

    E_i = Bisection_Search(E,E_list,n_E)
    If (E_i .EQ. 1) Then
        E1 = E_list(1)
        E2 = E_list(2)
        sig1 = sig_list(1)
        sig2 = sig_list(2)
    Else If (E_i .GT. n_E) Then
        E1 = E_list(n_E-1)
        E2 = E_list(n_E)
        sig1 = sig_list(n_E-1)
        sig2 = sig_list(n_E)
    Else
        E1 = E_list(E_i-1)
        E2 = E_list(E_i)
        sig1 = sig_list(E_i-1)
        sig2 = sig_list(E_i)
    End If
    sig = Linear_Interp(E,E1,E2,sig1,sig2)
End Function sig_res

Subroutine Broad_sig_start(E,M,T,v,gamma,vRmin,vRmax)
    Use Kinds, Only: dp
    Use Neutron_Utilities, Only: Neutron_Speed
    Implicit None
    Real(dp), Intent(In) :: E,M,T
    Real(dp), Intent(Out) :: v,gamma,vRmin,vRmax
    Real(dp) :: vT

    v = Neutron_Speed(E)
    gamma = Sqrt(M / (2._dp * k_b * T))
    vT = 4._dp / gamma  !same cutoff as SIGMA1 algorithm used by NJOY
    If (vT .LE. v) Then
        vRmin = v - vT
    Else
        !NOTE Use small VRmin instead of zero to prevent Log(E=0) from appearing in log interpolation cases... 1.E-16 equates to an energy of about 1.E-38 keV
        vRmin = 1.E-16_dp  !0._dp
    End If
    vRmax = v + vT
End Subroutine Broad_sig_start

Subroutine sig_T_A_broad(CS,E,T,sigT,sigA)
    Use Kinds, Only: dp
    Use Global, Only: SqrtPi
    Use Utilities, Only: Bisection_Search
    Use Neutron_Utilities, Only: Neutron_Energy
    Use Neutron_Utilities, Only: Neutron_Speed
    Implicit None
    Class(CS_Type), Intent(In) :: CS
    Real(dp), Intent(In) :: E ![kev]
    Real(dp), Intent(In) :: T ![K]
    Real(dp), Intent(Out) :: sigT,sigA
    Real(dp) :: v  ![km/s] neutron velocity
    Real(dp) :: gamma  ![km/s]
    Real(dp) :: vR_min,vR_max
    Integer :: iE_min,iE_max
    Real(dp) :: sig_T_A(1:2)
    Integer :: i
    Integer, Parameter :: total = 1
    Integer, Parameter :: absor = 2

    Call Broad_sig_start(E,CS%Mn,T,v,gamma,vR_min,vR_max)
    iE_min = Bisection_Search(Neutron_Energy(vR_min),CS%E_uni,CS%n_E_uni)
    If (Neutron_Energy(vR_max) .LE. CS%E_uni(iE_min)) Then  !only a single interval spanned in energy grid
        sig_T_A = Broad_Romberg_T_A(CS,iE_min,vR_min,vR_max,gamma,v) * gamma / (SqrtPi * v**2)
        sigT = sig_T_A(total)
        sigA = sig_T_A(absor)
        Return
    End If
    iE_max = Bisection_Search(Neutron_Energy(vR_max),CS%E_uni,CS%n_E_uni) - 1
    !first interval (from E_min to first indexed E) is partial
    sig_T_A = Broad_Romberg_T_A(CS,iE_min,vR_min,Neutron_Speed(CS%E_uni(iE_min)),gamma,v)
    !middle intervals
    Do i = iE_min,iE_max-1
        sig_T_A = sig_T_A + Broad_Romberg_T_A(CS,i+1,Neutron_Speed(CS%E_uni(i)),Neutron_Speed(CS%E_uni(i+1)),gamma,v)
    End Do
    !last interval (from second to last indexed E to Emax) is partial, also apply normalization
    sig_T_A = (sig_T_A + Broad_Romberg_T_A(CS,iE_max,Neutron_Speed(CS%E_uni(iE_max)),vR_max,gamma,v)) * gamma / (SqrtPi * v**2)
    sigT = sig_T_A(total)
    sigA = sig_T_A(absor)
End Subroutine sig_T_A_broad

Function Broad_Romberg_T_A(CS,iE,vR1,vR2,gamma,v) Result(sig_T_A)
    Use Kinds, Only: dp
    Use Neutron_Utilities, Only: Neutron_Energy
    Implicit None
    Real(dp) :: sig_T_A(1:2)
    Type(CS_type), Intent(In) :: CS
    Integer, Intent(In) :: iE
    Real(dp), Intent(In) :: vR1,vR2
    Real(dp), Intent(In) :: gamma,v
    Real(dp) :: sig_vR1(1:2),sig_vR2(1:2),sig_vR(1:2)
    Real(dp) :: s1(1:2),s2(1:2)
    Real(dp) :: vR
    Integer, Parameter :: total = 1
    Integer, Parameter :: absor = 2
    Integer, Parameter :: Tmax = 15  !maximum number of extrapolations in the table
    Real(dp) :: T(1:2,0:Tmax)  !Extrapolation table previous row
    Real(dp) :: Tk0(1:2),Tk(1:2)  !Extrapolation table current row values
    Integer :: i,j,k  !counters: i for table row, j for quadrature ordinates, k for table column
    Integer :: n      !number of intervals
    Real(dp) :: h0,h  !spacing between quadrature ordinates
    Real(dp) :: fk    !multiplier for extrapolation steps

    !Initial trapezoid estimate
    Call CS%sig_T_A(Neutron_Energy(vR1),sig_vR1(total),sig_vR1(absor),iE_put=iE)
    Call CS%sig_T_A(Neutron_Energy(vR2),sig_vR2(total),sig_vR2(absor),iE_put=iE)
    s1 = 0.5_dp * (Broad_Integrand( vR1,sig_vR1,gamma,v) + Broad_Integrand( vR2,sig_vR2,gamma,v))
    s2 = 0.5_dp * (Broad_Integrand(-vR1,sig_vR1,gamma,v) + Broad_Integrand(-vR2,sig_vR2,gamma,v))
    h0 = vR2 - vR1
    T(:,0) = h0 * (s1 - s2)
    n = 1
    Do i = 1,Tmax !up to Tmax rows in the table
        !Trapezoid estimate for the 0-th column of the i-th row of table
        n = n * 2
        h = h0 / Real(n,dp)
        Do j = 1,n-1,2  !only odd values of j, these are the NEW points at which to evaluate the integrand
            vR = vR1 + Real(j,dp)*h
            Call CS%sig_T_A(Neutron_Energy(vR),sig_vR(total),sig_vR(absor),iE_put=iE)
            s1 = s1 + Broad_Integrand( vR,sig_vR,gamma,v)
            s2 = s2 + Broad_Integrand(-vR,sig_vR,gamma,v)
        End Do
        Tk0 = h * (s1 - s2)
        !Fill i-th row, columns k = 1:i, with extrapolated estimates
        fk = 1._dp
        Do k = 1,i  !up to i columns this row
            fk = fk * 4._dp
            Tk = (fk * Tk0 - T(k-1)) / (fk - 1._dp)
            If (k .LT. i) Then
                T(k-1) = Tk0  !store Tk0 for next i
                Tk0 = Tk  !store Tk for next k
            End If !otherwise, skip storage steps if working final column
        End Do
        !check for convergence
        If ( All( Abs(T(:,i-1) - Tk) .LE. rTol * Abs(Tk) ) ) Then
            sig_T_A = Tk  !Tk is the highest precision converged value
            Return  !Normal exit
        Else  !prep for the next time though the loop
            !store Tk0 and Tk for next i
            T(i-1) = Tk0
            T(i) = Tk
        End If
    End Do
    !If we get this far, we did not converge
    Write(*,*)
    Write(*,'(A,I0,A)')        'ERROR:  n_Cross_sections: Broad_Romberg_T_A:  Failed to converge in ',Tmax,' extrapolations.'
    Write(*,'(A,2ES23.15)')    '        Final estimated value: ',Tk
    Write(*,'(A,2ES23.15)')    '        Prior estimated value: ',Tk0
    ERROR STOP
End Function Broad_Romberg_T_A

Function sig_T_broad(CS,E,T) Result(sigT)
    Use Kinds, Only: dp
    Use Global, Only: SqrtPi
    Use Utilities, Only: Bisection_Search
    Use Neutron_Utilities, Only: Neutron_Energy
    Use Neutron_Utilities, Only: Neutron_Speed
    Implicit None
    Real(dp) :: sigT
    Class(CS_Type), Intent(In) :: CS
    Real(dp), Intent(In) :: E ![kev]
    Real(dp), Intent(In) :: T ![K]
    Real(dp) :: v  ![km/s] neutron velocity
    Real(dp) :: gamma  ![km/s]
    Real(dp) :: vR_min,vR_max
    Integer :: iE_min,iE_max
    Integer :: i

    Call Broad_sig_start(E,CS%Mn,T,v,gamma,vR_min,vR_max)
    iE_min = Bisection_Search(Neutron_Energy(vR_min),CS%E_uni,CS%n_E_uni)
    If (Neutron_Energy(vR_max) .LE. CS%E_uni(iE_min)) Then  !only a single interval spanned in energy grid
        sigT = Broad_Romberg_T(CS,iE_min,vR_min,vR_max,gamma,v) *  gamma / (SqrtPi * v**2)
        Return
    End If
    iE_max = Bisection_Search(Neutron_Energy(vR_max),CS%E_uni,CS%n_E_uni) - 1
    !first interval (from E_min to first indexed E) is partial
    sigT = Broad_Romberg_T(CS,iE_min,vR_min,Neutron_Speed(CS%E_uni(iE_min)),gamma,v)
    !middle intervals
    Do i = iE_min,iE_max-1
        sigT = sigT + Broad_Romberg_T(CS,i+1,Neutron_Speed(CS%E_uni(i)),Neutron_Speed(CS%E_uni(i+1)),gamma,v)
    End Do
    !last interval (from second to last indexed E to Emax) is partial, also apply normalization
    sigT = (sigT + Broad_Romberg_T(CS,iE_max,Neutron_Speed(CS%E_uni(iE_max)),vR_max,gamma,v)) * gamma / (SqrtPi * v**2)
End Function sig_T_broad

Function Broad_Romberg_T(CS,iE,vR1,vR2,gamma,v) Result(sigT)
    Use Kinds, Only: dp
    Use Neutron_Utilities, Only: Neutron_Energy
    Implicit None
    Real(dp) :: sigT
    Type(CS_type), Intent(In) :: CS
    Integer, Intent(In) :: iE
    Real(dp), Intent(In) :: vR1,vR2
    Real(dp), Intent(In) :: gamma,v
    Real(dp) :: h
    Real(dp) :: sig_vR1,sig_vR2,sig_vR
    Real(dp) :: s1,s2
    Real(dp) :: vR
    Integer :: n,i,j
    Real(dp) :: R(0:10,0:10)

    sig_vR1 = CS%sig_T(Neutron_Energy(vR1),iE_put=iE)
    sig_vR2 = CS%sig_T(Neutron_Energy(vR2),iE_put=iE)
    h = vR2 - vR1
    s1 = 0.5_dp * (Broad_Integrand( vR1,sig_vR1,gamma,v) + Broad_Integrand( vR2,sig_vR2,gamma,v))
    s2 = 0.5_dp * (Broad_Integrand(-vR1,sig_vR1,gamma,v) + Broad_Integrand(-vR2,sig_vR2,gamma,v))
    R(0,0) = h * (s1 - s2)
    n = 1
    Do i = 1,10
        !compute trapezoid estimate for next row of table
        n = n*2
        h = (vR2 - vR1) / Real(n,dp)
        Do j = 1,n-1,2  !only odd values of j, these are the NEW points at which to evaluate the integrand
            vR = vR1 + Real(j,dp)*h
            sig_vR = CS%sig_T(Neutron_Energy(vR),iE_put=iE)
            s1 = s1 + Broad_Integrand( vR,sig_vR,gamma,v)
            s2 = s2 + Broad_Integrand(-vR,sig_vR,gamma,v)
        End Do
        R(0,i) = h * (s1 - s2)
        !fill out Romberg table row
        Do j = 1,i
            R(j,i) = Romb2(j) * (Romb1(j) * R(j-1,i) - R(j-1,i-1))
        End Do
        !check for convergence
        If ( Abs(R(i-1,i-1) - R(i,i)) .LE. rTol * Abs(R(i,i)) ) Then
             sigT = R(i,i)
            Return
        End If
    End Do
    !if we get this far, we failed to converge

End Function Broad_Romberg_T

Function sig_S_iso_broad(CS,iso,E,T) Result(sigS)
    Use Kinds, Only: dp
    Use Global, Only: SqrtPi
    Use Global, Only: neutron_mass
    Use Utilities, Only: Bisection_Search
    Use Neutron_Utilities, Only: Neutron_Energy
    Use Neutron_Utilities, Only: Neutron_Speed
    Implicit None
    Real(dp) :: sigS
    Class(CS_Type), Intent(In) :: CS
    Integer, Intent(In) :: iso
    Real(dp), Intent(In) :: E ![kev]
    Real(dp), Intent(In) :: T ![K]
    Real(dp) :: v  ![km/s] neutron velocity
    Real(dp) :: gamma  ![km/s]
    Real(dp) :: vR_min,vR_max
    Integer :: iE_min,iE_max
    Integer :: i

    Call Broad_sig_start(E,CS%An(iso)*neutron_mass,T,v,gamma,vR_min,vR_max)
    iE_min = Bisection_Search(Neutron_Energy(vR_min),CS%E_uni,CS%n_E_uni)
    If (Neutron_Energy(vR_max) .LE. CS%E_uni(iE_min)) Then  !only a single interval spanned in energy grid
        sigS = Broad_Romberg_S_iso(CS,iso,iE_min,vR_min,vR_max,gamma,v) *  gamma / (SqrtPi * v**2)
        Return
    End If
    iE_max = Bisection_Search(Neutron_Energy(vR_max),CS%E_uni,CS%n_E_uni) - 1
    !first interval (from E_min to first indexed E) is partial
    sigS = Broad_Romberg_S_iso(CS,iso,iE_min,vR_min,Neutron_Speed(CS%E_uni(iE_min)),gamma,v)
    !middle intervals
    Do i = iE_min,iE_max-1
        sigS = sigS + Broad_Romberg_S_iso(CS,iso,i+1,Neutron_Speed(CS%E_uni(i)),Neutron_Speed(CS%E_uni(i+1)),gamma,v)
    End Do
    !last interval (from second to last indexed E to Emax) is partial, also apply normalization
    sigS = (sigS + Broad_Romberg_S_iso(CS,iso,iE_max,Neutron_Speed(CS%E_uni(iE_max)),vR_max,gamma,v)) * gamma / (SqrtPi * v**2)
End Function sig_S_iso_broad

Function Broad_Romberg_S_iso(CS,iso,iE,vR1,vR2,gamma,v) Result(sigS)
    Use Kinds, Only: dp
    Use Neutron_Utilities, Only: Neutron_Energy
    Implicit None
    Real(dp) :: sigS
    Type(CS_type), Intent(In) :: CS
    Integer, Intent(In) :: iso
    Integer, Intent(In) :: iE
    Real(dp), Intent(In) :: vR1,vR2
    Real(dp), Intent(In) :: gamma,v
    Real(dp) :: h
    Real(dp) :: sig_vR1,sig_vR2,sig_vR
    Real(dp) :: s1,s2
    Real(dp) :: vR
    Integer :: n,i,j
    Real(dp) :: R(0:10,0:10)
    
    sig_vR1 = CS%sig_S_iso(iso,Neutron_Energy(vR1),iE_put=iE)
    sig_vR2 = CS%sig_S_iso(iso,Neutron_Energy(vR2),iE_put=iE)
    h = vR2 - vR1
    s1 = 0.5_dp * (Broad_Integrand( vR1,sig_vR1,gamma,v) + Broad_Integrand( vR2,sig_vR2,gamma,v))
    s2 = 0.5_dp * (Broad_Integrand(-vR1,sig_vR1,gamma,v) + Broad_Integrand(-vR2,sig_vR2,gamma,v))
    R(0,0) = h * (s1 - s2)
    n = 1
    Do i = 1,10
        n = n*2
        h = (vR2 - vR1) / Real(n,dp)
        Do j = 1,n-1,2  !only odd values of j, these are the NEW points at which to evaluate the integrand
            vR = vR1 + Real(j,dp)*h
            sig_vR = CS%sig_S_iso(iso,Neutron_Energy(vR),iE_put=iE)
            s1 = s1 + Broad_Integrand( vR,sig_vR,gamma,v)
            s2 = s2 + Broad_Integrand(-vR,sig_vR,gamma,v)
        End Do
        R(0,i) = h * (s1 - s2)
        Do j = 1,i
            R(j,i) = Romb2(j) * (Romb1(j) * R(j-1,i) - R(j-1,i-1))
        End Do
        If ( Abs(R(i-1,i-1) - R(i,i)) .LE. rTol * Abs(R(i,i)) ) Then  !check convergence with specified tolerances
             sigS = R(i,i)
            Return
        End If
    End Do
    !if we get this far, we failed to converge

End Function Broad_Romberg_S_iso

Elemental Function Broad_Integrand(vR,sig,gamma,v)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: Broad_Integrand
    Real(dp), Intent(In) :: vR,sig,gamma,v

    Broad_Integrand = vR**2 * sig * Exp( -((gamma * (v-vR))**2) )
End Function Broad_Integrand

Subroutine Write_Cross_Sections(CS,Broadened_CS,file_name)
    Use Kinds, Only: dp
    Use FileIO_Utilities, Only: half_dash_line
    Implicit None
    Type(CS_Type), Intent(In) :: CS
    Logical, Intent(In) :: Broadened_CS
    Character(*), Intent(In) :: file_name
    Integer :: unit,stat
    Integer :: i
    Real(dp) :: sT,sA,sTb,sAb

    Open(NEWUNIT = unit , FILE = file_name , STATUS = 'UNKNOWN' , ACTION = 'WRITE' , POSITION = 'APPEND' , IOSTAT = stat)
    If (stat .NE. 0) Then
        Print *,'ERROR:  Cross_Sections: Write_Cross_Sections:  File open error, '//file_name//', IOSTAT=',stat
        ERROR STOP
    End If
    Write(unit,'(A)') half_dash_line
    Write(unit,'(A)') 'CROSS SECTIONS INFORMATION'
    Write(unit,'(A)') half_dash_line
    Write(unit,'(A)') '  Microscopic cross sections for total atmosphere:'
    Write(unit,*)
    If (Broadened_CS) Then
        Write(unit,'(5A27)') 'Energy [keV]','Total CS [barns]','Total CS [barns]','Absorption CS [barns]','Absorption CS [barns]'
        Write(unit,'(5A27)') '','(0 deg K)','(273.15 deg K)','(0 deg K)','(273.15 deg K)'
        Write(unit,'(5A27)') '-----------------------','-----------------------','-----------------------','-----------------------','-----------------------'
        Do i = 1,CS%n_E_uni
            Call CS%sig_T_A(CS%E_uni(i),sT,sA)
            Call CS%sig_T_A_broad(CS%E_uni(i),273.15_dp,sTb,sAb)
            Write(unit,'(5ES27.16E3)') CS%E_uni(i),sT,sTb,sA,sAb
        End Do
    Else
        Write(unit,'(3A27)') 'E [keV]','Total CS [barns]','Absorption CS [barns]'
        Write(unit,'(3A27)') '-----------------------','-----------------------','-----------------------'
        Do i = 1,cs%n_E_uni
            Call cs%sig_T_A(CS%E_uni(i),sT,sA)
            Write(unit,'(3ES27.16E3)') CS%E_uni(i),sT,sA
        End Do
    End If
    Write(unit,*)
    Write(unit,*)
    Close(unit)
End Subroutine Write_Cross_Sections

End Module n_Cross_Sections
