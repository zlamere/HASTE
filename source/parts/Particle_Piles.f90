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
Module Particle_Piles

    Use Kinds, Only: dp
    Implicit None
    Private
    Public :: Pile_Type

    Type Pile_Particle_Type
        Real(dp) :: r(1:3)  !position
        Real(dp) :: OmegaHat(1:3)  !direction
        Real(dp) :: E  !energy
        Real(dp) :: w  !weight
        Real(dp) :: t  !time of emission
        Type(Pile_Particle_Type), Pointer :: up => Null()  !points to the previous space in the pile
        Type(Pile_Particle_Type), Pointer :: dn => Null()  !points to the next space in the pile
    Contains
        Procedure, Pass :: fill => Fill_Pile_Particle  !copies particle values into a space on the pile
        Procedure, Pass :: empty => Empty_Pile_Particle  !removes particle values from a space on the pile
    End Type

    Type Pile_Type
        Type(Pile_Particle_Type), Pointer :: head => Null()  !points to the first space on the pile
        Type(Pile_Particle_Type), Pointer :: tail => Null()  !points to the last space on the pile
        Type(Pile_Particle_Type), Pointer :: ptr => Null()  !points to the last particle currently on the pile
        Integer :: num_piled = 0  !tracks the current number of piled particles
        Integer :: max_piled = 0  !tracks the maximum extent of the pile
        Logical :: First_Time = .TRUE.  !indicates whether the file has been initialized
    Contains
        Procedure, Pass :: Put => Put_on_Pile  !puts a particle in the next empty space on the pile, adding space only when required
        Procedure, Pass :: Get => Get_from_Pile  !gets the last particle currently on the pile
    End Type

Contains

Function Setup_Particle_Pile() Result(p)
    Implicit None
    Type(Pile_Type) :: p

    Allocate(p%head)  !allocate the first space on the pile
    Nullify(p%head%up)  !no space is referenced before this one
    p%tail => p%head  !this space is currently also the last space on the pile
    Nullify(p%tail%dn)  !no space follows the last space on the pile
    Nullify(p%ptr)  !this space has not been filled, so do not yet reference it
    p%First_Time = .FALSE.  !mark the pile as set up
    p%num_piled = 0  !initialize piled counter
    p%max_piled = 1  !initialize space counter
End Function Setup_Particle_Pile

Subroutine Put_On_Pile(pile,r,OmegaHat,E,w,t)
    Use Kinds, Only: dp
    Implicit None
    Class(Pile_Type), Intent(InOut) :: pile  !this variable is automatically passed and need not be explicit in subroutine calls
    Real(dp), Intent(In) :: r(1:3)  !position
    Real(dp), Intent(In) :: OmegaHat(1:3)  !direction
    Real(dp), Intent(In) :: E  !energy
    Real(dp), Intent(In) :: w  !weight
    Real(dp), Intent(In) :: t  !time of emission

    If (.NOT. Associated(pile%ptr)) Then  ! There are no particles on the pile
        Call pile%head%fill(r,OmegaHat,E,w,t)  !Put the new particle in the space at the head of the pile
        pile%ptr => pile%head  !point to newly filled space
        pile%num_piled = pile%num_piled + 1  !increment the piled counter
        Return  !the structure of the pile has not changed, so there is nothing else to do
    Else If (.NOT. Associated(pile%ptr%dn)) Then  ! There is not pile space allocated for new values
        Allocate(pile%ptr%dn)  !Create new pile space
        pile%tail    => pile%ptr%dn  !the new space is now the last space in the pile
        pile%tail%up => pile%ptr  !the current pointer is the space before the new space
        Nullify(pile%tail%dn)  !no space follows the last space in the pile
        pile%max_piled = pile%max_piled + 1  !increment the space counter
    End If
    pile%ptr => pile%ptr%dn  !Put the new particle in the next empty space in the pile
    Call pile%ptr%fill(r,OmegaHat,E,w,t)  !copy the values into the particle
    pile%num_piled = pile%num_piled + 1  !increment the piled counter
End Subroutine Put_On_Pile

Subroutine Fill_Pile_Particle(p,r,OmegaHat,E,w,t)
    Use Kinds, Only: dp
    Implicit None
    Class(Pile_Particle_Type), Intent(InOut) :: p
    Real(dp), Intent(In) :: r(1:3)  !position
    Real(dp), Intent(In) :: OmegaHat(1:3)  !direction
    Real(dp), Intent(In) :: E  !energy
    Real(dp), Intent(In) :: w  !weight
    Real(dp), Intent(In) :: t  !time of emission

    p%r = r
    p%OmegaHat = OmegaHat
    p%E = E
    p%w = w
    p%t = t
End Subroutine Fill_Pile_Particle

Function Get_From_Pile(pile,r,OmegaHat,E,w,t) Result(particle_on_pile)
    Use Kinds, Only: dp
    Implicit None
    Logical :: particle_on_pile
    Class(Pile_Type), Intent(InOut) :: pile !this variable is automatically passed and need not be explicit in function calls
    Real(dp), Intent(Out) :: r(1:3)  !position
    Real(dp), Intent(Out) :: OmegaHat(1:3)  !direction
    Real(dp), Intent(Out) :: E  !energy
    Real(dp), Intent(Out) :: w  !weight
    Real(dp), Intent(Out) :: t  !time of emission

    If (Associated(pile%ptr)) Then  ! There is a particle on the pile, CONTINUE
        particle_on_pile = .TRUE.  !a particle will be returned by the function call
    Else  ! The pile is empty, RETURN
        particle_on_pile = .FALSE.  !a particle will NOT be returned by the function call
        Return
    End If
    Call pile%ptr%empty(r,OmegaHat,E,w,t)  !Get the particle from the bottom of the pile
    pile%ptr => pile%ptr%up  !Point to the new pile bottom above the particle just emptied
    pile%num_piled = pile%num_piled - 1  !decrement the piled counter
End Function Get_From_pile

Subroutine Empty_Pile_Particle(p,r,OmegaHat,E,w,t)
    Use Kinds, Only: dp
    Implicit None
    Class(Pile_Particle_Type), Intent(InOut) :: p
    Real(dp), Intent(Out) :: r(1:3)  !position
    Real(dp), Intent(Out) :: OmegaHat(1:3)  !direction
    Real(dp), Intent(Out) :: E  !energy
    Real(dp), Intent(Out) :: w  !weight
    Real(dp), Intent(Out) :: t  !time of emission

    !get values
    r = p%r
    OmegaHat = p%OmegaHat
    E = p%E
    w = p%w
    t = p%t
    !write nonsense vaues into the now empty particle
    p%r = 0._dp
    p%OmegaHat = 0._dp
    p%E = -1._dp
    p%w = 0._dp
    p%t = 0._dp
End Subroutine Empty_Pile_Particle

End Module Particle_Piles
