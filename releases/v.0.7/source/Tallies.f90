Module Tallies
    
    Use Kinds, Only: dp
    Implicit None
    Private
    Public :: Contrib_triplet
    Public :: Contrib_quadruplet
    Public :: Contrib_array
    Public :: Setup_Tallies
    
    Type :: Contrib_triplet
        Integer :: i1
        Integer :: i2
        Real(dp) :: f  !sum of contributions over all histories to the bin
    End Type
    
    Type, Extends(Contrib_triplet) :: Contrib_quadruplet
        Real(dp) :: f_sq  !sum of squared contributions over all histories to the bin
    End Type
    
    Type :: Contrib_array
        Real(dp) :: tot_f_sq  !sum of squared contributions for computing variance of total flux at the detector
        Real(dp), Allocatable :: f_sq_1(:)  !sum of squared contributions for computing variance of total flux in each 1st index bin at the detector
        Real(dp), Allocatable :: f_sq_2(:)  !sum of squared contributions for computing variance of total flux in each 2nd index bin at the detector
        Integer :: index  !number of bins in list with contributions, index is less than or equal to size
        Integer :: size  !extent of the list, size is greater than or equal to index
        Type(Contrib_quadruplet), Allocatable :: contribs(:)  !list of indexes and contributions to the time-energy grid
    Contains
        Procedure, Pass :: Tally_history
    End Type
    
Contains

Function Setup_Tallies(n_1_bins,n_2_bins) Result(t)
    Use Kinds, Only: dp
    Implicit None
    Type(Contrib_array) :: t
    Integer :: n_1_bins,n_2_bins
    
    t%tot_f_sq = 0._dp
    Allocate(t%f_sq_1(1:n_1_bins))
    t%f_sq_1 = 0._dp
    Allocate(t%f_sq_2(1:n_2_bins))
    t%f_sq_2 = 0._dp
    t%index = 0
    !use the smaller dimension of the grid as a starting length for the array
    t%size = Min(n_1_bins,n_2_bins)
    Allocate(t%contribs(1:t%size))
End Function

Subroutine Tally_History(t,n,contribs,n1,contribs1,n2,contribs2)
    Use Kinds, Only: dp
    Implicit None
    Class(Contrib_array), Intent(InOut) :: t
    Integer, Intent(InOut) :: n
    Type(Contrib_triplet), Intent(InOut) :: contribs(1:n)
    Integer, Intent(In) :: n1
    Real(dp), Intent(InOut) :: contribs1(1:n1)
    Integer, Intent(In) :: n2
    Real(dp), Intent(InOut) :: contribs2(1:n2)
    Integer :: i,j,k
    Real(dp) :: f_tot

    If (n .EQ. 0) Return  !no tallies for this neutron
    !resize the list if necessary
    If (t%size .LT. t%index + n) Call Resize_TE_Array(t,n)
    !check if the first contribution belongs at the head of the list
    If ( (t%index .EQ. 0) .OR. &  !this is the first set of tallies
       & (contribs(1)%i1 .LT. t%contribs(1)%i1) .OR. &  !this contribution is earlier in 1st index than the earliest contribution in the list
       & ((contribs(1)%i1 .EQ. t%contribs(1)%i1).AND.(contribs(1)%i2 .LT. t%contribs(1)%i2)) ) &  !this contribution is the same 1st index, but lower in 2nd index than the first contribution in the list
       & Then !add the bin to the head of the list
        !move the current entries down to make room at the top
        t%contribs(2:t%index+1) = t%contribs(1:t%index)
        t%index = t%index + 1
        !create an entry for the contribution at the top of the list
        t%contribs(1)%i1 = contribs(1)%i1
        t%contribs(1)%i2 = contribs(1)%i2
        t%contribs(1)%f = 0._dp
        t%contribs(1)%f_sq = 0._dp
    End If
    !now all contributions belong in or after the first bin in the list
    i = 1  !index in contribs
    j = 1  !index in t
    f_tot = 0._dp
    Do While (i .LE. n)  !i is index in contribs
        If (j .GT. t%index) Then  !the remaining contributions go at the end of t
            Do k = i,n
                t%index = t%index + 1
                !copy the contribution into the last entry of the list
                t%contribs(t%index)%i1 = contribs(k)%i1
                t%contribs(t%index)%i2 = contribs(k)%i2
                t%contribs(t%index)%f = contribs(k)%f
                t%contribs(t%index)%f_sq = contribs(k)%f**2
                f_tot = f_tot + contribs(k)%f
            End Do
            Exit
        End If
        !search for the right 1st index bin
        If (contribs(i)%i1 .LE. t%contribs(j)%i1) Then  !this might be the right 1st index bin
            If (contribs(i)%i1 .EQ. t%contribs(j)%i1) Then  !this is the right 1st index bin
                !search for the right 2nd index bin
                If (contribs(i)%i2 .LE. t%contribs(j)%i2) Then  !This might be the right 2nd index bin
                    If (contribs(i)%i2 .EQ. t%contribs(j)%i2) Then  !this is the right 2nd index bin
                        !add the contribution
                        t%contribs(j)%f = t%contribs(j)%f + contribs(i)%f
                        t%contribs(j)%f_sq = t%contribs(j)%f_sq + contribs(i)%f**2
                        f_tot = f_tot + contribs(i)%f
                    Else  !need to insert this contribution into the list
                        !move the current entries down to make room
                        t%contribs(j+1:t%index+1) = t%contribs(j:t%index)
                        t%index = t%index + 1
                        !copy the contribution into the j-th entry of the list
                        t%contribs(j)%i1 = contribs(i)%i1
                        t%contribs(j)%i2 = contribs(i)%i2
                        t%contribs(j)%f = contribs(i)%f
                        t%contribs(j)%f_sq = contribs(i)%f**2
                        f_tot = f_tot + contribs(i)%f
                    End If
                    i = i + 1
                    j = j + 1
                Else !this is not the right 2nd index bin, advance in the list and repeat the check
                    j = j + 1
                End If
            Else  !need to insert this contribution into the list
                !move the current entries down to make room
                t%contribs(j+1:t%index+1) = t%contribs(j:t%index)
                t%index = t%index + 1
                !copy the contribution into the j-th entry of the list
                t%contribs(j)%i1 = contribs(i)%i1
                t%contribs(j)%i2 = contribs(i)%i2
                t%contribs(j)%f = contribs(i)%f
                t%contribs(j)%f_sq = contribs(i)%f**2
                f_tot = f_tot + contribs(i)%f
                i = i + 1
                j = j + 1
            End If
        Else !this is not the right 1st index bin, advance
            j = j + 1
        End If
    End Do
    !update the summed squares for 1st and 2nd index bins, and reset the arrays for the next history after they are added
    t%f_sq_1 = t%f_sq_1 + contribs1**2
    contribs1 = 0._dp
    t%f_sq_2 = t%f_sq_2 + contribs2**2
    contribs2 = 0._dp
    !update the summed squares for the total estimate
    t%tot_f_sq = t%tot_f_sq + f_tot**2
    !reset the number of contributions to zero
    n = 0
End Subroutine Tally_History

Subroutine Resize_TE_Array(list,n)
    Use Kinds, Only: dp
    Implicit None
    Type(Contrib_array), Intent(InOut) :: list
    Integer, Intent(In) :: n  !number of elements to add
    Type(Contrib_quadruplet), Allocatable :: swap(:)
        
    Allocate(swap(1:list%size))
    swap = list%contribs
    Deallocate(list%contribs)
    Allocate(list%contribs(1:list%size+n))
    list%contribs(:)%i1 = -1
    list%contribs(:)%i2 = -1
    list%contribs(:)%f = 0._dp
    list%contribs(:)%f_sq = 0._dp
    list%contribs(1:list%size) = swap
    list%size = list%size + n
    Deallocate(swap)
End Subroutine Resize_TE_Array

End Module Tallies