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
Module Quadratures
    
    Use Kinds, Only: dp
    Implicit None
    Private
    Public :: Romberg_Quad
    Public :: Composite_Quad
    Public :: Composite_Trapezoid
    Public :: Composite_Simpson
    Public :: Adaptive_Quad
    Public :: Adaptive_Simpson
    Public :: MidPoint
    Public :: Trapezoid
    Public :: Simpson
    Public :: GaussLegendre5
    Public :: GaussLegendre16
    Public :: GaussLegendre96
    Public :: GaussLegendreN
    Public :: Progressive_GaussLegendre
    
    Real(dp), Parameter :: one_third = 1._dp / 3._dp
    Real(dp), Parameter :: one_sixth = 1._dp / 6._dp
    Real(dp), Parameter :: one_fifteenth = 1._dp / 15._dp
    
    Interface Romberg_Quad
        Module Procedure Romb_Quad
        Module Procedure Romb_Quad_ranges
    End Interface

Contains

Function Romb_Quad_ranges(f,ab,aTol,rTol,p) Result(q)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: q
    Interface
        Function f(x)    !the function to be integrated
            Use Kinds,Only: dp
            Implicit None
            Real(dp) :: f
            Real(dp), Intent(In) :: x
        End Function f
    End Interface
    Real(dp), Intent(In) :: ab(:)        !limits of integration
    Real(dp), Intent(In) :: rTol,aTol  !relative and absolute tolerances for convergence
    Real(dp), Intent(Out), Optional :: p  !estimated precision achived in the final estimate
    Real(dp) :: h,hi
    Integer :: n,i
    Real(dp) :: p_i

    q = 0._dp
    n = Size(ab)
    h = ab(n) - ab(1)
    Do i = 1,n-1
        hi = ab(i+1) - ab(i)
        If (Present(p)) Then
            q = q + Romb_Quad(f,ab(i),ab(i+1),(hi/h)*aTol,rTol,p_i)
            p = Min(p,p_i)
        Else
            q = q + Romb_Quad(f,ab(i),ab(i+1),(hi/h)*aTol,rTol)
        End If
    End Do
End Function Romb_Quad_ranges

Function Romb_Quad(f,a,b,aTol,rTol,p,n_ord,n_ext) Result(q)
    !Integrates f(x) on (a,b) by extrapolation on successive composite trapezoid
    Use Kinds, Only: dp
    Use Utilities, Only: Converged
    Use Utilities, Only: Prec
    Implicit None
    Real(dp):: q    !the result of the integration
    Interface
        Function f(x)    !the function to be integrated
            Use Kinds,Only: dp
            Implicit None
            Real(dp) :: f
            Real(dp), Intent(In) :: x
        End Function f
    End Interface
    Real(dp), Intent(In) :: a,b        !limits of integration
    Real(dp), Intent(In) :: rTol,aTol  !relative and absolute tolerances for convergence
    Real(dp), Intent(Out), Optional :: p  !estimated precision achived in the final estimate
    Integer, Intent(Out), Optional :: n_ord  !total number of function evaluations (ordinates)
    Integer, Intent(Out), Optional :: n_ext  !total number of extrapolation stages (rows in the table)
    Integer, Parameter :: Tmax = 15  !maximum number of extrapolations in the table
    Real(dp) :: T(0:Tmax)  !Extrapolation table previous row
    Real(dp) :: Tk0,Tk  !Extrapolation table current row values
    Integer :: i,j,k  !counters: i for table row, j for quadrature ordinates, k for table column
    Integer :: n      !number of intervals
    Real(dp) :: h0,h  !spacing between quadrature ordinates
    Real(dp) :: fk    !multiplier for extrapolation steps
    Real(dp) :: s     !sum of function values at quadrature ordinates
    Real(dp) :: Prec0 !precision of the previous extrapolation (used to monitor convergence)
    Integer :: ifin   !final number of extrapolations performed
    !    If the preprocessor directive ROMB_TABLES is set at compile-time, 
    ! then the tables of computed values are dumped to file in the working 
    ! directory as the routine runs.
#   if ROMB_TABLES
        Integer :: unit
#   endif

    !Initial trapezoid estimate
    n = 1
    s = 0.5_dp * (f(a) + f(b))
    h0 = b - a
    T(0) = h0 * s
    Prec0 = -1._dp
#   if ROMB_TABLES
        Open(NEWUNIT=unit,FILE='Romb_Tables.tst',ACTION='WRITE',STATUS='UNKNOWN',POSITION='APPEND')
        Do k = 0,Tmax
            Write(unit,'(I24)',ADVANCE='NO') k
        End Do
        Write(unit,'(ES24.15)') T(0)
        Close(unit)
#   endif
    Do i = 1,Tmax !up to Tmax rows in the table
        !Trapezoid estimate for the 0-th column of the i-th row of table
        n = n * 2
        h = h0 / Real(n,dp)
        Do j = 1,n-1,2  !Odd values of j are NEW points at which to evaluate f
            s = s + f(a + Real(j,dp)*h)
        End Do
        Tk0 = h * s
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
        !Check for convergence
        If (Converged(T(i-1),Tk,rTol,aTol)) Then
            q = Tk
#           if ROMB_TABLES
                T(i-1) = Tk0
                T(i) = Tk
                Open(NEWUNIT=unit,FILE='Romb_Tables.tst',ACTION='WRITE',STATUS='UNKNOWN',POSITION='APPEND')
                Do k = 0,i
                    Write(unit,'(ES24.15)',ADVANCE='NO') T(k)
                End Do
                Write(unit,'(A)') '*'
                Write(unit,*)
                Close(unit)
#           endif
            If (Present(p)) p = Prec(T(i-1),Tk)
            If (Present(n_ord)) n_ord = n
            If (Present(n_ext)) n_ext = i
            Return  !Normal exit
        Else  !Check for failures and prep for the next time though the loop
            !check for exit conditions other than convergence (failures)
            If (i .EQ. Tmax) Then  !maximum extrapolation has been reached without convergence
                ifin = Tmax
                Exit
            Else If (i.GT.Tmax/2) Then
                If (Prec(T(i-1),Tk).LT.Prec0) Then  !precision was LOST on this extrapolation, convergence will not occur
                    ifin = i
                    Exit
                End If
            End If
            !store Tk0 and Tk for next i, and update precision monitor
            Prec0 = Prec(T(i-1),Tk)
            T(i-1) = Tk0
            T(i) = Tk
#           if ROMB_TABLES
                Open(NEWUNIT=unit,FILE='Romb_Tables.tst',ACTION='WRITE',STATUS='UNKNOWN',POSITION='APPEND')
                Do k = 0,i
                    Write(unit,'(ES24.15)',ADVANCE='NO') T(k)
                End Do
                Write(unit,*)
                Close(unit)
#           endif
        End If
    End Do
    !If we get this far, we did not converge
    Write(*,*)
    Write(*,'(A,I0,A)')          'WARNING:  Quadratures: Romberg_Quad:  Failed to converge in ',ifin,' extrapolations.'
    If (ifin .LT. Tmax) Then
        Write(*,'(A10,A)')              '','Extrapolation was terminated for loss of precision.'
    End If
    Write(*,'(A10,A,ES23.15,A,F0.5,A)') '','Final estimated value: ',Tk,' (~',Prec(T(ifin-1),Tk),' good digits)'
    Write(*,'(A10,A,2(ES10.3,A))')      '','Final Extrapolation Error: ',Abs(Tk-T(ifin-1)),' (abs), ',Abs(Tk-T(ifin-1))/Tk,' (rel)'
    Write(*,'(A10,A,2(ES10.3,A))')      '','Convergence Criteria:      ',atol,             ' (abs), ',rtol,                ' (rel)'
    q = Tk
#   if ROMB_TABLES
        T(i-1) = Tk0
        T(i) = Tk
        Open(NEWUNIT=unit,FILE='Romb_Tables.tst',ACTION='WRITE',STATUS='UNKNOWN',POSITION='APPEND')
        Do k = 0,i
            Write(unit,'(ES24.15)',ADVANCE='NO') T(k)
        End Do
        Write(unit,'(A)') '*'
        Write(unit,*)
        Close(unit)
#   endif
    If (Present(p)) p = Prec(T(ifin-1),Tk)
    If (Present(n_ord)) n_ord = n
    If (Present(n_ext)) n_ext = ifin
    ! ERROR STOP
End Function Romb_Quad

Function Adaptive_Quad(f,a,b,Quad,rTol,aTol) Result(q)
    Use Kinds, Only: dp
    Use Utilities, Only: Converged
    Implicit None
    Real(dp) :: q    !the result of the integration
    Interface
        Function f(x)  !the function to be integrated
            Use Kinds,Only: dp
            Implicit None
            Real(dp) :: f
            Real(dp), Intent(In) :: x
        End Function f
        Function Quad(f,a,b)  !the quadrature rule to use
            Use Kinds,Only: dp
            Implicit None
            Real(dp) :: Quad  !the result of the integration
            Interface
                Function f(x)    !the function to be integrated
                    Use Kinds,Only: dp
                    Implicit None
                    Real(dp) :: f
                    Real(dp), Intent(In) :: x
                End Function f
            End Interface
            Real(dp),Intent(In) :: a,b  !limits of integration
        End Function Quad
    End Interface
    Real(dp), Intent(In) :: a,b    !limits of integration
    Real(dp), Intent(In) :: rTol,aTol  !relative and absolute tolerance to be used on each interval as criteria for convergence
    Real(dp) :: q0,q1,q2
    
    q0 = Quad(f,a,b)
    q1 = Quad(f,a,0.5_dp*(a+b))
    q2 = Quad(f,0.5_dp*(a+b),b)
    q = q1 + q2
    If (Converged(q,q0,rTol,aTol)) Return  !check convergence with specified tolerances
    !otherwise, recurse to refine grid
    q = Continue_Adaptive_Quad(f,a,0.5_dp*(a+b),Quad,q1,rTol,0.5_dp*aTol,level=2) + &
      & Continue_Adaptive_Quad(f,0.5_dp*(a+b),b,Quad,q2,rTol,0.5_dp*aTol,level=2)
End Function Adaptive_Quad

Recursive Function Continue_Adaptive_Quad(f,a,b,Quad,q0,rTol,aTol,level) Result(q)
    Use Kinds, Only: dp
    Use Utilities, Only: Converged
    Implicit None
    Real(dp) :: q    !the result of the integration
    Interface
        Function f(x)  !the function to be integrated
            Use Kinds,Only: dp
            Implicit None
            Real(dp) :: f
            Real(dp), Intent(In) :: x
        End Function f
        Function Quad(f,a,b)  !the quadrature rule to use
            Use Kinds,Only: dp
            Implicit None
            Real(dp) :: Quad  !the result of the integration
            Interface
                Function f(x)    !the function to be integrated
                    Use Kinds,Only: dp
                    Implicit None
                    Real(dp) :: f
                    Real(dp), Intent(In) :: x
                End Function f
            End Interface
            Real(dp),Intent(In) :: a,b  !limits of integration
        End Function Quad
    End Interface
    Real(dp), Intent(In) :: a,b    !limits of integration
    Real(dp), Intent(In) :: q0
    Real(dp), Intent(In) :: rTol,aTol  !relative and absolute tolerance to be used on each interval as criteria for convergence
    Integer, Intent(In) :: level
    Real(dp) :: q1,q2
    
    q1 = Quad(f,a,0.5_dp*(a+b))
    q2 = Quad(f,0.5_dp*(a+b),b)
    q = q1 + q2
    If (Converged(q,q0,rTol,aTol)) Return  !check convergence with specified tolerances
    If (level .GT. 30) Then !max allowed recursion depth, interval has been split 31 times...
        Write(*,'(A)') 'ERROR:  Quadratures: Continue_Adaptive_Quad:  Failed to converge before reaching max recursion depth.'
        ERROR STOP
    End If
    !otherwise, recurse to refine grid
    q = Continue_Adaptive_Quad(f,a,0.5_dp*(a+b),Quad,q1,rTol,0.5_dp*aTol,level+1) + &
      & Continue_Adaptive_Quad(f,0.5_dp*(a+b),b,Quad,q2,rTol,0.5_dp*aTol,level+1)
End Function Continue_Adaptive_Quad

Function Adaptive_Simpson(f,a,b,rTol,aTol) Result(q)
    Use Kinds,Only: dp
    Implicit None
    Real(dp) :: q    !the result of the integration
    Interface
        Function f(x)  !the function to be integrated
            Use Kinds,Only: dp
            Implicit None
            Real(dp) :: f
            Real(dp), Intent(In) :: x
        End Function f
    End Interface
    Real(dp), Intent(In) :: a,b    !limits of integration
    Real(dp), Intent(In) :: rTol,aTol  !relative and absolute tolerance to be used on each interval as criteria for convergence
    Real(dp) :: x2,x3,x4
    Real(dp) :: y1,y2,y3,y4,y5
    Real(dp) :: q0,q1,q2
    Real(dp) :: err_est
    
    x3 = 0.5_dp * (a + b)
    x2 = 0.5_dp * (a + x3)
    x4 = 0.5_dp * (x3 + b)
    y1 = f(a)
    y2 = f(x2)
    y3 = f(x3)
    y4 = f(x4)
    y5 = f(b)
    q0 = (b-a) * (y1 + 4._dp*y3 + y5) * one_sixth
    q1 = (x3-a) * (y1 + 4._dp*y2 + y3) * one_sixth
    q2 = (b-x3) * (y3 + 4._dp*y4 + y5) * one_sixth
    err_est = (q1 + q2 - q0) * one_fifteenth
    If (Abs(err_est).LE.aTol .OR. Abs(err_est).LE.rTol*Min(Abs(q1+q2),Abs(q0))) Then  !check convergence with specified tolerances
        q = q1 + q2 + err_est
        Return
    End If
    !otherwise, recurse to refine grid
    q = Continue_Adaptive_Simpson(f,a,x2,x3,y1,y2,y3,q1,rTol,0.5_dp*aTol,level=2) + &
      & Continue_Adaptive_Simpson(f,x3,x4,b,y3,y4,y5,q2,rTol,0.5_dp*aTol,level=2)
End Function Adaptive_Simpson

Recursive Function Continue_Adaptive_Simpson(f,x1,x3,x5,y1,y3,y5,q0,rTol,aTol,level) Result(q)
    Use Kinds, Only: dp
    Implicit none
    Real(dp) :: q
    Interface
        Function f(x)  !the function to be integrated
            Use Kinds,Only: dp
            Implicit None
            Real(dp) :: f
            Real(dp), Intent(In) :: x
        End Function f
    End Interface
    Real(dp), Intent(In) :: x1,x3,x5
    Real(dp), Intent(In) :: y1,y3,y5
    Real(dp), Intent(In) :: q0
    Real(dp), Intent(In) :: rTol,aTol
    Integer, Intent(In) :: level
    Real(dp) :: x2,x4
    Real(dp) :: y2,y4
    Real(dp) :: q1,q2
    Real(dp) :: err_est
    
    x2 = 0.5_dp * (x1 + x3)
    x4 = 0.5_dp * (x3 + x5)
    y2 = f(x2)
    y4 = f(x4)
    q1 = (x3-x1) * (y1 + 4._dp*y2 + y3) * one_sixth
    q2 = (x5-x3) * (y3 + 4._dp*y4 + y5) * one_sixth
    err_est = (q1 + q2 - q0) * one_fifteenth
    If (Abs(err_est).LE.aTol .OR. Abs(err_est).LE.rTol*Min(Abs(q1+q2),Abs(q0))) Then  !check convergence with specified tolerances
        q = q1 + q2 + err_est
        Return
    End If
    If (level .GT. 30) Then !max allowed recursion depth, interval has been split 31 times...
        Write(*,'(A)') 'ERROR:  Quadratures: Continue_Adaptive_Simpson:  Failed to converge before reaching max recursion depth.'
        ERROR STOP
    End If
    !Otherwise, recurse to refine grid
    q = Continue_Adaptive_Simpson(f,x1,x2,x3,y1,y2,y3,q1,rTol,0.5_dp*aTol,level+1) + &
      & Continue_Adaptive_Simpson(f,x3,x4,x5,y3,y4,y5,q2,rTol,0.5_dp*aTol,level+1)
End Function Continue_Adaptive_Simpson

Function Composite_Quad(f,a,b,Quad,aTol,rTol) Result(q)
    Use Kinds, Only: dp
    Use Utilities, Only: Converged
    Implicit None
    Real(dp) :: q    !the result of the integration
    Interface
        Function f(x)  !the function to be integrated
            Use Kinds, Only: dp
            Implicit None
            Real(dp) :: f
            Real(dp), Intent(In) :: x
        End Function f
        Function Quad(f,a,b)  !the quadrature rule to use
            Use Kinds, Only: dp
            Implicit None
            Real(dp) :: Quad  !the result of the integration
            Interface
                Function f(x)    !the function to be integrated
                    Use Kinds ,Only: dp
                    Implicit None
                    Real(dp) :: f
                    Real(dp), Intent(In) :: x
                End Function f
            End Interface
            Real(dp), Intent(In) :: a,b  !limits of integration
        End Function Quad
    End Interface
    Real(dp), Intent(In) :: a,b    !limits of integration
    Real(dp), Intent(In) :: aTol,rTol
    Integer :: i,j,n
    Real(dp) :: h,qOld
    
    n = 1
    qOld = Quad(f,a,b)
    Do i = 1,31  !max intervals is bounded by splitting interval 31 times
        q = 0._dp
        n = n*2
        h = (b-a) / Real(n,dp)
        Do j = 1,n
            q = q + Quad(f,a+(Real((j-1),dp)*h),a+(Real(j,dp)*h))
        End Do
        If (Converged(q,qOld,rTol,aTol)) Return  !check convergence with specified tolerances
        qOld = q
    End Do
    !if we get this far, we failed to converge
    Write(*,'(A,I0,A)') 'ERROR:  Quadratures: Composite_Quad:  Failed to converge on ',n,' intervals.'
    ERROR STOP
End Function Composite_Quad

Function MidPoint(f,a,b) Result(q)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: q    !the result of the integration
    Interface
        Function f(x)    !the function to be integrated
            Use Kinds, Only: dp
            Implicit None
            Real(dp) :: f
            Real(dp), Intent(In) :: x
        End Function f
    End Interface
    Real(dp), Intent(In) :: a,b    !limits of integration

    q = (b-a) * f(0.5_dp*(a+b))
End Function MidPoint

Function Trapezoid(f,a,b) Result(q)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: q    !the result of the integration
    Interface
        Function f(x)    !the function to be integrated
            Use Kinds, Only: dp
            Implicit None
            Real(dp) :: f
            Real(dp), Intent(In) :: x
        End Function f
    End Interface
    Real(dp), Intent(In) :: a,b    !limits of integration

    q = 0.5_dp * (b-a) * (f(a) + f(b))
End Function Trapezoid

Function Simpson(f,a,b) Result(q)
    Use Kinds, Only: dp
    Implicit None
    Real(dp):: q    !the result of the integration
    Interface
        Function f(x)    !the function to be integrated
            Use Kinds, Only: dp
            Implicit None
            Real(dp) :: f
            Real(dp), Intent(In) :: x
        End Function f
    End Interface
    Real(dp), Intent(In) :: a,b    !limits of integration
    
    q = (b-a) * (f(a) + 4._dp*f(0.5_dp*(a+b)) + f(b)) * one_sixth
End Function Simpson

Function Composite_Trapezoid(f,a,b,aTol,rTol) Result(q)
    Use Kinds, Only: dp
    Use Utilities, Only: Converged
    Implicit None
    Real(dp) :: q    !the result of the integration
    Interface
        Function f(x)  !the function to be integrated
            Use Kinds, Only: dp
            Implicit None
            Real(dp) :: f
            Real(dp), Intent(In) :: x
        End Function f
    End Interface
    Real(dp), Intent(In) :: a,b    !limits of integration
    Real(dp), Intent(In) :: aTol,rTol
    Integer :: i,j,n
    Real(dp) :: qOld,s,h
    
    n = 1
    s = 0.5_dp * (f(a) + f(b))
    qOld = (b - a) * s
    Do i = 1,31
        n = n * 2
        h = (b - a) / Real(n,dp)
        Do j = 1,n-1,2  !only odd values of j, these are the NEW points at which to evaluate f
            s = s + f(a + Real(j,dp)*h)
        End Do
        q = h * s
        If (Converged(q,qOld,rTol,aTol)) Return  !check convergence with specified tolerances
        qOld = q
    End Do
    !if we get this far, we failed to converge
    Write(*,'(A,I0,A)') 'ERROR:  Quadratures: Composite_Trapezoid:  Failed to converge on ',n,' intervals.'
    ERROR STOP
End Function Composite_Trapezoid

Function Composite_Simpson(f,a,b,aTol,rTol) Result(q)
    Use Kinds, Only: dp
    Use Utilities, Only: Converged
    Implicit None
    Real(dp) :: q    !the result of the integration
    Interface
        Function f(x)  !the function to be integrated
            Use Kinds, Only: dp
            Implicit None
            Real(dp) :: f
            Real(dp), Intent(In) :: x
        End Function f
    End Interface
    Real(dp), Intent(In) :: a,b    !limits of integration
    Real(dp), Intent(In) :: aTol,rTol
    Integer :: i,j,n
    Real(dp) :: qOld,s1,s2,s3,h
    
    n = 2
    s1 = f(a) + f(b)
    s2 = 0._dp
    s3 = f(0.5_dp * (a + b))
    qOld = (b - a) * (s1 + 4._dp*s3) * one_sixth
    Do i = 1,31
        n = n * 2
        h = (b - a) / Real(n,dp)
        s2 = s2 + s3
        s3 = 0._dp
        Do j = 1,n-1,2  !only odd values of j, these are the NEW points at which to evaluate f
            s3 = s3 + f(a + Real(j,dp)*h)
        End Do
        q = h * (s1 + 2._dp*s2 + 4._dp*s3) * one_third
        If (Converged(q,qOld,rTol,aTol)) Return  !check convergence with specified tolerances
        qOld = q
    End Do
    !if we get this far, we failed to converge
    Write(*,'(A)') 'ERROR:  Quadratures: Composite_Simpson:  Failed to converge on ',n,' intervals.'
    ERROR STOP
End Function Composite_Simpson

Function GaussLegendre(f,a,b,n,wi,xi) Result(q)
    Use Kinds, Only: dp
    Implicit None
    Real(dp) :: q
    Interface
        Function f(x)  !the function to integrate
            Use Kinds, Only: dp    
            Implicit None
            Real(dp) :: f
            Real(dp), Intent(In) :: x
        End Function f
    End Interface
    Real(dp), Intent(In) :: a,b  !limits of integration
    Integer, Intent(In) :: n
    Real(dp), Intent(In) :: wi(1:n),xi(1:n)
    Real(dp) :: c1,c2  !changes limits of integration from (a,b) to (-1,1)
    Integer :: i
    Real(dp) :: fi(1:n)
    
    c1 = 0.5_dp * (b-a)
    c2 = 0.5_dp * (b+a)
    Do i = 1,n
        fi(i) = f(c1 * xi(i) + c2)
    End Do
    q = c1 * Dot_Product(wi,fi)
End Function GaussLegendre

Function GaussLegendre5(f,a,b) Result(q)
    Use Kinds, Only: dp
    Use GaussLegendre_w_a, Only: n => n005
    Use GaussLegendre_w_a, Only: wi => w005
    Use GaussLegendre_w_a, Only: xi => a005
    Implicit None
    Real(dp) :: q  !result of the integration
    Interface
        Function f(x)  !the function to integrate
            Use Kinds, Only: dp    
            Implicit None
            Real(dp) :: f
            Real(dp), Intent(In) :: x
        End Function f
    End Interface
    Real(dp), Intent(In) :: a,b  !limits of integration
    
    q = GaussLegendre(f,a,b,n,wi,xi)
End Function GaussLegendre5

Function GaussLegendre16(f,a,b) Result(q)
    Use Kinds, Only: dp
    Use GaussLegendre_w_a, Only: n => n016
    Use GaussLegendre_w_a, Only: wi => w016
    Use GaussLegendre_w_a, Only: xi => a016
    Implicit None
    Real(dp) :: q  !result of the integration
    Interface
        Function f(x)  !the function to integrate
            Use Kinds, Only: dp    
            Implicit None
            Real(dp) :: f
            Real(dp), Intent(In) :: x
        End Function f
    End Interface
    Real(dp), Intent(In) :: a,b  !limits of integration
    
    q = GaussLegendre(f,a,b,n,wi,xi)
End Function GaussLegendre16

Function GaussLegendre96(f,a,b) Result(q)
    Use Kinds, Only: dp
    Use GaussLegendre_w_a, Only: n => n096
    Use GaussLegendre_w_a, Only: wi => w096
    Use GaussLegendre_w_a, Only: xi => a096
    Implicit None
    Real(dp) :: q  !result of the integration
    Interface
        Function f(x)  !the function to integrate
            Use Kinds, Only: dp    
            Implicit None
            Real(dp) :: f
            Real(dp), Intent(In) :: x
        End Function f
    End Interface
    Real(dp), Intent(In) :: a,b  !limits of integration
    
    q = GaussLegendre(f,a,b,n,wi,xi)
End Function GaussLegendre96

Function GaussLegendreN(N,f,a,b) Result(q)
    Use Kinds, Only: dp
    Use GaussLegendre_w_a
    Implicit None
    Real(dp) :: q  !result of the integration
    Integer, Intent(In) :: N  !number of quadrature points
    Interface
        Function f(x)  !the function to integrate
            Use Kinds, Only: dp    
            Implicit None
            Real(dp) :: f
            Real(dp), Intent(In) :: x
        End Function f
    End Interface
    Real(dp), Intent(In) :: a,b  !limits of integration
    Real(dp), Allocatable :: wi(:),ai(:)  !weights and abscissa
    
    Allocate(wi(1:n))
    Allocate(ai(1:n))
    Select Case (N)
        Case (1)
            wi = w001
            ai = a001
        Case (2)
            wi = w002
            ai = a002
        Case (3)
            wi = w003
            ai = a003
        Case (4)
            wi = w004
            ai = a004
        Case (5)
            wi = w005
            ai = a005
        Case (6)
            wi = w006
            ai = a006
        Case (7)
            wi = w007
            ai = a007
        Case (8)
            wi = w008
            ai = a008
        Case (9)
            wi = w009
            ai = a009
        Case (10)
            wi = w010
            ai = a010
        Case (11)
            wi = w011
            ai = a011
        Case (12)
            wi = w012
            ai = a012
        Case (13)
            wi = w013
            ai = a013
        Case (14)
            wi = w014
            ai = a014
        Case (15)
            wi = w015
            ai = a015
        Case (16)
            wi = w016
            ai = a016
        Case (17)
            wi = w017
            ai = a017
        Case (18)
            wi = w018
            ai = a018
        Case (19)
            wi = w019
            ai = a019
        Case (20)
            wi = w020
            ai = a020
        Case (21)
            wi = w021
            ai = a021
        Case (22)
            wi = w022
            ai = a022
        Case (23)
            wi = w023
            ai = a023
        Case (24)
            wi = w024
            ai = a024
        Case (25)
            wi = w025
            ai = a025
        Case (26)
            wi = w026
            ai = a026
        Case (27)
            wi = w027
            ai = a027
        Case (28)
            wi = w028
            ai = a028
        Case (29)
            wi = w029
            ai = a029
        Case (30)
            wi = w030
            ai = a030
        Case (31)
            wi = w031
            ai = a031
        Case (32)
            wi = w032
            ai = a032
        Case (33)
            wi = w033
            ai = a033
        Case (34)
            wi = w034
            ai = a034
        Case (35)
            wi = w035
            ai = a035
        Case (36)
            wi = w036
            ai = a036
        Case (37)
            wi = w037
            ai = a037
        Case (38)
            wi = w038
            ai = a038
        Case (39)
            wi = w039
            ai = a039
        Case (40)
            wi = w040
            ai = a040
        Case (41)
            wi = w041
            ai = a041
        Case (42)
            wi = w042
            ai = a042
        Case (43)
            wi = w043
            ai = a043
        Case (44)
            wi = w044
            ai = a044
        Case (45)
            wi = w045
            ai = a045
        Case (46)
            wi = w046
            ai = a046
        Case (47)
            wi = w047
            ai = a047
        Case (48)
            wi = w048
            ai = a048
        Case (49)
            wi = w049
            ai = a049
        Case (50)
            wi = w050
            ai = a050
        Case (51)
            wi = w051
            ai = a051
        Case (52)
            wi = w052
            ai = a052
        Case (53)
            wi = w053
            ai = a053
        Case (54)
            wi = w054
            ai = a054
        Case (55)
            wi = w055
            ai = a055
        Case (56)
            wi = w056
            ai = a056
        Case (57)
            wi = w057
            ai = a057
        Case (58)
            wi = w058
            ai = a058
        Case (59)
            wi = w059
            ai = a059
        Case (60)
            wi = w060
            ai = a060
        Case (61)
            wi = w061
            ai = a061
        Case (62)
            wi = w062
            ai = a062
        Case (63)
            wi = w063
            ai = a063
        Case (64)
            wi = w064
            ai = a064
        Case (65)
            wi = w065
            ai = a065
        Case (66)
            wi = w066
            ai = a066
        Case (67)
            wi = w067
            ai = a067
        Case (68)
            wi = w068
            ai = a068
        Case (69)
            wi = w069
            ai = a069
        Case (70)
            wi = w070
            ai = a070
        Case (71)
            wi = w071
            ai = a071
        Case (72)
            wi = w072
            ai = a072
        Case (73)
            wi = w073
            ai = a073
        Case (74)
            wi = w074
            ai = a074
        Case (75)
            wi = w075
            ai = a075
        Case (76)
            wi = w076
            ai = a076
        Case (77)
            wi = w077
            ai = a077
        Case (78)
            wi = w078
            ai = a078
        Case (79)
            wi = w079
            ai = a079
        Case (80)
            wi = w080
            ai = a080
        Case (81)
            wi = w081
            ai = a081
        Case (82)
            wi = w082
            ai = a082
        Case (83)
            wi = w083
            ai = a083
        Case (84)
            wi = w084
            ai = a084
        Case (85)
            wi = w085
            ai = a085
        Case (86)
            wi = w086
            ai = a086
        Case (87)
            wi = w087
            ai = a087
        Case (88)
            wi = w088
            ai = a088
        Case (89)
            wi = w089
            ai = a089
        Case (90)
            wi = w090
            ai = a090
        Case (91)
            wi = w091
            ai = a091
        Case (92)
            wi = w092
            ai = a092
        Case (93)
            wi = w093
            ai = a093
        Case (94)
            wi = w094
            ai = a094
        Case (95)
            wi = w095
            ai = a095
        Case (96)
            wi = w096
            ai = a096
        Case (97)
            wi = w097
            ai = a097
        Case (98)
            wi = w098
            ai = a098
        Case (99)
            wi = w099
            ai = a099
        Case (100)
            wi = w100
            ai = a100
        Case (101)
            wi = w101
            ai = a101
        Case (102)
            wi = w102
            ai = a102
        Case (103)
            wi = w103
            ai = a103
        Case (104)
            wi = w104
            ai = a104
        Case (105)
            wi = w105
            ai = a105
        Case (106)
            wi = w106
            ai = a106
        Case (107)
            wi = w107
            ai = a107
        Case (108)
            wi = w108
            ai = a108
        Case (109)
            wi = w109
            ai = a109
        Case (110)
            wi = w110
            ai = a110
        Case (111)
            wi = w111
            ai = a111
        Case (112)
            wi = w112
            ai = a112
        Case (113)
            wi = w113
            ai = a113
        Case (114)
            wi = w114
            ai = a114
        Case (115)
            wi = w115
            ai = a115
        Case (116)
            wi = w116
            ai = a116
        Case (117)
            wi = w117
            ai = a117
        Case (118)
            wi = w118
            ai = a118
        Case (119)
            wi = w119
            ai = a119
        Case (120)
            wi = w120
            ai = a120
        Case (121)
            wi = w121
            ai = a121
        Case (122)
            wi = w122
            ai = a122
        Case (123)
            wi = w123
            ai = a123
        Case (124)
            wi = w124
            ai = a124
        Case (125)
            wi = w125
            ai = a125
        Case (126)
            wi = w126
            ai = a126
        Case (127)
            wi = w127
            ai = a127
        Case (128)
            wi = w128
            ai = a128
        Case (129)
            wi = w129
            ai = a129
        Case (130)
            wi = w130
            ai = a130
        Case (131)
            wi = w131
            ai = a131
        Case (132)
            wi = w132
            ai = a132
        Case (133)
            wi = w133
            ai = a133
        Case (134)
            wi = w134
            ai = a134
        Case (135)
            wi = w135
            ai = a135
        Case (136)
            wi = w136
            ai = a136
        Case (137)
            wi = w137
            ai = a137
        Case (138)
            wi = w138
            ai = a138
        Case (139)
            wi = w139
            ai = a139
        Case (140)
            wi = w140
            ai = a140
        Case (141)
            wi = w141
            ai = a141
        Case (142)
            wi = w142
            ai = a142
        Case (143)
            wi = w143
            ai = a143
        Case (144)
            wi = w144
            ai = a144
        Case (145)
            wi = w145
            ai = a145
        Case (146)
            wi = w146
            ai = a146
        Case (147)
            wi = w147
            ai = a147
        Case (148)
            wi = w148
            ai = a148
        Case (149)
            wi = w149
            ai = a149
        Case (150)
            wi = w150
            ai = a150
        Case Default
            Write(*,'(A,I0)') 'ERROR:  Quadratures: GaussLegendreN:  Unsupported N, n = ',N
            ERROR STOP
    End Select
    q = GaussLegendre(f,a,b,n,wi,ai)
End Function GaussLegendreN

Function Progressive_GaussLegendre(f,a,b,rtol,atol,n_start,n_stride,n_done) Result(q)
    Use Kinds, Only: dp
    Use Utilities, Only: Converged
    Implicit None
    Real(dp) :: q
    Interface
        Function f(x)
            Use Kinds, Only: dp
            Implicit None
            Real(dp) :: f
            Real(dp), Intent(In) :: x
        End Function
    End Interface
    Real(dp), Intent(In) :: a,b
    Real(dp), Intent(In) :: rtol,atol
    Integer, Intent(In), Optional :: n_start,n_stride
    Integer, Intent(Out), Optional :: n_done
    Real(dp) :: q_old
    Integer :: n,dn
    Integer :: i
    
    If (Present(n_start)) Then
        n = n_start
    Else
        n = 1
    End If
    If (Present(n_stride)) Then
        dn = n_stride
    Else
        dn = 1
    End If
    q_old = GaussLegendreN(n,f,a,b)
    Do i = n+dn,150,dn !Gauss-Legendre weights and abscissa are available up to 150 points
        q = GaussLegendreN(i,f,a,b)
        If (Converged(q,q_old,rtol,atol)) Then
            If (Present(n_done)) n_done = i
            RETURN !NORMAL EXIT
        End If
        q_old = q
    End Do
    !If we get this far, we failed to converge on 150-point Gauss
    If (Present(n_done)) n_done = -1
    Write(*,'(A,I0,A)') 'ERROR:  Quadratures: Progressive_GaussLegendre:  Failed to converge with ',i-dn,' Gauss-points.'
    ERROR STOP
End Function Progressive_GaussLegendre

End Module Quadratures
