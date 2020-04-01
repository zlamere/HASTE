Program testKinds

Use Kinds

Implicit None

Write(*,*)
Write(*,'(A,I0,A,F10.6,A)') 'SINGLE precision (',sp,' bytes) gives ',-Log10(2._sp*Spacing(1._sp)),' digits of precision near 1.'
Write(*,'(A,I0,A,F10.6,A)') 'DOUBLE precision (',dp,' bytes) gives ',-Log10(2._dp*Spacing(1._dp)),' digits of precision near 1.'
!Write(*,'(A,I0,A,F10.6,A)') 'QUAD precision   (',qp,' bytes) gives ',-Log10(2._qp*Spacing(1._qp)),' digits of precision near 1.'

End Program testKinds