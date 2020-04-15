Program testKinds

Use Kinds

Implicit None

Write(*,*)
Write(*,'(A,I2,A,F10.6,A)') 'SINGLE precision (',sp,' bytes) gives ',-Log10(2._sp*Spacing(1._sp)),' digits of precision near 1'
Write(*,'(A,I2,A,F10.6,A)') 'DOUBLE precision (',dp,' bytes) gives ',-Log10(2._dp*Spacing(1._dp)),' digits of precision near 1'
# if EXTP
    Write(*,'(A,I2,A,F10.6,A)') 'QUAD precision   (',qp,' bytes) gives ',-Log10(2._qp*Spacing(1._qp)),' digits of precision near 1'
# endif
Write(*,*)
Write(*,'(A,I2,A,ES11.3E4,A,ES11.3E4)') 'SINGLE precision (',sp,' bytes) gives HUGE = ',HUGE(1._sp),' and TINY = ',TINY(1._sp)
Write(*,'(A,I2,A,ES11.3E4,A,ES11.3E4)') 'DOUBLE precision (',dp,' bytes) gives HUGE = ',HUGE(1._dp),' and TINY = ',TINY(1._dp)
# if EXTP
    Write(*,'(A,I2,A,ES11.3E4,A,ES11.3E4)') 'QUAD precision   (',qp,' bytes) gives HUGE = ',HUGE(1._qp),' and TINY = ',TINY(1._qp)
# endif
Write(*,*)
Write(*,'(A,I2,A,I0)') 'SHORT integer  (',is,' bytes) gives range +/- ',HUGE(1_is)
Write(*,'(A,I2,A,I0)') 'LONG integer   (',il,' bytes) gives range +/- ',HUGE(1_il)
Write(*,'(A,I2,A,I0)') 'DOUBLE integer (',id,' bytes) gives range +/- ',HUGE(1_id)

End Program testKinds