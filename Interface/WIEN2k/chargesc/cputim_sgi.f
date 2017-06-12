      SUBROUTINE CPUTIM(TIME)
      real*8 time                                                  
      time=second()
      end
      function second(dummy)
      real tarray(2)
      second=etime(tarray)
      end
      SUBROUTINE WALLTIM(DSEC)
      REAL*8 DSEC
      DSEC=0.0D0
      RETURN
      END

