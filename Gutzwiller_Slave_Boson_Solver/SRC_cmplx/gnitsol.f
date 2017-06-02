!**************************************************
! http://users.wpi.edu/~walker/NITSOL/readme.txt
!**************************************************
      subroutine gnitsol(f,n,x,ftol,stptol,input)
      implicit none
      integer n,input(10)
      real(8) x(n),ftol,stptol
      external f, dummy_jacv
! local
      integer info(6), ipar(1), iterm
      real(8) rpar(1)
!
      call nitsol(n, x, f, dummy_jacv, ftol, stptol, input, info, rpar, ipar, iterm)
      if (iterm.ne.0) then
        write(6,'(" WARNING: ERROR IN gnitsol with iterm = ",I2)') iterm
      endif
      write(6,'(" number of J*v evaluations      = ", I0)') info(2)
      write(6,'(" number of linear iterations    = ", I0)') info(4)
      write(6,'(" number of nonlinear iterations = ", I0)') info(5)
      write(6,'(" number of backtracks           = ", I0)') info(6)
      return
!
      end subroutine gnitsol
!
!**************************************************
      subroutine dummy_jacv()
      end subroutine dummy_jacv
