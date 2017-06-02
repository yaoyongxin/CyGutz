subroutine gcubatr_1(numfun, integrand, xrng, values, epsrel, epsabs, maxpts, iu)
  USE gprec; USE cui
  implicit none
  INTERFACE
    FUNCTION integrand(NUMFUN,X) RESULT(Value)
      USE gprec
      INTEGER, INTENT(IN) :: NUMFUN
      REAL(gq), DIMENSION(:), INTENT(IN) :: X
      REAL(gq), DIMENSION(NUMFUN) :: Value
    END FUNCTION integrand
  END INTERFACE
  integer, intent(in) :: numfun, maxpts, iu
  real(gq), intent(in) :: xrng(2), epsrel, epsabs
  real(gq), intent(out) :: values(numfun)
! local
  integer :: rgtype(1), neval
  real(gq) :: abserr(numfun), vertices(1, 0:1, 1)
!
  vertices(1, :, 1) = xrng
  rgtype = 1
  call cubatr(1, numfun, integrand, 1, vertices, rgtype, values, abserr, &
       &epsrel = epsrel, epsabs = epsabs, neval = neval, job = 2, maxpts = maxpts)
  if (iu >=0 )then
    write(iu, "( "" ESTIMATE OF ABSOLUTE ERROR: "")")
    write(iu, "(2X, 10es9.2)") abserr
    write(iu, "( "" NUMBER OF FUNCTION EVALATIONS = "", I10)") neval
  endif
  return
!
end subroutine gcubatr_1
