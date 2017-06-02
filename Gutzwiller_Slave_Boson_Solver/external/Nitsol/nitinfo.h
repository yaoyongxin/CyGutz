
      integer instep, newstep, krystat
      double precision avrate, fcurnrm
      common /nitinfo/ avrate, fcurnrm, instep, newstep, krystat

c
c If information on the current state of the nonlinear iteration is
c desired in a user-supplied subroutine (for example, deciding 
c whether to update a preconditioner), include this common block
c in the subroutine. The variables are as follows: 
c
c     instep - inexact Newton step number. 
c
c    newstep - set to 0 at the beginning of an inexact Newton step.
c              This may be checked in a user-supplied jacv to decide
c              whether to update the preconditioner.  If you test on
c              newstep .eq. 0 to determine whether to take some 
c              special action at the beginning of a nonlinear iteration, 
c              you must also set newstep to some nonzero value to
c              subsequently avoid taking that action unnecessarily. 
c
c    krystat - status of the Krylov iteration; same as itrmks (see 
c              the nitsol documentation). 
c
c    avrate  - average rate of convergence of the Krylov solver during
c              the previous inexact Newton step.  This may be checked
c              in a user-supplied jacv to decide when to update the
c              preconditioner.
c
c    fcurnrm - ||f(xcur)||. 
c
