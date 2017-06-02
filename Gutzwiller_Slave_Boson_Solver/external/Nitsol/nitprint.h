
      integer iplvl, ipunit
      common /nitprint/ iplvl, ipunit

c
c If diagnostic information is desired, include this common block in the 
c main program and set iplvl and ipunit according to the following: 
c
c     iplvl = 0 => no printout
c           = 1 => iteration numbers and F-norms
c           = 2 => ... + some stats, step norms, and linear model norms
c           = 3 => ... + some Krylov solver and backtrack information
c           = 4 => ... + more Krylov solver and backtrack information
c
c     ipunit = printout unit number, e.g., ipunit = 6 => standard output. 
c              NOTE: If ipunit = 0 on input, then it is set to 6.
c
