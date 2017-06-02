FAQ
===

* How to check the convergence of the Gutzwiller-slave-boson (CyGutz) solver?

The main output file of CyGutz at the current electron density iteration 
is `GUTZ.LOG`. The convergence can be checked by typing::

    grep -i maxerr GUTZ.LOG

The main output file of of CyGutz at the previous electron density iteration 
is `GL_LOG.OUT`. 

* Where do I quickly check the renormalization factor matrix 
  :math:`Z=R^\dagger R`, which is an indicator of electron correlation strength?

Open file `GUTZ.LOG` or `GL_LOG.OUT` and search for 
`eigen values of r\\\\daggerr`. The last entry is the final result. 

* What should I do if CyGutz fails to converge?

CyGutz usually converges for systems away from Mott transition 
by using the default initial guess of the solutions, 
which corresponds to the solution of zero interaction. 
For the `LDA+Gutzwiller slave boson` calculations, 
the initial guess of CyGutz at the electron density iteration `n (n>1)` 
takes the solution of CyGutz at the electron density iteration `n-1` 
by default. 
This usually works well. 
However, there are cases that such process would fail. 
In other words, the `LDA+Gutzwiller slave boson` calculation can converge 
in the first few electron density iterations, 
and suddenly it fails in the next iteration. 
This could happen for systems with low symmetries or 
with large dimension of the solution vector. 
In that case, one should check whether the system is 
close to Mott transition, 
by looking the eigen-values of :math:`Z`-matrix. 
If all the eigen-values are much larger than 0, 
say 0.2 or bigger, there is a good chance that CyGutz can converge 
with different initial guess. 
One good try is to use the default initial guess by typing::

    rm WH_RLNEF.INP

which removes the file storing the previous CyGutz solutions. If the system is indeed close to Mott transition, one should investigate the solutions of the Mott localized phases by typing::

    init_mott.py

and choose the appropriate options accordingly. 
