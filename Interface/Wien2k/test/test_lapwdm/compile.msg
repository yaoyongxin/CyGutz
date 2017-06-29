rm -f *.o _tmp_.* *.P .real .complex *.mod
if [ -f .complex ]; then \
	   make clean; \
	fi
touch .real
make TYPE='REAL' TYPE_COMMENT='!_REAL' ./lapwdm
make[1]: Entering directory `/home/ykent/OPT/apps/WIEN2K14.2/SRC_lapwdm'
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c reallocate.f
modules.frc: REAL version extracted
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c modules_tmp.f
mv modules_tmp.o modules.o
rm modules_tmp.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c abc.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c angle.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c atpar.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c dvbes1.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c dmat.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c errclr.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c errflg.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c euler.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c diracout.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c dergl.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c gtfnam2.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c harmon.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c inispl.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c inverss.f
l2main.frc: REAL version extracted
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c l2main_tmp.f
mv l2main_tmp.o l2main.o
rm l2main_tmp.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c lmtmat.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c lapwdm.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c latgen.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c outerr.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c outwin.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c orb.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c ph.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c rint13.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c inouh.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c inth.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c rotate.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c rotdef.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c sphbes.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c symop.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c t3j.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c SymmRot.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c symmetrization.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c ylm.f
gcc    -c -o cputim.o cputim.c
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c mknam.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c xsplt.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c locmat.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c timeinv.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c sym.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c readw.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c radint.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c rout.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c lomain.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c spin.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c couplx.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c printx.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c xoper.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c xqproduct.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c symmso.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c addtinvso.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c addtinv.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c output.f
mpiifort -o ./lapwdm reallocate.o modules.o abc.o angle.o atpar.o dvbes1.o dmat.o errclr.o errflg.o euler.o diracout.o dergl.o gtfnam2.o harmon.o inispl.o inverss.o l2main.o lmtmat.o lapwdm.o  latgen.o outerr.o outwin.o orb.o ph.o rint13.o inouh.o inth.o rotate.o rotdef.o sphbes.o  symop.o t3j.o SymmRot.o symmetrization.o ylm.o cputim.o mknam.o xsplt.o locmat.o timeinv.o sym.o readw.o radint.o rout.o lomain.o spin.o couplx.o printx.o xoper.o xqproduct.o symmso.o addtinvso.o addtinv.o output.o   -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -L/opt/intel/composer_xe_2015.2.164/mkl/lib/ -pthread -lfftw3 -L/home/ykent/OPT/lib/fftw-3.3.3-intel/lib -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread
make[1]: Leaving directory `/home/ykent/OPT/apps/WIEN2K14.2/SRC_lapwdm'
if [ -f .real ]; then \
	   make clean; \
	fi
make[1]: Entering directory `/home/ykent/OPT/apps/WIEN2K14.2/SRC_lapwdm'
rm -f *.o _tmp_.* *.P .real .complex *.mod
make[1]: Leaving directory `/home/ykent/OPT/apps/WIEN2K14.2/SRC_lapwdm'
touch .complex
make TYPE='COMPLEX' TYPE_COMMENT='!_COMPLEX' ./lapwdmc
make[1]: Entering directory `/home/ykent/OPT/apps/WIEN2K14.2/SRC_lapwdm'
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c reallocate.f
modules.frc: COMPLEX version extracted
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c modules_tmp.f
mv modules_tmp.o modules.o
rm modules_tmp.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c abc.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c angle.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c atpar.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c dvbes1.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c dmat.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c errclr.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c errflg.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c euler.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c diracout.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c dergl.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c gtfnam2.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c harmon.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c inispl.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c inverss.f
l2main.frc: COMPLEX version extracted
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c l2main_tmp.f
mv l2main_tmp.o l2main.o
rm l2main_tmp.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c lmtmat.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c lapwdm.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c latgen.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c outerr.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c outwin.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c orb.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c ph.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c rint13.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c inouh.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c inth.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c rotate.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c rotdef.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c sphbes.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c symop.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c t3j.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c SymmRot.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c symmetrization.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c ylm.f
gcc    -c -o cputim.o cputim.c
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c mknam.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c xsplt.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c locmat.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c timeinv.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c sym.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c readw.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c radint.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c rout.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c lomain.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c spin.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c couplx.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c printx.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c xoper.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c xqproduct.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c symmso.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c addtinvso.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c addtinv.f
mpiifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -c output.f
mpiifort -o ./lapwdmc reallocate.o modules.o abc.o angle.o atpar.o dvbes1.o dmat.o errclr.o errflg.o euler.o diracout.o dergl.o gtfnam2.o harmon.o inispl.o inverss.o l2main.o lmtmat.o lapwdm.o  latgen.o outerr.o outwin.o orb.o ph.o rint13.o inouh.o inth.o rotate.o rotdef.o sphbes.o  symop.o t3j.o SymmRot.o symmetrization.o ylm.o cputim.o mknam.o xsplt.o locmat.o timeinv.o sym.o readw.o radint.o rout.o lomain.o spin.o couplx.o printx.o xoper.o xqproduct.o symmso.o addtinvso.o addtinv.o output.o  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -assume buffered_io -I/opt/intel/mkl/include -DFFTW3 -I/home/ykent/OPT/lib/fftw-3.3.3-intel/include -L/opt/intel/composer_xe_2015.2.164/mkl/lib/ -pthread -lfftw3 -L/home/ykent/OPT/lib/fftw-3.3.3-intel/lib -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread
make[1]: Leaving directory `/home/ykent/OPT/apps/WIEN2K14.2/SRC_lapwdm'
