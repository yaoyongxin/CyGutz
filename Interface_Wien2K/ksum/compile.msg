rm -f *.o *_tmp_.* *.P .sequential .parallel *.mod
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c reallocate.f
reallocate.f(70): (col. 6) remark: LOOP WAS VECTORIZED.
reallocate.f(70): (col. 6) remark: LOOP WAS VECTORIZED.
reallocate.f(60): (col. 6) remark: LOOP WAS VECTORIZED.
reallocate.f(60): (col. 6) remark: LOOP WAS VECTORIZED.
reallocate.f(32): (col. 6) remark: LOOP WAS VECTORIZED.
reallocate.f(32): (col. 6) remark: LOOP WAS VECTORIZED.
reallocate.f(21): (col. 6) remark: LOOP WAS VECTORIZED.
reallocate.f(21): (col. 6) remark: LOOP WAS VECTORIZED.
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c modules.f
modules.f(63): (col. 5) remark: LOOP WAS VECTORIZED.
modules.f(130): (col. 4) remark: LOOP WAS VECTORIZED.
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c angle.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c abcd.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c atpar.f
atpar.f(29): (col. 10) remark: LOOP WAS VECTORIZED.
atpar.f(56): (col. 9) remark: LOOP WAS VECTORIZED.
atpar.f(106): (col. 10) remark: LOOP WAS VECTORIZED.
atpar.f(120): (col. 10) remark: LOOP WAS VECTORIZED.
atpar.f(134): (col. 10) remark: LOOP WAS VECTORIZED.
atpar.f(143): (col. 10) remark: LOOP WAS VECTORIZED.
atpar.f(147): (col. 10) remark: LOOP WAS VECTORIZED.
atpar.f(178): (col. 19) remark: LOOP WAS VECTORIZED.
atpar.f(193): (col. 12) remark: LOOP WAS VECTORIZED.
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c dvbes1.f
dvbes1.f(21): (col. 10) remark: LOOP WAS VECTORIZED.
dvbes1.f(14): (col. 10) remark: LOOP WAS VECTORIZED.
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c errclr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c errflg.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c gtfnam2.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c harmon.f
harmon.f(15): (col. 10) remark: LOOP WAS VECTORIZED.
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c inverss.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c l2main.f
l2main.f(47): (col. 2) remark: LOOP WAS VECTORIZED.
l2main.f(173): (col. 12) remark: LOOP WAS VECTORIZED.
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c qtlmain.f
qtlmain.f(214): (col. 12) remark: LOOP WAS VECTORIZED.
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c latgen.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c outerr.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c outwin.f
outwin.f(31): (col. 10) remark: LOOP WAS VECTORIZED.
outwin.f(60): (col. 11) remark: LOOP WAS VECTORIZED.
outwin.f(101): (col. 10) remark: LOOP WAS VECTORIZED.
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c ph.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c psplit.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c rint13.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c rotate.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c sphbes.f
sphbes.f(86): (col. 10) remark: LOOP WAS VECTORIZED.
sphbes.f(80): (col. 10) remark: LOOP WAS VECTORIZED.
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c symop.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c ylm.f
gcc  -c cputim.c
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c mknam.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c diracout.f
diracout.f(72): (col. 10) remark: LOOP WAS VECTORIZED.
diracout.f(87): (col. 7) remark: LOOP WAS VECTORIZED.
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c dergl.f
dergl.f(21): (col. 10) remark: PARTIAL LOOP WAS VECTORIZED.
dergl.f(21): (col. 10) remark: PARTIAL LOOP WAS VECTORIZED.
dergl.f(52): (col. 16) remark: LOOP WAS VECTORIZED.
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c xsplt.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c sym.f
sym.f(107): (col. 16) remark: LOOP WAS VECTORIZED.
sym.f(43): (col. 15) remark: LOOP WAS VECTORIZED.
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c outp.f
outp.f(29): (col. 9) remark: LOOP WAS VECTORIZED.
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c qtl.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c qtltext.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c radint.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c rout.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c lomain.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c readw.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c xsplt0.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c outdmat.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c inth.f
inth.f(22): (col. 10) remark: LOOP WAS VECTORIZED.
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c inouh.f
inouh.f(27): (col. 7) remark: LOOP WAS VECTORIZED.
inouh.f(56): (col. 7) remark: LOOP WAS VECTORIZED.
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c rotdef.f
rotdef.f(22): (col. 18) remark: LOOP WAS VECTORIZED.
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c transf.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c trafo.f
ifort  -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -c trafoso.f
trafoso.f(25): (col. 8) remark: PERMUTED LOOP WAS VECTORIZED.
ifort -o ./qtl reallocate.o modules.o angle.o abcd.o atpar.o dvbes1.o errclr.o errflg.o gtfnam2.o harmon.o inverss.o l2main.o qtlmain.o  latgen.o outerr.o outwin.o ph.o psplit.o rint13.o rotate.o sphbes.o symop.o ylm.o cputim.o mknam.o diracout.o dergl.o xsplt.o sym.o outp.o qtl.o qtltext.o radint.o rout.o lomain.o readw.o xsplt0.o outdmat.o inth.o inouh.o rotdef.o transf.o trafo.o trafoso.o -FR -mp1 -w -prec_div -pc80 -pad -ip -DINTEL_VML -traceback -L/opt/intel/mkl/10.0.011/lib/32 -lpthread -lmkl_lapack -lmkl -lguide -pthread 
