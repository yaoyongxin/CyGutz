PF90 = mpiifort

FFLAGS = -O0 -check bound -check pointer -check arg_temp_created \
		-traceback


FREE_FORM = -free


obj = t3j.o  Vcalc.o


all: ${obj} main.o
	${PF90} ${FREE_FORM} ${FFLAGS} -o main.out main.o ${obj} 


%.o: ./%.f
	${PF90} ${FREE_FORM} ${FFLAGS} -c $< 

clean:
	rm -f *o *out 
