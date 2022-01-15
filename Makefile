CC=gcc

all:
	$(CC) main.c aproksymator_na_bazie.c czebyszew_spl.c interpolator.c points.c prosta.c splines.c

aprox: main.o splines.o points.o aproksymator_na_bazie.o gaus/libge.a
	$(CC) -o aprox  main.o splines.o points.o aproksymator_na_bazie.o -L gaus -l ge

intrp: main.o splines.o points.o interpolator.o gaus/libge.a
	$(CC) -o intrp  main.o splines.o points.o interpolator.o -L gaus -l ge

prosta: main.o splines.o points.o prosta.o
	$(CC) -o prosta  main.o splines.o points.o prosta.o

czebyszew: main.o splines.o points.o czebyszew_spl.o gaus/libge.a
	$(CC) -o intrp  main.o splines.o points.o czebyszew_spl.o -L gaus -l ge

aproksymator_na_bazie.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c aproksymator_na_bazie.c

interpolator.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c interpolator.c

czebyszew_spl.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -I gaus -c czebyszew_spl.c

.PHONY: clean

clean:
	-rm *.o aprox intrp prosta
