CC=gcc


aprox: main.o splines.o points.o aproksymator_na_bazie.o gaus/libge.a
	$(CC) -ggdb -o aprox  main.o splines.o points.o aproksymator_na_bazie.o -L gaus -l ge

intrp: main.o splines.o points.o interpolator.o gaus/libge.a
	$(CC) -ggdb -o intrp  main.o splines.o points.o interpolator.o -L gaus -l ge

prosta: main.o splines.o points.o prosta.o
	$(CC) -ggdb  -o prosta  main.o splines.o points.o prosta.o

czebysz: main.o splines.o points.o czebysz.o gaus/libge.a
	$(CC) -ggdb -o czebysz  main.o splines.o points.o czebysz.o -L gaus -l ge

aproksymator_na_bazie.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -ggdb -I gaus -c aproksymator_na_bazie.c

interpolator.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -ggdb -I gaus -c interpolator.c

czebysz.o: makespl.h points.h gaus/piv_ge_solver.h
	$(CC) -ggdb -I gaus -c czebysz.c

.PHONY: clean


clean:
	-rm *.o aprox intrp prosta czebysz main


