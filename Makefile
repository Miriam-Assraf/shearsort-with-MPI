build:
	mpicxx -g -o ex2 ex2.c -lm

clean:
	rm -f *.o ./ex2

run:
	mpiexec -np 16 ./ex2

