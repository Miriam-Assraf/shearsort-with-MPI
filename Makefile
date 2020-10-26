build:
	mpicxx -fopenmp -c main.c -o main.o 
	mpicxx -fopenmp -c sorting_funcs.c -o sorting_funcs.o -lm
	mpicxx -fopenmp -c input_output_funcs.c -o input_output_funcs.o 
	mpicxx -fopenmp -o prog main.o sorting_funcs.o input_output_funcs.o

clean:
	rm -f *.o ./prog

run:
	mpiexec -np 16 ./prog

