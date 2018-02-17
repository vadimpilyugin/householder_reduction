CC = mpicc
EXECUTABLES = main

all: $(EXECUTABLES)

# Object files
main.o: main.c
	$(CC) -c -o main.o main.c

matrix.o: matrix.c
	$(CC) -c -o matrix.o matrix.c

# Executables
main: main.o matrix.o
	$(CC) -o $@ main.o matrix.o

run: main
	# =================
	mpiexec -n 7 ./main

clean:
	rm -f *.o $(EXECUTABLES)
