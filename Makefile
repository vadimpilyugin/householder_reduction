CC = mpicc
EXECUTABLES = main tests

all: $(EXECUTABLES)

# Object files
main.o: main.c
	$(CC) -c -o main.o main.c

matrix.o: matrix.c
	$(CC) -c -o matrix.o matrix.c

tests.o: tests.c
	$(CC) -c -o tests.o tests.c

# Executables
main: main.o matrix.o
	$(CC) -o $@ main.o matrix.o

tests: tests.o matrix.o
	$(CC) -o $@ tests.o matrix.o

# Phony targets
run: main
	# =================
	mpiexec -n 7 ./main

test: tests
	# =================
	mpiexec -n 7 ./tests

clean:
	rm -f *.o $(EXECUTABLES)
