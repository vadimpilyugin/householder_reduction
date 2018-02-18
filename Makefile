CC = mpicc
EXECUTABLES = build/main build/tests
FLAGS = -Wall -O3

all: mkdir $(EXECUTABLES)

mkdir: 
	mkdir -p build

# Object files
build/main.o: main.c
	$(CC) $(FLAGS) -c -o build/main.o main.c

build/matrix.o: matrix.c
	$(CC) $(FLAGS) -c -o build/matrix.o matrix.c

build/tests.o: tests.c
	$(CC) $(FLAGS) -c -o build/tests.o tests.c

# Executables
build/main: build/main.o build/matrix.o
	$(CC) $(FLAGS) -o $@ build/main.o build/matrix.o -lm

build/tests: build/tests.o build/matrix.o
	$(CC) $(FLAGS) -o $@ build/tests.o build/matrix.o -lm

# Phony targets
run: build/main
	# =================
	mpiexec -n 7 -oversubscribe ./build/main 1000

test: build/tests
	# =================
	mpiexec -n 7 -oversubscribe ./build/tests

clean:
	rm -rf *.o build $(EXECUTABLES)
