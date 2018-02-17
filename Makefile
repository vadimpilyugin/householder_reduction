EXECUTABLES = main

all: $(EXECUTABLES)

main: main.c
	mpicc -o $@ main.c

run: main
	# =================
	./main

clean:
	rm -f $(EXECUTABLES)