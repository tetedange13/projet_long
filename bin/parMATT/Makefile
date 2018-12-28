CC=mpicc
BIN=bin/parMatt

CFLAGS = -pthread\
	-ggdb\
	-Ithird_party\
	-Iinclude\
	-std=gnu99\
	-fopenmp\
	-O3

MATT = $(wildcard src/matt/*.c) 

MPI = $(wildcard src/mpi/*.c) 

obj = $(MATT:.c=.o) $(MPI:.c=.o)

$(shell mkdir -p bin)

parmatt: $(obj)
	$(CC) -o $(BIN) $^ $(CFLAGS) -lm


.PHONY: clean

clean:
	rm -f $(obj) bin/parMatt



