CFLAGS += -Wall -Wextra -pedantic -O3
N ?= 13
SEED ?= 42
NP ?= 36


all: pvc-seq pvc-par

seq: pvc-seq
	@./pvc-seq $(N) $(SEED)

par: pvc-par
	@mpirun -np $(NP) --hostfile host.txt pvc-par $(N) $(SEED)

pvc-seq: pvc-seq.c
	@$(CC) pvc-seq.c -o pvc-seq $(CFLAGS)

pvc-par: pvc-par.c
	@mpicc pvc-par.c -o pvc-par $(CFLAGS)

clean:
	@rm -f pvc-seq pvc-par


.PHONY: seq par clean
