CFLAGS += -Wall -Wextra -pedantic -O3
N ?= 13
SEED ?= 42


all:
	@$(CC) pvc-seq.c -o pvc-seq $(CFLAGS)
	@mpicc pvc-par.c -o pvc-par $(CFLAGS)

seq:
	@./pvc-seq $(N) $(SEED)

par:
	@mpirun --hostfile host.txt pvc-par $(N) $(SEED)

clean:
	@rm -f pvc-seq pvc-par


.PHONY: seq par clean
