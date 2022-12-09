CFLAGS += -Wall -Wextra -pedantic -O3
PARAM = 13 42


all:
	@$(CC) pvc-seq.c -o pvc-seq $(CFLAGS)
	@mpicc pvc-par.c -o pvc-par $(CFLAGS)

seq:
	@./pvc-seq $(PARAM)

par:
	@mpirun --hostfile host.txt pvc-par $(PARAM)

clean:
	@rm -f pvc-seq pvc-par


.PHONY: seq par clean
