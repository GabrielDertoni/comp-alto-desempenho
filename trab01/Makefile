CC := gcc
OMPFLAGS := -fopenmp
CFLAGS := -Wpedantic -Wall
FASTFLAGS := -Ofast -mavx2
DBGFLAGS := -fsanitize=address -g -DDEBUG
LDFLAGS := -lm

all: studentsseq studentspar

studentsseq: studentsseq.c
	$(CC) $(FASTFLAGS) $(OMPFLAGS) $(CFLAGS) $^ -o $@ $(LDFLAGS)

studentspar: studentspar.c
	$(CC) $(FASTFLAGS) $(OMPFLAGS) $(CFLAGS) $^ -o $@ $(LDFLAGS)

# seq $(nproc) | xargs -n 1 -I{} sh -c 'echo 1000 1000 1000 | OMP_NUM_THREADS={} ./studentspar | tail -n 1'

.PHONY: perf
