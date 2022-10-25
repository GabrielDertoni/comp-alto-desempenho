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

.PHONY: perf
