# CC := gcc
CC := clang
OMPFLAGS := -fopenmp
CFLAGS := -Wpedantic -Wall
FASTFLAGS := -Ofast -mavx2
DBGFLAGS := -fsanitize=address -g -DDEBUG
LDFLAGS := -lm

all: studentsseq studentspar

studentsseq: studentsseq.c
	# $(CC) -DPERF $(FASTFLAGS) $(CFLAGS) $^ -o $@ $(LDFLAGS)
	$(CC) -DPERF $(FASTFLAGS) $(OMPFLAGS) $(CFLAGS) $^ -o $@ $(LDFLAGS)

studentspar: studentspar.c
	$(CC) -DPERF -g $(FASTFLAGS) $(OMPFLAGS) $(CFLAGS) $^ -o $@ $(LDFLAGS)

studentsseq_v2: studentsseq_v2.c
	$(CC) -DPERF -g $(FASTFLAGS) $(CFLAGS) $^ -o $@ $(LDFLAGS)


.PHONY: perf
