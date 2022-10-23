OMPFLAGS := -fopenmp
CFLAGS := -Wpedantic -Wall
FASTFLAGS := -Ofast
DBGFLAGS := -fsanitize=address -g -DDEBUG
LDFLAGS := -lm

all: studentsseq_v2

# WARNING: REMOVER -DNO_OUT ANTES DE ENTREGAR O TRABALHO
studentsseq: studentsseq.c
	# gcc $(CFLAGS) $(FASTFLAGS) $^ -o $@ $(LDFLAGS)
	# gcc $(CFLAGS) $(DBGFLAGS)  $^ -o $@ $(LDFLAGS)
	gcc -DPERF -Ofast $(CFLAGS) -g -gdwarf-3 $^ -o $@ $(LDFLAGS)

studentsseq_v2: studentsseq_v2.c
	# gcc $(CFLAGS) $(FASTFLAGS) $^ -o $@ $(LDFLAGS)
	# gcc $(CFLAGS) $(DBGFLAGS)  $^ -o $@ $(LDFLAGS)
	gcc -DPERF -g -Ofast -mavx2 $(CFLAGS) $^ -o $@ $(LDFLAGS)
	# gcc -Ofast -mavx2 $(CFLAGS) $^ -o $@ $(LDFLAGS)


.PHONY: perf
