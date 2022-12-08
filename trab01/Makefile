CC := gcc
RM := rm -f
LDFLAGS := -lm
CFLAGS := -Ofast -fopenmp -Wall -Werror -pedantic


all:
	@$(CC) studentsseq.c -o studentsseq $(CFLAGS) $(LDFLAGS)
	@$(CC) studentspar.c -o studentspar $(CFLAGS) $(LDFLAGS)

seq:
	@./studentsseq

par:
	@./studentspar

clean:
	@$(RM) studentsseq studentspar


.PHONY: all clean
