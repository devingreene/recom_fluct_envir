exec: tree.c
	cc -Wall -Wextra -O2 -o exec tree.c -lm -lgsl -lgslcblas
