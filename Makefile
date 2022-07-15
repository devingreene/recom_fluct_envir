platform := $(shell uname)

ifeq ($(platform),Darwin)
	MACOSX := 1
else
	MACOSX := 0
endif

exec: tree.c
	cc -DMACOSX=$(MACOSX) $(CFLAGS) -Wall -Wextra -O2 -o \
		exec tree.c parse_mut.c -lm -lgsl -lgslcblas
