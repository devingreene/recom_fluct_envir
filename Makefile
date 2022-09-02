platform := $(shell uname)

ifeq ($(platform),Darwin)
	MACOSX := 1
else
	MACOSX := 0
endif

exec: parse_mut.c tree.c main.c this.h
	cc -DMACOSX=$(MACOSX) $(CFLAGS) -Wall -Wextra -O2 -o \
		exec parse_mut.c tree.c main.c -lm -lgsl -lgslcblas

unittests: parse_mut.c tree.c unittests.c this.h
	cc -DMACOSX=$(MACOSX) $(CFLAGS) -Wall -Wextra -O2 -o \
		unittests parse_mut.c tree.c unittests.c -lm -lgsl -lgslcblas
