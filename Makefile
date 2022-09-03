platform := $(shell uname)

ifeq ($(platform),Darwin)
	MACOSX := 1
else
	MACOSX := 0
endif

exec: tree.c parse_mut.c main.c this.h
	cc -DMACOSX=$(MACOSX) $(CFLAGS) -Wall -Wextra -O2 -o \
		exec tree.c parse_mut.c main.c -lm -lgsl -lgslcblas

unittests: tree.c parse_mut.c unittests.c this.h
	cc -DMACOSX=$(MACOSX) $(CFLAGS) -Wall -Wextra -O2 -o \
		unittests tree.c parse_mut.c unittests.c -lm -lgsl -lgslcblas
