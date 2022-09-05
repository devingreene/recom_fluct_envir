platform := $(shell uname)

ifeq ($(platform),Darwin)
	MACOSX := 1
else
	MACOSX := 0
endif

all: override CFLAGS += -O2
all: exec unittests

.PHONY: clean
clean:
	rm exec unittests

exec: parse_mut.c tree.c main.c recom_fluct_envir_many.h
	$(CC) -DMACOSX=$(MACOSX) $(CFLAGS) -Wall -Wextra -o \
		exec parse_mut.c tree.c main.c -lm -lgsl -lgslcblas

unittests: parse_mut.c tree.c unittests.c recom_fluct_envir_many.h
	$(CC) -DMACOSX=$(MACOSX) $(CFLAGS) -Wall -Wextra -o \
		unittests parse_mut.c tree.c unittests.c -lm -lgsl -lgslcblas

gdb: override CFLAGS += -ggdb
gdb: exec unittests
