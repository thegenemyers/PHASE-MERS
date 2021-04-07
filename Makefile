DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

CC = clang

ALL = Phasemer

all: $(ALL)

libfastk.c : gene_core.c
libfastk.h : gene_core.h

Phasemer: Phasemer.c libfastk.c libfastk.h
	$(CC) $(CFLAGS) -o Phasemer Phasemer.c libfastk.c -lpthread -lm

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f Phasing.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf Phasing.tar.gz LICENSE README.md Makefile *.h *.c
