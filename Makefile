DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

CC = clang

ALL = Phasemer Phasemap

all: $(ALL)

libfastk.c : gene_core.c
libfastk.h : gene_core.h

Phasemer: Phasemer.c libfastk.c libfastk.h
	$(CC) $(CFLAGS) -o Phasemer Phasemer.c libfastk.c -lpthread -lm

Phasemap: Phasemap.c libfastk.c libfastk.h reader.c sam.c sam.h DB.c DB.h QV.c QV.h
	$(CC) $(CFLAGS) -o Phasemap Phasemap.c libfastk.c reader.c sam.c DB.c QV.c -lpthread -lm -lz

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f Phasing.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf Phasing.tar.gz LICENSE README.md Makefile *.h *.c
