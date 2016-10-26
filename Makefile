CC=g++
CFLAGS=-c -m64 -O3 -Wall -static
all: entess

entess: EnTess.o apstring.o Tetran.o Atom.o Bond.o Molecule.o dTess.o Dataset.o Descriptor.o
	$(CC) EnTess.o apstring.o Atom.o Bond.o Molecule.o Dataset.o Tetran.o dTess.o Descriptor.o -m64 -o entess

apstring.o: apstring.cpp
	$(CC) $(CFLAGS) apstring.cpp

Tetran.o: Tetran.cpp
	$(CC) $(CFLAGS) Tetran.cpp

Atom.o: Atom.cpp
	$(CC) $(CFLAGS) Atom.cpp

Bond.o: Bond.cpp
	$(CC) $(CFLAGS) Bond.cpp

Molecule.o: Molecule.cpp
	$(CC) $(CFLAGS) Molecule.cpp

dTess.o: dTess.cpp
	$(CC) $(CFLAGS) dTess.cpp

Dataset.o: Dataset.cpp
	$(CC) $(CFLAGS) Dataset.cpp

Descriptor.o: Descriptor.cpp
	$(CC) $(CFLAGS) Descriptor.cpp

clean:
	rm -rf *.o
