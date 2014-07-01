CC = gcc
CFLAGS = -c -Wall -std=c99 -g
LDFLAGS = -lm
OBJECTS = main.o setup.o newdcdio.o pairdist.o covar.o pdbio.o


all: analysis

analysis: $(OBJECTS)
	$(CC) $(OBJECTS) -o analysis  $(LDFLAGS)

main.o: main.c 
	$(CC) $(CFLAGS) main.c -o main.o

setup.o: setup.c 
	$(CC) $(CFLAGS) setup.c -o setup.o

newdcdio.o: newdcdio.c 
	$(CC) $(CFLAGS) newdcdio.c -o newdcdio.o

pairdist.o: pairdist.c 
	$(CC) $(CFLAGS) pairdist.c -o pairdist.o

covar.o: covar.c
	$(CC) $(CFLAGS) covar.c -o covar.o

pdbio.o: pdbio.c
	$(CC) $(CFLAGS) pdbio.c -o pdbio.o

clean:
	rm -f *.o analysis
