CC=gcc
CFLAGS = -c -pg -O3 -Wall -msse -msse2 -g
LDFLAGS = -lz -lm  -pg
SOURCES = baseFAST.c CommandLineParser.c Common.c HashTable.c MrFAST.c Output.c Reads.c RefGenome.c Bitvectors.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = mrfast


all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS) 

.c.o:
	$(CC) $(CFLAGS) $< -o $@ 
clean:
	rm -f *.o *~ \#* mrfast
