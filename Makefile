CC=gcc
CFLAGS=-O3 -Wall 
LDFLAGS=

SOURCES=main.c sim.c init.c terminate.c double_vector.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=dpd_sim

all: $(SOURCES) $(EXECUTABLE)

# $@ is the name of file being generated
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@
	
clean:
	rm *.o $(EXECUTABLE)

.PHONY:
	clean