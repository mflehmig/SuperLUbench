include ../../config.inc

SOURCES=../dreadMM.c ../Util.c pdlinsolx.c 
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=pdlinsolx

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LD_SUPERLUMT)

.c.o:
	$(CC) $(CFLAGS) $(I_SUPERLUMT) $< -o $@

.PHONY: clean
clean:
	rm -f *.o $(EXECUTABLE)
