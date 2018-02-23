include ../../config.inc

SOURCES=../dreadMM.c ../Util.c dlinsolx.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=dlinsolx

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LD_SUPERLU) $(OBJECTS) -o $@

.c.o:
	$(CC) $(CFLAGS) $(I_SUPERLU) $< -o $@

.PHONY: clean
clean:
	rm -f *.o $(EXECUTABLE)
