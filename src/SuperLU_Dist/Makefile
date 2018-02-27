include ../../config.inc

SOURCES=../Util.c pddrive_ABglobal.c 
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=pddrive_ABglobal

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(MPICC) $(LD_SUPERLUDIST) $(OBJECTS) -o $@

.c.o:
	$(MPICC) $(CFLAGS) $(I_SUPERLUDIST) $< -o $@

.PHONY: clean
clean:
	rm -f *.o $(EXECUTABLE)
