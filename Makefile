
.PHONY: all

SUBDIRS = $(shell ls -d src/SuperLU*)
all:
	echo "Subdirs = $(SUBDIRS)..."
	for dir in $(SUBDIRS) ; do \
		make -C  $$dir ; \
	done
	
clean:
	echo "Subdirs = $(SUBDIRS)..."
	for dir in $(SUBDIRS) ; do \
		make -C  $$dir clean; \
	done

