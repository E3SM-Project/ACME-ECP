RRTMGP_ROOT = ../../..
BUILD_DIR = $(RRTMGP_ROOT)/build
include $(BUILD_DIR)/Makefile.conf
include $(BUILD_DIR)/Makefile.rules

VPATH = ./:$(BUILD_DIR):$(RRTMGP_ROOT)/test/util/src/

#
# RRTMGP library, module files
#
LDFLAGS   += -L$(BUILD_DIR)
LIBS      += -lrrtmgp -lrte
FCINCLUDE += -I$(BUILD_DIR)

#
# Extra sources -- extensions to RRTMGP classes, shared infrastructure, local sources
#
ADDITIONS = mo_test_files_io.o

#
# Targets
#
test_sources: $(ADDITIONS) test_sources.o

test_sources.o: $(ADDITIONS) test_sources.F90

clean:
	-rm ../test_sources *.o *.optrpt ../*.optrpt *.mod
