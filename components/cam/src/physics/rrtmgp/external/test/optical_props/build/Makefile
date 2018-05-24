BUILD_DIR = ../../../build
include $(BUILD_DIR)/Makefile.conf

VPATH = ../../../src:../../../build:../../../test/util/io/

#
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
test_optical_props: $(ADDITIONS) test_optical_props.o

test_optical_props.o: $(ADDITIONS) test_optical_props.F90

clean:
	-rm ../test_optical_props *.o *.optrpt ../*.optrpt *.mod
