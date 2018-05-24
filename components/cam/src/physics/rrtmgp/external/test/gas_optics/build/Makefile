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
ADDITIONS = mo_test_files_io.o mo_load_coefficients.o

all: test_gas_optics test_gas_optics_options

test_gas_optics:   $(ADDITIONS) test_gas_optics.o
test_gas_optics.o: $(ADDITIONS) test_gas_optics.F90

test_gas_optics_options:   $(ADDITIONS) test_gas_optics_options.o
test_gas_optics_options.o: $(ADDITIONS) test_gas_optics_options.F90

clean:
	-rm ../test_gas_optics *.o *.optrpt ../*.optrpt *.mod
