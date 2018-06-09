RRTMGP_ROOT = ../../..
BUILD_DIR = $(RRTMGP_ROOT)/build
include $(BUILD_DIR)/Makefile.conf
include $(BUILD_DIR)/Makefile.rules

VPATH = ./:$(BUILD_DIR):$(RRTMGP_ROOT)/test/util/src/:$(RRTMGP_ROOT)/extensions

#
# RRTMGP library, module files
#
LDFLAGS   += -L$(BUILD_DIR)
LIBS      += -lrrtmgp -lrte
FCINCLUDE += -I$(BUILD_DIR)

#
# Extra sources -- extensions to RRTMGP classes, shared infrastructure, local sources
#
ADDITIONS = mo_test_files_io.o mo_load_coefficients.o mo_fluxes_bygpoint.o

#
# Targets
#

test_lw_solver: $(ADDITIONS) test_lw_solver.o

test_lw_solver.o: $(ADDITIONS) test_lw_solver.F90


clean:
	-rm ../test_lw_solver . *.o *.optrpt ../*.optrpt *.mod
