#
# Meant to be called as a sub-make; relies on variables set in ../../Makefile
#
RRTMGP_ROOT = ../../..
BUILD_DIR = $(RRTMGP_ROOT)/build
include $(BUILD_DIR)/Makefile.rules

VPATH = ./:$(BUILD_DIR):$(RRTMGP_ROOT)/extensions/:$(RRTMGP_ROOT)/test/util/src/

# RRTMGP library, module files
LDFLAGS   += -L$(BUILD_DIR)
LIBS      += -lrrtmgp -lrte
FCINCLUDE += -I$(BUILD_DIR)

# Extra sources -- extensions to RRTMGP classes,
# shared infrastructure, local sources
ADDITIONS = mo_load_coefficients.o mo_test_files_io.o mo_fluxes_byband.o mo_fluxes_byband_kernels.o

# Targets
test_flux_reduce: $(ADDITIONS) test_flux_reduce.o
test_flux_reduce.o: $(ADDITIONS) test_flux_reduce.F90

mo_fluxes_byband.o: mo_fluxes_byband_kernels.o mo_fluxes_byband.F90

clean:
	-rm ../test_flux_reduce *.o *.optrpt *.mod ../*.optrpt
