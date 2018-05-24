RRTMGP_ROOT = ../../..
BUILD_DIR = $(RRTMGP_ROOT)/build
include $(BUILD_DIR)/Makefile.conf
include $(BUILD_DIR)/Makefile.rules

VPATH = ./:$(BUILD_DIR):$(RRTMGP_ROOT)/extensions/:$(RRTMGP_ROOT)/test/util/src/

#
# RRTMGP library, module files
#
LDFLAGS   += -L$(BUILD_DIR)
LIBS      += -lrrtmgp -lrte
FCINCLUDE += -I$(BUILD_DIR)

#
# Extra sources -- extensions to RRTMGP classes, shared infrastructure, local sources
#
ADDITIONS = mo_test_files_io.o mo_load_coefficients.o mo_heating_rates.o mo_fluxes_byband.o mo_fluxes_byband_kernels.o mo_rrtmgp_clr_all_sky.o

#
# Targets
#
all: test_flux_compute test_flux_compute_orig_interface test_flux_compute_from_optics

test_flux_compute: $(ADDITIONS) test_flux_compute.o

test_flux_compute.o: $(ADDITIONS) test_flux_compute.F90

test_flux_compute_orig_interface:   $(ADDITIONS) test_flux_compute_orig_interface.o

test_flux_compute_orig_interface.o: $(ADDITIONS) test_flux_compute_orig_interface.F90

test_flux_compute_from_optics:   $(ADDITIONS) test_flux_compute_from_optics.o

test_flux_compute_from_optics.o: $(ADDITIONS) test_flux_compute_from_optics.F90

mo_fluxes_byband.o: mo_fluxes_byband.F90 mo_fluxes_byband_kernels.o

clean:
	-rm ../test_flux_compute *.o *.optrpt ../*.optrpt *.mod
