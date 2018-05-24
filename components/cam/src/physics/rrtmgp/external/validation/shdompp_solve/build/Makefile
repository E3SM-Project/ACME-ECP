RRTMGP_ROOT = ../../..
BUILD_DIR = $(RRTMGP_ROOT)/build
include $(BUILD_DIR)/Makefile.conf
include $(BUILD_DIR)/Makefile.rules

# Override default rules
%: %.o
	$(ECOMPILE) -o ../$@ $^ $(LDFLAGS) $(LIBS)
%.o: %.f
	$(FC) $(F77FLAGS) $(FCINCLUDE) -c  $<

VPATH = ./:$(BUILD_DIR):$(RRTMGP_ROOT)/test/util/src:$(RRTMGP_ROOT)/extensions

#
# RRTMGP library, module files
#
LDFLAGS   += -L$(BUILD_DIR)
LIBS      += -lrrtmgp -lrte
FCINCLUDE += -I$(BUILD_DIR)

#
# Extra sources -- extensions to RRTMGP classes, shared infrastructure, local sources
#
ADDITIONS = mo_test_files_io.o mo_load_coefficients.o fftpack.o mo_fluxes_byband_kernels.o mo_fluxes_byband.o

F77FLAGS += -O2

shdompp_solve: shdompp_solve.o shdompp_rrtm.o $(ADDITIONS)

shdompp_solve.o: shdompp_rrtm.o $(ADDITIONS)

mo_fluxes_byband.o: mo_fluxes_byband_kernels.o mo_fluxes_byband.F90

clean:
	-rm ../shdompp_solve ../change_albedo_emis_sza *.o *.mod *.optrpt ../*.optrpt
