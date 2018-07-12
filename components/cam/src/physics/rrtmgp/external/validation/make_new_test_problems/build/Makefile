#
# Meant to be called as a sub-make; relies on variables set in ../../Makefile
#
RRTMGP_ROOT = ../../..
BUILD_DIR = $(RRTMGP_ROOT)/build
include $(BUILD_DIR)/Makefile.rules
%.o: %.f
	$(FC) $(F77FLAGS) $(FCINCLUDE) -c  $<

VPATH = ./:$(BUILD_DIR):$(RRTMGP_ROOT)/test/util/src:$(RRTMGP_ROOT)/validation/shdompp_solve/build

#
# RRTMGP library, module files
#
LDFLAGS   += -L$(BUILD_DIR)
LIBS      += -lrrtmgp -lrte
FCINCLUDE += -I$(BUILD_DIR)

#
# Extra sources -- extensions to RRTMGP classes, shared infrastructure, local sources
#
ADDITIONS = mo_test_files_io.o mo_load_coefficients.o fftpack.o mie_gamma_water.o
F77FLAGS += -O2

all: add_cloud_optics change_albedo_emis_sza

add_cloud_optics: add_cloud_optics.o $(ADDITIONS)
add_cloud_optics.o: $(ADDITIONS) add_cloud_optics.F90

change_albedo_emis_sza: change_albedo_emis_sza.o $(ADDITIONS)
change_albedo_emis_sza.o: $(ADDITIONS) change_albedo_emis_sza.F90

clean:
	-rm ../add_cloud_optics ../change_albedo_emis_sza *.o *.mod *.optrpt ../*.optrpt
