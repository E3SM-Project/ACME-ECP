RRTMGP_ROOT = ../../..
BUILD_DIR = $(RRTMGP_ROOT)/build
include $(BUILD_DIR)/Makefile.rules

ECRAD_ROOT = /users/robert/PowerFolders/Other-peoples-codes/ecrad-1.0.1-orig
FCINCLUDE += -I$(ECRAD_ROOT)/mod

# Override default rules to put executables in directory above
%: %.o
	$(ECOMPILE) -o ../$@ $^ $(LDFLAGS) $(LIBS)

RTE_DIR = $(RRTMGP_ROOT)/rte
RTE_KERNEL_DIR = $(RRTMGP_ROOT)/rte/kernels
RRTMGP_DIR = $(RRTMGP_ROOT)/rrtmgp
RRTMGP_KERNEL_DIR = $(RRTMGP_ROOT)/rrtmgp/kernels
VPATH = $(RTE_DIR):$(RTE_KERNEL_DIR):$(RRTMGP_DIR):$(RRTMGP_KERNEL_DIR)
VPATH += $(RRTMGP_ROOT)/extensions/:$(RRTMGP_ROOT)/test/util/src/

RTE_SRC = \
  mo_rte_kind.o \
  mo_optical_props_kernels.o \
  mo_optical_props.o \
	mo_source_functions.o \
	mo_fluxes_broadband_kernels.o \
  mo_fluxes.o \
	mo_rte_lw_solver.o \
  mo_rte_lw.o \
  mo_rte_sw_solver.o \
  mo_rte_sw.o

#
# Extra sources -- extensions to RRTMGP classes, shared infrastructure, local sources
#
ADDITIONS = mo_test_files_io.o  mo_heating_rates.o mo_fluxes_byband_kernels.o mo_fluxes_byband.o mo_gas_concentrations.o mo_reorder_kernels.o mo_util_reorder.o mo_rrtmgp_constants.o
#
# radiation_io from ecrad; needed by radiation_flux.
#
ADDITIONS += radiation_io.o easy_netcdf.o
VPATH += $(ECRAD_ROOT)/utilities:$(ECRAD_ROOT)/radiation:$(ECRAD_ROOT)/ifsaux
#
# RRTMGP library, module files
#
LDFLAGS += -L$(ECRAD_ROOT)/lib -L.
LIBS    += -lrte -lradiation -lifsaux
#
# Targets
#
all: test_flux_compute_from_optics test_two_stream
#
# test_flux_compute_from_optics builds a local version of librte
#   which calls the homogenous solvers from ecrad.
# Build a local librte and link to existing ecrad libraries
#
test_flux_compute_from_optics:   test_flux_compute_from_optics.o   librte.a $(ADDITIONS)
test_flux_compute_from_optics.o: test_flux_compute_from_optics.F90 librte.a $(ADDITIONS)

clean:
	-rm *.o *.optrpt ../*.optrpt *.mod ../test_flux_compute* librte.a

.PHONY: radiation_two_stream.F90
##################################
#
# test_two_stream uses local bits of ecrad code
#   radiation_two_stream is modifed to allow mu0 to be specified for each element
#   (since ecrad is vectorized over g-point and rte vectorized over column)
#
test_two_stream:   test_two_stream.o   radiation_two_stream.o parkind1.o yomhook.o $(ADDITIONS)
test_two_stream.o: test_two_stream.F90 radiation_two_stream.o $(ADDITIONS)

# Small bits of ecrad source for testing two-stream
yomhook.o:                 parkind1.o yomhook.F90
radiation_two_stream.o: parkind1.o yomhook.o mo_optical_props.o mo_fluxes.o radiation_two_stream.F90
radiation_two_stream.F90:
	cp $(ECRAD_ROOT)/radiation/radiation_two_stream.F90 .
	patch radiation_two_stream.F90 < radiation_two_stream.F90.patch

##################################
# Build a local, bare-minimum copy of librte
librte.a: $(RTE_SRC)
	ar -rvs librte.a $(RTE_SRC)

mo_rte_kind.o:         mo_rte_kind.F90
mo_rrtmgp_constants.o: mo_rte_kind.o mo_rrtmgp_constants.F90
#
# Optical properties
#
mo_optical_props_kernels.o: mo_rte_kind.o mo_optical_props_kernels.F90
mo_optical_props.o:         mo_rte_kind.o mo_optical_props_kernels.o   mo_optical_props.F90

mo_source_functions.o:      mo_rte_kind.o mo_optical_props.o mo_source_functions.F90
#
# Flux reduction
#
mo_fluxes_broadband_kernels.o: mo_rte_kind.o mo_fluxes_broadband_kernels.F90
mo_fluxes.o:        mo_rte_kind.o mo_fluxes_broadband_kernels.o mo_fluxes.F90
#
# Radiative transfer
#
mo_rte_lw_solver.o: mo_rte_kind.o \
	mo_optical_props.o \
	mo_rte_lw_solver.F90

mo_rte_sw_solver.o: mo_rte_kind.o \
	mo_optical_props.o \
	mo_rte_sw_solver.F90

mo_rte_lw.o: mo_rte_kind.o \
	mo_optical_props.o  \
	mo_source_functions.o \
	mo_fluxes.o \
	mo_rte_lw_solver.o \
	mo_rte_lw.F90

mo_rte_sw.o: mo_rte_kind.o \
	mo_optical_props.o  \
	mo_fluxes.o \
	mo_rte_sw_solver.o \
	mo_rte_sw.F90

##################################
# Test file infrastructure; this relies on bits of rrtmgp
mo_fluxes_byband.o: mo_fluxes_byband_kernels.o librte.a

mo_test_files_io.o: mo_rte_kind.o \
	mo_optical_props.o  \
	mo_gas_concentrations.o \
	mo_util_reorder.o \
	mo_test_files_io.F90

mo_util_reorder.o:  mo_rte_kind.o \
	mo_reorder_kernels.o \
	mo_util_reorder.F90

mo_gas_contentrations.o:  mo_rte_kind.o \
	mo_gas_contentrations.F90

mo_heating_rates.o:  mo_rte_kind.o \
	mo_rrtmgp_constants.o \
	mo_heating_rates.F90

mo_rrtmgp_constants.o:  mo_rte_kind.o \
	mo_rrtmgp_constants.F90
