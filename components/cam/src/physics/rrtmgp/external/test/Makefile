#!/bin/make

include ../build/Makefile.conf ../test/Makefile.libs

# compile tests (all)
all:  test_adding \
	test_two_stream \
	test_sources \
	test_flux_reduce \
	test_gas_optics \
	test_flux_compute \
	test_lw_solver

# tests adding
test_adding:
	@echo --- building test adding ---
	@echo $(FCINCLUDE) $(LDFLAGS) $(LIBS)
	( cd ../test/adding/build ; make )

# tests two_stream
test_two_stream:
	@echo --- building test test_two_stream ---
	( cd ../test/two_stream/build ; make )

	# test source
test_sources:
	@echo --- building test test_sources ---
	( cd ../test/sources/build ; make )

# tests flux_reduce
test_flux_reduce:
	@echo --- building test flux_reduce ---
	( cd ../test/flux_reduce/build ; make )

# tests gas_optics
test_gas_optics:
	@echo --- building test gas_optics ---
	( cd ../test/gas_optics/build ; make )

# tests flux_compute
test_flux_compute:
	@echo --- building test flux_compute ---
	( cd ../test/flux_compute/build ; make )

# tests lw_solver
test_lw_solver:
	@echo --- building test lw_solver ---
	( cd ../test/lw_solver/build ; make )

# clean
clean:
	( cd ../test/adding/build ; make clean)
	( cd ../test/two_stream/build ; make clean )
	( cd ../test/sources/build ; make clean )
	( cd ../test/flux_reduce/build ; make clean )
	( cd ../test/gas_optics/build ; make clean )
	( cd ../test/flux_compute/build ; make clean )
	( cd ../test/lw_solver/build ; make clean )
