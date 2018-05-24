#!/bin/make

include ../build/Makefile.conf

# compile tests (all)
all:  make_new_test_problems shdompp_solve

# tests adding
make_new_test_problems:
	@echo --- building make_new_test_problems ---
	( cd ../validation/make_new_test_problems/build ; make )

shdompp_solve:
	@echo --- building shdompp_solve ---
	( cd ../validation/shdompp_solve/build ; make )

ecrad:
	@echo --- building ecrad ---
	( cd ../validation/ecrad/build ; make)

# clean
clean:
	( cd ./make_new_test_problems/build ; make clean)
	( cd ./shdompp_solve/build ; make clean)
	( cd ./ecrad/build ; make clean)
	( rm *.pdf; cd data; rm solver* shdompp* rrtmgp* ecrad* one-clear* garand*)

.PHONY: make_new_test_problems shdompp_solve ecrad
