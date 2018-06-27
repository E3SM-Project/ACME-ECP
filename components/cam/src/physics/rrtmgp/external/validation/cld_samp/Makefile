#
# Meant to be called as a sub-make; relies on variables set in ../../Makefile
#
BUILD_DIR = ../../build
VPATH = ../../test/util/io/

#
# RRTMGP library, module files
#
LDFLAGS   += -L$(BUILD_DIR)
LIBS      += -lrrtmgp
FCINCLUDE += -I$(BUILD_DIR)

COMPILE = $(FC) $(FCFLAGS) $(FCINCLUDE) -c
LINK    = $(FC) $(FCFLAGS) $(LDFLAGS) $(LIBS)

%.o: %.F90
	$(COMPILE) $<

%: %.o
	$(LINK) -o $@ $^

test_cld_samp: test_cld_samp.o
