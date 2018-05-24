RRTMGP_ROOT = ../../..
BUILD_DIR = $(RRTMGP_ROOT)/build
include $(BUILD_DIR)/Makefile.conf
include $(BUILD_DIR)/Makefile.rules

VPATH = ./:$(BUILD_DIR):$(RRTMGP_ROOT)/test/util/src/

# RRTMGP library, module files
LDFLAGS   += -L$(BUILD_DIR)
LIBS      += -lrrtmgp -lrte
FCINCLUDE += -I$(BUILD_DIR)

# Extra sources -- extensions to RRTMGP classes,
# shared infrastructure, local sources
ADDITIONS = mo_test_files_io.o

# Targets
all: test_adding
test_adding: $(ADDITIONS) test_adding.o
test_adding.o: $(ADDITIONS) test_adding.F90

clean:
	-rm ../test_adding *.o *.optrpt *.mod ../*.optrpt
