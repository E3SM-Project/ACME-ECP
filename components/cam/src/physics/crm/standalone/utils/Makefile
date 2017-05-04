
include make.inc

FC := ftn
CC := cc
FREEFLAGS := -Mextend
FFLAGS := -O2 ${CPPDEFS}
CFLAGS := -O2 ${CPPDEFS}

SOURCES := $(shell cat Srcfiles)
OBJS    := $(addsuffix .o, $(basename $(SOURCES)))

all: libcrm

libcrm: crm_module.o
	ar rcs libcrm.a *.o

clean:
	rm -f *.mod *.o

.SUFFIXES:
.SUFFIXES: .F .F90 .f90 .c .s .o

%.o: %.mod

.F.o:
	$(FC) -c $(FIXEDFLAGS) $(FFLAGS) $<

.F90.o:
	$(FC) -c $(FREEFLAGS) $(FFLAGS) $<

.f90.o:
	$(FC) -c $(FREEFLAGS) $(FFLAGS) $<

.c.o:
	$(CC) -c $(CFLAGS) $<

.s.o:
	$(AS) -m  $<

include Depends




