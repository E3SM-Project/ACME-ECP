
include make.inc

FC := ${FC}
CC := ${CC}
FREEFLAGS := ${FREEFLAGS}
FIXEDFLAGS := ${FIXEDFLAGS}
FFLAGS := ${FFLAGS} ${CPPDEFS}
CFLAGS := ${CFLAGS} ${CPPDEFS}

SOURCES := $(shell cat Srcfiles)
OBJS    := $(addsuffix .o, $(basename $(SOURCES)))

all: standalone

standalone: crm_standalone.o task_init.o crmsurface.o kurant.o abcoefs.o zero.o buoyancy.o forcing.o damping.o ice_fall.o boundaries.o coriolis.o pressure.o adams.o periodic.o task_util_NOMPI.o uvw.o diagnose.o sat.o bound_exchange.o utils.o press_rhs.o fft.o bound_duvdt.o press_grad.o diffuse_mom.o diffuse_mom2D.o diffuse_mom3D.o diffuse_scalar.o diffuse_scalar2D.o diffuse_scalar3D.o tke_full.o shear_prod2D.o advect_mom.o advect2_mom_z.o advect2_mom_xy.o advect_all_scalars.o advect_scalar.o advect_scalar2D.o
	${FC} -o crm_standalone *.o ${LDFLAGS}

clean:
	rm -f *.mod *.o crm_standalone

.SUFFIXES:
.SUFFIXES: .F .F90 .f90 .c .o

%.o: %.mod

.F.o:
	$(FC) ${INCLUDE} -c $(FIXEDFLAGS) $(FFLAGS) $<

.F90.o:
	$(FC) ${INCLUDE} -c $(FREEFLAGS) $(FFLAGS) $<

.f90.o:
	$(FC) ${INCLUDE} -c $(FREEFLAGS) $(FFLAGS) $<

.c.o:
	$(CC) ${INCLUDE} -c $(CFLAGS) $<

include Depends




