geopk.o:geopk.F90
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS) $(FREEFLAGS) -fcray-pointer $<

micro_mg_data.o:micro_mg_data.F90
	$(FC) -c $(INCLDIR) $(INCS) $(FFLAGS_NOOPT) $(FREEFLAGS) -fcray-pointer $<
