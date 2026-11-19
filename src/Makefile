# Makefile for EarthCARE Driver

BINDIR = .
BIN  = $(BINDIR)/ecare_driver

OBJS = parameters.o             \
       numbers.o                \
	   pdfs.o               \
       input.o                  \
	   netcdf_writer.o      \
	   summary_processor.o  \
       output.o                 \
       ecare_driver.o

FC     = ifort

### Compile Options:
LDFLAGS = -I/usr/include -I/usr/local/include -L/usr/local/lib
LDLIBS = -lz -lm -ldf -lsz -ljpeg -lc -lpthread -lhdf5_hl -lhdf5_fortran -lhdf5  -lnetcdff -lnetcdf
FFLAGS = -mcmodel=large -O2 -g -traceback
HDF_INC = /usr/local/include
HDF_LIB = /usr/local/lib
NF_INC = /usr/local/include
NF_LIB = /usr/local/lib
INC = /usr/include
LIB = /usr/lib
HDF5_LDFLAGS = -I$(INC) -I$(HDF_INC) -L$(LIB)
NETCDF_LDFLAGS = -I${NF_INC} -L${NF_LIB}

### Dependency:
$(BIN): $(OBJS)
	$(FC) -o $@ $(OBJS) $(FFLAGS) $(LDFLAGS) $(LDLIBS)

parameters.o: parameters.f90
	$(FC) -c $< $(FFLAGS) 

numbers.o: numbers.f90
	$(FC) -c $< $(FFLAGS) 

pdfs.o: pdfs.f90
	$(FC) -c $< $(FFLAGS)

input.o: input.f90 parameters.o
	$(FC) -c $< $(FFLAGS) $(HDF5_LDFLAGS)

netcdf_writer.o: netcdf_writer.f90 parameters.o
	$(FC) -c $< $(FFLAGS) $(NETCDF_LDFLAGS)

summary_processor.o: summary_processor.f90 parameters.o netcdf_writer.o
	$(FC) -c $< $(FFLAGS) $(NETCDF_LDFLAGS)

output.o: output.f90 parameters.o numbers.o pdfs.o summary_processor.o
	$(FC) -c $< $(FFLAGS) $(NETCDF_LDFLAGS)

ecare_driver.o: ecare_driver.f90 parameters.o numbers.o input.o output.o
	$(FC) -c $< $(FFLAGS) $(HDF5_LDFLAGS) $(NETCDF_LDFLAGS)

clean :
	rm -f $(BIN) core *.o *.mod
