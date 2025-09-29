include ../configure

PROGS   = obsope obssim obsmake #dec_prepbufr

OBJS    = obsope_tools.o

OBJS_PAWR = radar_obs.o radar_tools.o dec_pawr_nml.o \
            read_toshiba$(TOSHIBA_OPT).o read_toshiba$(TOSHIBA_OPT)_f.o

ifeq ($(USE_PAWR_DECODER),T)
PROGS     += dec_pawr
endif

COMMON_DIR = ../common

COMMON_OBJS = $(COMMON_DIR)/SFMT.o           \
              $(COMMON_DIR)/common.o         \
              $(COMMON_DIR)/common_mpi.o     \
              $(COMMON_DIR)/common_ncio.o    \
              $(COMMON_DIR)/common_scale.o     \
              $(COMMON_DIR)/common_obs_scale.o \
              $(COMMON_DIR)/common_mpi_scale.o \
              $(COMMON_DIR)/common_nml.o

INCS    = $(SCALE_INC) $(NETCDF_INC) -I$(SCALE_RM_OBJDIR)

LIBS    = $(SCALE_LIB) $(NETCDF_LIB) $(LAPACK_LIB)

.PHONY:	all clean FORCE

all:	$(PROGS)

obsope:	obsope.o $(OBJS) $(COMMON_DIR)
	$(FC) $(FOPTS) $(FMCMODEL) -o obsope obsope.o \
	$(COMMON_OBJS) $(OBJS) $(SCALE_RM_OBJS) $(LIBS) -I$(SCALE_RM_OBJDIR)
obsmake:	obsmake.o $(OBS) $(COMMON_DIR)
	$(FC) $(FOPTS) $(FMCMODEL) -o obsmake obsmake.o \
	$(COMMON_OBJS) $(OBJS) $(SCALE_RM_OBJS) $(LIBS)
obssim:	obssim.o $(OBS) $(COMMON_DIR)
	$(FC) $(FOPTS) $(FMCMODEL) -o obssim obssim.o \
	$(COMMON_OBJS) $(OBJS) $(SCALE_RM_OBJS) $(LIBS)
dec_prepbufr:	dec_prepbufr.o $(COMMON_DIR)
	$(FC) $(FOPTS) -o dec_prepbufr dec_prepbufr.o $(COMMON_OBJS) \
	$(SCALE_RM_OBJS) $(LIBS) $(BUFR_LIB)
dec_pawr: dec_pawr.o $(COMMON_DIR)
	$(FC) $(FOPTS) $(FMCMODEL) -o dec_pawr dec_pawr.o \
	$(COMMON_OBJS) $(OBJS_PAWR) $(SCALE_RM_OBJS) $(LIBS)

obsope.o:	obsope.f90 $(OBJS) $(COMMON_DIR)
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) $(INCS)
obsmake.o:	obsmake.f90 $(OBJS) $(COMMON_DIR)
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) $(INCS)
obssim.o:	obssim.f90 $(OBJS) $(COMMON_DIR)
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) $(INCS)
dec_prepbufr.o:	dec_prepbufr.f90 $(COMMON_DIR)
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) $(INCS)

obsope_tools.o:	obsope_tools.f90 $(COMMON_DIR)
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) $(INCS)

dec_pawr.o:   pawr/dec_pawr.f90 radar_obs.o radar_tools.o dec_pawr_nml.o read_toshiba$(TOSHIBA_OPT)_f.o $(COMMON_DIR)
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) $(INCS)
dec_pawr_nml.o:  pawr/dec_pawr_nml.f90 
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) $(INCS)


radar_obs.o:  pawr/radar_obs.f90 radar_tools.o read_toshiba$(TOSHIBA_OPT)_f.o dec_pawr_nml.o  $(COMMON_DIR)
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) $(INCS)
radar_tools.o:  pawr/radar_tools.f90 dec_pawr_nml.o  $(COMMON_DIR)
	$(FC) $(FOPTS) -c $< -I$(COMMON_DIR) $(INCS)

read_toshiba$(TOSHIBA_OPT).o:   pawr/read_toshiba$(TOSHIBA_OPT).c pawr/read_toshiba$(TOSHIBA_OPT).h	
	$(CC) $(COPTS) -c $<
read_toshiba$(TOSHIBA_OPT)_f.o: pawr/read_toshiba$(TOSHIBA_OPT)_f.f90 read_toshiba$(TOSHIBA_OPT).o
	$(FC) $(FOPTS) -c $<

clean:
	rm -f *.o *.mod *.lst $(PROGS) obsope

FORCE:
