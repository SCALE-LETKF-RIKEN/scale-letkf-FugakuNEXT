################################################################################
#
# Makefile for library unit test
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../..)
SYSDEP_DIR  = $(TOPDIR)/sysdep
BUILD_DIR  := $(PWD)/.libs
BINDIR      = $(PWD)

SCALE_DEBUG=T
export SCALE_DEBUG
SCALE_USE_AGRESSIVEOPT=F
export SCALE_USE_AGRESSIVEOPT

include $(SYSDEP_DIR)/Makedef.$(SCALE_SYS)
include $(TOPDIR)/Mkinclude

VPATH = $(BUILD_DIR)

BINNAME = unit
BINNAME_ORDER = unit_comm_order
BINNAME_CHECK = unit_comm_check
RUBY ?= ruby
MPIRUN ?= mpirun
OMP_NUM_THREADS ?= 1
VALIDATION_TIMEOUT ?= 120
VALIDATION_RUN = OMP_NUM_THREADS=$(OMP_NUM_THREADS) $(RUBY) $(TOPDIR)/scalelib/test/unit/run_with_timeout.rb $(VALIDATION_TIMEOUT) $(MPIRUN)
# Shell commands (for a recipe line) of an MPI job that must be rejected: it must stop with a
# non-zero status before the timeout (124 and 137 are the timeout), print the diagnostic in the
# shell variable expected, and have no failed assertion. The log is removed if it is rejected.
# $(1): log file, $(2): arguments of $(VALIDATION_RUN) (ranks, program, arguments), $(3): label of the job
EXPECT_REJECTION = status=0; $(VALIDATION_RUN) $(2) > "$(1)" 2>&1 || status=$$?; \
  if [ "$$status" -eq 0 ] || [ "$$status" -eq 124 ] || [ "$$status" -eq 137 ] || ! grep -Fq "$$expected" "$(1)" || grep -Fq '*** Error [Assert' "$(1)"; then \
    cat "$(1)"; echo "Expected $(3) to be rejected with its diagnostic"; exit 1; \
  fi; \
  rm -f "$(1)"

LIBS = $(LIBDIR)/libscale$(POSTFIX).a $(LIBDIR)/libdcutils.a

OBJS = \
	test_comm.o             \
	test_atmos_phy_mp_common.o \
	test_atmos_phy_tb_smg.o \
	test_atmos_dyn.o

all:
	$(MAKE) makedir
	$(MAKE) makelib
	$(MAKE) makebin
	$(MAKE) run

makedir:
	mkdir -p $(LIBDIR)
	mkdir -p $(BUILD_DIR)

makelib:
	$(MAKE) -C $(DCUTILSDIR)
	$(MAKE) -C $(SCALELIBDIR)/src

makebin: $(BINDIR)/$(BINNAME) $(BINDIR)/$(BINNAME_CHECK) $(BINDIR)/$(BINNAME_ORDER)

run: run-comm-check run-comm-order
	$(MPIRUN) -n 4 ./$(BINNAME) run.conf


# The checks of the variable ID in COMM_cartesC abort the run: each mode is a separate MPI job,
# which must fail with the diagnostic of the check (a short part of the message is searched).
# The valid uses of the IDs are tested in ./unit (test_comm). wait_other, other_grid,
# foreign_array, pc_same_vid, pc_dir, pc_released and pc_rel_nowait need persistent IDs: they are
# run with comm_check_pc.conf, which sets two-sided persistent communication explicitly (also where
# one-sided is the default).
# The checks are done in the public routines before the dispatch to the communication mode (and
# the check at COMM_finalize before the windows are freed), so they are not run with one-sided
# communication again.
# The checks of the groups of persistent variables (grp_*, issue #3594) need persistent IDs too:
# all are run with comm_check_pc.conf (the group is checked, its members exchanged alone), and the
# runtime checks also with comm_check_batch.conf (the group in batch). batch_nol2 is rejected in
# COMM_setup by the configuration (COMM_USE_BATCH without COMM_USE_PC_LAYOUT2), also without
# persistent communication.
# The modes of the IDs not bound to arrays (free_*, issue #3624) are run with comm_check_pc.conf, and
# the runtime checks of their group also with comm_check_batch.conf (in batch: COMM_USE_NONPC_BATCH).
# npb_nobatch is rejected in COMM_setup by the configuration (COMM_USE_NONPC_BATCH without
# COMM_USE_BATCH, comm_check_nonpc_batch_nobatch.conf).
# COMM_USE_PC_LAYOUT2 and COMM_USE_BATCH are .true. by default with OpenACC (issue #3615): there
# the runs with run.conf and the other files that do not set them (here, in run-average and in
# run) use the layout 2 and the batch, and the files named comm_*.conf set both explicitly.
COMM_CHECK_MODES = packbuf wait_nopost post_nowait wait_rank wait_array
COMM_CHECK_PC_MODES = wait_other other_grid foreign_array pc_same_vid pc_dir pc_released pc_rel_nowait
COMM_CHECK_GRP_RUN_MODES = grp_wait_mem grp_repost grp_rel_round grp_in_flight
COMM_CHECK_FREE_MODES = free_dir free_2d free_same_vid
COMM_CHECK_FREE_RUN_MODES = free_wait_arr free_wait_mem
COMM_CHECK_GRP_DEF_MODES = grp_same_pid grp_dup grp_4dir grp_range grp_rank_n grp_rank_id
COMM_CHECK_RUNS = $(foreach m,$(COMM_CHECK_MODES),run.conf:$(m)) $(foreach m,$(COMM_CHECK_PC_MODES),comm_check_pc.conf:$(m)) \
  $(foreach m,$(COMM_CHECK_GRP_RUN_MODES) $(COMM_CHECK_GRP_DEF_MODES),comm_check_pc.conf:$(m)) \
  $(foreach m,$(COMM_CHECK_GRP_RUN_MODES),comm_check_batch.conf:$(m)) \
  $(foreach m,$(COMM_CHECK_FREE_MODES) $(COMM_CHECK_FREE_RUN_MODES),comm_check_pc.conf:$(m)) \
  $(foreach m,$(COMM_CHECK_FREE_RUN_MODES),comm_check_batch.conf:$(m)) \
  comm_check_batch_nolayout2.conf:batch_nol2 comm_check_batch_nopc.conf:batch_nol2 \
  comm_check_nonpc_batch_nobatch.conf:npb_nobatch

run-comm-check:
	@set -eu; for run in $(COMM_CHECK_RUNS); do \
	  conf=$${run%%:*}; mode=$${run#*:}; \
	  case "$$mode" in \
	    packbuf|pc_same_vid) expected='already in use' ;; \
	    wait_nopost) expected='without a matching post' ;; \
	    wait_other)  expected='in use by another ID' ;; \
	    post_nowait) expected='posted but never waited' ;; \
	    wait_rank)   expected='rank of the array waited differs' ;; \
	    wait_array)  expected='not the array posted' ;; \
	    pc_dir)      expected='registered for the other direction' ;; \
	    pc_released) expected='vid is invalid: not initialized' ;; \
	    pc_rel_nowait) expected='still active at the release' ;; \
	    free_dir)    expected='registered for the other direction' ;; \
	    free_2d)     expected='vid is invalid' ;; \
	    free_wait_arr) expected='not the array posted' ;; \
	    free_same_vid) expected='two members share the packing buffer' ;; \
	    free_wait_mem) expected='waited before all its members are posted' ;; \
	    grp_wait_mem) expected='waited before all its members are posted' ;; \
	    grp_repost)  expected='posted again before all the members are waited' ;; \
	    grp_rel_round) expected='a round of the group is in progress at the release' ;; \
	    grp_in_flight) expected='a member is in an exchange in progress' ;; \
	    grp_same_pid) expected='two members share the packing buffer' ;; \
	    grp_dup)     expected='already in another group' ;; \
	    grp_4dir)    expected='registered for the 4-direction exchange' ;; \
	    grp_range)   expected='vid of a member is invalid' ;; \
	    grp_rank_n)  expected='the group differs among the processes' ;; \
	    grp_rank_id) expected='IDs of the members differ among the processes' ;; \
	    npb_nobatch) expected='COMM_USE_NONPC_BATCH = .true. requires COMM_USE_BATCH' ;; \
	    batch_nol2)  expected='COMM_USE_BATCH = .true. requires COMM_USE_PC_LAYOUT2' ;; \
	    *)           expected='not the one registered' ;; \
	  esac; \
	  $(call EXPECT_REJECTION,test_comm_check_$$mode.log,-n 4 ./$(BINNAME_CHECK) $$conf $$mode,$$mode ($$conf)); \
	  echo "COMM check $$mode rejected ($$conf)"; \
	done

# The order of the posts and waits of the exchanges in progress at the same time differs among the
# processes (issue #3594). Each run is layout:number of processes; all the cases are run for each of
# the communication configurations in COMM_ORDER_CONFS in turn in one job (they must use persistent
# communication), and the results with the first one are compared with those with the others on
# every cell; then all of them are set up again twice in turn with a single exchange.
# comm_order_pc_extend.conf starts with a capacity of two persistent variables, which is extended.
# comm_order_pc_layout2.conf uses the layout 2 of the persistent 8-direction exchange
# (COMM_USE_PC_LAYOUT2), with the same small capacity; the number of the persistent requests of
# each variable is checked against the layout that each configuration selects.
# comm_order_nonpc_layout2.conf uses the layout 2 of the non-persistent 8-direction exchange of 3-D
# variables (COMM_USE_NONPC_LAYOUT2, issue #3624) only, and comm_order_pc_batch.conf uses it too.
# comm_order_pc_batch.conf adds the batch of the groups of persistent variables (COMM_USE_BATCH),
# which is not used in 2-D (run_comm_2d.conf, run_comm_2d_edge.conf: the groups are checked only).
# In each configuration, the cases are run again with groups declared (COMM_vars_group_init),
# whose state (checked, or in batch) is checked, and their results are compared with those kept.
# run_comm_halo_x.conf has IHALO = 3 and JHALO = 2 (the other 3-D layouts have IHALO = JHALO), so
# that the halo widths in x and y taken one for the other are detected; it is periodic in x, so that
# every process has both the W and the E halos (and non-periodic in y).
# run_comm_barrier.conf adds the barriers of PROF_MPI_BARRIER, and skips the case whose posts and
# waits differ among the processes (it hangs with the barriers of the current code).
# A wrong order of the messages may hang, so each job has a timeout.
# In run_comm_large.conf (KA = 64) the messages exceed the eager limit of Open MPI (4 kB by default
# for shared memory), so they are sent with the rendezvous protocol: in single precision, a message
# of the 8-direction exchange is 12,288 B per row in y (N/S) and 24,576 B in x (W/E), one of the
# 4-direction exchange in y is 26,624 B, and a corner is 512 B (twice as large in double precision).
COMM_ORDER_TIMEOUT ?= 120
COMM_ORDER_RUN = OMP_NUM_THREADS=$(OMP_NUM_THREADS) $(RUBY) $(TOPDIR)/scalelib/test/unit/run_with_timeout.rb $(COMM_ORDER_TIMEOUT) $(MPIRUN)
COMM_ORDER_CONFS = comm_order_pc.conf comm_order_pc_extend.conf comm_order_pc_layout2.conf comm_order_nonpc_layout2.conf comm_order_pc_batch.conf
# run_comm_large.conf takes long, so it is run without comm_order_pc_extend.conf (the capacity is
# extended in comm_order_pc_layout2.conf too)
COMM_ORDER_CONFS_LARGE = comm_order_pc.conf comm_order_pc_layout2.conf comm_order_nonpc_layout2.conf comm_order_pc_batch.conf
COMM_ORDER_LAYOUTS = periodic:4 corner:4 corner_xy:4 empty_we:2 empty_all:1 2d:2 2d_edge:3 3columns:3 self_we:2 self_all:1 halo_x:4 barrier:4 large:4

run-comm-order:
	@set -eu; for run in $(COMM_ORDER_LAYOUTS); do \
	  layout=$${run%%:*}; ranks=$${run#*:}; \
	  if [ "$$layout" = large ]; then confs="$(COMM_ORDER_CONFS_LARGE)"; else confs="$(COMM_ORDER_CONFS)"; fi; \
	  $(COMM_ORDER_RUN) -n $$ranks ./$(BINNAME_ORDER) run_comm_$$layout.conf $$confs; \
	  echo "COMM order passed (run_comm_$$layout.conf, $$ranks processes)"; \
	done

$(BINDIR)/$(BINNAME)  : $(BUILD_DIR)/$(BINNAME).o  $(patsubst %,$(BUILD_DIR)/%,$(OBJS)) $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS) $(CONTRIB_LIBS)

$(BINDIR)/$(BINNAME_CHECK): $(BUILD_DIR)/test_comm_check.o $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS) $(CONTRIB_LIBS)

$(BUILD_DIR)/test_comm_check.o: test_comm_check.f90 $(LIBS)

$(BINDIR)/$(BINNAME_ORDER) : $(BUILD_DIR)/test_comm_order.o $(LIBS)
	$(LD) $(LDFLAGS) $(ADDITIONAL_FFLAGS) -o $@ $^ $(SCALE_NETCDF_LIBS) $(SCALE_MATHLIB_LIBS) $(SCALE_PAPI_LIBS) $(CONTRIB_LIBS)

$(BUILD_DIR)/test_comm_order.o: test_comm_order.f90 $(LIBS)

$(BUILD_DIR)/$(BINNAME).o  : $(BINNAME).f90  $(LIBS) $(patsubst %,$(BUILD_DIR)/%,$(OBJS))

$(BUILD_DIR)/test_comm.o:             test_comm.f90
$(BUILD_DIR)/test_atmos_phy_mp_common.o: test_atmos_phy_mp_common.f90
$(BUILD_DIR)/test_atmos_phy_tb_smg.o: test_atmos_phy_tb_smg.f90
$(BUILD_DIR)/test_atmos_dyn.o:        test_atmos_dyn.f90
$(BUILD_DIR)/unit.o:                  unit.f90 $(patsubst %,$(BUILD_DIR)/%,$(MODS))



allclean: distclean
	$(MAKE) -C $(DCUTILSDIR)      allclean
	$(MAKE) -C $(SCALELIBDIR)/src allclean
	rm -f $(BINDIR)/$(BINNAME) $(BINDIR)/$(BINNAME_CHECK) $(BINDIR)/$(BINNAME_ORDER)

distclean: clean
	rm -f $(BUILD_DIR)/$(BINNAME)
	rm -f LOG.pe000000 $(BINDIR)/test_comm_check_*.log

clean:
	rm -rf $(BUILD_DIR)/*.dSYM
	rm -f  $(BUILD_DIR)/*.*



.SUFFIXES:
.SUFFIXES: .o .f90 .mod

%.mod: %.f90
	make $(patsubst %.f90,%.o,$<)

$(BUILD_DIR)/%.o : %.f90
	$(FC) $(FFLAGS) $(ADDITIONAL_FFLAGS) -DVERSION_MACRO=\"$(VERSION)\" -I$(BUILD_DIR) -I$(DCUTILSDIR) -I$(MODDIR) -I$(SCALELIBDIR)/include -I$(SCALERMDIR)/include $(MODDIROPT) $(BUILD_DIR) -o $@ -c $<

.PHONY : clean distclean allclean dcutils scalelib

include $(TOPDIR)/utils/make/Make_environments
