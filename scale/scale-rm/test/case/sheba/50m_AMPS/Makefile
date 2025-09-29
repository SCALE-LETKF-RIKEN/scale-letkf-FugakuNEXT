################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..


# user-defined source files
CODE_DIR    = ..
ORG_SRCS    = mod_user.F90

# parameters for run
INITCONF    = init.conf
RUNCONF     = run.conf
SNOCONF     = sno.conf
TPROC       = 10
TPROC_SNO   = 1

# required data (parameters,distributed files)
DATPARAM    = "[ ${SCALE_DB}/AMPS        AMPS_DATA ] \
               [ ../AMPSTASK.F           AMPSTASK.F ] \
	       [ ../env.txt              env.txt ] \
	       [ ../largeScaledForce.txt largeScaledForce.txt ] \
	       [ ../mipas_sheba.txt      mipas_sheba.txt ] \
	       [ $(TOPDIR)/data/mp/AMPS/seedini_0m seedini_0m ] \
	       [ ${TOPDIR}/data/rad      rad ]"
DATDISTS    =



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common
