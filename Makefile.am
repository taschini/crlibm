## Process this file with automake to produce Makefile.in

lib_LIBRARIES = libcrlibm.a

# shared libraries need libtool. We'll do it later.
#lib_LTLIBRARIES = libcrlibm.la


# removed as of 0.10
#	log_fast.c log_fast.h log_accurate.h log_accurate.c 
#	exp_fast.h exp_fast.c exp_accurate.h exp_accurate.c 
#	log10_accurate.h log10_accurate.c 
#	log2_accurate.h  log2_accurate.c 

# IMPORTANT: be sure that the extra_dist below includes all files for both alternatives

if USE_HARDWARE_DE
libcrlibm_a_SOURCES =  \
	crlibm.h crlibm_private.h crlibm_private.c double-extended.h\
	exp-td.h exp-td.c \
	log-de.c log-de.h  \
	log10-td.h log10-td.c \
	log2-td.h  log2-td.c \
	rem_pio2_accurate.h rem_pio2_accurate.c \
	trigo_fast.c trigo_fast.h trigo_accurate.c trigo_accurate.h \
	atan_fast.c atan_fast.h atan_accurate.h atan_accurate.c \
	csh_fast.h csh_fast.c exp_accurate.h exp_accurate.c \
	asin-td.h asin-td.c expm1.h expm1.c log1p.h log1p.c
else 
libcrlibm_a_SOURCES =  \
	crlibm.h crlibm_private.h crlibm_private.c triple-double.h\
	exp-td.h exp-td.c \
	log-td.c log-td.h  \
	log10-td.h log10-td.c \
	log2-td.h  log2-td.c \
	rem_pio2_accurate.h rem_pio2_accurate.c \
	trigo_fast.c trigo_fast.h trigo_accurate.c trigo_accurate.h \
	atan_fast.c atan_fast.h atan_accurate.h atan_accurate.c \
	csh_fast.h csh_fast.c exp_accurate.h exp_accurate.c \
	asin-td.h asin-td.c expm1.h expm1.c log1p.h log1p.c
endif


# add all the scs_lib useful functions
libcrlibm_a_LIBADD = scs_lib/scs_private.o scs_lib/addition_scs.o  \
	scs_lib/division_scs.o scs_lib/print_scs.o\
	scs_lib/double2scs.o scs_lib/zero_scs.o\
	scs_lib/multiplication_scs.o scs_lib/scs2double.o

include_HEADERS = crlibm.h

SUBDIRS = scs_lib . tests

# Add your compiler here. Default is assumed to be gcc.
if COMPILER_ICC
# mp means maintain precision. Should be useless with std=c99, this is an ICC bug
AM_CFLAGS =  -std=c99 -mp -Qoption,cpp,--extended_float_types -IPF_fp_speculationsafe
else
AM_CFLAGS =  -std=c99 -fPIC -Wall -Wshadow -Wpointer-arith -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline 
endif



EXTRA_DIST = VERSION TODO\
	docs/crlibm.pdf \
	maple/common-procedures.mpl \
	maple/double-extended.mpl maple/triple-double.mpl\
	maple/exp-td.mpl maple/csh.mpl \
	maple/log-td.mpl  maple/log-de.mpl \
	maple/log2-td.mpl maple/log10-td.mpl \
	maple/trigo.mpl maple/atan.mpl \
	maple/asin-td.mpl \
	maple/sqrt.mpl \
	maple/expm1.mpl \
	maple/log1p.mpl \
	gappa/trigoSinCosCase3.gappa gappa/trigoTanCase2.gappa \
	gappa/trigoTanCase2.gappa \
	gappa/trigoSinCosCase3.gappa gappa/trigoTanCase2.gappa \
	gappa/log-de-acc-index0-E0.gappa    gappa/log-de-index0-E0.gappa \
	gappa/log-de-acc-index0-E1N.gappa   gappa/log-de-index0-E1N.gappa \
	gappa/log-de-acc-index1N-E0.gappa   gappa/log-de-index1N-E0N.gappa \
	gappa/log-de-acc-index1N-E1N.gappa  gappa/log-de.README \
	gappa/log-td.gappa gappa/log-td-E0.gappa gappa/log-td-E0-logir0.gappa \
	gappa/log-td-accurate.gappa gappa/log-td-accurate-E0.gappa gappa/log-td-accurate-E0-logir0.gappa \
	gappa/exp-td-accurate1.gappa gappa/exp-td-accurate2.gappa gappa/exp-td-accurate3.gappa gappa/exp-td-accurate4.gappa gappa/exp-td-proof.sh \
	docs/latex/crlibm.tex \
	docs/latex/0_intro.tex docs/latex/0_getting-started.tex\
	docs/latex/1_common.tex docs/latex/2_exp.tex docs/latex/3_log.tex\
	docs/latex/4_atan.tex docs/latex/5_csh.tex docs/latex/6_trigo.tex \
	docs/latex/fig_scs docs/latex/fig_exp\
	docs/latex/elem-fun.bib\
	scs_lib/tests/tbx_timing.h \
	log-de.c log-de.h double-extended.h log-td.c log-td.h triple-double.h atan-pentium.c atan-itanium.c \
	maple/asin-td.mpl expm1.h expm1.c log1p.h log1p.c
 
