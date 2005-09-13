## Process this file with automake to produce Makefile.in

lib_LIBRARIES = libcrlibm.a

# shared libraries need libtool. We'll do it later.
#lib_LTLIBRARIES = libcrlibm.la


# removed as of 0.10
#	log_fast.c log_fast.h log_accurate.h log_accurate.c \
#	exp_fast.h exp_fast.c exp_accurate.h exp_accurate.c \
#	log10_accurate.h log10_accurate.c \
#	log2_accurate.h  log2_accurate.c \


if HAS_HARDWARE_DE
libcrlibm_a_SOURCES =  \
	crlibm.h crlibm_private.h crlibm_private.c double-extended.h\
	exp-td.h exp-td.c \
	log-de.c log-de.h  \
	log10-td.h log10-td.c \
	log2-td.h  log2-td.c \
	rem_pio2_accurate.h rem_pio2_accurate.c \
	trigo_fast.c trigo_fast.h trigo_accurate.c trigo_accurate.h \
	atan_fast.c atan_fast.h atan_accurate.h atan_accurate.c \
	csh_fast.h csh_fast.c exp_accurate.h exp_accurate.c
else 
libcrlibm_a_SOURCES =  \
	crlibm.h crlibm_private.h crlibm_private.c\
	exp-td.h exp-td.c \
	log-td.c log-td.h  \
	log10-td.h log10-td.c \
	log2-td.h  log2-td.c \
	rem_pio2_accurate.h rem_pio2_accurate.c \
	trigo_fast.c trigo_fast.h trigo_accurate.c trigo_accurate.h \
	atan_fast.c atan_fast.h atan_accurate.h atan_accurate.c \
	csh_fast.h csh_fast.c exp_accurate.h exp_accurate.c
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
AM_CFLAGS =  -std=c99 -Wall  -mp
else
AM_CFLAGS =  -std=c99 -fPIC -Wall -Wshadow -Wpointer-arith -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline 
endif

EXTRA_DIST = VERSION\
	maple/common-procedures.mpl \
	maple/log.mpl maple/trigo.mpl maple/atan.mpl maple/csh.mpl \
	maple/trigoSinCosCase3.gappa maple/trigoTanCase2.gappa \
	docs/crlibm.pdf \
	docs/latex/crlibm.tex \
	docs/latex/0_intro.tex docs/latex/0_getting-started.tex\
	docs/latex/1_common.tex docs/latex/2_exp.tex docs/latex/3_log.tex\
	docs/latex/4_atan.tex docs/latex/5_csh.tex docs/latex/6_trigo.tex \
	docs/latex/fig_scs docs/latex/fig_exp\
	docs/latex/elem-fun.bib
 
