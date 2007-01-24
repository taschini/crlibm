## Process this file with automake to produce Makefile.in


lib_LIBRARIES = libscs.a
INCLUDES=-I$(top_srcdir)


libscs_a_SOURCES = scs.h scs_private.h scs_private.c \
	double2scs.c print_scs.c \
	division_scs.c\
	addition_scs.c multiplication_scs.c\
	scs2double.c zero_scs.c \
	scs2mpfr.c 


include_HEADERS = scs.h ../crlibm_config.h

