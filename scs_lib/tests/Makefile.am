## Process this file with automake to produce Makefile.in

bin_PROGRAMS =  scs_test_timing scs_test_accuracy

LDADD= $(top_builddir)/libscs.a
INCLUDES=-I$(top_srcdir)

../libscs.a:
	cd $(top_srcdir); make libscs.a 

scs_test_timing_SOURCES = test_timing.c tbx_timing.h

scs_test_accuracy_SOURCES = test_accuracy.c