## Process this file with automake to produce Makefile.in

bin_PROGRAMS =  crlibm_testval crlibm_soaktest crlibm_testperf

LDADD= $(top_builddir)/libcrlibm.a
INCLUDES=-I$(top_srcdir)
CFLAGS = -g -O0
../libcrlibm.a:
	cd $(top_srcdir); make libcrlibm.a 

crlibm_testval_SOURCES = test_val.c test_common.c test_common.h
crlibm_soaktest_SOURCES = soak_test.c test_common.c test_common.h
crlibm_testperf_SOURCES = test_perf.c test_common.c test_common.h
