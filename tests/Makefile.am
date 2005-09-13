## Process this file with automake to produce Makefile.in

bin_PROGRAMS =  crlibm_testval crlibm_soaktest crlibm_testperf crlibm_blindtest crlibm_generate_test_vectors

LDADD= $(top_builddir)/libcrlibm.a
INCLUDES=-I$(top_srcdir)
AM_CFLAGS = -g -O0 -std=c99
../libcrlibm.a:
	cd $(top_srcdir); make libcrlibm.a 

crlibm_testval_SOURCES = test_val.c test_common.c test_common.h
crlibm_soaktest_SOURCES = soak_test.c test_common.c test_common.h
crlibm_testperf_SOURCES = test_perf.c test_common.c test_common.h
crlibm_blindtest_SOURCES = blind_test.c test_common.c test_common.h
crlibm_generate_test_vectors_SOURCES = generate_test_vectors.c test_common.c test_common.h

check_exp: crlibm_blindtest
	echo ./crlibm_blindtest exp.testdata>check_exp; chmod a+rx check_exp
check_log: crlibm_blindtest
	echo ./crlibm_blindtest log.testdata>check_log; chmod a+rx check_log
check_log2: crlibm_blindtest
	echo ./crlibm_blindtest log2.testdata>check_log2; chmod a+rx check_log2
check_log10: crlibm_blindtest
	echo ./crlibm_blindtest log10.testdata>check_log10; chmod a+rx check_log10
check_sin: crlibm_blindtest
	echo ./crlibm_blindtest sin.testdata>check_sin; chmod a+rx check_sin
check_cos: crlibm_blindtest
	echo ./crlibm_blindtest cos.testdata>check_cos; chmod a+rx check_cos
check_tan: crlibm_blindtest
	echo ./crlibm_blindtest tan.testdata>check_tan; chmod a+rx check_tan
check_atan: crlibm_blindtest
	echo ./crlibm_blindtest atan.testdata>check_atan; chmod a+rx check_atan
check_sinh: crlibm_blindtest
	echo ./crlibm_blindtest sinh.testdata>check_sinh; chmod a+rx check_sinh
check_cosh: crlibm_blindtest
	echo ./crlibm_blindtest cosh.testdata>check_cosh; chmod a+rx check_cosh

CRLIBM_CHECKLIST = check_exp check_log check_log2 check_log10  check_sin check_cos check_tan check_atan check_sinh check_cosh
TESTS = $(CRLIBM_CHECKLIST)

EXTRA_DIST = exp.testdata log.testdata log2.testdata log10.testdata sin.testdata cos.testdata tan.testdata atan.testdata sinh.testdata cosh.testdata

CLEANFILES = $(CRLIBM_CHECKLIST)