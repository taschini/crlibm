## Process this file with automake to produce Makefile.in

AM_CFLAGS = -O0 -std=gnu99 # produces a warning, but I really need -O0 and I don't know how to get it.

if USE_INTERVAL_FUNCTIONS
bin_PROGRAMS =  crlibm_testval crlibm_soaktest crlibm_soaktest_interval crlibm_testperf crlibm_testperf_interval crlibm_blindtest crlibm_generate_test_vectors
else
bin_PROGRAMS =  crlibm_testval crlibm_soaktest crlibm_testperf  crlibm_blindtest crlibm_generate_test_vectors
endif

LDADD= $(top_builddir)/libcrlibm.a
INCLUDES=-I$(top_srcdir)
../libcrlibm.a:
	cd $(top_srcdir); make libcrlibm.a 

crlibm_testval_SOURCES = test_val.c test_common.c test_common.h
crlibm_soaktest_SOURCES = soak_test.c test_common.c test_common.h
crlibm_soaktest_interval_SOURCES = soak_test-interval.c test_common.c test_common.h
crlibm_testperf_SOURCES = test_perf.c test_common.c test_common.h  powmidpoint.h

crlibm_testperf_interval_SOURCES = test_perf-interval.c test_common.c test_common.h
crlibm_blindtest_SOURCES = blind_test.c test_common.c test_common.h
crlibm_generate_test_vectors_SOURCES = generate_test_vectors.c test_common.c test_common.h

check_exp: crlibm_blindtest
	echo ./crlibm_blindtest exp.testdata>check_exp; chmod a+rx check_exp
check_expm1: crlibm_blindtest
	echo ./crlibm_blindtest expm1.testdata>check_expm1; chmod a+rx check_expm1
check_log: crlibm_blindtest
	echo ./crlibm_blindtest log.testdata>check_log; chmod a+rx check_log
check_log2: crlibm_blindtest
	echo ./crlibm_blindtest log2.testdata>check_log2; chmod a+rx check_log2
check_log10: crlibm_blindtest
	echo ./crlibm_blindtest log10.testdata>check_log10; chmod a+rx check_log10
check_log1p: crlibm_blindtest
	echo ./crlibm_blindtest log1p.testdata>check_log1p; chmod a+rx check_log1p
check_sin: crlibm_blindtest
	echo ./crlibm_blindtest sin.testdata>check_sin; chmod a+rx check_sin
check_cos: crlibm_blindtest
	echo ./crlibm_blindtest cos.testdata>check_cos; chmod a+rx check_cos
check_tan: crlibm_blindtest
	echo ./crlibm_blindtest tan.testdata>check_tan; chmod a+rx check_tan
check_asin: crlibm_blindtest
	echo ./crlibm_blindtest asin.testdata>check_asin; chmod a+rx check_asin
check_acos: crlibm_blindtest
	echo ./crlibm_blindtest acos.testdata>check_acos; chmod a+rx check_acos
check_atan: crlibm_blindtest
	echo ./crlibm_blindtest atan.testdata>check_atan; chmod a+rx check_atan
check_asinpi: crlibm_blindtest
	echo ./crlibm_blindtest asinpi.testdata>check_asinpi; chmod a+rx check_asinpi
check_acospi: crlibm_blindtest
	echo ./crlibm_blindtest acospi.testdata>check_acospi; chmod a+rx check_acospi
check_atanpi: crlibm_blindtest
	echo ./crlibm_blindtest atanpi.testdata>check_atanpi; chmod a+rx check_atanpi
check_sinh: crlibm_blindtest
	echo ./crlibm_blindtest sinh.testdata>check_sinh; chmod a+rx check_sinh
check_cosh: crlibm_blindtest
	echo ./crlibm_blindtest cosh.testdata>check_cosh; chmod a+rx check_cosh
check_sinpi: crlibm_blindtest
	echo ./crlibm_blindtest sinpi.testdata>check_sinpi; chmod a+rx check_sinpi
check_cospi: crlibm_blindtest
	echo ./crlibm_blindtest cospi.testdata>check_cospi; chmod a+rx check_cospi
check_tanpi: crlibm_blindtest
	echo ./crlibm_blindtest tanpi.testdata>check_tanpi; chmod a+rx check_tanpi
check_pow: crlibm_blindtest
	echo ./crlibm_blindtest pow.testdata>check_pow; chmod a+rx check_pow

CRLIBM_CHECKLIST = check_exp check_expm1 check_log check_log2 check_log10 check_log1p \
check_sin check_cos check_tan check_asin check_acos check_atan  check_sinh check_cosh  \
check_sinpi check_cospi check_tanpi check_asinpi check_acospi check_atanpi check_pow

TESTS = $(CRLIBM_CHECKLIST)

EXTRA_DIST = exp.testdata expm1.testdata log.testdata log2.testdata log10.testdata log1p.testdata  \
sin.testdata cos.testdata tan.testdata asin.testdata acos.testdata atan.testdata sinh.testdata cosh.testdata \
sinpi.testdata cospi.testdata tanpi.testdata asinpi.testdata acospi.testdata atanpi.testdata pow.testdata

CLEANFILES = $(CRLIBM_CHECKLIST)