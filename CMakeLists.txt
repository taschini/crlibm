PROJECT(CRLIBM C)
MESSAGE("Trying to build crlibm on a " ${CMAKE_SYSTEM_PROCESSOR} " under " ${CMAKE_SYSTEM_NAME} ", compiler is "${CMAKE_C_COMPILER})

ADD_LIBRARY(crlibm SHARED
	crlibm.h crlibm_private.h crlibm_private.c triple-double.h
	exp-td.h exp-td.c exp-td-standalone.c
	expm1-standalone.c exp_accurate.h exp_accurate.c 
	expm1.h expm1.c 
	log-td.c log-td.h 
	log1p.c 
	log10-td.h log10-td.c 
	log2-td.h  log2-td.c 
	rem_pio2_accurate.h rem_pio2_accurate.c 
	trigo_fast.c trigo_fast.h trigo_accurate.c trigo_accurate.h 
	asin-td.h asin-td.c 
	acos-td.h acos-td.c 
	atan_fast.c atan_fast.h atan_accurate.h atan_accurate.c 
	csh_fast.h csh_fast.c 
	scs_lib/scs_private.c scs_lib/addition_scs.c
	scs_lib/division_scs.c scs_lib/print_scs.c
	scs_lib/double2scs.c scs_lib/zero_scs.c
	scs_lib/multiplication_scs.c scs_lib/scs2double.c
	scs_lib/tests/tbx_timing.h
)

ADD_DEFINITIONS(-DSCS_NB_WORDS=8 -DSCS_NB_BITS=30)


ADD_SUBDIRECTORY(tests)

ENABLE_TESTING()

ADD_TEST(exp tests/blind_test tests/exp.testdata)
ADD_TEST(log tests/blind_test tests/log.testdata)
ADD_TEST(sin tests/blind_test tests/sin.testdata)

ADD_TEST(cos tests/blind_test tests/cos.testdata)
ADD_TEST(tan tests/blind_test tests/tan.testdata)
ADD_TEST(atan tests/blind_test tests/atan.testdata)
ADD_TEST(sinh tests/blind_test tests/sinh.testdata)
ADD_TEST(cosh tests/blind_test tests/cosh.testdata)



# Define all the variables that were defined by the configure scripts
# We define them by -D flags, not by config.h

# First check the host processor

IF(CMAKE_SYSTEM_PROCESSOR MATCHES i*86|k6|k7|pentium*|athlon*)
  ADD_DEFINITIONS(-DCRLIBM_TYPECPU_X86)
ENDIF(CMAKE_SYSTEM_PROCESSOR MATCHES i*86|k6|k7|pentium*|athlon*)

IF(CMAKE_SYSTEM_PROCESSOR MATCHES alpha)
  ADD_DEFINITIONS(-DCRLIBM_TYPECPU_ALPHA)
ENDIF(CMAKE_SYSTEM_PROCESSOR MATCHES alpha)

IF(CMAKE_SYSTEM_PROCESSOR MATCHES powerpc*)
  ADD_DEFINITIONS(-DCRLIBM_TYPECPU_POWERPC)
ENDIF(CMAKE_SYSTEM_PROCESSOR MATCHES powerpc*)

IF(CMAKE_SYSTEM_PROCESSOR MATCHES sparc)
  ADD_DEFINITIONS(-DCRLIBM_TYPECPU_SPARC)
ENDIF(CMAKE_SYSTEM_PROCESSOR MATCHES sparc)

IF(CMAKE_SYSTEM_PROCESSOR MATCHES x86_64|amd64)
  ADD_DEFINITIONS(-DCRLIBM_TYPECPU_AMD64)
ENDIF(CMAKE_SYSTEM_PROCESSOR MATCHES x86_64|amd64)

IF(CMAKE_SYSTEM_PROCESSOR MATCHES ia64)
  ADD_DEFINITIONS(-DCRLIBM_TYPECPU_ITANIUM)
ENDIF(CMAKE_SYSTEM_PROCESSOR MATCHES ia64)


FIND_PATH(INTTYPES_H
	  inttypes.h
	  /usr/include/ /usr/local/include/
	  DOC "C99 file that defines int32_t, int64_t, etc"
)

IF(INTTYPES_H)
  ADD_DEFINITIONS(-DHAVE_INTTYPES_H)
ENDIF(INTTYPES_H)


# Now look for libraries

#   MPFR

FIND_PATH(MPFR_H
	  mpfr.h
	  /usr/include/ /usr/local/include/
	  DOC "Path of mpfr.h, the include file for GNU MPFR library, useful only for testing and developing crlibm"
)

FIND_LIBRARY(MPFR_LIB 
             NAMES mpfr 
             PATH /usr/lib/ /usr/local/lib/  
             DOC "Directory of the GNU MPFR library, useful for testing and developing crlibm")

# for compatibility with configure
IF(MPFR_H AND MPFR_LIB)
  ADD_DEFINITIONS(-DHAVE_MPFR_H)
  INCLUDE_DIRECTORIES(${MPFR_H})   # add include files
  LINK_DIRECTORIES(${MPFR_LIB})    # add link dir and link flag
  SET(LIBMPFR mpfr)
  SET(LIBGMP gmp)         # variables used in tests/CMakeList.txt
ENDIF(MPFR_H AND MPFR_LIB)


# IBM's Libultim

FIND_PATH(IBM_H
	  MathLib.h
	  /usr/include/ /usr/local/include/
	  DOC "Path of Mathlib.h, the include file for IBM's libultim, for comparison purpose"
)

FIND_LIBRARY(IBM_LIB 
             NAMES ultim 
             PATH /usr/lib/ /usr/local/lib/  
             DOC "Directory of IBM's libultim, for comparison purpose")

# for compatibility with configure
IF(IBM_H AND IBM_LIB)
  ADD_DEFINITIONS(-DHAVE_MATHLIB_H)
  INCLUDE_DIRECTORIES(${IBM_H})   # add include files
  LINK_DIRECTORIES(${IBM_LIB})    # add link dir and link flag
  SET(LIBIBM ultim)
ENDIF(IBM_H AND IBM_LIB)

  
# Sun's libmcr

FIND_PATH(MCR_H
	  mcr.h
	  /usr/include/ /usr/local/include/
	  DOC "Path of mcr.h, the include file for Sun's libmcr, for comparison purpose"
)

FIND_LIBRARY(MCR_LIB 
             NAMES mcr 
             PATH /usr/lib/ /usr/local/lib/  
             DOC "Directory of Sun's libmcr, for comparison purpose")

# for compatibility with configure
IF(MCR_H AND MCR_LIB)
  ADD_DEFINITIONS(-DHAVE_MCR_H)
  INCLUDE_DIRECTORIES(${MCR_H})   # add include files
  LINK_DIRECTORIES(${MCR_LIB})    # add link dir and link flag
  SET(LIBMCR mcr)
ENDIF(MCR_H AND MCR_LIB)

  
