dnl Process this file with autoconf to produce a configure script.
AC_INIT

crlibm_version=`cat $srcdir/VERSION`
dnl For Automake
AM_INIT_AUTOMAKE(crlibm, $crlibm_version)

AC_CONFIG_HEADERS([crlibm_config.h])


AC_SUBST(VERSION)

ISODATE=`date +%Y-%m-%d`
AC_SUBST(ISODATE)

dnl Check for system information
AC_CANONICAL_HOST
case $host_os in
hpux*)
	AC_DEFINE_UNQUOTED(CRLIBM_TYPEOS_HPUX,1,[OS type]);;
*cygwin* | *mingw* | *uwin* | *djgpp | *emx*)
	AC_DEFINE_UNQUOTED(CRLIBM_TYPEOS_CYGWIN,1,[OS type]);;
freebsd*|netbsd*|darwin*)
	AC_DEFINE_UNQUOTED(CRLIBM_TYPEOS_BSD,1,[OS type]);;
esac;


case $host_cpu in
alpha)     
	AC_DEFINE_UNQUOTED(CRLIBM_TYPECPU_ALPHA,1,[Processor type]) ;;
sparc)     
	AC_DEFINE_UNQUOTED(CRLIBM_TYPECPU_SPARC,1,[Processor type]) ;;
powerpc*)  
	AC_DEFINE_UNQUOTED(CRLIBM_TYPECPU_POWERPC,1,[Processor type]);;
i*86|k6|k7|pentium*|athlon*)
	AC_DEFINE_UNQUOTED(CRLIBM_TYPECPU_X86,1,[Processor type])
	has_ia32_de=true;;
x86_64|amd64)
        AC_DEFINE_UNQUOTED(CRLIBM_TYPECPU_AMD64,1,[Processor type])
	has_ia32_de=true;;
ia64)
	AC_DEFINE_UNQUOTED(CRLIBM_TYPECPU_ITANIUM,1,[Processor type])
	has_ia64_de=true;;
esac;


AC_CONFIG_SRCDIR([crlibm.h])
AC_PROG_MAKE_SET

dnl Check for compiler

AC_PROG_CC(gcc icc cc)

# Try to set compiler flags to get floating point arithmetics 
# to work. This is an old heritage and should become useless within 30 years
#
if test "${GCC}" != "yes" ; then
    case "${build}" in
        *-*-hpux* )
	    AM_CFLAGS="${AM_CFLAGS} -Wp,-H30000"
            ;;
        *-dec-osf* )
            AM_CFLAGS="${AM_CFLAGS} -ieee"
            ;;
	*-apple-darwin* )
	    AM_CFLAGS="${AM_CFLAGS} -fno-common"
	    ;;
	*-sgi-irix* )
	    ;;
    esac
else
    case "${CC}" in
	*icc* )
	    # -Wall causes relocation errors with icc. 
	    AM_CFLAGS="${AM_CFLAGS} -wd269"
	    ;;
        * )
            case "${build}" in
	        *-dec-osf* )
	            AM_CFLAGS="${AM_CFLAGS} -mieee -Wall"
	            ;;
                *alpha*-*-linux* )
                    AM_CFLAGS="${AM_CFLAGS} -mieee -Wall"
                    ;;
                *alpha*-*-freebsd* )
                    AM_CFLAGS="${AM_CFLAGS} -mieee -Wall"
                    ;;
	        *-apple-darwin* )
	            AM_CFLAGS="${AM_CFLAGS} -fno-common -Wall"
	            ;;
	        * )
	            AM_CFLAGS="${AM_CFLAGS} -Wall"
	            ;;
            esac
	    ;;
    esac
fi






dnl send this information to automake to set up compilation flags
dnl (add your compiler here)
AM_CONDITIONAL(COMPILER_ICC, test x$CC = xicc)

dnl Checks for programs.
AC_PROG_INSTALL
AC_PROG_RANLIB
AC_PROG_LN_S

dnl Checks for standard C stuff.
AC_HEADER_STDC
AC_C_CONST
AC_C_INLINE


AC_CHECK_FUNCS([sqrt])
AC_CHECK_HEADERS([fenv.h])
AC_CHECK_HEADERS([float.h])
AC_CHECK_HEADER(fpu_control.h, has_fpu_control=yes, has_fpu_control=no)

dnl Send this information to automake
AM_CONDITIONAL(HAS_FPU_CONTROL, 	test x$has_fpu_control = xyes)

if test x$has_fpu_control = xyes ; then
AC_DEFINE_UNQUOTED(CRLIBM_HAS_FPU_CONTROL,1,[architecture- and system-specific FPU control header file])
fi

dnl adds -lm
AC_CHECK_LIB([m], [log])

dnl Other checks
AC_C_BIGENDIAN(,,)
dnl defines WORDS_BIGENDIAN, or not 

AC_C_LONG_DOUBLE
dnl defines HAVE_LONG_DOUBLE, or not


AC_CHECK_SIZEOF(int)
dnl defines SIZEOF_INT to be the size in bytes of int





dnl All the rest is related to MPFR and to Ziv's lib, useful for testing.

dnl Sometimes the following function is provided by aclocal, sometimes not.
dnl Copying it here is probably a bad way of managing this problem
AC_DEFUN([AC_ADD_PATH_LIBS],
[
if ` test "$1" `
then  
  AC_MSG_CHECKING($2 library)
	if  test -r "$1/lib$2.a" -o -r "$1/lib$2.so"
	then
	  LDFLAGS="$LDFLAGS -L$1"
          LIBS   ="$LIBS -l$2" 	
	else
	   AC_MSG_ERROR($2 not found)
	fi
  AC_MSG_RESULT(yes)
else
  AC_CHECK_LIB($2, main, , AC_MSG_ERROR($2 not found))
fi
]
)





dnl Command-line arguments to ./configure:


dnl switch for SSE2
AC_ARG_ENABLE(sse2,
  [  --enable-sse2    forces the use of sse2 arith  [[default=no]]],
  [case $enableval in
    yes)
      AC_DEFINE_UNQUOTED(HAS_SSE2,1,[Use the sse2 operators ]);;
    no) ;;
    *) AC_MSG_ERROR([bad value $enableval for --enable-sse2, need yes or no]) ;;
  esac],
  [enable_sse2=no])
dnl Send this information to automake 
AM_CONDITIONAL(USE_SSE2, test x$enable_sse2 = xyes  )




dnl switch for compiling the interval functions
AC_ARG_ENABLE(interval-functions,
  [  --enable-interval-functions    enables filib-compatible interval functions (DOES NOT WORK YET)  [[default=no]]],
  [case $enableval in
    yes)
      AC_DEFINE_UNQUOTED(BUILD_INTERVAL_FUNCTIONS,1,[Compile the filib-compatible interval functions]);;
    no) ;;
    *) AC_MSG_ERROR([bad value $enableval for --enable-interval-functions, need yes or no]) ;;
  esac],
  [enable_interval_functions=no])
dnl Send this information to automake 
AM_CONDITIONAL(USE_INTERVAL_FUNCTIONS, test x$enable_interval_functions = xyes  )



dnl A switch to disable use of double-extended optimised version
AC_ARG_ENABLE(double-extended,        
  [  --enable-double-extended       enable double-extended optimised versions of functions if hardware supports it  [[default=yes]]],
  [case $enableval in
    yes|no) ;;
    *) AC_MSG_ERROR([bad value $enableval for --enable-double-extended, need yes or no]) ;;
  esac],
  [enable_double_extended=yes])

dnl following line needs a little bit more work on log-de.c and double-extended.h
dnl AM_CONDITIONAL(USE_HARDWARE_DE, test x$has_ia32_de = xtrue -a x$enable_double_extended = xyes -a x$has_fpu_control=xyes   -o  x$has_ia64_de = xtrue -a x$CC = xicc -a x$enable_double_extended = xyes  )

AM_CONDITIONAL(USE_HARDWARE_DE, test x$has_ia32_de = xtrue  -a x$has_fpu_control = xyes  -a x$enable_double_extended = xyes -a x$enable_interval_functions = xno    )
dnl  Note that at the moment, enabling intervals disables double-extended








dnl If at some point we wish to include SPARC optimization back,
dnl  the default should take the value 10 and 25 in this case 
dnl  for around 200 bits of precision

AC_ARG_ENABLE(nbdigits,
  [  --enable-nbdigits=N            set the number of digits in the SCS structure to N [[default=8]] ],
  [case $enableval in
	      *[[0-9]])  AC_DEFINE_UNQUOTED(SCS_NB_WORDS, $enableval, [Number of digits in the SCS structure]);;
	      *) AC_MSG_ERROR([bad value $enableval for --enable-nbdigits, number expected]) ;;
	 esac],
  [AC_DEFINE_UNQUOTED(SCS_NB_WORDS, 8, [Number of digits in the SCS structure])]
)

AC_ARG_ENABLE(digitsize,
  [  --enable-digitsize=N           set the size of a digit to N bits [[default=30]] ],
  [case $enableval in
	      *[[0-9]])  AC_DEFINE_UNQUOTED(SCS_NB_BITS, $enableval, [Size of an SCS digit]);;
	      *) AC_MSG_ERROR([bad value $enableval for --enable-digitsize, number expected]) ;;
	 esac],
  [AC_DEFINE_UNQUOTED(SCS_NB_BITS, 30, [Size of an SCS digit])]
)


dnl  ... for  MPFR
AC_ARG_ENABLE(mpfr,       
  [  --enable-mpfr                  enable MPFR function [[default=no]]],
  [case $enableval in 
     yes|no) ;;
     *) AC_MSG_ERROR([bad value $enableval for --enable-mpfr, need yes or no]) ;;
   esac],
  [enable_mpfr=no])


dnl  ... for GMP
AC_ARG_ENABLE(gmp,
  [  --enable-gmp                   enable GMP function [[default=no]]],
  [case $enableval in
     yes|no) ;;
     *) AC_MSG_ERROR([bad value $enableval for --enable-gmp, need yes or no]) ;;
   esac],
  [enable_gmp=no])


dnl  ... for IBM Accurate Portable Library
AC_ARG_ENABLE(ultim,
  [  --enable-ultim                 enable comparison with libultim, IBM's Accurate Portable Math Library [[default=no]]],
  [case $enableval in
     yes|no) ;;
     *) AC_MSG_ERROR([bad value $enableval for --enable-ultim, need yes or no]) ;;
   esac],
   [enable_ultim=no])


dnl  ... for Sun's libmcr
AC_ARG_ENABLE(mcr,
  [  --enable-mcr                   enable comparison with Sun's LIBMCR  [[default=no]]],
  [case $enableval in
     yes|no) ;;
     *) AC_MSG_ERROR([bad value $enableval for --enable-mcr, need yes or no]) ;;
   esac],
  [enable_mcr=no])


dnl Now consider the result of the previous to set directories:


dnl ... for GMP
if test "$enable_gmp" = "yes" -o "$enable_mpfr" = "yes"; then

AC_ARG_WITH(gmp_include, [  --with-gmp-include=DIR    GMP include directory ],
 with_gmp_include=$withval)
AC_ARG_WITH(gmp_lib,     [  --with-gmp-lib=DIR        GMP lib directory ],
 with_gmp_lib=$withval)
AC_ARG_WITH(gmp,         [  --with-gmp=DIR            GMP directory ],
 with_gmp_include=$withval/include with_gmp_lib=$withval/lib)

if test -d "$with_gmp_include"; then
  CPPFLAGS="$CPPFLAGS -I$with_gmp_include"
else
  with_gmp_include=
fi
AC_ADD_PATH_LIBS($with_gmp_lib, gmp)

dnl If gmp.h is found, #defines HAVE_GMP_H
AC_CHECK_HEADERS([gmp.h])

dnl If libgmp is found, adds -lgmp to LIBS
dnl correct in theory:  AC_CHECK_LIB([gmp], [__gmp_fprintf])
AC_CHECK_LIB([gmp], [main])

dnl end of test for gmp
fi



dnl  and for MPFR 

if test "$enable_mpfr" = "yes"; then
AC_ARG_WITH(mpfr_include, [  --with-mpfr-include=DIR   MPFR include directory ], 
with_mpfr_include=$withval)
AC_ARG_WITH(mpfr_lib,     [  --with-mpfr-lib=DIR       MPFR lib directory ], 
with_mpfr_lib=$withval)
AC_ARG_WITH(mpfr,         [  --with-mpfr=DIR           MPFR directory ],
with_mpfr_include=$withval/include with_mpfr_lib=$withval/lib)

if test -d "$with_mpfr_include"; then
  CPPFLAGS="$CPPFLAGS -I$with_mpfr_include"
else
  with_mpfr_include=
fi
AC_ADD_PATH_LIBS($with_mpfr_lib, mpfr)


dnl If mpfr.h is found, #defines HAVE_MPFR_H
dnl this line should be
dnl AC_CHECK_HEADERS([mpfr.h], [], [], [#if HAVE_GMP_H #include<gmp.h> #endif])
dnl but it doesn't work properly so we assume that gmp is installed on the system
dnl is mpfr is !!!
AC_CHECK_HEADERS([mpfr.h], [], [], [#include<gmp.h>])

dnl If MPFR is found, adds -lmpfr to LIBS
dnl correct in theory:  AC_CHECK_LIB([mpfr], [mpfr_init2])
AC_CHECK_LIB([mpfr], [main])

dnl end of test for mpfr
fi



dnl ... for IBM
if test "$enable_ultim" = "yes"; then

AC_ARG_WITH(ultim_include, [  --with-ultim-include=DIR  ULTIM include directory ],
 with_ultim_include=$withval)
AC_ARG_WITH(ultim_lib,     [  --with-ultim-lib=DIR      ULTIM lib directory ],
 with_ultim_lib=$withval)
AC_ARG_WITH(ultim,         [  --with-ultim=DIR          ULTIM directory ],
 with_ultim_include=$withval/include with_ultim_lib=$withval/lib)

if test -d "$with_ultim_include"; then
  CPPFLAGS="$CPPFLAGS -I$with_ultim_include"
else
  with_ultim_include=
fi
AC_ADD_PATH_LIBS($with_ultim_lib, ultim)

AC_CHECK_HEADERS([MathLib.h])

AC_CHECK_LIB([ultim], [main])

dnl end of test for IBM libultim
fi



dnl ... for SUN
if test "$enable_mcr" = "yes"; then

AC_ARG_WITH(mcr_include, [  --with-mcr-include=DIR    LIBMCR include directory ],
 with_mcr_include=$withval)
AC_ARG_WITH(mcr_lib,     [  --with-mcr-lib=DIR        LIBMCR lib directory ],
 with_mcr_lib=$withval)
AC_ARG_WITH(mcr,         [  --with-mcr=DIR            LIBMCR directory ],
 with_mcr_include=$withval/include with_mcr_lib=$withval/lib)

if test -d "$with_mcr_include"; then
  CPPFLAGS="$CPPFLAGS -I$with_mcr_include"
else
  with_mcr_include=
fi
AC_ADD_PATH_LIBS($with_mcr_lib, mcr)

AC_CHECK_HEADERS([libmcr.h])

AC_CHECK_LIB([mcr], [main])

dnl end of test for IBM libultim
fi


#AC_CONFIG_SUBDIRS([scs_lib])



AC_CONFIG_FILES([Makefile tests/Makefile scs_lib/Makefile])
AC_OUTPUT


