dnl Process this file with autoconf to produce a configure script.
AC_INIT

crlibm_version=`cat $srcdir/VERSION`
dnl For Automake
AM_INIT_AUTOMAKE(crlibm, $crlibm_version)

AM_CONFIG_HEADER(crlibm_config.h)


AC_SUBST(VERSION)

ISODATE=`date +%Y-%m-%d`
AC_SUBST(ISODATE)

dnl Check for system information
AC_CANONICAL_HOST

case $host_cpu in
alpha)     
	AC_DEFINE_UNQUOTED(CRLIBM_TYPECPU_ALPHA,1,[Processor type]) ;;
sparc)     
	AC_DEFINE_UNQUOTED(CRLIBM_TYPECPU_SPARC,1,[Processor type]) ;;
powerpc*)  
	AC_DEFINE_UNQUOTED(CRLIBM_TYPECPU_POWERPC,1,[Processor type]);;
i*86|k6|k7|pentium*|athlon*)
	AC_DEFINE_UNQUOTED(CRLIBM_TYPECPU_X86,1,[Processor type]);;
x86_64|amd64)
        AC_DEFINE_UNQUOTED(CRLIBM_TYPECPU_AMD64,1,[Processor type]);;
ia64)
	AC_DEFINE_UNQUOTED(CRLIBM_TYPECPU_ITANIUM,1,[Processor type]);;
esac;


AC_CONFIG_SRCDIR([exp.c])
AC_PROG_MAKE_SET

dnl Check for compiler

AC_PROG_CC(gcc icc cc)

dnl send this information to automake to set up compilation flags
dnl (add your compiler here)
AM_CONDITIONAL(COMPILER_ICC, test x$CC = xicc)

dnl Checks for programs.
AC_PROG_INSTALL
AC_PROG_RANLIB
AC_PROG_AWK
AC_PROG_LN_S

dnl Checks for standard C stuff.
AC_HEADER_STDC
AC_FUNC_MALLOC
AC_C_CONST
AC_C_INLINE


AC_CHECK_FUNCS([sqrt])
AC_CHECK_HEADERS([fenv.h])
AC_CHECK_HEADERS([float.h])

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
AC_DEFUN(AC_ADD_PATH_LIBS,
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

dnl  ... for  MPFR
AC_ARG_ENABLE(mpfr,       [  --enable-mpfr           enable MPFR function [[default=no]]],
[case $enableval in
yes|no) ;;
*) AC_MSG_ERROR([bad value $enableval for --enable-mpfr, need yes or no]) ;;
esac],
[enable_mpfr=no])


dnl  ... for GMP
AC_ARG_ENABLE(gmp,        [  --enable-gmp            enable GMP function [[default=no]]],
[case $enableval in
yes|no) ;;
*) AC_MSG_ERROR([bad value $enableval for --enable-gmp, need yes or no]) ;;
esac],
[enable_gmp=no])


dnl  ... for IBM Accurate Portable Library
AC_ARG_ENABLE(ultim,        [  --enable-ultim        enable libultim, IBM's Accurate Portable Math Library [[default=no]]],
[case $enableval in
yes|no) ;;
*) AC_MSG_ERROR([bad value $enableval for --enable-ultim, need yes or no]) ;;
esac],
[enable_ultim=no])


dnl  ... for Sun's libmcr
AC_ARG_ENABLE(mcr,        [  --enable-mcr        enable Sun's LIBMCR  [[default=no]]],
[case $enableval in
yes|no) ;;
*) AC_MSG_ERROR([bad value $enableval for --enable-mcr, need yes or no]) ;;
esac],
[enable_mcr=no])


dnl Now consider the result of the previous to set directories:


dnl ... for GMP
if test "$enable_gmp" = "yes" -o "$enable_mpfr" = "yes"; then

AC_ARG_WITH(gmp_include, [ --with-gmp-include=DIR   GMP include directory ],
 with_gmp_include=$withval)
AC_ARG_WITH(gmp_lib,     [ --with-gmp-lib=DIR       GMP lib directory ],
 with_gmp_lib=$withval)
AC_ARG_WITH(gmp,         [ --with-gmp=DIR           GMP directory ],
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
AC_ARG_WITH(mpfr_include, [ --with-mpfr-include=DIR  MPFR include directory ], 
with_mpfr_include=$withval)
AC_ARG_WITH(mpfr_lib,     [ --with-mpfr-lib=DIR      MPFR lib directory ], 
with_mpfr_lib=$withval)
AC_ARG_WITH(mpfr,         [ --with-mpfr=DIR          MPFR directory ],
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

AC_ARG_WITH(ultim_include, [ --with-ultim-include=DIR   ULTIM include directory ],
 with_ultim_include=$withval)
AC_ARG_WITH(ultim_lib,     [ --with-ultim-lib=DIR       ULTIM lib directory ],
 with_ultim_lib=$withval)
AC_ARG_WITH(ultim,         [ --with-ultim=DIR           ULTIM directory ],
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

AC_ARG_WITH(mcr_include, [ --with-mcr-include=DIR   LIBMCR include directory ],
 with_mcr_include=$withval)
AC_ARG_WITH(mcr_lib,     [ --with-mcr-lib=DIR   LIBMCR lib directory ],
 with_mcr_lib=$withval)
AC_ARG_WITH(mcr,         [ --with-mcr=DIR       LIBMCR directory ],
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


AC_CONFIG_SUBDIRS([scs_lib])



AC_OUTPUT(Makefile tests/Makefile)


