dnl Process this file with autoconf to produce a configure scsript.
AC_INIT

scs_version=`cat $srcdir/VERSION`
dnl For Automake
AM_INIT_AUTOMAKE(scslib, $scs_version)

AM_CONFIG_HEADER(scs_config.h)


AC_SUBST(VERSION)

ISODATE=`date +%Y-%m-%d`
AC_SUBST(ISODATE)


dnl Check for system information
AC_CANONICAL_HOST
case $host_os in
hpux*)
	AC_DEFINE_UNQUOTED(SCS_TYPEOS_HPUX,1,[OS type]);;
esac;


case $host_cpu in
alpha)     
	AC_DEFINE_UNQUOTED(SCS_TYPECPU_ALPHA,1,[Processor type]) ;;
sparc)     
	AC_DEFINE_UNQUOTED(SCS_TYPECPU_SPARC,1,[Processor type]) ;;
powerpc*)  
	AC_DEFINE_UNQUOTED(SCS_TYPECPU_POWERPC,1,[Processor type]);;
i*86|k6|k7|pentium*|athlon*)
	AC_DEFINE_UNQUOTED(SCS_TYPECPU_X86,1,[Processor type]);;
x86_64|amd64)
        AC_DEFINE_UNQUOTED(SCS_TYPECPU_AMD64,1,[Processor type]);;
ia64)
	AC_DEFINE_UNQUOTED(SCS_TYPECPU_ITANIUM,1,[Processor type]);;
esac;


AC_CONFIG_SRCDIR([print_scs.c])
AC_PROG_MAKE_SET

dnl Check for compiler


AC_PROG_CC(gcc-3.3 gcc-3.2 gcc icc cc)
# Try to set compiler flags to get floating point arithmetics 
# to work
#
if test "${GCC}" != "yes" ; then
    case "${build}" in
        *-*-hpux* )
	    CFLAGS="${CFLAGS} -Wp,-H30000"
            ;;
        *-dec-osf* )
            CFLAGS="${CFLAGS} -ieee"
            ;;
	*-apple-darwin* )
	    CFLAGS="${CFLAGS} -fno-common"
	    ;;
	*-sgi-irix* )
	    ;;
    esac
else
    case "${CC}" in
	*icc* )
	    # -Wall causes relocation errors with icc
	    CFLAGS="${CFLAGS} -wd269"
	    ;;
        * )
            case "${build}" in
	        *-dec-osf* )
	            CFLAGS="${CFLAGS} -mieee -Wall"
	            ;;
                *alpha*-*-linux* )
                    CFLAGS="${CFLAGS} -mieee -Wall"
                    ;;
                *alpha*-*-freebsd* )
                    CFLAGS="${CFLAGS} -mieee -Wall"
                    ;;
	        *-apple-darwin* )
	            CFLAGS="${CFLAGS} -fno-common -Wall"
	            ;;
	        * )
	            CFLAGS="${CFLAGS} -Wall"
	            ;;
            esac
	    ;;
    esac
fi



dnl send this information to automake
AM_CONDITIONAL(COMPILER_ICC, test x$CC = xicc)

dnl Checks for other programs.
AC_PROG_INSTALL
AC_PROG_RANLIB

dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS(limits.h)






dnl Optional compile-time arguments to set the precision.
dnl 

dnl If at some point we wish to include SPARC optimization back,
dnl  the default should take the value 10 and 25 in this case 
dnl  for around 200 bits of precision

AC_ARG_ENABLE(nbdigits,
	[  --enable-nbdigits=N      set the number of digits in the SCS structure to N [[default=8]] ],
	[case $enableval in
	      *[[0-9]])  AC_DEFINE_UNQUOTED(SCS_NB_WORDS, $enableval, [Number of digits in the SCS structure]);;
	      *) AC_MSG_ERROR([bad value $enableval for --enable-nbdigits, number expected]) ;;
	 esac],
	[AC_DEFINE_UNQUOTED(SCS_NB_WORDS, 8, [Number of digits in the SCS structure])]
)

AC_ARG_ENABLE(digitsize,
	[ --enable-digitsize=N      set the size of a digit to N bits [[default=30]] ],
	[case $enableval in
	      *[[0-9]])  AC_DEFINE_UNQUOTED(SCS_NB_BITS, $enableval, [Size of an SCS digit]);;
	      *) AC_MSG_ERROR([bad value $enableval for --enable-digitsize, number expected]) ;;
	 esac],
	[AC_DEFINE_UNQUOTED(SCS_NB_BITS, 30, [Size of an SCS digit])]
)



dnl Sometimes the following function is provided by aclocal, sometimes not.
dnl Copying it here is probably a bad way of managing this problem

AC_DEFUN(AC_ADD_PATH_LIBS,
[
if ` test "$1" `
then  
  AC_MSG_CHECKING($2 library)
	if  test -r "$1/lib$2.a"
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


dnl  *******************
dnl  ****** MPFR *******
dnl  *******************
AC_ARG_ENABLE(mpfr,       [  --enable-mpfr           enable MPFR function [[default=no]]],
[case $enableval in
yes|no) ;;
*) AC_MSG_ERROR([bad value $enableval for --enable-mpfr, need yes or no]) ;;
esac],
[enable_mpfr=no])


dnl  *******************
dnl  ******* GMP *******
dnl  *******************
AC_ARG_ENABLE(gmp,        [  --enable-gmp            enable GMP function [[default=no]]],
[case $enableval in
yes|no) ;;
*) AC_MSG_ERROR([bad value $enableval for --enable-gmp, need yes or no]) ;;
esac],
[enable_gmp=no])





dnl  *******************
dnl  ******* GMP *******
dnl  *******************
if test "$enable_gmp" = "yes" -o "$enable_mpfr" = "yes"; then
dnl Directories setting for GMP
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



dnl  *******************
dnl  ****** MPFR *******
dnl  *******************
if test "$enable_mpfr" = "yes"; then
dnl Directories setting for MPFR
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




dnl adds -lm
AC_CHECK_LIB([m], [log])


dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE


dnl Other checks
AC_C_BIGENDIAN(,,)
dnl defines WORDS_BIGENDIAN, or not 

AC_C_LONG_DOUBLE
dnl defines HAVE_LONG_DOUBLE, or not

AC_CHECK_SIZEOF(int)
dnl defines SIZEOF_INT to be the size in bytes of int


AC_CONFIG_FILES([Makefile tests/Makefile])
AC_OUTPUT


