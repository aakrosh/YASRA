#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([YASRA], [2.0], [ratan@bx.psu.edu])
AC_CONFIG_SRCDIR([src/assembler.c])
AC_CONFIG_HEADERS([src/config.h])
AM_INIT_AUTOMAKE([-W -Wall -Werror])

AC_CANONICAL_HOST

case "${host}" in
	i?86*) my_inline="-minline-all-stringops -momit-leaf-frame-pointer" ;;
	x86_64*) my_inline="-minline-all-stringops -momit-leaf-frame-pointer" ;;
	* ) my_inline="" ;;
esac

AC_SUBST([ARCH],[$my_inline])

## debug
AC_ARG_ENABLE([debug],
[  --enable-debug      Compile for debugging],
[case "${enableval}" in
  yes)	debug=true ;;
  no)	debug==false ;;
  *)	AC_MSG_ERROR([bad value ${enableval} for -enable-debug]) ;;
esac],[debug=false])
AM_CONDITIONAL([DEBUG], [test x$debug = xtrue])

## profile
AC_ARG_ENABLE([profile],
[  --enable-profile    Compile for profiling],
[case "${enableval}" in
  yes)	profile=true ;;
  no)	profile==false ;;
  *)	AC_MSG_ERROR([bad value ${enableval} for -enable-profile]) ;;
esac],[profile=false])
AM_CONDITIONAL([PROFILE], [test x$profile = xtrue])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([inttypes.h limits.h stdint.h stdlib.h string.h unistd.h wchar.h])

# Checks for other programs
AC_CHECK_PROG(LASTZ, lastz, [yes], [no])
if test "$LASTZ" == "yes"; then
    AC_MSG_RESULT("yes")
else
    AC_MSG_RESULT("Please install LASTZ and add it to your PATH before installing YASRA")
    exit 1;
fi


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset pow strchr strdup strerror strrchr strstr strtol])

AC_CONFIG_FILES([
  Makefile
  src/Makefile
])
AC_OUTPUT
