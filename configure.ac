dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.63)
AC_INIT([lorg_tools], [git], [jleroux@computing.dcu.ie])

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])

LT_INIT

CXXFLAGS=""
CFLAGS=""

AC_PROG_CC_STDC

AC_PROG_CXX

AC_LANG_CPLUSPLUS


AX_CXX_COMPILE_STDCXX_11([noext])
AX_CHECK_COMPILE_FLAG([-Wall -Wextra -pedantic -Wno-unused-local-typedefs -Wno-unknown-warning-option -Wno-unknown-pragmas -march=native],\
[CXXFLAGS="$CXXFLAGS -Wall -Wextra -pedantic -Wno-unused-local-typedefs -Wno-unknown-warning-option -Wno-unknown-pragmas -march=native"])


AC_CONFIG_HEADERS([config.h])

AC_SUBST(CXXFLAGS)
AC_SUBST(LDFLAGS)

AC_PROG_SED

##not sure why we need that
#AC_PROG_LEX

AC_C_INLINE
AC_HEADER_STDBOOL
AC_CHECK_HEADER([string],,AC_MSG_ERROR(STL classes missing ?))
AC_CHECK_HEADER([list],,AC_MSG_ERROR(STL classes missing ?))
AC_CHECK_HEADER([vector],,AC_MSG_ERROR(STL classes missing ?))

AC_CHECK_FUNCS([memset],,,AC_MSG_ERROR(memset is missing))
AC_CHECK_FUNCS([strdup],,,AC_MSG_ERROR(strdup is missing))


#Fix a bug in Linux with tbb tick_count (need -ltr)
case "$host_os" in
    linux*)
        AM_CONDITIONAL(WANT_LTR, true);;
    *)
        AM_CONDITIONAL(WANT_LTR, false);;
esac


AC_ARG_WITH(tcmalloc, [message],
                      [
                      if test "$withval" = "no"; then
                        want_tcmalloc="no"
                      elif test "$withval" = "yes"; then
                        want_tcmalloc="yes"
                      else
                        want_tcmalloc="yes"
                      fi
                   ], [ want_tcmalloc="no" ])

AM_CONDITIONAL([WANT_TCMALLOC], [test "x$want_tcmalloc" != "xno"])


AX_BOOST_BASE([1.56])
AX_BOOST_REGEX
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_SYSTEM
AX_BOOST_FILESYSTEM
AX_BOOST_UNIT_TEST_FRAMEWORK
AX_BOOST_SERIALIZATION


dnl AC_ARG_WITH(tbb, [message],
dnl                  [
dnl                  if test "$withval" = "no"; then
dnl                     want_tbb="no"
dnl                  else
dnl                     AX_TBB([use_tbb=yes], AC_MSG_WARN([cannot find the Intel TBB library consider to give --with-tbb to link it. Will compile without thread (parser only)]))
dnl                  fi
dnl                  ]
dnl                  )


AX_TBB([use_tbb=yes], AC_MSG_WARN([cannot find the Intel TBB library consider to give --with-tbb to link it. Will compile without thread (parser only)]))


#echo $want_tbb

AM_CONDITIONAL([WANT_TBB],[test "x$use_tbb" != "xno"])

AM_CONDITIONAL([WANT_TBB], [test "x$want_tbb" != "xno"])



# Check for DYNET
AC_ARG_WITH(dynet,
	[AC_HELP_STRING([--with-dynet=DIR], [dynet in DIR])],
	[dynet_dir="${withval}"], [dynet_dir=''])
if test "x$dynet_dir" != "x"; then
  DYNET_CPPFLAGS="-I${dynet_dir}"
  DYNET_LDFLAGS="-L${dynet_dir}/build/dynet/ -ldynet"
else
  AS_ERROR("Must specify --with-dynet=[DYNET_DIR]")
fi

# Check if we have CUDA
AC_ARG_WITH(cuda,
	[AC_HELP_STRING([--with-cuda=DIR], [cuda in DIR])],
	[cuda_dir="${withval}"], [cuda_dir=''])
if test "x$cuda_dir" != "x"; then
  DYNET_CPPFLAGS="-I${dynet_dir} -I${cuda_dir}/include -DHAVE_CUDA"
  DYNET_LDFLAGS="-L${dynet_dir}/build/dynet -ldynet -ldynetcuda -L${cuda_dir}/lib64 -lcudart -lcublas"
fi
AC_SUBST(DYNET_CPPFLAGS)
AC_SUBST(DYNET_LDFLAGS)

# Check for Eigen
AC_ARG_WITH(eigen,
	[AC_HELP_STRING([--with-eigen=DIR], [eigen in DIR])],
	[eigen_dir="${withval}"], [eigen_dir=''])
if test "x$eigen_dir" != "x"; then
  EIGEN_CPPFLAGS="-I${eigen_dir}"
else
  AS_ERROR("Must specify --with-eigen=[EIGEN_DIR]")
fi
AC_SUBST(EIGEN_CPPFLAGS)



AC_CONFIG_FILES([
        Makefile
	src/Makefile
	src/training/Makefile
        test/Makefile
])
AC_OUTPUT
