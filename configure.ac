AC_INIT(Bio++ NumCalc, 1.1.0, Julien.Dutheil@univ-montp2.fr, Bpp-NumCalc)
AC_CONFIG_SRCDIR(NumCalc/NumTools.h)
AM_INIT_AUTOMAKE()
AC_LANG_CPLUSPLUS
AC_CHECK_LIB(bpp-utils, main, , [AC_MSG_ERROR([Library bpp-utils not found.])])
AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LIBTOOL
AC_MSG_CHECKING([Covariant return with virtual inheritance])
AC_ARG_WITH(
	[virtual-cov],
	[  --with-virtual-cov      Enable covariant return type with
	                  virtual inheritance (GCC < 3.4)],
  ,
  with_virtual_cov="yes"
)
AC_MSG_RESULT($with_virtual_cov);
if test "$with_virtual_cov" == "no"
then
	AC_MSG_WARN([
+--------------------------------------------------------+
|Covariant return type with virtual inheritance disabled.|
+--------------------------------------------------------+
]);
else
	AC_DEFINE(VIRTUAL_COV, yes)
fi
AC_OUTPUT(Makefile NumCalc/Makefile)
