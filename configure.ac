
AC_INIT(RubberBand, 1.9.0, chris.cannam@breakfastquay.com)

AC_ARG_ENABLE(programs,
  [  --disable-programs      disable compilation of programs],
  [case "${enableval}" in
     yes | no ) WITH_PROGRAMS="${enableval}" ;;
     *) AC_MSG_ERROR(bad value ${enableval} for --disable-programs) ;;
   esac],
  [WITH_PROGRAMS="yes"]
)
AC_SUBST(WITH_PROGRAMS)

AC_ARG_ENABLE(vamp,
  [  --disable-vamp          disable compilation of Vamp plugin],
  [case "${enableval}" in
     yes | no ) WITH_VAMP="${enableval}" ;;
     *) AC_MSG_ERROR(bad value ${enableval} for --disable-vamp) ;;
   esac],
  [WITH_VAMP="yes"]
)
AC_SUBST(WITH_VAMP)

AC_ARG_ENABLE(ladspa,
  [  --disable-ladspa        disable compilation of LADSPA plugin],
  [case "${enableval}" in
     yes | no ) WITH_LADSPA="${enableval}" ;;
     *) AC_MSG_ERROR(bad value ${enableval} for --disable-ladspa) ;;
   esac],
  [WITH_LADSPA="yes"]
)
AC_SUBST(WITH_LADSPA)

AC_CONFIG_SRCDIR(src/StretcherImpl.h)
AC_PROG_CXX
AC_PROG_CC
AC_HEADER_STDC
AC_C_BIGENDIAN

PKG_CHECK_MODULES([SRC],[samplerate])
AC_SUBST(SRC_CFLAGS)
AC_SUBST(SRC_LIBS)

if test "x$WITH_PROGRAMS" = "xyes"; then
  PKG_CHECK_MODULES([SNDFILE],[sndfile])
fi
AC_SUBST(SNDFILE_CFLAGS)
AC_SUBST(SNDFILE_LIBS)

PKG_CHECK_MODULES([FFTW],[fftw3])
AC_SUBST(FFTW_CFLAGS)
AC_SUBST(FFTW_LIBS)

if test "x$WITH_LADSPA" = "xyes"; then
  AC_CHECK_HEADERS(ladspa.h)
fi
AC_CHECK_HEADERS(pthread.h)

if test "x$WITH_VAMP" = "xyes"; then
  PKG_CHECK_MODULES([Vamp],[vamp-sdk])
fi
AC_SUBST(Vamp_CFLAGS)
AC_SUBST(Vamp_LIBS)

changequote(,)dnl
if test "x$GCC" = "xyes"; then
  case " $CFLAGS " in
    *[\ \	]-fPIC\ -Wall[\ \	]*) ;;
    *) CFLAGS="$CFLAGS -fPIC -Wall" ;;
  esac
  case " $CXXFLAGS " in
    *[\ \	]-fPIC\ -Wall[\ \	]*) ;;
    *) CXXFLAGS="$CXXFLAGS -fPIC -Wall" ;;
  esac
fi
changequote([,])dnl

AC_OUTPUT([Makefile])

