AC_INIT([resample_image], [1.0.0], [glotfeltyk@si.edu])
AM_INIT_AUTOMAKE([foreign -Wall])
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CXX
AC_PROG_SED

# -- Hack to provide missing & corrupt pkg-config files

## Missing .pc files: ascdm, region, and readline(?)
for p in pkgconfig/*pc.in
do
  o=`echo $p | sed 's,\.in$,,'`  
  sed s,@prefix@,$prefix,g < $p > $o
done

## caldb in ciao4.11 is bad; fixed in 4.12
cat $prefix/lib/pkgconfig/caldb4.pc | \
  sed 's,^Libs:.*,Libs: -L${libdir} -lcaldb4,' | \
  sed 's,^Cflags:.*,Cflags: -I${includedir},' > pkgconfig/caldb4.pc 

## Some OTS have hard coded /proj stuff
if test -e $prefix/ots/lib/pkgconfig
then
    for p in cfitsio fftw3 fftw3f gsl xpa
    do
        sed "s,^prefix=.*,prefix=$prefix/ots," < $prefix/ots/lib/pkgconfig/${p}.pc > pkgconfig/${p}.pc
    done
fi

# -- Readline -- there has got to be a better way to do this ?

AC_CHECK_FUNC(tgetent,
   [READLINE_SUPPLIB=],
     [AC_CHECK_LIB(tinfo, tgetent,
     [READLINE_SUPPLIB=-ltinfo],
       [AC_CHECK_LIB(curses, tgetent,
         [READLINE_SUPPLIB=-lcurses],
         [AC_CHECK_LIB(ncurses, tgetent,
            [READLINE_SUPPLIB=-lncurses],
            [AC_CHECK_LIB(termcap, tgetent,
               [READLINE_SUPPLIB=-ltermcap],
               [AC_MSG_WARN(No termcap nor curses library found)])
            ])
         ])
      ])
   ])

AC_CHECK_LIB(readline,readline,[READLINELIB="-lreadline -lhistory -lncurses"],
	AC_MSG_ERROR([
You need the GNU readline library(ftp://ftp.gnu.org/gnu/readline/ ) to build
this program!
]), [$READLINE_SUPPLIB])

LIBS="$LIBS $PROCSTATLIB $READLINELIB $READLINE_SUPPLIB"

# ---------- CIAO PKG_CONFIG stuff
PKG_CONFIG_PATH=$prefix/lib/pkgconfig:$prefix/ots/lib/pkgconfig:$PKG_CONFIG_PATH
PKG_CONFIG_PATH=./pkgconfig:${PKG_CONFIG_PATH}
export PKG_CONFIG_PATH
PKG_CHECK_MODULES( [CIAO], [
  ascdm
  cxcparam
  stk
  err
  ds
  dmimgio
])

# -------- dmimgio.h ------------
foo=$CPPFLAGS
CPPFLAGS=$CIAO_CFLAGS
AC_CHECK_HEADER( dmimgio.h, [], 
  AC_MSG_ERROR([Cannot locate dmimgio.h header file], [1])
)
CPPFLAGS=$foo


# ------------

AM_CONDITIONAL([LINUX], [test x`uname` == xLinux])



# -------------- Which files to configure ----------

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  wrapper/Makefile
  test/Makefile
])

AC_CONFIG_FILES([wrapper/resample_image:wrapper/wrapper.sh],[chmod +x wrapper/resample_image])

# -- Go! -----------

AC_OUTPUT
  
