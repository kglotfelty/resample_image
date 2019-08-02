AC_INIT([resample_image], [1.0.0], [glotfeltyk@si.edu])
AM_INIT_AUTOMAKE([foreign -Wall])
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CXX

# -- Hack to provide missing pkg-config files

sed s,@prefix@,$prefix,g < pkgconfig/ascdm.pc.in > pkgconfig/ascdm.pc
sed s,@prefix@,$prefix,g < pkgconfig/readline.pc.in > pkgconfig/readline.pc
sed s,@prefix@,$prefix,g < pkgconfig/region.pc.in > pkgconfig/region.pc

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
export PKG_CONFIG_PATH=$prefix/lib/pkgconfig:$prefix/ots/lib/pkgconfig:$PKG_CONFIG_PATH
PKG_CONFIG_PATH=./pkgconfig:${PKG_CONFIG_PATH}
PKG_CHECK_MODULES( [CIAO], [
  ascdm
  cxcparam
  stk
  err
  ds
])



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
  
