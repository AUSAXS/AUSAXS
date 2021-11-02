AC_DEFUN([DSR_PARSE_CGAL_MAKEFILE],
[AC_MSG_CHECKING(CGAL_MAKEFILE)
rm -rf include/dsrpdb/config.h
if test \! -z "$CGAL_MAKEFILE" -a -e "$CGAL_MAKEFILE"; then
    tname=`mktemp /tmp/cgal_makefile_dsrXXXXXX`
    
cat > $tname << _ACEOF
include $CGAL_MAKEFILE

cppflags:
	@echo \$(CGAL_CXXFLAGS)

cxxflags:
	@echo
ldflags:
	@echo \$(CGAL_LDFLAGS)
ldlibs:
	@echo -lCGAL  -lgmp -lm 
_ACEOF
    CPPFLAGS="`make -s -f $tname cppflags` $CPPFLAGS"
    CXXFLAGS="`make -s -f $tname cxxflags` $CXXFLAGS"
    LDFLAGS="`make -s -f $tname ldflags` `make -s -f $tname ldlibs` $LDFLAGS"
dnl    LDLIBS="`make -s -f $tname ldlibs` $LDLIBS"
    rm -f $tname
	echo "#define PDB_USE_CGAL" | cat > include/dsrpdb/config.h
    AC_MSG_RESULT(yes)
else
  	AC_MSG_RESULT(no)
	AC_CHECK_HEADER(CGAL/Exact_predicates_inexact_constructions_kernel.h, pdb_have_cgal=yes, pdb_have_cgal=no)
	AC_MSG_CHECKING(CGAL)
	if test "$pdb_have_cgal" == "yes"; then
		AC_MSG_RESULT(yes)
		LDFLAGS="$LDFLAGS -lCGAL "
		echo "#define PDB_USE_CGAL" | cat > include/dsrpdb/config.h
	else
		AC_MSG_RESULT(no)
	fi 
  
fi]) # AC_PROG_CC_G

AC_DEFUN([DSR_CHECK_PROGRAM_OPTIONS],
[AC_CHECK_HEADER(boost/program_options.hpp, pdb_prog_opts=yes, pdb_prog_opts=no)
AC_CHECK_LIB(boost_program_options, main, pdb_prog_opts_lib=yes, pdb_prog_opts_lib=no)
if test "$pdb_prog_opts" == "yes" -a  "$pdb_prog_opts_lib" == "yes"; then 
echo "#define PDB_USE_BOOST_PROGRAM_OPTIONS" | cat >> include/dsrpdb/config.h
LDFLAGS="$LDFLAGS -lboost_program_options"
pdb_have_program_options=yes
else
pdb_have_program_options=no
fi])


AC_DEFUN([DSR_CHECK_IMAGEMAGICK],
[AC_CHECK_LIB([Magick++], main, pdb_libmagick=yes, pdb_libmagick=no)
AC_CHECK_HEADER([Magick++.h], pdb_magickh=yes, pdb_magickh=no)
if test "$pdb_libmagick" == "yes" -a "$pdb_magickh" == "yes"; then
pdb_have_magick=yes
LDFLAGS="$LDFLAGS -lMagick++"
echo "#define PDB_USE_MAGICK" | cat >> include/dsrpdb/config.h
else
pdb_have_magick=no
fi

])

AC_DEFUN([DSR_CHECK_MULTIARRAY],
[AC_CHECK_HEADER([boost/multi_array.hpp], pdb_multiarray=yes, pdb_multiarray=no)
if test "$pdb_multiarray" == "yes"; then
pdb_have_multiarray=yes
echo "#define PDB_USE_BOOST_MULTI_ARRAY" | cat >> include/dsrpdb/config.h
else
pdb_have_multiarray=no
fi

])
