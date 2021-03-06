#!/bin/sh

LIBS="-L/home/xiaohan/research/program/NjwFem/LinearAlgebra/Pastix/src/../install -lpastix -lgfortran -lm -lrt -L/home/xiaohan/research/program/NjwFem/LinearAlgebra/Scotch/lib -lscotch -lscotcherrexit -lpthread";
INC="-I/home/xiaohan/research/program/NjwFem/LinearAlgebra/Pastix/src/../install";
CC="gcc -Wall";
CXX=""
CCOPT="-O3  -DFORCE_NOMPI  -I/home/xiaohan/research/program/NjwFem/LinearAlgebra/Scotch/include -DWITH_SCOTCH  -DVERSION='4030' -DX_ARCHi686_pc_linux -DDOF_CONSTANT";
CL="gcc -Wall";
FC="gfortran -ffree-form -x f95-cpp-input";
FCOPT="";
FL="gfortran ";
OPTS=" -DFORCE_NOMPI  -DFORCE_DOUBLE -DPREC_DOUBLE  -I/home/xiaohan/research/program/NjwFem/LinearAlgebra/Scotch/include -DWITH_SCOTCH";
BLAS="-lblas";
VERSION="4030";
LIBS_MURGE=`echo $LIBS | sed -e 's/-lpastix/-lpastix_murge -lpastix/g'`

echo $OPTS | grep \\\-DDISTRIBUTED >/dev/null

if [ $? -ne 0 ]
then
    LIBS_MURGE="-DMURGE_MESSAGE=\"To use Murge interface"
    LIBS_MURGE="$LIBS_MURGE compile with -DDISTRIBUTED\"";
fi


usage="usage : $0 [options] - Shows PaStiX libs, includes and compiler\n";
usage="$usage    options : \n";
usage="$usage        --libs               - prints librairies\n";
usage="$usage        --libs_murge         - prints librairies\n";
usage="$usage        --incs               - prints includes\n";
usage="$usage        --cc                 - prints C compiler\n";
usage="$usage        --cxx                - prints C++ compiler\n";
usage="$usage        --ccopts             - prints C compiler options\n";
usage="$usage        --cl                 - prints C linker\n";
usage="$usage        --fc                 - prints fortran compiler\n";
usage="$usage        --fcopts             - prints fortran compiler options\n";
usage="$usage        --fl                 - prints fortran linker\n";
usage="$usage        --opts               - prints PaStiX compiling options\n";
usage="$usage        --vers               - prints PaStiX version\n";
usage="$usage        --blas               - prints blas choosen in config.in\n";

if [ $# = 0 ]
then
    echo "Librairies               : $LIBS" ;
    echo "Librairies with murge    : $LIBS_MURGE";
    echo "Incs                     : $INC" ;
    echo "C Compiler               : $CC" ;
    echo "C++ Compiler             : $CXX" ;
    echo "C Compiler options       : $CCOPT" ;
    echo "Fortran Compiler         : $FC" ;
    echo "Fortran Compiler options : $FCOPT" ;
    echo "C Linker                 : $CL" ;
    echo "Fortran Linker           : $FL" ;
    echo "Options                  : $OPTS" ;
    echo "Version                  : $VERSION" ;
    echo "Blas                     : $BLAS" ;
elif [ $# = 1 ]
then
    case $1 in
        --libs)
            echo $LIBS;;
        --libs_murge)
            echo $LIBS_MURGE;;
        --incs)
            echo $INC;;
        --cc)
            echo $CC;;
        --cxx)
            echo $CXX;;
        --ccopts)
            echo $CCOPT;;
        --fc)
            echo $FC;;
        --fcopts)
            echo $FCOPT;;
        --cl)
            echo $CL;;
        --fl)
            echo $FL;;
        --opts)
            echo $OPTS;;
        --blas)
            echo $BLAS;;
        --vers)
            echo $VERSION;;

        *)
            echo -e $usage
    esac;
else
    echo -e $usage
fi;
