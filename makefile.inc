CFLAGS = -O3 -Wall

#HOME = /home/xiaohan/Planet/Work/PileUps/Codes
NJWFEM_DIR = $(NJW_DIR)

BOOST_INC=
#BOOST_INC=-I/usr/include/boost-1_33_1/

#PETSC 2.3
#include ${PETSC_DIR}/bmake/common/base
#PETSC 3.0
#include ${PETSC_DIR}/conf/base
#PETSC 3.5
#include ${PETSC_DIR}/conf/variables
#include ${PETSC_DIR}/conf/rules

SPOOLES_DIR = $(NJWFEM_DIR)/LinearAlgebra/Spooles
SPOOLES_I = -I$(SPOOLES_DIR)/LinSol
SPOOLES_A = $(SPOOLES_DIR)/LinSol/srcMT/BridgeMT.a $(SPOOLES_DIR)/LinSol/srcST/Bridge.a  $(SPOOLES_DIR)/MT/src/spoolesMT.a $(SPOOLES_DIR)/spooles.a

MUMPS_DIR =  $(NJWFEM_DIR)/LinearAlgebra/Mumps
MUMPS_I = -I$(MUMPS_DIR)/include -I$(MUMPS_DIR)/libseq
MUMPS_A = $(MUMPS_DIR)/lib/libdmumps.a $(MUMPS_DIR)/lib/libmumps_common.a
MUMPS_L = -L$(MUMPS_DIR)/lib/ -lpord  -L$(MUMPS_DIR)/libseq -lmpiseq -lblas

COMMON_INC = $(DIR3_INC)  $(SPOOLES_I) $(MUMPS_I) $(SUPERLU_I) $(BOOST_INC) $(PETSC_INCLUDE)

COMMON_LIB = $(DIR3_LIB)  $(SPOOLES_A) $(SUPERLU_A) $(MUMPS_A) $(MUMPS_L) ${PETSC_KSP_LIB} $(PARDISO_LIB) 
#-llapack -lpthread -lm -lgfortran

COMMON_O = $(NJWFEM_DIR)/FiniteElement/lagrangeBasisFns.o

PASTIX_DIR =  $(NJWFEM_DIR)/LinearAlgebra/Pastix/install
COMMON_INC := $(COMMON_INC) -I$(PASTIX_DIR)
COMMON_LIB := $(COMMON_LIB) $(PASTIX_DIR)/libpastix.a
COMMON_LIB := $(COMMON_LIB) -L$(PASTIX_DIR) -L$(NJWFEM_DIR)/LinearAlgebra/Scotch/lib -lpastix -lscotch -lscotcherrexit -lrt -lblas
