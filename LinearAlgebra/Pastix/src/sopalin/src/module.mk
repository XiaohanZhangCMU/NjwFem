MODULE:=sopalin

MY_SRC:= bordi.c		\
	 sopalin_thread.c


MY_SRC_MULT_ARCH :=	sopalin_init.c		\
			sopalin_option.c	\
			pastix.c		\
			ooc.c			\
			debug_dump.c		\
			csc_intern_io.c	\
			csc_intern_compute.c	\
			csc_intern_build.c	\
			csc_intern_solve.c	\
			csc_intern_updown.c	\
			tools.c		\
			coefinit.c		\
			sopalin_init.c		\
			debug_dump.c		\
			csc_intern_io.c	\
			csc_intern_build.c	\
			csc_intern_solve.c	\
			csc_intern_updown.c	\
			ooc.c			\
			sopalin3d.c		\
			pastix.c		\
			pastix_fortran.c	\
			cscd_utils.c		\
			csc_utils.c		\
			cscd_utils_fortran.c    \
			sparse_gemm_cpu.c
ifeq (1, $(words $(findstring -DWITH_MAGMABLAS, ${CCOPT})))
MY_SRC_MULT_ARCH := ${MY_SRC_MULT_ARCH}	zgetrf_stapiv_gpu.c
MY_SRC_MULT_ARCH := ${MY_SRC_MULT_ARCH}	zpotrf_stapiv_gpu.c
MY_SRC_MULT_ARCH := ${MY_SRC_MULT_ARCH}	zhetrf_stapiv_gpu.c
MY_SRC_MULT_ARCH := ${MY_SRC_MULT_ARCH}	zsytrf_stapiv_gpu.c
endif

ifeq (1, $(words $(findstring -DWITH_STARPU, ${CCOPT})))
ifeq (0, $(words $(findstring -DFORCE_NO_CUDA, ${CCOPT})))
MY_SRC_MULT_ARCH:= $(MY_SRC_MULT_ARCH) sparse_gemm.cu
MY_SRC_MULT_ARCH:= ${MY_SRC_MULT_ARCH} geadd_cuda.cu getra_cuda.cu
ifeq (1, $(words $(findstring -DCUDA_SM_VERSION=20, ${CCOPT})))
MY_SRC_MULT_ARCH:= $(MY_SRC_MULT_ARCH) sparse_gemm_fermi.cu
MY_SRC_MULT_ARCH:= $(MY_SRC_MULT_ARCH) sparse_gemdm_fermi.cu
endif
MY_SRC_MULT_ARCH:= $(MY_SRC_MULT_ARCH) sparse_gemdm.cu
endif
endif

MY_SRC_FACT_TYPE :=	sopalin3d.c		\
			starpu_submit_tasks.c	\
			csc_intern_compute.c	\
			starpu_updo.c

MY_SRC_MURGE		:= murge_fortran.c murge.c
MY_SRC_MURGE_MULT_ARCH  :=


OJTDIR:=$(MODULE)/obj/$(HOSTARCH)
DEPDIR:=$(MODULE)/dep
SRCDIR:=$(MODULE)/src

MY_SRC		 := $(patsubst %,$(SRCDIR)/%,$(MY_SRC))
MY_SRC_MULT_ARCH := $(patsubst %,$(SRCDIR)/%,$(MY_SRC_MULT_ARCH))
MY_SRC_FACT_TYPE := $(patsubst %,$(SRCDIR)/%,$(MY_SRC_FACT_TYPE))

MY_SRC_MULT_ARCH := $(MY_SRC_MULT_ARCH)		\
	$(MY_SRC_FACT_TYPE)				\
	$(patsubst %.c,%_ge.c,$(MY_SRC_FACT_TYPE))	\
	$(patsubst %.c,%_sy.c,$(MY_SRC_FACT_TYPE))      \
	$(patsubst %.c,%_he.c,$(MY_SRC_FACT_TYPE))

MY_OBJ :=  $(patsubst %.c,%.o,$(filter %.c,$(MY_SRC)))			\
	$(patsubst %.cu,%.o,$(filter %.cu,$(MY_SRC)))			\
	$(patsubst %.c,%.o,$(filter %.c,$(MY_SRC_MULT_ARCH)))		\
	$(patsubst %.cu,%.o,$(filter %.cu,$(MY_SRC_MULT_ARCH)))	\
	$(patsubst %.c,%_s.o,$(filter %.c,$(MY_SRC_MULT_ARCH)))	\
	$(patsubst %.cu,%_s.o,$(filter %.cu,$(MY_SRC_MULT_ARCH)))	\
	$(patsubst %.c,%_d.o,$(filter %.c,$(MY_SRC_MULT_ARCH)))	\
	$(patsubst %.cu,%_d.o,$(filter %.cu,$(MY_SRC_MULT_ARCH)))	\
	$(patsubst %.c,%_c.o,$(filter %.c,$(MY_SRC_MULT_ARCH)))	\
	$(patsubst %.cu,%_c.o,$(filter %.cu,$(MY_SRC_MULT_ARCH)))	\
	$(patsubst %.c,%_z.o,$(filter %.c,$(MY_SRC_MULT_ARCH)))	\
	$(patsubst %.cu,%_z.o,$(filter %.cu,$(MY_SRC_MULT_ARCH)))

ifeq (0, $(words $(findstring -DCUDA_SM_VERSION=20, ${CCOPT})))
ifeq (2, $(words $(findstring -DTYPE_COMPLEX, ${CCOPT})	\
		 $(findstring -DPREC_DOUBLE, ${CCOPT})))
MY_OBJ := $(filter-out %/sparse_gemm.o, ${MY_OBJ})
MY_OBJ := $(filter-out %/sparse_gemdm.o, ${MY_OBJ})
endif
MY_OBJ := $(filter-out %/sparse_gemm_z.o, ${MY_OBJ})
MY_OBJ := $(filter-out %/sparse_gemdm_z.o, ${MY_OBJ})
endif

MY_OBJ := $(patsubst $(SRCDIR)/%,$(OJTDIR)/%,$(MY_OBJ))

OBJ_LIB := $(OBJ_LIB) $(MY_OBJ)
DEP_LIB := $(DEP_LIB) $(patsubst $(OJTDIR)/%.o,$(DEPDIR)/%.d,$(MY_OBJ))

MY_SRC_MURGE		:= $(patsubst %,$(SRCDIR)/%,$(MY_SRC_MURGE))
MY_SRC_MURGE_MULT_ARCH  := $(patsubst %,$(SRCDIR)/%,$(MY_SRC_MURGE_MULT_ARCH))

MY_OBJ := $(patsubst %.c,%.o,$(filter %.c,$(MY_SRC_MURGE)))		\
	$(patsubst %.c,%.o,$(filter %.c,$(MY_SRC_MURGE_MULT_ARCH)))	\
	$(patsubst %.c,%_s.o,$(filter %.c,$(MY_SRC_MURGE_MULT_ARCH)))	\
	$(patsubst %.c,%_d.o,$(filter %.c,$(MY_SRC_MURGE_MULT_ARCH)))	\
	$(patsubst %.c,%_c.o,$(filter %.c,$(MY_SRC_MURGE_MULT_ARCH)))	\
	$(patsubst %.c,%_z.o,$(filter %.c,$(MY_SRC_MURGE_MULT_ARCH)))
MY_OBJ := $(patsubst $(SRCDIR)/%,$(OJTDIR)/%,$(MY_OBJ))

OBJ_MURGE := $(OBJ_MURGE) $(MY_OBJ)
DEP_MURGE := $(DEP_MURGE) $(patsubst $(OJTDIR)/%.o,$(DEPDIR)/%.d,$(MY_OBJ))


ifeq (1, $(words $(findstring debug, ${MAKECMDGOALS})))
NVCCOPT := ${NVCCOPT} -G
endif

include all_mp_rules.mk


MURGE_FORTRAN_INC = common_pastix.h

$(SRCDIR)/murge_fortran.c:	murge/include/murge.h			\
				$(SRCDIR)/murge_pastix_fortran.c	\
					murge/scripts/geninterface.pl
	-$(RM) $@
	for f in $(MURGE_FORTRAN_INC) ; do			\
		echo "#include \""$${f}"\"" >> $@;		\
	done
	echo "#ifdef FORCE_NOMPI"   >> $@
	echo "#include \"nompi.h\"" >> $@
	echo "#else		  " >> $@
	echo "#include <mpi.h>	  " >> $@
	echo "#endif              " >> $@
	murge/scripts/geninterface.pl -f $< >> $@
	cat sopalin/src/murge_pastix_fortran.c  >> $@

NODESIZE2 = 16
TEST_LIBS = -DNODESIZE=8 -L${BUILD_LIBDIR} -lpastix ${EXTRALIB} ${BLASLIB} $

test: test_sparse_gemm test_geadd test_getra test_getrf
test_sparse_build: 	${BUILD_TSTDIR}/test_sparse_gemm 	\
			${BUILD_TSTDIR}/test_sparse_gemm_s	\
			${BUILD_TSTDIR}/test_sparse_gemm_d 	\
			${BUILD_TSTDIR}/test_sparse_gemm_c 	\
			${BUILD_TSTDIR}/test_sparse_gemm_z
.PHONY: test_sparse_gemm
test_sparse_gemm: test_sparse_build
	${CUDA_MEMCHECK} ${BUILD_TSTDIR}/test_sparse_gemm   1000 100 100
	${CUDA_MEMCHECK} ${BUILD_TSTDIR}/test_sparse_gemm_s 1000 100 100
	${CUDA_MEMCHECK} ${BUILD_TSTDIR}/test_sparse_gemm_d 1000 100 100
	${CUDA_MEMCHECK} ${BUILD_TSTDIR}/test_sparse_gemm_c 1000 100 100
	${CUDA_MEMCHECK} ${BUILD_TSTDIR}/test_sparse_gemm_z 1000 100 100

test_geadd_build: 	${BUILD_TSTDIR}/test_geadd 	\
			${BUILD_TSTDIR}/test_geadd_s	\
			${BUILD_TSTDIR}/test_geadd_d 	\
			${BUILD_TSTDIR}/test_geadd_c 	\
			${BUILD_TSTDIR}/test_geadd_z
.PHONY: test_geadd
test_geadd: test_geadd_build
	${CUDA_MEMCHECK} ${BUILD_TSTDIR}/test_geadd   1000 10000
	${CUDA_MEMCHECK} ${BUILD_TSTDIR}/test_geadd_s 1000 10000
	${CUDA_MEMCHECK} ${BUILD_TSTDIR}/test_geadd_d 1000 10000
	${CUDA_MEMCHECK} ${BUILD_TSTDIR}/test_geadd_c 1000 10000
	${CUDA_MEMCHECK} ${BUILD_TSTDIR}/test_geadd_z 1000 10000

test_getra_build: 	${BUILD_TSTDIR}/test_getra 	\
			${BUILD_TSTDIR}/test_getra_s	\
			${BUILD_TSTDIR}/test_getra_d 	\
			${BUILD_TSTDIR}/test_getra_c 	\
			${BUILD_TSTDIR}/test_getra_z
.PHONY: test_getra
test_getra: test_getra_build
	${CUDA_MEMCHECK} ${BUILD_TSTDIR}/test_getra   1000 10000
	${CUDA_MEMCHECK} ${BUILD_TSTDIR}/test_getra_s 1000 10000
	${CUDA_MEMCHECK} ${BUILD_TSTDIR}/test_getra_d 1000 10000
	${CUDA_MEMCHECK} ${BUILD_TSTDIR}/test_getra_c 1000 10000
	${CUDA_MEMCHECK} ${BUILD_TSTDIR}/test_getra_z 1000 10000

test_getrf_build: 	${BUILD_TSTDIR}/test_getrf 	\
			${BUILD_TSTDIR}/test_getrf_s	\
			${BUILD_TSTDIR}/test_getrf_d 	\
			${BUILD_TSTDIR}/test_getrf_c 	\
			${BUILD_TSTDIR}/test_getrf_z
.PHONY: test_getrf
test_getrf: test_getrf_build
	${CUDA_MEMCHECK} ${BUILD_TSTDIR}/test_getrf   1000 10000
	${CUDA_MEMCHECK} ${BUILD_TSTDIR}/test_getrf_s 1000 10000
	${CUDA_MEMCHECK} ${BUILD_TSTDIR}/test_getrf_d 1000 10000
	${CUDA_MEMCHECK} ${BUILD_TSTDIR}/test_getrf_c 1000 10000
	${CUDA_MEMCHECK} ${BUILD_TSTDIR}/test_getrf_z 1000 10000

${BUILD_TSTDIR}/% :	$(SRCDIR)/%.c ${BUILD_LIBDIR}/libpastix.a
	${MPCCPRG} ${CCOPT} $< -o $@ ${CCINC} ${TEST_LIBS} $(CCTYPESFLT);
${BUILD_TSTDIR}/% :	$(OJTDIR)/%.c ${BUILD_LIBDIR}/libpastix.a
	${MPCCPRG} ${CCOPT} $< -o $@ ${CCINC} ${TEST_LIBS} $(CCTYPESFLT);

test_gemm: 	${BUILD_TSTDIR}/test_gemm \
		${BUILD_TSTDIR}/test_gemm_s \
		${BUILD_TSTDIR}/test_gemm_d \
		${BUILD_TSTDIR}/test_gemm_c \
		${BUILD_TSTDIR}/test_gemm_z
	${BUILD_TSTDIR}/test_gemm_d

test_trsm: 	${BUILD_TSTDIR}/test_trsm \
		${BUILD_TSTDIR}/test_trsm_s \
		${BUILD_TSTDIR}/test_trsm_d \
		${BUILD_TSTDIR}/test_trsm_c \
		${BUILD_TSTDIR}/test_trsm_z
	${BUILD_TSTDIR}/test_trsm_d

test_geam: 	${BUILD_TSTDIR}/test_geam \
		${BUILD_TSTDIR}/test_geam_s \
		${BUILD_TSTDIR}/test_geam_d \
		${BUILD_TSTDIR}/test_geam_c \
		${BUILD_TSTDIR}/test_geam_z
	${BUILD_TSTDIR}/test_geam_d


test_pof: 	${BUILD_TSTDIR}/test_pof \
		${BUILD_TSTDIR}/test_pof_s \
		${BUILD_TSTDIR}/test_pof_d \
		${BUILD_TSTDIR}/test_pof_c \
		${BUILD_TSTDIR}/test_pof_z
	${BUILD_TSTDIR}/test_pof_d

test_axpy: 	${BUILD_TSTDIR}/test_axpy \
		${BUILD_TSTDIR}/test_axpy_s \
		${BUILD_TSTDIR}/test_axpy_d \
		${BUILD_TSTDIR}/test_axpy_c \
		${BUILD_TSTDIR}/test_axpy_z
	${BUILD_TSTDIR}/test_axpy_d

test_copy: 	${BUILD_TSTDIR}/test_copy \
		${BUILD_TSTDIR}/test_copy_s \
		${BUILD_TSTDIR}/test_copy_d \
		${BUILD_TSTDIR}/test_copy_c \
		${BUILD_TSTDIR}/test_copy_z
	${BUILD_TSTDIR}/test_copy_d


test_scal: 	${BUILD_TSTDIR}/test_scal \
		${BUILD_TSTDIR}/test_scal_s \
		${BUILD_TSTDIR}/test_scal_d \
		${BUILD_TSTDIR}/test_scal_c \
		${BUILD_TSTDIR}/test_scal_z
	${BUILD_TSTDIR}/test_scal_d

test_gemv: 	${BUILD_TSTDIR}/test_gemv \
		${BUILD_TSTDIR}/test_gemv_s \
		${BUILD_TSTDIR}/test_gemv_d \
		${BUILD_TSTDIR}/test_gemv_c \
		${BUILD_TSTDIR}/test_gemv_z
	${BUILD_TSTDIR}/test_gemv_d

test_trsv: 	${BUILD_TSTDIR}/test_trsv \
		${BUILD_TSTDIR}/test_trsv_s \
		${BUILD_TSTDIR}/test_trsv_d \
		${BUILD_TSTDIR}/test_trsv_c \
		${BUILD_TSTDIR}/test_trsv_z
	${BUILD_TSTDIR}/test_trsv_d

pingpong_1: 	${BUILD_TSTDIR}/pingpong_1
	mpirun -np 2 ${BUILD_TSTDIR}/pingpong_1

pingpong_1e: 	${BUILD_TSTDIR}/pingpong_1e
	mpirun -np ${NODESIZE2} ${BUILD_TSTDIR}/pingpong_1e

pingpong_2: 	${BUILD_TSTDIR}/pingpong_2
	mpirun -np 4 ${BUILD_TSTDIR}/pingpong_2

pingpong_2e: 	${BUILD_TSTDIR}/pingpong_2e
	mpirun -np ${NODESIZE2} ${BUILD_TSTDIR}/pingpong_2e

pingpong_4: 	${BUILD_TSTDIR}/pingpong_4
	mpirun -np 8 ${BUILD_TSTDIR}/pingpong_4

pingpong_4e: 	${BUILD_TSTDIR}/pingpong_4e
	mpirun -np ${NODESIZE2} ${BUILD_TSTDIR}/pingpong_4e

pingpong_8: 	${BUILD_TSTDIR}/pingpong_8
	mpirun -np 16 ${BUILD_TSTDIR}/pingpong_8

pingpong_8e: 	${BUILD_TSTDIR}/pingpong_8e
	mpirun -np ${NODESIZE2} ${BUILD_TSTDIR}/pingpong_8e

pingpong_n1: 	${BUILD_TSTDIR}/pingpong_n1
	mpirun -np 2 ${BUILD_TSTDIR}/pingpong_n1

pingpong_n1e: 	${BUILD_TSTDIR}/pingpong_n1e
	mpirun -np ${NODESIZE2} ${BUILD_TSTDIR}/pingpong_n1e

pingpong_n2: 	${BUILD_TSTDIR}/pingpong_n2
	mpirun -np 4 ${BUILD_TSTDIR}/pingpong_n2

pingpong_n2e: 	${BUILD_TSTDIR}/pingpong_n2e
	mpirun -np ${NODESIZE2} ${BUILD_TSTDIR}/pingpong_n2e

pingpong_n4: 	${BUILD_TSTDIR}/pingpong_n4
	mpirun -np 8 ${BUILD_TSTDIR}/pingpong_n4

pingpong_n4e: 	${BUILD_TSTDIR}/pingpong_n4e
	mpirun -np ${NODESIZE2} ${BUILD_TSTDIR}/pingpong_n4e

pingpong_n8: 	${BUILD_TSTDIR}/pingpong_n8
	mpirun -np 16 ${BUILD_TSTDIR}/pingpong_n8

pingpong_n8e: 	${BUILD_TSTDIR}/pingpong_n8e
	mpirun -np ${NODESIZE2} ${BUILD_TSTDIR}/pingpong_n8e
