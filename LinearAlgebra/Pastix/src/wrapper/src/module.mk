MODULE=wrapper

MY_SRC:=
MY_SRC_MULT_ARCH :=
MY_SRC_FACT_TYPE :=

MY_SRC_MURGE 		:= 
MY_SRC_MURGE_MULT_ARCH  := 
include all_modules.mk
include all_rules.mk

$(OJTDIR) : 
	mkdir -p $@

$(OJTDIR)/pastix.i : $(SRCDIR)/pastix.i $(BUILD_INCDIR)/pastix.h $(filter-out $(wildcard $(OJTDIR)), $(OJTDIR))

$(OJTDIR)/%.i : $(SRCDIR)/%.i
	cat $< > $@.tmp2
	cat $(BUILD_INCDIR)/pastix.h >> $@.tmp2

ifeq (0, $(DBL_FLG))
ifeq (0, $(CPLX_FLG))
	sed -e 's/pastix_float_t/float         /g' $@.tmp2 > $@.tmp
else
	sed -e 's/pastix_float_t/float complex /g' $@.tmp2 > $@.tmp
endif
else
ifeq (0, $(CPLX_FLG))
	sed -e 's/pastix_float_t/double         /g' $@.tmp2 > $@.tmp
else
	sed -e 's/pastix_float_t/double complex /g' $@.tmp2 > $@.tmp
endif
endif

ifeq (1, $(INT_FLG))
	sed -e 's/pastix_int_t/int64_t     /g' $@.tmp > $@
else
ifeq (2, $(INT_FLG))
	sed -e 's/pastix_int_t/int32_t     /g' $@.tmp > $@
else
ifeq (3, $(INT_FLG))
	sed -e 's/pastix_int_t/long        /g' $@.tmp > $@
else
	sed -e 's/pastix_int_t/int         /g' $@.tmp > $@
endif
endif
endif
	rm $@.tmp
	rm $@.tmp2

$(OJTDIR)/%.c :$(filter-out $(wildcard $(OJTDIR)), $(OJTDIR))
$(OJTDIR)/%.o :$(filter-out $(wildcard $(OJTDIR)), $(OJTDIR))

$(OJTDIR)/pastix_wrap_python.c : $(OJTDIR)/pastix.i $(BUILD_INCDIR)/pastix.h 
	swig -python -I$(MPI4PY_INC)/mpi4py -I$(BUILD_INCDIR) -o $@ $<
	mv $(patsubst %.i,%.py, $<) $(LIBDIR)

$(OJTDIR)/pastix_wrap_octave.c : $(OJTDIR)/pastix.i $(BUILD_INCDIR)/pastix.h 
	swig -octave -I$(MPI4PY_INC)/mpi4py -I$(BUILD_INCDIR)  -o $@ $<

$(OJTDIR)/pastix_wrap_python.o : $(OJTDIR)/pastix_wrap_python.c 
	$(MPCCPRG) -I$(MPI4PY_INC) $(shell ${BUILD_BINDIR}/pastix-conf --incs) -c -fpic $< -I${PYTHON_INC} -o $@


$(LIBDIR)/_pastix.so : $(BUILD_LIBDIR)/libpastix$(LIB) ${BUILD_BINDIR}/pastix-conf 
$(LIBDIR)/_pastix.so : $(filter-out $(wildcard $(LIBDIR)), $(LIBDIR))
$(LIBDIR)/_pastix.so : $(OJTDIR)/pastix_wrap_python.o
	$(MPCCPRG) -shared $< $(shell ${BUILD_BINDIR}/pastix-conf --libs) $(shell ${BUILD_BINDIR}/pastix-conf --blas)  -o $@
