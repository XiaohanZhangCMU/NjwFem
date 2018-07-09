MODULE:=blend
MY_SRC = task.c				\
	fanboth2.c			\
	splitfunc.c			\
	blend_symbol_cost.c		\
	distribPart.c			\
	simu.c				\
	extendVector.c 			\
	cost.c	 			\
	costfunc.c 			\
	write_ps.c			\
	blendctrl.c			\
	param_blend.c			\
	partbuild.c			\
	extrastruct.c 			\
	splitpart.c			\
	splitpartlocal.c		\
	blend.c 			\
	solver_check.c			\
	elimin.c 			\
	solverRealloc.c 		\
	eliminfunc.c 			\
	assemblyGener.c 		\
	bulles.c			\
	solver_io.c			\
	solverMatrixGen.c		\
	queue.c				\
	symbolrand.c

MY_SRC_MULT_ARCH = 

OJTDIR:=blend/obj/$(HOSTARCH)
SRCDIR:=blend/src
DEPDIR:=blend/dep

include all_modules.mk
include all_rules.mk
