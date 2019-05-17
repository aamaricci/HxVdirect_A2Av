FC=mpif90

HERE   = `pwd`
PEXE1= test_MPI_LANCZOS
PEXE2= test_MPI_ARPACK

DIREXE = ~/.bin

##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90

OBJS     = COMMON_VARS.o SETUP.o MATVEC_PRODUCT.o  SF_SP_LINALG.o 

FFLAG = -O3 -ffast-math -march=native -funroll-all-loops
DFLAG = -O0 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none


INCS=$(shell pkg-config --cflags dmft_tools scifor)
ARGS=$(shell pkg-config --libs   dmft_tools scifor)


all: $(OBJS)
	${FC} ${FFLAG} -cpp -D_MPI ${OBJS} ${PEXE1}.f90 -o ${DIREXE}/${PEXE1} ${INCS} ${ARGS}
	${FC} ${FFLAG} -cpp -D_MPI ${OBJS} ${PEXE2}.f90 -o ${DIREXE}/${PEXE2} ${INCS} ${ARGS}
	@echo ""


debug: FFLAG=${DFLAG} 
debug: $(OBJS)
	${FC} ${FFLAG} -cpp -D_MPI ${OBJS} ${PEXE1}.f90 -o ${DIREXE}/${PEXE1} ${INCS} ${ARGS}
	${FC} ${FFLAG} -cpp -D_MPI ${OBJS} ${PEXE2}.f90 -o ${DIREXE}/${PEXE2} ${INCS} ${ARGS}
	@echo ""


objs: ${OBJS}

clean: 
	@echo "Cleaning:"
	@rm -f *.mod
	@rm -f *.o
	@rm -f *~
	@echo ""

.f90.o:	
	${FC} ${FFLAG} -cpp -D_MPI -c $< $(INCS)


