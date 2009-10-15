
OPT			+=	-DTREE
#OPT			+=	-DLINKLIST

CC			=	gcc
OPTIMIZE	=	-O2	-W -Wall
GSL_INCL	=	-I/usr/local/include
GSL_LIBS	=	-L/usr/local/lib
FFTW_INCL	=	-I/usr/local/include
FFTW_LIBS	=	-L/usr/local/lib

OPTIONS		=	$(OPTIMIZE) $(OPT)

EXEC		=	cuco

OBJS		=	allvars.o begrun.o init.o read_ic.o \
				main.o driftfac.o longrange.o pm_periodic.o \
				run.o predict.o accel.o io.o timestep.o

ifeq (TREE, $(findstring TREE, $(OPT)))
	OBJS	+=	gravtree.o
endif

ifeq (LINKLIST, $(findstring LINKLIST, $(OPT)))
	OBJS	+=	gravlist.o
endif

INCL		=	allvars.h proto.h Makefile

CFLAGS		=	$(OPTIONS) $(GSL_INCL) $(FFTW_INCL)

LIBS		=	$(GSL_LIBS) -lgsl -lgslcblas -lm \
				$(FFTW_LIBS) -lrfftw -lfftw

$(EXEC)	:	$(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)

$(OBJS)	:	$(INCL)

clean:
	rm -f $(OBJS) $(EXEC)
