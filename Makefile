#--------------------------------------- Interpolation methods
OPT += -DNGP
OPT += -DDECAY_DARK_MATTER

#--------------------------------------- Others
#OPT += -DNOFILE_OUTPUT
#--------------------------------------- Select target computer
#SYSTYPE="Workstation"
#SYSTYPE="dlcheng"
#SYSTYPE="ITSC"
SYSTYPE="Mac"
#SYSTYPE="uv2000"

#--------------------------------------- Adjust settings for target computer
ifeq ($(SYSTYPE),"uv2000")
CC       =   gcc
OPTIMIZE =   -O3 -Wall
GSL_INCL =  -I/home/dlcheng/lib/gsl-1.6/include
GSL_LIBS =  -L/home/dlcheng/lib/gsl-1.6/lib
endif

ifeq ($(SYSTYPE),"Mac")
CC       =   /usr/local/Cellar/gcc/5.2.0/bin/gcc-5  
OPTIMIZE =   -O3 -Wall
GSL_LIBS=   -L/usr/local/Cellar/gsl/1.16/lib 
GSL_INCL =  -I/usr/local/Cellar/gsl/1.16/include
endif

ifeq ($(SYSTYPE),"Workstation")
CC       =   gcc   
OPTIMIZE =   -O3 -Wall 
GSL_INCL =  -I/home/dlcheng/Install/gsl/include
GSL_LIBS =  -L/home/dlcheng/Install/gsl/lib
endif

ifeq ($(SYSTYPE),"dlcheng")
CC       =   gcc   
OPTIMIZE =   -O3 -Wall 
GSL_INCL =  -I/home/dalong/Install/gsl/include
GSL_LIBS =  -L/home/dalong/Install/gsl/lib
endif

ifeq ($(SYSTYPE),"ITSC")
CC       =   gcc   
OPTIMIZE =   -O3 -Wall
GSL_LIBS=   -L/users/s0902248/Lib/gsl-1.9/lib  -Wl,"-R /users/s0902248/Lib/gsl-1.9/lib"
GSL_INCL =  -I/users/s0902248/Lib/gsl-1.9/include
endif



OPTIONS =  $(OPTIMIZE) $(OPT) 

EXEC   = Gpower

OBJS   = allvars.o assignment.o fft.o init.o load_gadget.o main.o memory_control.o out_put.o \
         power_bin.o warn_end.o

INCL   = allvars.h proto.h define.h Makefile


CFLAGS = $(OPTIONS) $(GSL_INCL) -fopenmp


LIBS   = $(GSL_LIBS) -lgsl -lgslcblas -lm -fopenmp

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC) *.gch

