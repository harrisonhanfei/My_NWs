#  Declaration
CXXFLAGS = -Wall -v -O3 -std=c++0x
PARFLAGS = -fopenmp
CC	=	gcc
C++	=	g++
F77    	=	f77

.SUFFIXES:

.SUFFIXES: .f .F .cpp .c .o .h

%.o: %.f
	${F77} ${FFLAGS} $< -c
	
%.o: %.F
	${F77} ${FFLAGS} $< -c

%.o: %.cpp %.h
	${C++} ${CXXFLAGS} $< -c

%.o: %.c %.h
	${CC} ${CFLAGS} $< -c


objects =  App_Network_2D.o Background_vectors.o Cutoff_Wins.o Fem_3D.o\
			  Gauss.o GenNetwork_2D.o Geometry_3D.o \
           	  Hns.o Input_Reader.o MathMatrix.o Tecplot_Export.o MainPro.o \
	                   
necn : $(objects)        
	${C++} ${PARFLAGS} -o necn $(objects)
	              
MainPro.o : MainPro.cpp
	${C++} ${CXXFLAGS} MainPro.cpp -c

# αĿ���ļ�
.PHONY : clean
clean :
	-rm necn $(objects)
