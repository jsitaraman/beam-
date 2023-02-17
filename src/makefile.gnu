MODULENAME=beam
F90= gfortran
CC = gcc
CXX= g++
AR = ar -rvs
FFLAGS =  -fdefault-real-8 -fconvert=big-endian -O2
INCLUDES =  
CFLAGS =  -g 
OBJF90 = 
OBJECTS = beam.o EulerBernoulliBeam.o main.o
LDFLAGS = 

default: $(OBJECTS) $(OBJF90) *.h
	$(CXX) $(CFLAGS) $(OBJECTS) $(OBJF90) $(OBJEXEC) $(LDFLAGS) -lm -o $(MODULENAME).exe

clean : 
	rm -r *.o *.mod $(MODULENAME).exe

%.o:%.cu
	$(CUC)  $(CFLAGS) -c $< -o $*.o
%.o:%.C
	$(CXX) $(CFLAGS) -c $< -o $*.o
%.o:%.cpp
	$(CXX) $(INCLUDES) $(CFLAGS) -c $< -o $*.o
%.o:%.c
	$(CC) $(INCLUDES) $(CFLAGS) -c $< -o $*.o
%.o:%.F90
	$(F90) $(FFLAGS) -c $< -o $*.o
%.o:%.f90
	$(F90) $(FFLAGS) -c $< -o $*.o
%.o:%.f
	$(F90) $(FFLAGS) -c $< -o $*.o
