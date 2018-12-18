CXX = icpc
FLAGS = -c -O2

OBJECTS = 3D_RegStokes_case_1.o

ctreecode.exe : $(OBJECTS)
        $(CXX) -o ctreecode *.o -lm

3D_RegStokes_case_1.o : 3D_RegStokes_case_1.cpp
        $(CXX) $(FLAGS) 3D_RegStokes_case_1.cpp

clean :
        \rm -f .o ~ ctreecode