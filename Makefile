OFILES = 3D_RegStokes_case_1.o
TARGET = 3D_RegStokes_case_1.exe

OFILESA = 3D_RegStokes_case_2.o
TARGETA = 3D_RegStokes_case_2.exe

OFILESB = 3D_SingStokes_Taylor_case_1.o
TARGETB = 3D_SingStokes_Taylor_case_1.exe

OFILESC = 3D_SingStokes_Taylor_case_2.o
TARGETC = 3D_SingStokes_Taylor_case_2.exe


CXX = icpc
CXXFLAGS = -O2


all: $(TARGET) $(TARGETA) $(TARGETB) $(TARGETC)

$(TARGET): $(OFILES)
	$(CXX) $(CXXFLAGS) $(OFILES) -o $@ -lm

$(TARGETA): $(OFILESA)
	$(CXX) $(CXXFLAGS) $(OFILESA) -o $@ -lm

$(TARGETB): $(OFILESB)
	$(CXX) $(CXXFLAGS) $(OFILESB) -o $@ -lm

$(TARGETC): $(OFILESC)
	$(CXX) $(CXXFLAGS) $(OFILESC) -o $@ -lm

clean:
	rm -rf $(OFILES) $(TARGET) $(OFILESA) $(TARGETA) 
	rm -rf $(OFILESB) $(TARGETB) $(OFILESC) $(TARGETC)


3D_RegStokes_case_1.o : 3D_RegStokes_case_1.cpp
3D_RegStokes_case_2.o : 3D_RegStokes_case_2.cpp
3D_SingStokes_Taylor_case_1.o : 3D_SingStokes_Taylor_case_1.cpp
3D_SingStokes_Taylor_case_2.o : 3D_SingStokes_Taylor_case_2.cpp
