# MakeFile

#CXX = icpc
CXX = g++
#CFLAGS = -O3 -std=c++14 -xhost -qopenmp -I../BaseLibs -ipo
CFLAGS = -O3 -std=c++14 -mavx2 -fopenmp -I../BaseLibs
#OPTS = -mkl
OPTS = -lgomp
LIBS = -L../BaseLibs -lBaseFuncs -lTestObjFuncs
OBJDIR = ./obj

SOURCES  = $(wildcard *.cpp) $(wildcard Generators/*.cpp) $(wildcard Optimizers/*.cpp) 
OBJECTS  = $(addprefix $(OBJDIR)/, $(SOURCES:.cpp=.o))

#target file
TARGET = OptimizerTest

$(TARGET): $(OBJECTS) 
	$(CXX) -o $@ $^ $(LDFLAGS) $(LIBS) $(OPTS)
$(OBJDIR)/%.o: %.cpp
	@if [ ! -d $(OBJDIR) ]; \
		then echo "mkdir -p $(OBJDIR) $(OBJDIR)/Generators  $(OBJDIR)/Optimizers"; mkdir -p $(OBJDIR); mkdir -p $(OBJDIR)/Generators;  mkdir -p $(OBJDIR)/Optimizers;\
	fi
	@[ -d $(OBJDIR) ]
	$(CXX) -c $(CFLAGS) $(INDS) -o $@ $<

clean:
	rm $(OBJDIR)/*.o
	rm $(OBJDIR)/Generators/*.o
	rm $(OBJDIR)/Optimizers/*.o
	rm $(TARGET)
depend:
	makedepend -- -Y -- $(SOURCES)


# DO NOT DELETE
