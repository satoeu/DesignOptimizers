# MakeFile

#CXX = icpc
#ARX = xiar
#CFLAGS = -O3 -std=c++14 -xhost -qopenmp -I../ -ipo
#OPTS = -mkl
#LIBS = -L../

CXX = g++
ARX = ar
CFLAGS = -O3 -std=c++14 -I../ -mavx2 -fopenmp
OPTS = -lgomp
LIBS = -L../

OBJDIR = ./obj

SOURCES  = $(wildcard *.cpp)
OBJECTS  = $(addprefix $(OBJDIR)/, $(SOURCES:.cpp=.o))

#target file
TARGET = ../libTestObjFuncs.a

$(TARGET): $(OBJECTS) 
	$(ARX) rcs -o $@ $^
$(OBJDIR)/%.o: %.cpp
	@if [ ! -d $(OBJDIR) ]; \
		then echo "mkdir -p $(OBJDIR)"; mkdir -p $(OBJDIR); \
	fi
	@[ -d $(OBJDIR) ]
	$(CXX) -c $(CFLAGS) $(INDS) -o $@ $<

clean:
	rm -r $(OBJDIR)
	rm $(TARGET)
depend:
	makedepend -- -Y -- $(SOURCES)


# DO NOT DELETE

