CXX       := g++
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

CXXFLAGS := -O3 -fPIC -std=c++14 $(ROOTCFLAGS)

SHLIB      := libTMinFit.so
EXEC       := genSystErr

SHSRC      := TMinFit.cpp
SHOBJ      := $(SHSRC:.cpp=.o)

EXESRC     := genSystErr.cpp SimpleCardReader.cpp
EXEOBJ     := $(EXESRC:.cpp=.o)

all: $(SHLIB) $(EXEC)

#
# Shared library
#
$(SHLIB): $(SHOBJ)
	$(CXX) -shared -o $@ $^ $(ROOTLIBS) -lMinuit

#
# Executable
#
$(EXEC): $(EXEOBJ) $(SHLIB)
	$(CXX) -o $@ $(EXEOBJ) \
	    -L. -lTMinFit \
	    $(ROOTLIBS) -lMinuit \
	    -Wl,-rpath,'$$ORIGIN'

#
# Object files
#
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o *.so $(EXEC)
