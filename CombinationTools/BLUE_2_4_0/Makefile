ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLDFLAGS   = $(shell $(ROOTSYS)/bin/root-config --ldflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)

# Linux
CXX           = g++
CXXFLAGS      = -O4 -fPIC -I. 
SOFLAGS       = -shared $(ROOTLDFLAGS)

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS)

.SUFFIXES: .cxx
.PHONY: clean

code:= Blue.cxx BlueDict.cxx
objs:= $(foreach file, $(code), $(basename $(file)).o)

libBlue.so: $(objs)
	$(CXX) -o $@ $(SOFLAGS) $^ $(LIBS)

clean:
	rm -f $(objs) BlueDict.cxx BlueDict_rdict.pcm libBlue.so

BlueDict.cxx: Blue.h BlueLinkDef.h
	rootcling -f $@ -c $^

.cxx.o:
	$(CXX) $(CXXFLAGS) -o $@ -c $<
