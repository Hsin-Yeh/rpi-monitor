SRCDIR= src
BUILDDIR= obj_file
TARGET= makePlots

SRCTXT= cc
SOURCES= $(shell find $(SRCDIR) -type f -name *.$(SRCTXT))
OBJECTS= $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCTXT)=.o))
INC=-I include
INCDIR= include

CXXFLAGS=-g -m64 -O2 -Wall -Wno-unused-variable -Wno-unused-command-line-argument -std=c++0x $(INC)
ROOTFLAGS=$(shell root-config --libs --cflags --glibs)
LDFLAGS+=-lutil -lboost_iostreams -lboost_system -lboost_filesystem

$(TARGET): $(OBJECTS)
	g++ $^ -o $@ $(CXXFLAGS) $(ROOTFLAGS) $(LDFLAGS)

$(BUILDDIR)/main.o: $(SRCDIR)/main.cc
	g++ -c $(CXXFLAGS) $(ROOTFLAGS) $(LDFLAGS) $< -o $@

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCTXT) $(INCDIR)/%.h
	g++ -c $(CXXFLAGS) $(ROOTFLAGS) $(LDFLAGS) $< -o $@

clean:
	rm -f $(TARGET) $(BUILDDIR)/*.o include/*~ $(SRCDIR)/*~ ./*~

