#    Copyright (C) 2015 Andrew D. Smith
#
#    Authors: Jenny Qu,Andrew D. Smith
#
#    This code is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

PROGS = treetool

ifndef SMITHLAB_CPP
$(error SMITHLAB_CPP variable undefined)
endif

## Macro below (careful of whitespace when calling)
not-containing = $(foreach v,$2,$(if $(findstring $1,$v),,$v))

SOURCES = $(wildcard *.cpp)
OBJSOURCES = $(call not-containing,$(PROGS),$(SOURCES))
OBJECTS = $(patsubst %.cpp,%.o,$(OBJSOURCES))

CXX = g++
CXXFLAGS = -Wall
OPTFLAGS = -O2
DEBUGFLAGS = -g

INCLUDEDIRS :=  $(SMITHLAB_CPP) $(METHPIPE_ROOT)
INCLUDEARGS := $(addprefix -I,$(INCLUDEDIRS))

# Flags for Mavericks
ifeq "$(shell uname)" "Darwin"
CXXFLAGS += -arch x86_64
INCLUDEARGS += -I/usr/local/include
LIBS += -L/usr/local/lib
ifeq "$(shell if [ `sysctl -n kern.osrelease | cut -d . -f 1` -ge 13 ];\
then echo 'true'; fi)" "true"
CXXFLAGS += -stdlib=libstdc++
endif
endif

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif
ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

all: $(PROGS)

$(PROGS): PhyloTree.o \
	$(addprefix $(SMITHLAB_CPP)/, \
	OptionParser.o smithlab_utils.o smithlab_os.o)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

clean:
	@-rm -f *.o *~ $(PROGS)
.PHONY: clean
