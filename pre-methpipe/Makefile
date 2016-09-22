#    Copyright (C) 2014 Andrew D. Smith
#
#    Authors: Andrew D. Smith
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

ifndef SMITHLAB_CPP
$(error Must define SMITHLAB_CPP variable)
endif

PROGS =

SOURCES = $(wildcard *.cpp)
INCLUDEDIRS = $(SMITHLAB_CPP)
LIBS = -lgsl -lgslcblas # -lefence

ifdef METHPIPE_ROOT
PROGS += methcounts2 bsrate2 roimethstat2 hairpinreads binmeth levels2 select-sites
INCLUDEDIRS += $(METHPIPE_ROOT)/src/common
endif

INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -fmessage-length=72
OPTFLAGS = -O2
DEBUGFLAGS = -g

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

all: $(PROGS)

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o)

roimethstat2 methcounts2 binmeth levels2: \
	$(addprefix $(METHPIPE_ROOT)/src/common/, MethpipeFiles.o)

select-sites: \
	$(addprefix $(METHPIPE_ROOT)/src/common/, MethpipeSite.o)

mapscores bsrate2 methcounts2: $(addprefix $(SMITHLAB_CPP)/, MappedRead.o)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
