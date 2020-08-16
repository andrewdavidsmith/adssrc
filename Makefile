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

PROGS = collapsebed countoverlaps baseliftover bedgraph_liftover expand_chains \
	adjust_intervals collapse_bedgraph dataframe-euclidean-dist varyingrows \
	readliftover collapse3col transpose_data_frame

SOURCES = $(wildcard *.cpp)
INCLUDEDIRS = $(SMITHLAB_CPP)
LIBS = -lgsl -lgslcblas

ifdef METHPIPE_ROOT
PROGS += tsscpgplot smoothmeth binmeth autocorr
INCLUDEDIRS += $(METHPIPE_ROOT)/src/common
endif

INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

CXX = g++
CXXFLAGS = -Wall -Wextra -fmessage-length=72 -std=c++11
OPTFLAGS = -O2
DEBUGFLAGS = -g

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

all: $(PROGS)

chain_file_utils.o: chain_file_utils.cpp chain_file_utils.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o)

tsscpgplot smoothmeth: \
	$(addprefix $(METHPIPE_ROOT)/src/common/, MethpipeFiles.o)

smoothmeth: $(addprefix $(METHPIPE_ROOT)/src/common/, TwoStateHMM.o)

reorder:    $(addprefix $(SMITHLAB_CPP)/, MappedRead.o)

autocorr: $(addprefix $(METHPIPE_ROOT)/src/common/, MethpipeSite.o)

expand_chains: chain_file_utils.o

# No rule to make object bsutils.o in Methpipe, so this won't compile.
#majormethstate: $(addprefix $(METHPIPE_ROOT)/src/common/, bsutils.o)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
