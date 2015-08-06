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

SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))
CXX = g++
CXXFLAGS = -Wall
OPTFLAGS = -O2
DEBUGFLAGS = -g
COMMON_DIR = $(SMITHLAB_CPP)

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

all: $(OBJECTS)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -I$(COMMON_DIR) -c -o $@ $<
clean:
	@-rm -f *.o *~

.PHONY: clean
