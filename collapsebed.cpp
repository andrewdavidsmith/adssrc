/*    collapsebed: report contiguous genomic intervals covered at
 *    least the specified number of times by given intervals.
 *
 *    Copyright (C) 2014 Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;

struct end_point {
  end_point(const string c, const size_t s, const bool isf) : 
    chr(c), start(s), is_first(isf) {}
  bool operator<(const end_point &other) const {
    return (chr < other.chr || 
	    (chr == other.chr && 
	     (start < other.start ||
	      (start == other.start && 
	       is_first < other.is_first))));
  }
  string chr;
  size_t start;
  bool is_first;
};

static end_point
first_end_point(const GenomicRegion &r) {
  return end_point(r.get_chrom(), r.get_start(), true);
}

static end_point
second_end_point(const GenomicRegion &r) {
  return end_point(r.get_chrom(), r.get_end(), false);
}

int
main(int argc, const char **argv) {

  try {

    /* FILES */
    string outfile;
    bool VERBOSE = false;
    size_t cutoff = 1;

    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(strip_path(argv[0]), "", "<interval-files>");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("cutoff", 'c', "report intervals covered this many times",
		      false , cutoff);
    // opt_parse.add_opt("verbose", 'v', "print more run info", 
    // 		      false , VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() < 1) { 
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const vector<string> interval_files(leftover_args);
    /**********************************************************************/
    
    vector<end_point> end_points;
    for (size_t i = 0; i < interval_files.size(); ++i) {
      if (VERBOSE)
	cerr << interval_files[i] << endl;
      vector<GenomicRegion> intervals;
      ReadBEDFile(interval_files[i], intervals);
      for (size_t j = 0; j < intervals.size(); ++j) {
	end_points.push_back(first_end_point(intervals[j]));
	end_points.push_back(second_end_point(intervals[j]));
      }
    }
    sort(end_points.begin(), end_points.end());
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    
    size_t count = 0;
    GenomicRegion region;
    for (size_t i = 0; i < end_points.size(); ++i) {
      if (end_points[i].is_first) {
	++count;
	if (count == cutoff) {
	  region.set_chrom(end_points[i].chr);
	  region.set_start(end_points[i].start);
	}
      }
      else {
	if (count == cutoff) {
	  region.set_end(end_points[i].start);
	  out << region << endl;
	}
	--count;
      }
    }
    
  }
  catch (SMITHLABException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
