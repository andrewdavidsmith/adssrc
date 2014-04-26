/*    countoverlaps: for each target interval, count the number of
 *    others that overlap it
 *
 *    Copyright (C) 2009-2014 Andrew D. Smith
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

#include <fstream>
#include <vector>
#include <tr1/unordered_set>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::pair;
using std::tr1::unordered_set;

struct end_pt {
  end_pt(const string c, const size_t p) : chr(c), pos(p) {}
  bool operator<(const end_pt &other) const {
    return chr < other.chr || (chr == other.chr && pos < other.pos);
  }
  bool operator==(const end_pt &other) const {
    return chr == other.chr && pos == other.pos;
  }
  string chr;
  size_t pos;
};

static end_pt
get_left(const GenomicRegion &r) {return end_pt(r.get_chrom(), r.get_start());}

static end_pt
get_right(const GenomicRegion &r) {return end_pt(r.get_chrom(), r.get_end());}


static void
make_segment_finder(const vector<GenomicRegion> &intervals, 
		    vector<pair<end_pt, vector<size_t> > > &segments) {
  
  // make a set of sorted end-points for all intervals
  for (size_t i = 0; i < intervals.size(); ++i) {
    segments.push_back(make_pair(get_left(intervals[i]), vector<size_t>(1, i)));
    segments.push_back(make_pair(get_right(intervals[i]), vector<size_t>(1, i)));
  }
  sort(segments.begin(), segments.end());
  
  // merge common end-points
  size_t j = 0;
  for (size_t i = 1; i < segments.size(); ++i)
    if (segments[j].first == segments[i].first)
      segments[j].second.push_back(segments[i].second.front());
    else segments[++j] = segments[i];
  segments.erase(segments.begin() + j + 1, segments.end());
  
  // ensure every sub-interval contains all ids of overlapping
  // intervals in the segment finder
  unordered_set<size_t> ids;
  for (size_t i = 0; i < segments.size(); ++i) {
    for (size_t j = 0; j < segments[i].second.size(); ++j)
      if (ids.find(segments[i].second[j]) == ids.end())
	ids.insert(segments[i].second[j]);
      else ids.erase(segments[i].second[j]);
    segments[i].second.clear();
    copy(ids.begin(), ids.end(), back_inserter(segments[i].second));
  }
}


int
main(int argc, const char **argv) {

  /* FILES */
  string outfile;
  bool VERBOSE = false;
  bool PRINT_IDS = false;

  /****************** GET COMMAND LINE ARGUMENTS ***************************/
  OptionParser opt_parse(strip_path(argv[0]), "count regions overlapping "
			 "each of a set of target regions",
			 "<targets> <to-count>");
  opt_parse.add_opt("output", 'o', "output file (default: stdout)", 
		    false , outfile);
  opt_parse.add_opt("ids", 'I', "report ids", false , PRINT_IDS);
  opt_parse.add_opt("verbose", 'v', "print more run info", false , VERBOSE);
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
  if (leftover_args.size() != 2) {
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }
  const string targets_file = leftover_args.front();
  const string to_count_file = leftover_args.back();
  /**********************************************************************/
  
  try {
    
    vector<GenomicRegion> to_count;
    ReadBEDFile(to_count_file, to_count);
    
    vector<pair<end_pt, vector<size_t> > > segment_finder;
    make_segment_finder(to_count, segment_finder);
    
    vector<end_pt> keys;
    vector<vector<size_t> > values(segment_finder.size());
    for (size_t i = 0; i < segment_finder.size(); ++i) {
      keys.push_back(segment_finder[i].first);
      values[i].swap(segment_finder[i].second);
    }
    vector<pair<end_pt, vector<size_t> > >().swap(segment_finder);
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    std::ifstream in(targets_file.c_str());
    
    GenomicRegion target;
    while (in >> target) {
      size_t left = lower_bound(keys.begin(), keys.end(), 
				get_left(target)) - keys.begin();
      const size_t right = lower_bound(keys.begin(), keys.end(),
				       get_right(target)) - keys.begin();
      
      // required because ids are stored at starts of segments
      if (left > 0 && !(keys[left] == get_left(target))) --left;
      
      unordered_set<size_t> ids;
      for (; left < right; ++left)
	copy(values[left].begin(), values[left].end(), 
	     std::inserter(ids, ids.end()));
      target.set_score(ids.size());
      
      out << target;
      if (PRINT_IDS)
	for (unordered_set<size_t>::const_iterator 
	       j(ids.begin()); j != ids.end(); ++j) 
	  out << '\t' << to_count[*j].get_name();
      out << endl;
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
