/*    collapse_bedgraph:
 *
 *    Copyright (C) 2018 Andrew D. Smith
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

#include <exception>
#include <iostream>
#include <fstream>

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::runtime_error;

struct BedGraphInterval {
  string chr;
  size_t start;
  size_t end;
  double score;
  BedGraphInterval(const string &c, const size_t s,
                   const size_t e, const double sc) :
    chr(c), start(s), end(e), score(sc) {}
  BedGraphInterval() {}
};

std::istream &
operator>>(std::istream &in, BedGraphInterval &bgi) {
  return in >> bgi.chr >> bgi.start >> bgi.end >> bgi.score;
}

std::ostream &
operator<<(std::ostream &out, const BedGraphInterval &bgi) {
  return out << bgi.chr << '\t'
             << bgi.start << '\t' << bgi.end << '\t' << bgi.score;
}

struct end_point_score {
  end_point_score(const string c, const size_t s, const bool isf,
                  const double sc) :
    chr(c), start(s), count(1), total(sc), is_first(isf) {}
  bool operator<(const end_point_score &other) const {
    return (chr < other.chr ||
            (chr == other.chr &&
             (start < other.start ||
              (start == other.start &&
               is_first < other.is_first))));
  }
  string chr;
  size_t start;
  size_t count;
  double total;
  bool is_first;
};

static end_point_score
first_end_point_score(const BedGraphInterval &r) {
  return end_point_score(r.chr, r.start, true, r.score);
}

static end_point_score
second_end_point_score(const BedGraphInterval &r) {
  return end_point_score(r.chr, r.end, false, r.score);
}

static bool
equal_end_points_score(const end_point_score &a, const end_point_score &b) {
  return a.chr == b.chr && a.start == b.start && a.is_first == b.is_first;
}


static void
collapse_end_points(vector<end_point_score> &ep) {
  size_t j = 0;
  for (size_t i = 1; i < ep.size(); ++i) {
    if (equal_end_points_score(ep[i], ep[j])) {
      ep[j].count++;
      ep[j].total += ep[i].total;
    }
    else ep[++j] = ep[i];
  }
  ep.erase(ep.begin() + j + 1, ep.end());
}


static void
collapse_bed(vector<end_point_score> &end_points,
             vector<BedGraphInterval> &collapsed) {

  collapse_end_points(end_points);

  size_t count = 0;
  double total = 0.0;
  for (size_t i = 0; i < end_points.size() - 1; ++i) {
    if (end_points[i].is_first) {
      count += end_points[i].count;
      total += end_points[i].total;
    }
    else {
      count -= end_points[i].count;
      total -= end_points[i].total;
    }
    if (count > 0 && end_points[i].start < end_points[i + 1].start) {
      collapsed.push_back(BedGraphInterval(end_points[i].chr,
                                           end_points[i].start,
                                           end_points[i + 1].start,
                                           total/count));
    }
  }
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

int
main(int argc, const char **argv) {

  try {

    /* FILES */
    string outfile;
    bool VERBOSE = false;

    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(strip_path(argv[0]), "", "<interval-files>");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false , VERBOSE);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string bedgraph_file(leftover_args.front());
    /**********************************************************************/

    if (VERBOSE)
      cerr << "reading bedgraph file: " << bedgraph_file << endl;

    std::ifstream in(bedgraph_file);
    if (!in)
      throw runtime_error("cannot open file: " + bedgraph_file);

    vector<BedGraphInterval> intervals;
    BedGraphInterval b;
    while (in >> b)
      intervals.push_back(b);

    if (VERBOSE)
      cerr << "separating endpoints" << endl;
    vector<end_point_score> end_points;
    for (size_t j = 0; j < intervals.size(); ++j) {
      end_points.push_back(first_end_point_score(intervals[j]));
      end_points.push_back(second_end_point_score(intervals[j]));
    }

    if (VERBOSE)
      cerr << "sorting endpoints" << endl;
    sort(end_points.begin(), end_points.end());

    if (VERBOSE)
      cerr << "collapsing overlapping intervals" << endl;
    vector<BedGraphInterval> collapsed;
    collapse_bed(end_points, collapsed);

    if (VERBOSE)
      cerr << "writing output" << endl;
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    copy(collapsed.begin(), collapsed.end(),
         std::ostream_iterator<BedGraphInterval>(out, "\n"));

  }
  catch (const runtime_error &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (const std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
