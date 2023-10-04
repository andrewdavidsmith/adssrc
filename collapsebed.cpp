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
first_end_point_score(const GenomicRegion &r) {
  return end_point_score(r.get_chrom(), r.get_start(), true, r.get_score());
}

static end_point_score
second_end_point_score(const GenomicRegion &r) {
  return end_point_score(r.get_chrom(), r.get_end(), false, r.get_score());
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
collapse_bed(const size_t cutoff, vector<end_point_score> &end_points,
             vector<GenomicRegion> &collapsed) {

  collapse_end_points(end_points);

  size_t count = 0;
  double total = 0.0;
  GenomicRegion region;
  for (size_t i = 0; i < end_points.size() - 1; ++i) {
    if (end_points[i].is_first) {
      count += end_points[i].count;
      total += end_points[i].total;
    }
    else {
      count -= end_points[i].count;
      total -= end_points[i].total;
    }
    if (count >= cutoff) {
      region.set_chrom(end_points[i].chr);
      region.set_start(end_points[i].start);
      region.set_end(end_points[i + 1].start);
      region.set_score(total/count);
      collapsed.push_back(region);
    }
  }
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


struct end_point {
  end_point(const string c, const size_t s, const bool isf) :
    chr(c), start(s), count(1), is_first(isf) {}
  bool operator<(const end_point &other) const {
    return (chr < other.chr ||
            (chr == other.chr &&
             (start < other.start ||
              (start == other.start &&
               is_first < other.is_first))));
  }
  string chr;
  size_t start;
  size_t count;
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

static bool
equal_end_points(const end_point &a, const end_point &b) {
  return a.chr == b.chr && a.start == b.start && a.is_first == b.is_first;
}


static void
collapse_end_points(vector<end_point> &ep) {
  size_t j = 0;
  for (size_t i = 1; i < ep.size(); ++i) {
    if (equal_end_points(ep[i], ep[j]))
      ep[j].count++;
    else ep[++j] = ep[i];
  }
  ep.erase(ep.begin() + j + 1, ep.end());
}


static void
collapse_bed(const size_t cutoff, vector<end_point> &end_points,
             vector<GenomicRegion> &collapsed) {

  collapse_end_points(end_points);

  size_t count = 0;
  GenomicRegion region;
  for (size_t i = 0; i < end_points.size() - 1; ++i) {
    if (end_points[i].is_first)
      count += end_points[i].count;
    else
      count -= end_points[i].count;
    if (count >= cutoff) {
      region.set_chrom(end_points[i].chr);
      region.set_start(end_points[i].start);
      region.set_end(end_points[i + 1].start);
      region.set_score(count);
      collapsed.push_back(region);
    }
  }
}



static void
collapse_bed_merged(const size_t cutoff, const vector<end_point> &end_points,
                    vector<GenomicRegion> &collapsed) {
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
        collapsed.push_back(region);
      }
      --count;
    }
  }
}


int
main(int argc, const char **argv) {

  try {

    /* FILES */
    string outfile;
    bool VERBOSE = false;
    bool FRAGMENTS = false;
    size_t cutoff = 1;

    bool AVERAGE_SCORES = false;

    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(strip_path(argv[0]), "", "<interval-files>");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("frag", 'f', "keep fragments separate",
                      false , FRAGMENTS);
    opt_parse.add_opt("cutoff", 'c', "report intervals covered this many times",
                      false , cutoff);
    opt_parse.add_opt("average", 'a', "average scores (implies -f)",
                      false, AVERAGE_SCORES);
    // opt_parse.add_opt("verbose", 'v', "print more run info",
    //                false , VERBOSE);
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

    vector<GenomicRegion> collapsed;
    if (AVERAGE_SCORES) {

      vector<end_point_score> end_points;
      for (size_t i = 0; i < interval_files.size(); ++i) {
        if (VERBOSE)
          cerr << interval_files[i] << endl;
        vector<GenomicRegion> intervals;
        ReadBEDFile(interval_files[i], intervals);
        for (size_t j = 0; j < intervals.size(); ++j) {
          end_points.push_back(first_end_point_score(intervals[j]));
          end_points.push_back(second_end_point_score(intervals[j]));
        }
      }
      sort(end_points.begin(), end_points.end());

      collapse_bed(cutoff, end_points, collapsed);
    }
    else {

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

      if (FRAGMENTS)
        collapse_bed(cutoff, end_points, collapsed);
      else
        collapse_bed_merged(cutoff, end_points, collapsed);
    }

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    copy(collapsed.begin(), collapsed.end(),
         std::ostream_iterator<GenomicRegion>(out, "\n"));

  }
  catch (std::exception &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
