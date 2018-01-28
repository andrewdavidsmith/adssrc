/*    syncbins: take a set of non-overlapping intervals with
 *    associated scores (bedgraph format) and position the
 *    corresponding information relative to a set of bins.
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
#include "GenomicRegion.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::unordered_map;
using std::pair;

struct bedgraph {
  string chr;
  size_t start;
  size_t end;
  double score;
};

static std::istream&
operator>>(std::istream &in, bedgraph &bg) {
  string line;
  if (getline(in, line)) {
    std::istringstream iss(line);
    if (!(iss >> bg.chr >> bg.start >> bg.end >> bg.score))
      throw std::runtime_error("bad line: " + line);
  }
  return in;
}

static bool
ends_before(bedgraph &bg, const string &chr, const size_t start) {
  return bg.chr < chr || (bg.chr == chr && bg.end <= start);
}

static bool
doesnt_extend_past(bedgraph &bg, const string &chr, const size_t end) {
  return bg.chr < chr || (bg.chr == chr && (bg.end <= end));
}

static bool
overlaps(bedgraph &bg, const string &chr, const size_t start,
         const size_t end) {
  return bg.chr == chr &&
    (std::max(bg.start, start) < std::min(bg.end, end));
}

static void
update_bin(const size_t start, const size_t end,
           const bedgraph &bg, double &total, size_t &count) {
  const size_t overlap_size =
    std::min(bg.end, end) - std::max(bg.start, start);
  total += bg.score*overlap_size;
  count += overlap_size;
}

int
main(int argc, const char **argv) {

  try {

    /* FILES */
    string outfile;
    bool VERBOSE = false;
    bool require_data = false;
    size_t bin_size = 100;

    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(strip_path(argv[0]), "", "<chroms> <bedgraph>");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("bin-size", 'b', "bin size (default: 100)",
                      false , bin_size);
    opt_parse.add_opt("with-data", 'd', "only output bins with data",
                      false , require_data);
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
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string chroms_file(leftover_args.front());
    const string bg_file(leftover_args.back());
    /**********************************************************************/

    std::ifstream chroms_in(chroms_file.c_str());
    if (!chroms_in)
      throw std::runtime_error("bad file: " + chroms_file);

    string line;
    vector<pair<string, size_t> > chroms;
    while (getline(chroms_in, line)) {
      std::istringstream iss(line);
      string chr;
      size_t sz = 0;
      if (!(iss >> chr >> sz))
        throw std::runtime_error("bad line: " + line);
      chroms.push_back(make_pair(chr, sz));
    }
    sort(chroms.begin(), chroms.end());
    if (VERBOSE)
      for (auto &i : chroms)
        cerr << i.first << '\t' << i.second << endl;

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    bedgraph curr;
    std::ifstream in(bg_file.c_str());
    in >> curr;

    for (size_t i = 0; i < chroms.size(); ++i) {

      const string chr(chroms[i].first);
      const size_t chr_size = chroms[i].second;
      for (size_t start = 0; start < chr_size; start += bin_size) {
        const size_t end = std::min(chr_size, start + bin_size);

        double total = 0.0; // sum of (bases * scores) for covered bases in bin
        size_t count = 0; // sum of (bases) == number of covered bases

        // move past all earlier intervals; might not be needed as a
        // loop, since bins cover entire genome
        if (in && ends_before(curr, chr, start))
          in >> curr;

        // process overlapping intervals not relevant for the next bin
        while (in && doesnt_extend_past(curr, chr, end)) {
          update_bin(start, end, curr, total, count);
          in >> curr;
        }
        // process any overlapping intervals relevant for the next bin
        if (in && overlaps(curr, chr, start, end))
          update_bin(start, end, curr, total, count);

        if (!require_data || count > 0)
          out << chr << '\t'
              << start << '\t' << end << '\t'
              << (count > 0 ? total/count : 0) << '\n';
      }
    }
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  catch (std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
