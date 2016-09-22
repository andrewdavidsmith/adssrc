/*  select-sites: get all sites from methcounts output in (or outside)
 *  a set of BED intervals
 *
 *  Copyright (C) 2016 University of Southern California
 *                     Andrew D. Smith
 *
 *  Authors: Andrew D. Smith
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "MethpipeSite.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::pair;
using std::ios_base;


static void
move_to_start_of_line(std::istream &in) {
  char next;
  while (in.good() && in.get(next) && next != '\n') {
    in.seekg(-2, ios_base::cur);
  }
  if (in.bad())
    // hope this only happens when hitting the start of the file
    in.clear();
}


static void
find_start_line(const string &chr, const size_t idx, std::istream &in) {

  in.seekg(0, ios_base::beg);
  const size_t begin_pos = in.tellg();
  in.seekg(0, ios_base::end);
  const size_t end_pos = in.tellg();

  if (end_pos - begin_pos < 2)
    throw SMITHLABException("empty meth file");

  size_t step_size = (end_pos - begin_pos)/2;

  in.seekg(0, ios_base::beg);
  string low_chr;
  size_t low_idx = 0;
  in >> low_chr >> low_idx;

  // MAGIC: need the -2 here to get past the EOF and possibly a '\n'
  in.seekg(-2, ios_base::end);
  move_to_start_of_line(in);
  string high_chr;
  size_t high_idx;
  in >> high_chr >> high_idx;

  size_t pos = step_size;
  in.seekg(pos, ios_base::beg);
  move_to_start_of_line(in);

  while (step_size > 0) {
    string mid_chr;
    size_t mid_idx = 0;
    in >> mid_chr >> mid_idx;
    step_size /= 2;
    if (chr < mid_chr || (chr == mid_chr && idx <= mid_idx)) {
      std::swap(mid_chr, high_chr);
      std::swap(mid_idx, high_idx);
      pos -= step_size;
    }
    else {
      std::swap(mid_chr, low_chr);
      std::swap(mid_idx, low_idx);
      pos += step_size;
    }
    in.seekg(pos, ios_base::beg);
    move_to_start_of_line(in);
  }
}


static void
process_interval(std::istream &in,
                 const GenomicRegion &region, std::ostream &out) {

  const string chrom(region.get_chrom());
  const size_t start_pos = region.get_start();
  const size_t end_pos = region.get_end();
  find_start_line(chrom, start_pos, in);

  MSite curr_site;
  while ((in >> curr_site) &&
         (curr_site.chrom == chrom &&
          curr_site.pos < end_pos)) {
    if (start_pos <= curr_site.pos)
      out << curr_site << endl;
  }
  in.clear();
}


int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    string outfile;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "Compute average CpG "
                           "methylation in each of a set of genomic intervals",
                           "<intervals-bed> <sites-meth-file>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
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
    const string regions_file = leftover_args.front();
    const string sites_file = leftover_args.back();
    /****************** END COMMAND LINE OPTIONS *****************/

    vector<GenomicRegion> regions;
    ReadBEDFile(regions_file, regions);
    if (!check_sorted(regions))
      throw SMITHLABException("regions not sorted: " + regions_file);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    std::ifstream in(sites_file.c_str());
    if (!in)
      throw SMITHLABException("bad sites file: " + sites_file);

    for (size_t i = 0; i < regions.size(); ++i)
      process_interval(in, regions[i], out);

  }
  catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
