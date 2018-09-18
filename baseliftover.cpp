/* baseliftover: lift over individual bases. Currently this has
 * problems and I can't remember what the base_liftover program does
 * differently.
 *
 * Copyright (C) 2018 Andrew D. Smith
 *
 * Authors: Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <stdexcept>
#include <algorithm>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "GenomicRegion.hpp"

using std::string;
using std::ios_base;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::runtime_error;


struct site {
  site() {}
  site(const string &c, const size_t p) : chrom(c), pos(p) {}
  string chrom;
  size_t pos;
};


struct block {
  string r_chrom;
  size_t r_start;
  size_t r_end;
  string q_chrom;
  size_t q_start;
  size_t q_end;
  char q_strand;
  bool operator<(const block &o) const {
    // not checking the r_end because these should not overlap
    return r_chrom < o.r_chrom || (r_chrom == o.r_chrom && r_start < o.r_start);
  }
};


std::ostream &
operator<<(std::ostream &os, const block &b) {
  return os << b.r_chrom << '\t'
            << b.r_start << '\t'
            << b.r_end << '\t'
            << b.q_chrom << '\t'
            << b.q_start << '\t'
            << b.q_end << '\t'
            << b.q_strand;
}


std::istream &
operator>>(std::istream &is, block &b) {
  return (is
          >> b.r_chrom
          >> b.r_start
          >> b.r_end
          >> b.q_chrom
          >> b.q_start
          >> b.q_end
          >> b.q_strand);
}


static bool
precedes(const block &b, const site &s) {
  return b.r_chrom < s.chrom || (b.r_chrom == s.chrom && b.r_end <= s.pos);
}


static bool
contains(const block &b, const site &s) {
  return b.r_start <= s.pos && s.pos < b.r_end && b.r_chrom == s.chrom;
}


static void
lift_site(const block &b, const site &work, GenomicRegion &lifted) {

  const size_t diff = work.pos - b.r_start;

  lifted.set_chrom(b.q_chrom);
  if (b.q_strand == '+')
    lifted.set_start(b.q_start + diff);
  else {
    lifted.set_start(b.q_end - diff - 1);
    lifted.set_strand(lifted.pos_strand() ? '-' : '+');
  }
  lifted.set_end(lifted.get_start());
}


static bool
lift_pos_strand(std::ostream &out,
                const vector<block> &the_blocks, size_t curr_block,
                site work, const GenomicRegion &region) {

  const size_t n_blocks = the_blocks.size();

  // make a temporary site based on 5' end of read to locate
  const size_t w = region.get_width();

  bool lifted_position = false;

  for (size_t i = 0; i < w; ++i) {
    // find relevant block
    while (curr_block < n_blocks && precedes(the_blocks[curr_block], work))
      ++curr_block;

    // if we have a relevant block, do the lifting
    if (curr_block < n_blocks && contains(the_blocks[curr_block], work)) {
      GenomicRegion lifted;
      lift_site(the_blocks[curr_block], work, lifted);
      out << lifted << endl;
      lifted_position = true;
    }
    ++work.pos;
  }
  return lifted_position;
}


static bool
lift_neg_strand(std::ostream &out,
                const vector<block> &the_blocks, size_t curr_block,
                site work, const GenomicRegion &region) {

  const size_t n_blocks = the_blocks.size();

  // make a temporary site based on 5' end of read to locate
  const size_t w = region.get_width();

  bool lifted_position = false;

  for (size_t i = 0; i < w; ++i) {
    // find relevant block
    while (curr_block < n_blocks && precedes(the_blocks[curr_block], work))
      ++curr_block;

    // if we have a relevant block, do the lifting
    if (curr_block < n_blocks && contains(the_blocks[curr_block], work)) {
      GenomicRegion lifted;
      lift_site(the_blocks[curr_block], work, lifted);
      out << lifted << endl;
      lifted_position = true;
    }
    ++work.pos;
  }
  return lifted_position;
}


int
main(int argc, const char **argv) {
  try{
    string outfile;
    string fails_file;

    bool VERBOSE = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "liftover each base possible in a set of intervals",
                           "<chain> <intervals>");
    opt_parse.add_opt("outfile", 'o', "output file", false, outfile);
    opt_parse.add_opt("fails", 'f', "failed liftovers",
                      false, fails_file);
    opt_parse.add_opt("verbose", 'v', "(optional) Print more information",
                      false, VERBOSE);
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
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    const string chain_file(leftover_args.front());
    const string intervals_file(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    /* READ IN AND VERIFY THE CHAIN BLOCKS FOR LIFTING OVER
     */
    if (VERBOSE)
      cerr << "[reading input blocks]" << endl;
    std::ifstream in_blocks(chain_file.c_str());
    if (!in_blocks)
      throw runtime_error("bad input file: " + chain_file);

    vector<block> the_blocks;
    block tmp_block;
    while (in_blocks >> tmp_block)
      the_blocks.push_back(tmp_block);

    if (!check_sorted(the_blocks))
      throw runtime_error("blocks not sorted in: " + chain_file);

    const size_t n_blocks = the_blocks.size();
    if (VERBOSE)
      cerr << "[total blocks = " << n_blocks << "]" << endl;

    /* PREPARE THE OUTPUT FILES (FOR THE LIFTED AND THE FAILS)
     */
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    std::ofstream ff;
    if (!fails_file.empty()) ff.open(fails_file.c_str());
    std::ostream fails(ff.rdbuf());

    /* NOW PROCESS THE SITES
     */
    std::ifstream in_intervals(intervals_file.c_str());
    if (!in_intervals)
      throw runtime_error("could not open file: " + intervals_file);

    if (VERBOSE)
      cerr << "[processing sites]" << endl;

    size_t curr_block = 0;
    GenomicRegion tmp_gr;
    size_t lines_processed = 0;

    while (in_intervals >> tmp_gr) {
      ++lines_processed;

      site work(tmp_gr.get_chrom(), tmp_gr.get_start());

      // find relevant block
      while (curr_block < n_blocks && precedes(the_blocks[curr_block], work))
        ++curr_block;

      size_t local_block = curr_block;

      GenomicRegion lifted_gr(tmp_gr);
      const bool lift_success = tmp_gr.pos_strand() ?
        lift_pos_strand(out, the_blocks, local_block, work, lifted_gr) :
        lift_neg_strand(out, the_blocks, local_block, work, lifted_gr);

      if (!lift_success && !fails_file.empty())
        fails << tmp_gr << '\n';
    }
    if (VERBOSE)
      cerr << "[lines processed = " << lines_processed << "]" << endl;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
