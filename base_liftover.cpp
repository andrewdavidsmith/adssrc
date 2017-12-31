/*    base_liftover: lift over sites using index file
 *
 *    Copyright (C) 2017 University of Southern California and
 *                       Andrew D. Smith
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

struct site {
  string chrom;
  size_t position;
  double data;
};

std::ostream &
operator<<(std::ostream &os, const site &s) {
  return os << s.chrom << '\t' << s.position << '\t' << s.data;
}

std::istream &
operator>>(std::istream &in, site &s) {
  size_t dummy;
  in >> s.chrom >> s.position >> dummy >> s.data;
  return in;
}

struct block {
  SimpleGenomicRegion ref;
  SimpleGenomicRegion query;
  char query_strand;
  bool operator<(const block &other) const {
    return ref < other.ref;
  }
};

std::ostream &
operator<<(std::ostream &os, const block &b) {
  os << b.ref << '\t' << b.query << '\t' << b.query_strand;
  return os;
}

std::istream &
operator>>(std::istream &is, block &b) {

  string tmp_chrom;
  is >> tmp_chrom;
  b.ref.set_chrom(tmp_chrom);

  size_t tmp_pos;
  is >> tmp_pos;
  b.ref.set_start(tmp_pos);
  is >> tmp_pos;
  b.ref.set_end(tmp_pos);

  is >> tmp_chrom;
  b.query.set_chrom(tmp_chrom);

  is >> tmp_pos;
  b.query.set_start(tmp_pos);
  is >> tmp_pos;
  b.query.set_end(tmp_pos);

  is >> b.query_strand;

  return is;
}

static bool
precedes(const block &b, const site &s) {
  return (b.ref.get_chrom() < s.chrom ||
          (b.ref.get_chrom() == s.chrom &&
           b.ref.get_end() <= s.position));
}

static bool
contains(const block &b, const site &s) {
  return b.ref.get_chrom() == s.chrom &&
    b.ref.get_start() <= s.position &&
    s.position < b.ref.get_end();
}

static void
lift_site(const block &b, const site &orig, site &lifted) {
  const size_t diff = orig.position - b.ref.get_start();
  lifted.chrom = b.query.get_chrom();
  lifted.data = orig.data;
  if (b.query_strand == '+')
    lifted.position = b.query.get_start() + diff;
  else lifted.position = b.query.get_end() - diff - 1;
}

int
main(int argc, const char **argv) {
  try{
    string outfile;
    string fails_file;

    bool VERBOSE = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "liftover individual sites sorted",
                           "<chain> <sites>");
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
    const string sites_file(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "[reading input chains]" << endl;
    vector<block> the_blocks;
    block tmp_block;
    std::ifstream in_blocks(chain_file.c_str());
    while (in_blocks >> tmp_block) {
      the_blocks.push_back(tmp_block);
    }
    const size_t n_blocks = the_blocks.size();

    if (VERBOSE)
      cerr << "[sorting input blocks]" << endl;
    sort(the_blocks.begin(), the_blocks.end());

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    std::ofstream ff;
    if (!fails_file.empty()) ff.open(fails_file.c_str());
    std::ostream fails(ff.rdbuf());

    if (VERBOSE)
      cerr << "[processing sites]" << endl;

    std::ifstream in_sites(sites_file.c_str());

    size_t curr_block = 0;
    site tmp_site;

    while (in_sites >> tmp_site) {

      // find relevant block
      while (curr_block < n_blocks &&
             precedes(the_blocks[curr_block], tmp_site))
        ++curr_block;

      if (curr_block < n_blocks &&
          contains(the_blocks[curr_block], tmp_site)) {
        // the site is contained in the block
        site lifted;
        lift_site(the_blocks[curr_block], tmp_site, lifted);
        out << lifted << endl; // "\t" << tmp_site << endl;
      }
      else {
        if (!fails_file.empty())
          fails << tmp_site << endl;
      }
    }
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
