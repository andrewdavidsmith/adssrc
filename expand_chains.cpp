/* Copyright (C) 2017 University of Southern California and
 *                    Andrew D. Smith
 *
 * Authors: Andrew D. Smith
 *
 * expand_chains: expand chains to intervals. The chains should be
 * in available from the UCSC Genome Browser downloads and in a
 * directory called liftOver for a particular species. The names
 * should look like mm10ToHg19.over.chain.gz. An example can be
 * found here:
 *
 * http://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToHg19.over.chain.gz
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <exception>

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


struct aln_dat {
  size_t size; /* the size of the ungapped alignment */
  size_t dt;   /* the difference between the end of this block and the
                  beginning of the next block (reference sequence) */
  size_t dq;   /* the difference between the end of this block and the
                  beginning of the next block (query sequence) */
};

std::ostream &
operator<<(std::ostream &os, const aln_dat &ad) {
  os << ad.size;
  if (ad.dt > 0 || ad.dq > 0)
    os << '\t' << ad.dt << '\t' << ad.dq;
  return os;
}

std::istream &
operator>>(std::istream &is, aln_dat &ad) {
  is >> ad.size;
  if (is >> ad.dt)
    is >> ad.dq;
  else {
    ad.dt = 0;
    ad.dq = 0;
    is.clear();
  }
  return is;
}

struct chain_head {
  size_t score; // chain score
  string tName; // chromosome (reference sequence)
  size_t tSize; // chromosome size (reference sequence)
  char tStrand; // strand (reference sequence)
  size_t tStart; // alignment start position (reference sequence)
  size_t tEnd; // alignment end position (reference sequence)
  string qName; // chromosome (query sequence)
  size_t qSize; // chromosome size (query sequence)
  char qStrand; // strand (query sequence)
  size_t qStart; // alignment start position (query sequence)
  size_t qEnd; // alignment end position (query sequence)
  size_t id; // chain ID
};

std::ostream &
operator<<(std::ostream &os, const chain_head &ch) {
  return os << "chain" << ' '
     << ch.score << ' '
     << ch.tName << ' '
     << ch.tSize << ' '
     << ch.tStrand << ' '
     << ch.tStart << ' '
     << ch.tEnd << ' '
     << ch.qName << ' '
     << ch.qSize << ' '
     << ch.qStrand << ' '
     << ch.qStart << ' '
     << ch.qEnd << ' '
     << ch.id;
}

std::istream &
operator>>(std::istream &is, chain_head &ch) {
  string dummy;
  is >> dummy; // this should always equal "chain";
  is >> ch.score;
  is >> ch.tName;
  is >> ch.tSize;
  is >> ch.tStrand;
  is >> ch.tStart;
  is >> ch.tEnd;
  is >> ch.qName;
  is >> ch.qSize;
  is >> ch.qStrand;
  is >> ch.qStart;
  is >> ch.qEnd;
  is >> ch.id;

  char should_end_line = is.get(); // warnings... I know...

  return is;
};

struct chain {
  chain_head head;
  vector<aln_dat> alns; // alignment data lines
};

std::ostream &
operator<<(std::ostream &os, const chain &c) {
  const string ref_chrom(c.head.tName);
  const string query_chrom(c.head.qName);
  size_t ref_position = c.head.tStart;
  size_t query_position = c.head.qStart;
  if (c.head.qStrand == '+') {
    for (size_t i = 0; i < c.alns.size(); ++i) {
      os << ref_chrom << '\t'
         << ref_position << '\t'
         << ref_position + c.alns[i].size << '\t';
      os << query_chrom << '\t'
         << query_position << '\t'
         << query_position + c.alns[i].size << '\t'
         << '+' << endl;
      ref_position += (c.alns[i].size + c.alns[i].dt);
      query_position += (c.alns[i].size + c.alns[i].dq);
    }
  }
  else {
    const size_t query_chrom_size = c.head.qSize;
    for (size_t i = 0; i < c.alns.size(); ++i) {
      os << ref_chrom << '\t'
         << ref_position << '\t'
         << ref_position + c.alns[i].size << '\t';
      os << query_chrom << '\t'
         << (query_chrom_size - (query_position + c.alns[i].size)) << '\t'
         << (query_chrom_size - query_position) << '\t'
         << '-' << endl;
      ref_position += (c.alns[i].size + c.alns[i].dt);
      query_position += (c.alns[i].size + c.alns[i].dq);
    }
  }
  return os;
}

std::istream &
operator>>(std::istream &in, chain &c) {
  if (in >> c.head) {
    string line;
    aln_dat ad_in;
    vector<aln_dat> tmp_alns;
    while (getline(in, line) && !line.empty()) {
      std::istringstream iss(line);
      if (iss >> ad_in)
        tmp_alns.push_back(ad_in);
    }
    c.alns.swap(tmp_alns);
  }
  return in;
};

int
main(int argc, const char **argv) {
  try{
    bool VERBOSE = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "expand a chain format file",
                           "<input> <output>");
    opt_parse.add_opt("verbose", 'v', "print more information",
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
    const string outfile(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "[reading input chains]" << endl;
    vector<chain> the_chains;
    chain tmp_chain;
    std::ifstream in_chains(chain_file.c_str());
    while (in_chains >> tmp_chain)
      the_chains.push_back(tmp_chain);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    if (VERBOSE)
      cerr << "[writing expanded chains]" << endl;

    for (size_t i = 0; i < the_chains.size(); ++i)
      out << the_chains[i];

  }
  catch (const std::runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
