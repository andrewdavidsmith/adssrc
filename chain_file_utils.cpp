/* Copyright (C) 2018 University of Southern California and
 *                    Andrew D. Smith
 *
 * Authors: Andrew D. Smith
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

#include "chain_file_utils.hpp"

#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

using std::string;
using std::ostream;
using std::istream;
using std::vector;
using std::endl;

ostream &
operator<<(ostream &os, const aln_dat &ad) {
  os << ad.size;
  if (ad.dt > 0 || ad.dq > 0)
    os << '\t' << ad.dt << '\t' << ad.dq;
  return os;
}

istream &
operator>>(istream &is, aln_dat &ad) {
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

ostream &
operator<<(ostream &os, const chain_head &ch) {
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

istream &
operator>>(istream &is, chain_head &ch) {
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

ostream &
operator<<(ostream &os, const chain &c) {
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

istream &
operator>>(istream &in, chain &c) {
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

void
parse_species_from_chain_file_name(const string &chain_file,
                                   string &spec1, string &spec2) {
  const size_t first_dot = chain_file.find('.');
  spec1 = chain_file.substr(0, first_dot);
  const size_t second_dot = chain_file.find('.', first_dot + 1);
  spec2 = chain_file.substr(first_dot + 1, second_dot - first_dot - 1);
}

std::ostream &
operator<<(std::ostream &os, const chain_block &b) {
  os << b.ref << '\t' << b.query << '\t' << b.query_strand;
  return os;
}

std::istream &
operator>>(std::istream &in, chain_block &b) {

  string r_chrom;
  size_t r_start = 0, r_end = 0;

  string q_chrom;
  size_t q_start = 0, q_end = 0;
  char q_strand;

  if (in
      >> r_chrom >> r_start >> r_end
      >> q_chrom >> q_start >> q_end >> q_strand) {
    b.ref = SimpleGenomicRegion(r_chrom, r_start, r_end);
    b.query = SimpleGenomicRegion(q_chrom, q_start, q_end);
    b.query_strand = q_strand;
  }
  return in;
}

