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

#ifndef CHAIN_FILE_UTILS_HPP
#define CHAIN_FILE_UTILS_HPP

#include <string>
#include <vector>
#include <iostream>

struct aln_dat {
  size_t size; /* the size of the ungapped alignment */
  size_t dt;   /* the difference between the end of this block and the
                  beginning of the next block (reference sequence) */
  size_t dq;   /* the difference between the end of this block and the
                  beginning of the next block (query sequence) */
};

std::ostream &
operator<<(std::ostream &os, const aln_dat &ad);

std::istream &
operator>>(std::istream &is, aln_dat &ad);

struct chain_head {
  size_t score; // chain score
  std::string tName; // chromosome (reference sequence)
  size_t tSize; // chromosome size (reference sequence)
  char tStrand; // strand (reference sequence)
  size_t tStart; // alignment start position (reference sequence)
  size_t tEnd; // alignment end position (reference sequence)
  std::string qName; // chromosome (query sequence)
  size_t qSize; // chromosome size (query sequence)
  char qStrand; // strand (query sequence)
  size_t qStart; // alignment start position (query sequence)
  size_t qEnd; // alignment end position (query sequence)
  size_t id; // chain ID
};

std::ostream &
operator<<(std::ostream &os, const chain_head &ch);

std::istream &
operator>>(std::istream &is, chain_head &ch);

struct chain {
  chain_head head;
  std::vector<aln_dat> alns; // alignment data lines
};

std::ostream &
operator<<(std::ostream &os, const chain &c);

std::istream &
operator>>(std::istream &in, chain &c);

#include <GenomicRegion.hpp>

struct chain_block {
  SimpleGenomicRegion ref;
  SimpleGenomicRegion query;
  char query_strand;
  bool operator<(const chain_block &other) const {return ref < other.ref;}
};

std::ostream &
operator<<(std::ostream &os, const chain_block &b);

std::istream &
operator>>(std::istream &in, chain_block &b);

void
parse_species_from_chain_file_name(const std::string &chain_file,
                                   std::string &spec1, std::string &spec2);

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// void
// add_compact_site_pairs_from_blocks(const AssemblyData &asd,
//                                    const std::string &spec1,
//                                    const std::string &spec2,
//                                    const chain_block &b,
//                                    std::vector<std::pair<compact_site,
//                                                compact_site> > &cs);


#endif
