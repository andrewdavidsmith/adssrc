/*    binmeth: a program to compute methylation levels in genomic bins
 *    (fixed width intervals that partition the genome).
 *
 *    Copyright (C) 2015  University of Southern California and
 *                        Andrew D. Smith
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
#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>
#include <tr1/cmath>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "MethpipeFiles.hpp"

#include "bsutils.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::tr1::unordered_map;
using std::pair;
using std::make_pair;

struct Site {
  string chrom;
  size_t pos;
  string strand;
  string context;
  double meth;
  size_t n_reads;

  size_t n_meth() const {return std::tr1::round(meth*n_reads);}

  void add(const Site &other) {
    if (!is_mutated() && other.is_mutated())
      context += 'x';
    // ADS: order matters below as n_reads update invalidates n_meth()
    // function until meth has been updated
    const size_t total_c_reads = n_meth() + other.n_meth();
    n_reads += other.n_reads;
    meth = static_cast<double>(total_c_reads)/n_reads;
  }

  // ADS: function below has redundant check for is_cpg, which is
  // expensive and might be ok to remove
  bool is_mate_of(const Site &first) {
    return (first.pos + 1 == pos && first.is_cpg() && is_cpg() &&
            first.strand == "+" && strand == "-");
  }
  ////////////////////////////////////////////////////////////////////////
  /////  Functions below test the type of site. These are CpG, CHH, and
  /////  CHG divided into two kinds: CCG and CXG, the former including a
  /////  CpG within. Also included is a function that tests if a site
  /////  has a mutation.
  /////  WARNING: None of these functions test for the length of their
  /////  argument string, which could cause problems.
  ////////////////////////////////////////////////////////////////////////
  bool is_cpg() const {
    return (context[0] == 'C' && context[1] == 'p' && context[2] == 'G');
  }
  bool is_chh() const {
    return (context[0] == 'C' && context[1] == 'H' && context[2] == 'H');
  }
  bool is_ccg() const {
    return (context[0] == 'C' && context[1] == 'C' && context[2] == 'G');
  }
  bool is_cxg() const {
    return (context[0] == 'C' && context[1] == 'X' && context[2] == 'G');
  }
  bool is_mutated() const {
    return context[3] == 'x';
  }
};


struct CountSet {
  size_t total_sites;
  size_t sites_covered;
  size_t max_coverage;
  size_t mutations;
  size_t total_c, total_t;
  size_t called_meth, called_unmeth;
  double mean_agg;
  CountSet() : total_sites(0), sites_covered(0), max_coverage(0),
    mutations(0), total_c(0), total_t(0),
    called_meth(0), called_unmeth(0),
    mean_agg(0.0) {}

  static void
  clear_counts(CountSet &cpg, CountSet &cpg_symm, CountSet &chh,
               CountSet &cxg, CountSet &ccg, CountSet &all_c);

  void update(const Site &s) {
    if (s.is_mutated()) {
      ++mutations;
    }
    else if (s.n_reads > 0) {
      ++sites_covered;
      max_coverage = std::max(max_coverage, s.n_reads);
      total_c += s.n_meth();
      total_t += s.n_reads - s.n_meth();
      mean_agg += s.meth;
      double lower = 0.0, upper = 0.0;
      wilson_ci_for_binomial(alpha, s.n_reads, s.meth, lower, upper);
      called_meth += (lower > 0.5);
      called_unmeth += (upper < 0.5);
    }
    ++total_sites;
  }

  size_t coverage() const {return total_c + total_t;}
  size_t total_called() const {return called_meth + called_unmeth;}

  double weighted_mean_meth() const {
    return static_cast<double>(total_c)/coverage();
  }
  double fractional_meth() const {
    return static_cast<double>(called_meth)/total_called();
  }
  double mean_meth() const {
    return mean_agg/sites_covered;
  }

  string format_summary(const string &context) const {
    std::ostringstream oss;
    const bool good = (sites_covered != 0);
    oss << "METHYLATION LEVELS (" + context + " CONTEXT):\n"
        << '\t' << "sites" << '\t' << total_sites << '\n'
        << '\t' << "sites_covered" << '\t' << sites_covered << '\n'
        << '\t' << "fraction_covered" << '\t'
        << static_cast<double>(sites_covered)/total_sites << '\n'
        << '\t' << "mean_depth" << '\t'
        << static_cast<double>(coverage())/total_sites << '\n'
        << '\t' << "mean_depth_covered" << '\t'
        << static_cast<double>(coverage())/sites_covered << '\n'
        << '\t' << "max_depth" << '\t' << max_coverage << '\n'
        << '\t' << "mutations" << '\t' << mutations << '\n'
        << '\t' << "mean_meth" << '\t'
        << (good ? toa(mean_meth()) : "N/A")  << '\n'
        << '\t' << "w_mean_meth" << '\t'
        << (good ? toa(weighted_mean_meth()) : "N/A") << '\n'
        << '\t' << "frac_meth" << '\t'
        << (good ? toa(fractional_meth()) : "N/A");
    return oss.str();
  }

  static double alpha;
};

double CountSet::alpha = 0.95;



void
CountSet::clear_counts(CountSet &cpg, CountSet &cpg_symm, CountSet &chh,
                       CountSet &cxg, CountSet &ccg, CountSet &all_c) {
  cpg = CountSet();
  cpg_symm = CountSet();
  chh = CountSet();
  cxg = CountSet();
  ccg = CountSet();
  all_c = CountSet();
}



static bool
get_meth_unmeth(std::ifstream &in, Site &site) {
  return methpipe::read_site(in, site.chrom, site.pos, site.strand,
                             site.context, site.meth, site.n_reads);
}



static void
write_interval(const string &chrom_name,
               const size_t start_pos, const size_t bin_size,
               const CountSet &cpg, const CountSet &cpg_symm,
               const CountSet &chh, const CountSet &cxg,
               const CountSet &ccg, const CountSet &all_c,
               std::ostream &out) {

  out << chrom_name << '\t'
      << start_pos << '\t'
      << start_pos + bin_size << '\t'
      << cpg.total_sites << ':'
      << cpg.sites_covered << ':'
      << cpg.total_c << ':'
      << cpg.coverage() << '\t'
      << cpg.weighted_mean_meth() << '\t'
      << '+' << endl;
}



static bool
bin_precedes_site(const string &chrom_name,
                  const size_t end_pos, const Site &site) {
  return chrom_name != site.chrom || end_pos <= site.pos;
}



static bool
load_chrom_sizes(std::ifstream &in,
                 vector<string> &chrom_names,
                 vector<size_t> &chrom_sizes) {
  string name;
  size_t size = 0;
  vector<pair<string, size_t> > ns;
  while (in >> name >> size)
    ns.push_back(make_pair(name, size));

  sort(ns.begin(), ns.end());

  for (size_t i = 0; i < ns.size(); ++i) {
    chrom_names.push_back(ns[i].first);
    chrom_sizes.push_back(ns[i].second);
  }
  return true;
}



static bool
increment_bin(const size_t bin_size,
              const vector<string> &chrom_names,
              const vector<size_t> &chrom_sizes,
              size_t &chrom_index,
              string &chrom_name,
              size_t &chrom_size,
              size_t &start_pos) {

  if (start_pos + bin_size < chrom_size)
    start_pos += bin_size;

  else {
    start_pos = 0;
    ++chrom_index;
    if (chrom_index >= chrom_sizes.size())
      return false;
    chrom_size = chrom_sizes[chrom_index];
    chrom_name = chrom_names[chrom_index];
  }
  return true;
}



int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    size_t bin_size = 1000;
    string outfile;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "methylation levels in bins",
                           "<chrom-sizes> <methcounts-file>");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("binsize", 'b', "size of bins",
                      false, bin_size);
    opt_parse.add_opt("alpha", 'a', "alpha for confidence interval",
                      false, CountSet::alpha);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (opt_parse.help_requested()) {
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
    const string chrom_sizes_file = leftover_args.front();
    const string meth_file = leftover_args.back();
    /****************** END COMMAND LINE OPTIONS *****************/

    vector<string> chrom_names;
    vector<size_t> chrom_sizes;
    std::ifstream chrom_sizes_in(chrom_sizes_file.c_str());
    if (!chrom_sizes_in ||
        !load_chrom_sizes(chrom_sizes_in, chrom_names, chrom_sizes))
      throw SMITHLABException("bad chrom sizes file: " + chrom_sizes_file);
    size_t chrom_index = 0;
    string chrom_name = chrom_names[chrom_index];
    size_t chrom_size = chrom_sizes[chrom_index];

    std::ifstream in(meth_file.c_str());
    if (!in)
      throw SMITHLABException("bad input file: " + meth_file);

    CountSet cpg, cpg_symm, chh, cxg, ccg, all_c;
    Site site, prev_site;

    size_t start_pos = 0ul;

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    while (get_meth_unmeth(in, site)) {

      // need to make sure chrom variable has right chrom for printing
      while (bin_precedes_site(chrom_name,
                               start_pos + bin_size, site)) {
        write_interval(chrom_name, start_pos, bin_size,
                       cpg, cpg_symm, chh, cxg, ccg, all_c, out);

        CountSet::clear_counts(cpg, cpg_symm, chh, cxg, ccg, all_c);
        increment_bin(bin_size, chrom_names, chrom_sizes,
                      chrom_index, chrom_name, chrom_size, start_pos);
      }

      if (site.chrom != prev_site.chrom)
        if (VERBOSE)
          cerr << "PROCESSING:\t" << site.chrom << "\n";

      if (site.is_cpg()) {
        cpg.update(site);
        if (site.is_mate_of(prev_site)) {
          site.add(prev_site);
          cpg_symm.update(site);
        }
      }
      else if (site.is_chh())
        chh.update(site);
      else if (site.is_ccg())
        ccg.update(site);
      else if (site.is_cxg())
        cxg.update(site);
      else
        throw SMITHLABException("bad site context: " + site.context);

      all_c.update(site);

      prev_site = site;
    }

    do {
      write_interval(chrom_name, start_pos, bin_size,
                     cpg, cpg_symm, chh, cxg, ccg, all_c, out);
      CountSet::clear_counts(cpg, cpg_symm, chh, cxg, ccg, all_c);
    }
    while (increment_bin(bin_size, chrom_names, chrom_sizes,
                         chrom_index, chrom_name, chrom_size, start_pos));
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
