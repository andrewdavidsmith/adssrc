/*    binmeth: a program to compute methylation levels in genomic bins
 *    (fixed width intervals that partition the genome).
 *
 *    Copyright (C) 2016  University of Southern California and
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
#include <cmath>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "MethpipeFiles.hpp"
#include "MethpipeSite.hpp"

#include "bsutils.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
// using std::unordered_map;
using std::pair;
using std::make_pair;


// function used to ensure input sites are sorted
static bool
precedes(const MSite &prev, const MSite &curr) {
  const int chrom_comp = prev.chrom.compare(curr.chrom);
  return (chrom_comp < 0 ||
          (chrom_comp == 0 && prev.pos < curr.pos));
}


// aggregation of data related to bins; intended use as singleton
struct BinInfo {

  string chrom_name;
  size_t chrom_size;
  size_t start_pos;
  size_t chrom_index;
  size_t bin_size;

  BinInfo(const size_t input_bin_size,
          const vector<string> &chrom_names,
          const vector<size_t> &chrom_sizes) {
    bin_size = input_bin_size;
    chrom_index = 0;
    start_pos = 0;
    chrom_name = chrom_names[chrom_index];
    chrom_size = chrom_sizes[chrom_index];
  }

  bool
  increment(const vector<string> &chrom_names,
            const vector<size_t> &chrom_sizes);
};


bool
BinInfo::increment(const vector<string> &chrom_names,
                   const vector<size_t> &chrom_sizes) {

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


static bool
bin_precedes_site(BinInfo &bin, const MSite &site) {
  const int chrom_comp = bin.chrom_name.compare(site.chrom);
  return (chrom_comp < 0 ||
          (chrom_comp == 0 &&
           (bin.start_pos + bin.bin_size <= site.pos)));
}


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

  void update(const MSite &s) {
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


static void
write_interval(const bool PRINT_ADDITIONAL_LEVELS,
               const BinInfo &bin,
               const CountSet &cpg,
               // const CountSet &cpg_symm,
               // const CountSet &chh, const CountSet &cxg,
               // const CountSet &ccg, const CountSet &all_c,
               std::ostream &out) {

  out << bin.chrom_name << '\t'
      << bin.start_pos << '\t'
      << bin.start_pos + bin.bin_size << '\t'
      << cpg.total_sites << ':'
      << cpg.sites_covered << ':'
      << cpg.total_c << ':'
      << cpg.coverage() << '\t'
      << cpg.weighted_mean_meth() << '\t'
      << '+';

  if (PRINT_ADDITIONAL_LEVELS)
    out << '\t'
        << cpg.fractional_meth() << '\t'
        << cpg.mean_meth();

  out << '\n';
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


static std::istream &
get_meth_unmeth(std::istream &in, MSite &site) {
  return methpipe::read_site(in, site.chrom, site.pos, site.strand,
                             site.context, site.meth, site.n_reads);
}


int
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    size_t bin_size = 1000;
    string outfile;
    bool PRINT_ADDITIONAL_LEVELS = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "methylation levels in bins",
                           "<chrom-sizes> <methcounts-file>");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("binsize", 'b', "size of bins",
                      false, bin_size);
    opt_parse.add_opt("alpha", 'a', "alpha for confidence interval",
                      false, CountSet::alpha);
    opt_parse.add_opt("more-levels", 'M', "print more meth level information",
                      false, PRINT_ADDITIONAL_LEVELS);
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

    std::ifstream in(meth_file.c_str());
    if (!in)
      throw SMITHLABException("bad input file: " + meth_file);

    CountSet cpg, cpg_symm, chh, cxg, ccg, all_c;
    MSite site, prev_site;

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    // BinInfo (singleton) to keep track of the current chromosome
    // while iterating over sites
    BinInfo bin(bin_size, chrom_names, chrom_sizes);

    // iterate over sites (lines in methcounts output file)
    while (get_meth_unmeth(in, site)) {

      if (!precedes(prev_site, site))
        throw SMITHLABException("sites not sorted in: " + meth_file);

      // need to make sure curr_chrom_name variable has right chrom
      // for printing
      while (bin_precedes_site(bin, site)) {
        write_interval(PRINT_ADDITIONAL_LEVELS,
                       bin, cpg, out); //, cpg_symm, chh, cxg, ccg, all_c, out);
        CountSet::clear_counts(cpg, cpg_symm, chh, cxg, ccg, all_c);
        bin.increment(chrom_names, chrom_sizes);
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
      else if (site.is_chh()) chh.update(site);
      else if (site.is_ccg()) ccg.update(site);
      else if (site.is_cxg()) cxg.update(site);
      else
        throw SMITHLABException("bad site context: " + site.context);

      all_c.update(site);

      prev_site = site;
    }

    do {
      write_interval(PRINT_ADDITIONAL_LEVELS,
                     bin, cpg, out); //, cpg_symm, chh, cxg, ccg, all_c, out);
      CountSet::clear_counts(cpg, cpg_symm, chh, cxg, ccg, all_c);
    }
    while (bin.increment(chrom_names, chrom_sizes));
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
