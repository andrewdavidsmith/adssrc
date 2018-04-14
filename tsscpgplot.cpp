/*    tsscpgplot: get data to plot methylation level around a TSS
 *
 *    Copyright (C) 2010-2016 Andrew D. Smith
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

#include <fstream>
#include <unordered_map>
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "MethpipeFiles.hpp"

using std::string;
using std::vector;
using std::ostream;
using std::endl;
using std::cerr;
using std::cout;
using std::pair;
using std::make_pair;
using std::sort;
using std::unordered_map;

static size_t
get_reads(const GenomicRegion &cpg) {
  const string region_name(cpg.get_name());
  const size_t colon_offset = region_name.find(":");
  return atoi(region_name.substr(colon_offset + 1).c_str());
}


static void
process_chrom(const size_t region_size,
              const vector<GenomicRegion> &cpg,
              const vector<GenomicRegion> &tss,
              vector<double> &totals,
              vector<size_t> &counts,
              vector<size_t> &density) {

  for (size_t i = 0; i < tss.size(); ++i) {

    GenomicRegion a(tss[i]);
    a.set_start(a.get_start() - region_size);
    a.set_end(a.get_start() + 1);

    GenomicRegion b(tss[i]);
    b.set_end(b.get_end() + region_size);
    b.set_start(b.get_end() - 1);

    const size_t left = a.get_start();
    const size_t right = left + 2*region_size;
    vector<GenomicRegion>::const_iterator start(find_closest(cpg, a));
    vector<GenomicRegion>::const_iterator end(find_closest(cpg, b));

    if ((*start < tss[i]) && (tss[i] < *end)) {

      const size_t idx = (start - cpg.begin());
      const size_t lim = (end - cpg.begin());

      if (tss[i].pos_strand()) {
        for (size_t k = idx; k <= lim; ++k)
          if (cpg[k].get_start() >= left && cpg[k].get_start() < right) {
            const size_t base_idx = cpg[k].get_start() - left;
            const size_t reads = get_reads(cpg[k]);
            totals[base_idx] += cpg[k].get_score()*reads;
            counts[base_idx] += reads;
            density[base_idx] += 1;
          }
      }
      else {
        for (size_t k = idx; k <= lim; ++k)
          if (cpg[k].get_start() >= left && cpg[k].get_start() < right) {
            const size_t base_idx = right - 1 - cpg[k].get_start();
            const size_t reads = get_reads(cpg[k]);
            totals[base_idx] += cpg[k].get_score()*reads;
            counts[base_idx] += reads;
            density[base_idx] += 1;
          }
      }
    }
  }
}


static void
collapse_bins(const size_t bin_size, vector<double> &totals,
              vector<size_t> &counts, vector<size_t> &density) {

  const size_t n_bins =
    std::ceil(static_cast<double>(totals.size())/bin_size);

  vector<double> t(n_bins, 0.0);
  vector<size_t> c(n_bins, 0ul);
  vector<size_t> d(n_bins,0ul);
  for (size_t i = 0; i < totals.size(); ++i) {
    t[i/bin_size] += totals[i];
    c[i/bin_size] += counts[i];
    d[i/bin_size] += density[i];
  }
  totals.swap(t);
  counts.swap(c);
  density.swap(d);
}



static void
extract_methpipe_format_cpg_region(const string &buffer,
                                   GenomicRegion &r) {
  std::istringstream is(buffer);
  string chrom, name;
  size_t pos = 0ul, coverage = 0ul;
  char strand;
  double meth_level;
  is >> chrom >> pos >> strand >> name >> meth_level >> coverage;
  r = GenomicRegion(chrom, pos, pos + 1, name + ":" + toa(coverage),
                    meth_level, strand);
}



int main(int argc, const char **argv) {

  try {
    /* FILES */
    string outfile;
    size_t region_size = 10000;
    bool VERBOSE = false;
    size_t bin_size = 50;

    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(strip_path(argv[0]), "", "<TSS_FILE> <CPG_FILE>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("size", 's', "region size",
                      false , region_size);
    opt_parse.add_opt("bin", 'b', "bin size",
                      false , bin_size);
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
    const string tss_file_name = leftover_args.front();
    const string cpg_file_name = leftover_args.back();
    /**********************************************************************/

    if (VERBOSE)
      cerr << "loading tss data" << endl;
    vector<GenomicRegion> tss_all;
    ReadBEDFile(tss_file_name, tss_all);
    assert(check_sorted(tss_all));

    if (VERBOSE)
      cerr << "loading cpg data" << endl;
    vector<GenomicRegion> cpg_all;
    if (methpipe::is_methpipe_file_single(cpg_file_name)) {
      if (VERBOSE)
        cerr << "format is methpipe" << endl;
      std::ifstream cpgin(cpg_file_name.c_str());
      if (!cpgin)
        throw SMITHLABException("bad file: " + cpg_file_name);
      GenomicRegion r;
      string buffer;
      while (getline(cpgin, buffer)) {
        extract_methpipe_format_cpg_region(buffer, r);
        cpg_all.push_back(r);
      }
      assert(check_sorted(cpg_all));
    }
    else {
      if (VERBOSE)
        cerr << "format is bed" << endl;
      std::ifstream cpgin(cpg_file_name.c_str());
      if (!cpgin)
        throw SMITHLABException("bad file: " + cpg_file_name);
      GenomicRegion r;
      while (cpgin >> r)
        cpg_all.push_back(r);
      assert(check_sorted(cpg_all));
    }
    if (VERBOSE)
      cerr << "number of CpG sites: " << cpg_all.size() << endl;

    vector<vector<GenomicRegion> > tss;
    separate_chromosomes(tss_all, tss);

    vector<vector<GenomicRegion> > cpg;
    separate_chromosomes(cpg_all, cpg);

    if (VERBOSE)
      cerr << "making cpg lookup table" << endl;
    unordered_map<string, size_t> cpg_lookup;
    for (size_t i = 0; i < cpg.size(); ++i)
      cpg_lookup[cpg[i].front().get_chrom()] = i;

    vector<double> totals(2*region_size, 0.0);
    vector<size_t> counts(2*region_size, 0ul);
    vector<size_t> density(2*region_size,0ul);

    size_t total_tss = 0;
    for (size_t i = 0; i < tss.size(); ++i) {
      const unordered_map<string, size_t>::const_iterator j =
        cpg_lookup.find(tss[i][0].get_chrom());
      if (j != cpg_lookup.end()) {
        total_tss += tss[i].size();
        process_chrom(region_size, cpg[j->second], tss[i], totals, counts, density);
      }
    }

    collapse_bins(bin_size, totals, counts, density);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    if (VERBOSE)
      cerr << "format: "
           << "position" << '\t'
           << "fraction_of_tss" << '\t'
           << "fraction_of_counts" << '\t'
           << "mean count per region" << '\t' 
           << "mean density per region" << endl;

    for (size_t i = 0; i < totals.size(); ++i)
      out << i << "\t"
          << counts[i]/double(total_tss) << "\t"
          << totals[i]/counts[i] << "\t"
          << counts[i]/double(total_tss) << "\t"
          << density[i]/double(total_tss) << endl;
  }
  catch (SMITHLABException &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
