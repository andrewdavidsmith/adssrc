/*    autocorr: get autocorrelation stats from a methcounts file
 *
 *    Copyright (C) 2015 University of Southern California
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or
 *    modify it under the terms of the GNU General Public License as
 *    published by the Free Software Foundation, either version 3 of
 *    the License, or (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public
 *    License along with this program.  If not, see
 *    <http://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <tr1/cmath>
#include <numeric>

#include "OptionParser.hpp"
#include "GenomicRegion.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "MethpipeSite.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::sort;

typedef float methType;


template <class T> static T
corr(const vector<T> &x, const vector<T> &y) {

  const T n = x.size();

  const T x_mean = accumulate(x.begin(), x.end(), 0.0)/n;
  const T y_mean = accumulate(y.begin(), y.end(), 0.0)/n;

  const T ip = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);

  const T x_ss = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
  const T y_ss = std::inner_product(y.begin(), y.end(), y.begin(), 0.0);

  return (ip - n*x_mean*y_mean)/
    (std::sqrt(x_ss - n*x_mean*x_mean)*(std::sqrt(y_ss - n*y_mean*y_mean)));

}


static bool
region_precedes_site(const vector<GenomicRegion> &regions,
                     const size_t region_idx, const MSite &s) {
  return (region_idx < regions.size() &&
          (regions[region_idx].get_chrom() < s.chrom ||
           (regions[region_idx].get_chrom() == s.chrom &&
            regions[region_idx].get_end() <= s.pos)));
}

static bool
region_contains_site(const vector<GenomicRegion> &regions,
                     const size_t region_idx, const MSite &s) {
  // !!! ASSUMES REGION DOES NOT PRECEDE SITE
  return (region_idx < regions.size() &&
          regions[region_idx].get_chrom() == s.chrom &&
          regions[region_idx].get_start() <= s.pos);
}


static size_t
boundary_next_to_site(const vector<GenomicRegion> &regions,
                      size_t &region_idx, const MSite &s) {
  size_t b = s.pos;
  while (region_precedes_site(regions, region_idx, s))
    ++region_idx;
  const bool contained = region_contains_site(regions, region_idx, s);
  if (contained) b = regions[region_idx].get_end();
  else
    if (regions[region_idx].get_chrom() == s.chrom)
      b = regions[region_idx].get_start();
    else b = std::numeric_limits<size_t>::max(); 

  return b;
}


static void
process_chrom(const bool not_span, const vector<GenomicRegion> &regions,
              size_t &region_idx,
              const size_t min_reads, const size_t max_dist,
              const vector<MSite> &sites,
              vector<vector<methType> > &a, vector<vector<methType> > &b) {
  
  for (size_t i = 0; i < sites.size(); ++i) {
    if (sites[i].n_reads >= min_reads) {
      size_t pos_limit = sites[i].pos + max_dist;
      if (not_span)
        pos_limit = std::min(pos_limit,
            boundary_next_to_site(regions, region_idx, sites[i-1]));
      size_t j = i + 1;
      while (j < sites.size() && sites[j].pos <= pos_limit) {
        if (sites[j].n_reads >= min_reads) {
          const size_t dist = sites[j].pos - sites[i].pos;
          a[dist].push_back(sites[i].meth);
          b[dist].push_back(sites[j].meth);
        }
        ++j;
      }
    }
  }
}


static void
process_chrom_nearest(const bool not_span,
                      const vector<GenomicRegion> &regions,
                      size_t &region_idx,
                      const size_t min_reads, const size_t max_dist,
                      const vector<MSite> &sites,
                      vector<vector<methType> > &a,
                      vector<vector<methType> > &b) {

  for (size_t i = 1; i < sites.size(); ++i)
    if (sites[i].n_reads >= min_reads &&
        sites[i - 1].n_reads >= min_reads) {
      const size_t dist = sites[i].pos - sites[i - 1].pos;
      const size_t r = not_span
                       ? boundary_next_to_site(regions, region_idx, sites[i-1])
                       : 0;
      if (dist <= max_dist && (!not_span || sites[i].pos < r)) {
        a[dist].push_back(sites[i-1].meth);
        b[dist].push_back(sites[i].meth);
      }
    }
}


static void
report_progress_to_stderr(std::ifstream &in,
                          const size_t report_frequency,
                          const size_t filesize) {
  if (in.tellg() % report_frequency == 0)
    cerr << smithlab::toa(percent(in.tellg(), filesize)) << "%\r";
}



static void
load_regions(const string &regions_file,
             vector<GenomicRegion> &regions) {

  std::ifstream in(regions_file.c_str());
  if (!in)
    throw SMITHLABException("bad regions file: " + regions_file);

  GenomicRegion r;
  while (in >> r)
    regions.push_back(r);

  if (!check_sorted(regions))
    throw SMITHLABException("regions file not sorted: " + regions_file);
}


static bool
site_allowed(const bool exclude_regions,
             const vector<GenomicRegion> &regions,
             const MSite &s, size_t &region_idx) {
  while (region_precedes_site(regions, region_idx, s))
    ++region_idx;
  const bool contained = region_contains_site(regions, region_idx, s);
  return (exclude_regions ? !contained : contained);
}


int main(int argc, const char **argv) {

  try {

    static const size_t report_frequency = 1024;

    /* FILES */
    string outfile;
    size_t max_dist = 4000;
    size_t min_reads = 10;
    size_t min_sites = 500;

    string regions_file;

    bool VERBOSE = false;
    bool PROGRESS = false;

    bool nearest_neighbors = false;
    bool not_span = false;
    bool exclude_regions = false;

    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(strip_path(argv[0]), "", "<meth-file>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("dist", 'd', "max dist",
                      false , max_dist);
    opt_parse.add_opt("sites", 's', "min sites for correlation",
                      false , min_sites);
    opt_parse.add_opt("progress", 'P', "report progress", false, PROGRESS);
    opt_parse.add_opt("nearest", 'N', "use only nearest neighbors",
                      false , nearest_neighbors);
    opt_parse.add_opt("notspan", 'S', "cancel corr spanning distinct regions",
                      false , not_span);
    opt_parse.add_opt("reads", 'r', "min reads", false, min_reads);
    opt_parse.add_opt("regions", '\0', "regions file", false, regions_file);
    opt_parse.add_opt("exclude", 'E', "exclude regions (default include)",
                      false, exclude_regions);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string cpg_file_name = leftover_args.front();
    /**********************************************************************/

    const size_t filesize = get_filesize(cpg_file_name);

    if (VERBOSE)
      cerr << "min_reads : " << min_reads << endl
           << "max_dist : " << max_dist << endl
           << "input_file : " << cpg_file_name << endl
           << "input_size : " << filesize << endl;

    std::ifstream in(cpg_file_name.c_str());
    if (!in)
      throw SMITHLABException("bad input file: " + cpg_file_name);

    vector<MSite> sites;
    MSite s;
    string prev_chrom;
    vector<vector<methType> > x(max_dist + 1), y(max_dist + 1);

    vector<GenomicRegion> regions;
    if (!regions_file.empty())
      load_regions(regions_file, regions);
    size_t region_idx = 0;
    size_t region_process_idx = 0;

    while (in >> s) {
      if (s.chrom != prev_chrom) {
        if (nearest_neighbors)
          process_chrom_nearest(not_span, regions, region_process_idx,
                                min_reads, max_dist, sites, x, y);
        else
          process_chrom(not_span, regions, region_process_idx,
                        min_reads, max_dist, sites, x, y);
        sites.clear();
      }

      if (regions_file.empty() ||
          site_allowed(exclude_regions, regions, s, region_idx))
        sites.push_back(s);
      prev_chrom.swap(s.chrom);
      
      if (PROGRESS)
        report_progress_to_stderr(in, report_frequency, filesize);
    }

    if (nearest_neighbors)
      process_chrom_nearest(not_span, regions, region_process_idx,
                            min_reads, max_dist, sites, x, y);
    else
      process_chrom(not_span, regions, region_process_idx,
                    min_reads, max_dist, sites, x, y);
    if (PROGRESS) cerr << "100%" << endl;

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    for (size_t i = 1; i < x.size(); ++i)
      if (x[i].size() >= min_sites)
        out << i << '\t' << corr(x[i], y[i])  << '\t' << x[i].size() << endl;

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
