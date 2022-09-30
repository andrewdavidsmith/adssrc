/* autocorr: get autocorrelation stats from a counts file
 *
 * Copyright (C) 2015-2022 University of Southern California
 *                         Andrew D. Smith
 *
 * Author: Andrew D. Smith
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
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <numeric>
#include <random>

#include "OptionParser.hpp"
#include "GenomicRegion.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "MSite.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::sort;
using std::inner_product;
using std::runtime_error;


template <class T>
static double
corr(const vector<T> &X, const vector<T> &Y, T &sdX, T &sdY, T &covXY) {
  // calculate Pearson correlation coefficient of two vectors X and Y

  auto X_beg = begin(X);
  auto X_end = end(X);
  auto Y_beg = begin(Y);
  auto Y_end = end(Y);

  const T N = X.size(); // length of vectors X and Y (hope they are equal...)

  const T X_mean = accumulate(X_beg, X_end, 0.0)/N; // mean(X)
  const T Y_mean = accumulate(Y_beg, Y_end, 0.0)/N; // mean(Y)

  const T ip = inner_product(X_beg, X_end, Y_beg, 0.0); // X.Y

  const T X_ss = inner_product(X_beg, X_end, X_beg, 0.0); // Sum of sequares X
  const T Y_ss = inner_product(Y_beg, Y_end, Y_beg, 0.0); // Sum of sequares Y

  // sqrt of coveriance;
  // Sum XY - N.mu(X).mu(Y)
  covXY = ip - N*X_mean*Y_mean;

  // sqrt[SSX - N.mu(X).mu(X)]
  sdX = std::sqrt(X_ss - N*X_mean*X_mean);
  // sqrt[SSY - N.mu(Y).mu(Y)]
  sdY = std::sqrt(Y_ss - N*Y_mean*Y_mean);

  // Pearson correlation
  const double rXY = covXY/(sdX * sdY);

  const double size_factor = 1.0/std::sqrt(N - 1.0);
  sdX *= size_factor;
  sdY *= size_factor;

  return rXY;
}


static bool
region_precedes_site(const SimpleGenomicRegion &region, const MSite &site) {
  // check if the region precedes a the site; [a, b) doesn't contain x
  return (region.get_chrom() < site.chrom ||
          (region.get_chrom() == site.chrom && region.get_end() <= site.pos));
}


static bool
region_contains_site(const SimpleGenomicRegion region, const MSite &site) {
  // check if a given site is contained in a given region location
  // Containment is for half open intervals [a, b)
  return (region.get_chrom() == site.chrom &&
          region.get_start() <= site.pos && site.pos < region.get_end());
}


// Find the position of the region boundary (inclusive counting of
// bases) to the right of the current site; **assume the region
// contains the site**. Index of the region is updated.
static size_t
boundary_position(const vector<SimpleGenomicRegion> &regions,
                  size_t &idx, const MSite &site) {
  // move index of regions so the region doesn't entirely precede it
  while (idx < regions.size() && region_precedes_site(regions[idx], site))
    ++idx;

  // check if the region contains the site
  if (idx < regions.size() && region_contains_site(regions[idx], site))
    // the right side (closed interval) of containing region; need to
    // return the site that is within the region, since it will be
    // used to know when to stop, and that must be inclusive
    return regions[idx].get_end()-1;
  // by default return a very far boundary position
  else return std::numeric_limits<size_t>::max();

}


static bool
strands_are_good(const int strand, const MSite &a, const MSite &b) {
  return ((strand == 0) ||
          (strand == 1 && a.strand == b.strand) ||
          (strand == -1 && a.strand != b.strand));
}


static std::pair<size_t, size_t>
get_window_endpoint_bins(const size_t window_size, const size_t the_distance,
                         const size_t max_distance) {
  return std::make_pair(
    (the_distance - window_size) > 1 ? (the_distance - window_size) : 1,
    ((the_distance + window_size) < max_distance ?
     (the_distance + window_size) : max_distance)
                        );
}


static void
process_chrom(const bool require_same_region,
              const vector<SimpleGenomicRegion> &regions,
              size_t &region_idx,
              const size_t min_reads, const size_t max_dist,
              const size_t window_size, const vector<MSite> &sites,
              const size_t the_neighbor, const int strand,
              vector<vector<double> > &X, vector<vector<double> > &Y) {

  // assign CpGs in the given vector to the distance tables a & b;
  // each row corresponds to a gap distance in bp
  for (size_t i = 0; i < sites.size(); ++i) {
    // check if current site is covered by enough reads
    if (sites[i].n_reads >= min_reads) {

      // determine the limit of sites to consider
      size_t position_limit = sites[i].pos + max_dist;
      if (require_same_region)
        position_limit =
          std::min(position_limit,
                   boundary_position(regions, region_idx, sites[i]));

      if (the_neighbor == 0) { // use any neighbor in range
        size_t j = i + 1; // start with 1st neighbor
        for (; j < sites.size() && sites[j].pos <= position_limit; ++j)

          // check other site has enough reads and right orientation
          if (sites[j].n_reads >= min_reads &&
              strands_are_good(strand, sites[i], sites[j])) {

            // get the distance between current site and neighbor.
            // This seems directional, assuming j is beyond i, but the
            // opposite direction is included in the same bins at a
            // different iteration on i
            const size_t curr_dist = sites[j].pos - sites[i].pos;

            // get range of posns current sites should contribute to
            const std::pair<size_t, size_t> wb =
              get_window_endpoint_bins(window_size, curr_dist, max_dist);

            // include data for current sites in each relevant bin
            for (size_t d = wb.first; d <= wb.second; ++d) {
              X[d].push_back(sites[i].meth);
              Y[d].push_back(sites[j].meth);
            }
          }
      }
      else {
        size_t offset = i + the_neighbor; // use specified neighbor
        // make sure other site is within range
        if (offset < sites.size() && sites[offset].pos <= position_limit)

          // check if other site has enough reads
          if (sites[offset].n_reads >= min_reads &&
              strands_are_good(strand, sites[i], sites[offset])) {

            // as above: distance between current site and neighbor
            const size_t curr_dist = sites[offset].pos - sites[i].pos;

            // get range of posns current sites should contribute to
            const std::pair<size_t, size_t> wb =
              get_window_endpoint_bins(window_size, curr_dist, max_dist);

            // include data for current sites in each relevant bin
            for (size_t d = wb.first; d <= wb.second; ++d) {
              X[d].push_back(sites[i].meth);
              Y[d].push_back(sites[offset].meth);
            }
          }
      }
    }
  }
}


static void
load_regions(const string &regions_file, vector<SimpleGenomicRegion> &regions) {
  // load genomic intervals (regions)
  std::ifstream in(regions_file);
  if (!in)
    throw std::runtime_error("bad regions file: " + regions_file);

  SimpleGenomicRegion r;
  while (in >> r)
    regions.push_back(r);

  if (!is_sorted(begin(regions), end(regions)))
    throw std::runtime_error("regions file not sorted: " + regions_file);
}


static bool
site_allowed(const vector<SimpleGenomicRegion> &regions,
             const MSite &site, size_t &idx) {
  // check if a site is allowed to be used for correlation calculation
  // depending on whether you want to exlude or include certain regions
  while (idx < regions.size() && region_precedes_site(regions[idx], site))
    ++idx;
  return idx < regions.size() && region_contains_site(regions[idx], site);
}


static void
report_values(const string &valfile, const size_t num_val,
              const vector<vector<double> > &x,
              const vector<vector<double> > &y) {
  std::ofstream out_val(valfile);
  if (!out_val)
    throw runtime_error("cannot open values output file: " + valfile);

  std::random_device rd;
  std::mt19937 generator(rd());
  for (size_t i = 1; i < x.size(); ++i)
    if (!x[i].empty()) {
      vector<size_t> idx(x[i].size());
      std::iota(begin(idx), end(idx), 0);
      shuffle(begin(idx), end(idx), generator);
      idx.resize(std::min(idx.size(), num_val));
      // ADS: below likely not needed
      sort(begin(idx), end(idx));
      out_val << i;
      for (size_t j = 0; j < idx.size(); j++)
        out_val << '\t' << x[i][idx[j]] << '\t' << y[i][idx[j]];
      out_val << endl;
    }
}


int main(int argc, const char **argv) {

  try {

    string description = R"""(

This program computes statistics on the autocorrelation of methylation
levels. The input file must be in "counts" format from dnmtools. The
output file has the format of one line per distance between
sites. These distances begin with 2 since we assume CpG sites. Each
line has the follwing values: distance, correlation, N, sdX, sdY,
covXY. Typically only the first two of those are of interest.

  )""";

    string outfile;
    size_t max_dist = 4000;
    size_t min_reads = 10;
    size_t min_sites = 500;
    size_t the_neighbor = 0; // 0=nearest
    size_t window_size = 0;
    int strand = 0; // code: same=1, different=-1 and any=0

    string valfile; // report vals for debug here
    size_t num_val = 1000; // report vals for debug, how many?

    string regions_file;

    bool VERBOSE = false;
    bool PROGRESS = false;

    bool require_same_region = false;

    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(strip_path(argv[0]), description, "<meth-file>");
    opt_parse.add_opt("output", 'o', "name of output file (default: stdout)",
                      false , outfile);
    opt_parse.add_opt("max-dist", 'd', "maximum distance for pairs of sites",
                      false , max_dist);
    opt_parse.add_opt("min-sites", 's', "minimum sites needed for correlation",
                      false , min_sites);
    opt_parse.add_opt("min-coverage", 'c', "minimum coverage needed at a site",
                      false, min_reads);
    opt_parse.add_opt("nth-neighbor", 'n',
                      "use only nth neighbor (0=any, 1=nearest, ...)",
                      false , the_neighbor);
    opt_parse.add_opt("window-size", 'w', "smoothing window size for distances",
                      false , window_size);
    opt_parse.add_opt("strand", 't', "strand: same=1, opposite=-1, any=0",
                      false , strand);
    opt_parse.add_opt("regions", 'R', "file of regions to process "
                      "(bed format)", false, regions_file);
    opt_parse.add_opt("same-region", 'S', "require both sites in same region",
                      false , require_same_region);
    opt_parse.add_opt("valout", '\0', "write values in this file "
                      "(default: none)",
                      false , valfile);
    opt_parse.add_opt("report-n-vals", '\0', "number of values to output",
                      false , num_val);
    opt_parse.add_opt("progress", 'P', "report progress", false, PROGRESS);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false , VERBOSE);
    opt_parse.set_show_defaults();
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      cerr << opt_parse.about_message() << endl;
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
    const string input_filename = leftover_args.front();
    /*************************************************************************/

    const size_t input_filesize = get_filesize(input_filename);

    if (VERBOSE)
      cerr << "min_reads : " << min_reads << endl
           << "max_dist : " << max_dist << endl
           << "input_filename : " << input_filename << endl
           << "input_filesize : " << input_filesize << endl;

    std::ifstream in(input_filename);
    if (!in)
      throw std::runtime_error("bad input file: " + input_filename);

    vector<MSite> sites;
    MSite the_site;
    string prev_chrom;

    // a pair of vectors for each distance in the range of interest;
    // each pair of vectors includes each pair of values at the
    // corresponding distance
    vector<vector<double> > X(max_dist + 1), Y(max_dist + 1);

    vector<SimpleGenomicRegion> regions;
    if (!regions_file.empty()) {
      load_regions(regions_file, regions);
      if (VERBOSE)
        cerr << "loading regions: " << regions.size() << endl;
    }
    size_t region_idx_add = 0;
    size_t region_idx_process = 0;

    while (in >> the_site) {
      if (the_site.chrom != prev_chrom) {
        process_chrom(require_same_region, regions, region_idx_process,
                      min_reads, max_dist, window_size, sites,
                      the_neighbor, strand, X, Y);
        if (PROGRESS)
          cerr << "processing: " << the_site.chrom << endl;
        sites.clear();
      }
      if (regions.empty() || site_allowed(regions, the_site, region_idx_add))
        sites.push_back(the_site);
      prev_chrom.swap(the_site.chrom);
    }
    process_chrom(require_same_region, regions, region_idx_process,
                  min_reads, max_dist, window_size, sites,
                  the_neighbor, strand, X, Y);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile);
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    for (size_t i = 1; i < X.size(); ++i)
      if (X[i].size() >= min_sites) {
        double sdX = 0.0, sdY = 0.0, covXY = 0.0;
        const double rXY = corr(X[i], Y[i], sdX, sdY, covXY);
        const size_t N = X[i].size();
        out << i << '\t'
            << rXY << '\t'
            << N << '\t'
            << sdX << '\t'
            << sdY << '\t'
            << covXY << endl;
      }

    if (!valfile.empty()) // report (some or all of) the values
      report_values(valfile, num_val, X, Y);
  }
  catch (std::exception &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
