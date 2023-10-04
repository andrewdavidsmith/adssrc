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

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "MSite.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::runtime_error;
using std::end;
using std::begin;

struct genomic_interval {
  string chrom;
  size_t start_pos;
  size_t end_pos;
  bool operator<(const genomic_interval &other) const {
    const int x = chrom.compare(other.chrom);
    return (x < 0 ||
            (x == 0 &&
             (start_pos < other.start_pos ||
              (start_pos == other.start_pos &&
               (end_pos < other.end_pos)))));
  }
};


std::istream &
operator>>(std::istream &in, genomic_interval &gi) {
  string line;
  if (getline(in, line)) {
    std::istringstream iss(line);
    if (!(iss >> gi.chrom >> gi.start_pos >> gi.end_pos))
      in.setstate(std::ios_base::failbit);
  }
  return in;
}


static bool
region_precedes_site(const genomic_interval &region, const MSite &site) {
  // check if the region precedes a the site; [a, b) doesn't contain x
  const int x = region.chrom.compare(site.chrom);
  return x < 0 || (x == 0 && region.end_pos <= site.pos);
}


static bool
region_contains_site(const genomic_interval region, const MSite &site) {
  // check if a given site is contained in a given region location
  // Containment is for half open intervals [a, b)
  return (region.chrom == site.chrom &&
          region.start_pos <= site.pos && site.pos < region.end_pos);
}


// Find the position of the region boundary (inclusive counting of
// bases) to the right of the current site; **assume the region
// contains the site**. Index of the region is updated.
static size_t
boundary_position(const vector<genomic_interval> &regions,
                  size_t &idx, const MSite &site) {
  // move index of regions so the region doesn't entirely precede it
  while (idx < regions.size() && region_precedes_site(regions[idx], site))
    ++idx;

  // check if the region contains the site
  if (idx < regions.size() && region_contains_site(regions[idx], site))
    // the right side (closed interval) of containing region; need to
    // return the site that is within the region, since it will be
    // used to know when to stop, and that must be inclusive
    return regions[idx].end_pos - 1;
  // by default return a very far boundary position
  return std::numeric_limits<size_t>::max();
}


static bool
strands_are_good(const int strand, const MSite &a, const MSite &b) {
  return strand == 0 ||
    (strand == 1 && a.strand == b.strand) ||
    (strand == -1 && a.strand != b.strand);
}


struct sum_stats {
  size_t N;      // counts
  double X, Y;   // sums
  double XY;     // cross product
  double XX, YY; // sums of squares

  void
  update(const double x, const double y) {
    ++N;
    X += x;
    Y += y;
    XY += x*y;
    XX += x*x;
    YY += y*y;
  }

  double
  get_values(double &sdX, double &sdY, double &covXY) const {

    // Sum XY - N.mu(X).mu(Y) = Sum XY - Sum(X)Sum(Y)/N
    covXY = XY - (X*Y)/N;
    // sqrt[SSX - N.mu(X).mu(X)]
    sdX = std::sqrt(XX - (X*X)/N);
    // sqrt[SSY - N.mu(Y).mu(Y)]
    sdY = std::sqrt(YY - (Y*Y)/N);

    // Pearson correlation
    const double rXY = covXY/(sdX*sdY);

    const double size_factor = 1.0/std::sqrt(N - 1.0);
    sdX *= size_factor;
    sdY *= size_factor;

    return rXY;
  }
};


static bool
site_allowed(const vector<genomic_interval> &regions,
             const MSite &site, size_t &idx) {
  // check if a site is allowed to be used for correlation calculation
  // depending on whether you want to exlude or include certain regions
  while (idx < regions.size() && region_precedes_site(regions[idx], site))
    ++idx;
  return idx < regions.size() && region_contains_site(regions[idx], site);
}


static void
process_chrom(const bool require_same_region,
              const vector<genomic_interval> &regions,
              size_t &region_idx, // region_idx changes along chrom
              const size_t min_reads, const size_t max_dist,
              const size_t window_size, const vector<MSite> &sites,
              const int strand, vector<sum_stats> &the_stats) {

  // early exit if there are not enough sites
  if (sites.size() <= 1)
    return;

  // assign CpGs in the given vector to the distance tables a & b;
  // each row corresponds to a gap distance in bp
  const auto j_lim = end(sites);
  const auto i_lim = j_lim - 1;
  for (auto i = begin(sites); i != i_lim; ++i) {

    // check if current site is covered by enough reads
    if (i->n_reads >= min_reads) {

      // determine the limit of sites to consider
      size_t position_limit = i->pos + max_dist;
      if (require_same_region)
        position_limit =
          std::min(position_limit,
                   boundary_position(regions, region_idx, *i));

      // start j at the 1st neighbor
      for (auto j = i + 1; j != j_lim && j->pos <= position_limit; ++j) {

        // check other site has enough reads and right orientation
        if (j->n_reads >= min_reads && strands_are_good(strand, *i, *j)) {

          // get range of posns current sites should contribute to
          const size_t the_dist = j->pos - i->pos;
          size_t d = (the_dist > window_size) ? the_dist - window_size : 1;
          const size_t d_lim = std::min(the_dist + window_size, max_dist);

          // include data for current sites in each relevant bin
          while (d <= d_lim)
            the_stats[d++].update(i->meth, j->meth);
        }
      }
    }
  }
}


static void
process_sites_all_neighbors(const bool report_progress,
                            const string &filename,
                            const bool require_same_region,
                            const vector<genomic_interval> &regions,
                            const size_t min_reads, const size_t max_dist,
                            const size_t window_size, const int strand,
                            vector<sum_stats> &the_stats) {

  std::ifstream in(filename);
  if (!in)
    throw std::runtime_error("bad input file: " + filename);

  vector<MSite> sites;
  MSite the_site;
  string prev_chrom;

  size_t region_idx_add = 0;
  size_t region_idx_process = 0;

  while (in >> the_site) {
    if (the_site.chrom != prev_chrom) {
      process_chrom(require_same_region, regions, region_idx_process,
                    min_reads, max_dist, window_size, sites, strand, the_stats);
      if (report_progress)
        cerr << "processing: " << the_site.chrom << endl;
      sites.clear();
    }
    if (regions.empty() || site_allowed(regions, the_site, region_idx_add))
      sites.push_back(the_site);
    prev_chrom = std::move(the_site.chrom);
  }
  process_chrom(require_same_region, regions, region_idx_process,
                min_reads, max_dist, window_size, sites, strand, the_stats);
}


static void
process_chrom_kth(const bool require_same_region,
                  const vector<genomic_interval> &regions,
                  size_t &region_idx, // region_idx changes along chrom
                  const size_t min_reads, const size_t max_dist,
                  const size_t window_size, const vector<MSite> &sites,
                  const int strand, const size_t the_neighbor,
                  vector<sum_stats> &the_stats) {

  // early exit if there are not enough sites for specified neighbor
  if (sites.size() <= the_neighbor)
    return;

  // limit of outer iteration to allow for k-th neighbors
  const auto j_lim = end(sites);
  const auto i_lim = j_lim - the_neighbor - 1;

  // iterate over "left" site to get all pairs
  for (auto i = begin(sites); i != i_lim; ++i) {

    // check if current site is covered by enough reads
    if (i->n_reads >= min_reads) {

      // determine the limit of "right" sites to consider
      size_t position_limit = i->pos + max_dist;
      if (require_same_region)
        position_limit =
          std::min(position_limit,
                   boundary_position(regions, region_idx, *i));

      const auto j = i + the_neighbor; // get "right" neighbor site
      // make sure right neighbor is in range
      if (j != j_lim && j->pos <= position_limit) {

        // check if other site has enough reads
        if (j->n_reads >= min_reads && strands_are_good(strand, *i, *j)) {

          // get range of posns current sites should contribute to
          const size_t the_dist = j->pos - i->pos;
          size_t d = (the_dist > window_size) ? the_dist - window_size : 1;
          const size_t d_lim = std::min(the_dist + window_size, max_dist);
          while (d <= d_lim)
            the_stats[d++].update(i->meth, j->meth);
        }
      }
    }
  }
}


static void
process_sites_kth_neighbor(const bool report_progress,
                           const string &filename,
                           const bool require_same_region,
                           const vector<genomic_interval> &regions,
                           const size_t min_reads, const size_t max_dist,
                           const size_t window_size, const int strand,
                           const size_t the_neighbor,
                           vector<sum_stats> &the_stats) {

  std::ifstream in(filename);
  if (!in)
    throw std::runtime_error("bad input file: " + filename);

  vector<MSite> sites;
  MSite the_site;
  string prev_chrom;

  size_t region_idx_add = 0;
  size_t region_idx_process = 0;

  while (in >> the_site) {
    if (the_site.chrom != prev_chrom) {
      process_chrom_kth(require_same_region, regions, region_idx_process,
                        min_reads, max_dist, window_size, sites,
                        strand, the_neighbor, the_stats);
      if (report_progress)
        cerr << "processing: " << the_site.chrom << endl;
      sites.clear();
    }
    if (regions.empty() || site_allowed(regions, the_site, region_idx_add))
      sites.push_back(the_site);
    prev_chrom = std::move(the_site.chrom);
  }
  process_chrom_kth(require_same_region, regions, region_idx_process,
                    min_reads, max_dist, window_size, sites,
                    strand, the_neighbor, the_stats);
}


static void
load_regions(const string &regions_file, vector<genomic_interval> &regions) {
  // load genomic intervals (regions)
  std::ifstream in(regions_file);
  if (!in)
    throw std::runtime_error("bad regions file: " + regions_file);

  genomic_interval r;
  while (in >> r)
    regions.push_back(r);

  if (!is_sorted(begin(regions), end(regions)))
    throw std::runtime_error("regions file not sorted: " + regions_file);
}


int main(int argc, const char **argv) {

  try {

    string description = R"""(

This program computes statistics on the autocorrelation of methylation
levels. The input file must be in "counts" format from dnmtools. The
output file has one line for each distance between sites. Sites are
assumed to be CpG sites, but need not be. Each line of output has the
follwing values: distance, correlation, N, sdX, sdY, covXY. The value
of N is the number of observations contributing.

  )""";

    string outfile;
    size_t max_dist = 4000;
    size_t min_reads = 10;
    size_t min_sites = 500;
    size_t the_neighbor = 0; // 0=nearest
    size_t window_size = 0;
    int strand = 0; // code: same=1, different=-1 and any=0

    const double megabytes = 1024*1024;

    string regions_file;

    bool verbose = false;
    bool report_progress = false;

    bool require_same_region = false;

    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(strip_path(argv[0]), description, "<counts-file>");
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
    opt_parse.add_opt("progress", 'P', "report progress", false,
                      report_progress);
    opt_parse.add_opt("verbose", 'v', "print more run info",
                      false , verbose);
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

    if (verbose)
      cerr << "min_reads: " << min_reads << endl
           << "max_dist: " << max_dist << endl
           << "min_sites: " << min_sites << endl
           << "the_neighbor: " << the_neighbor << endl
           << "input_filename: " << input_filename << endl
           << "input_filesize: "
           << std::fixed << std::setprecision(2)
           << input_filesize/megabytes << "MB" << endl;

    vector<sum_stats> the_stats(max_dist + 1);

    vector<genomic_interval> regions;
    if (!regions_file.empty()) {
      if (verbose)
        cerr << "loading regions: " << regions_file << endl;
      load_regions(regions_file, regions);
      if (verbose)
        cerr << "total regions: " << regions.size() << endl;
    }

    if (the_neighbor == 0)
      process_sites_all_neighbors(report_progress,
                                  input_filename, require_same_region, regions,
                                  min_reads, max_dist, window_size, strand,
                                  the_stats);
    else
      process_sites_kth_neighbor(report_progress,
                                 input_filename, require_same_region, regions,
                                 min_reads, max_dist, window_size, strand,
                                 the_neighbor, the_stats);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    for (size_t i = 1; i <= max_dist; ++i) {
      double sdX = 0.0, sdY = 0.0, covXY = 0.0;
      const size_t N = the_stats[i].N;
      const double rXY = the_stats[i].get_values(sdX, sdY, covXY);
      if (N >= min_sites)
        out << i << '\t' << rXY << '\t' << N << '\t'
            << sdX << '\t' << sdY << '\t' << covXY << endl;
    }
  }
  catch (std::exception &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
