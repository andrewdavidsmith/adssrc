/* smoothmeth: program to smooth methylation levels at CpG sites using
 * an HMM
 *
 * Copyright (C) 2014 Andrew D Smith
 *
 * Author: Andrew D. Smith
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
 * License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 */

#include <numeric>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "TwoStateHMM.hpp"
#include "MethpipeFiles.hpp"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;
using std::numeric_limits;
using std::max;
using std::min;
using std::pair;


template <class T, class U> static void
separate_regions(const bool VERBOSE, const size_t desert_size, 
		 vector<SimpleGenomicRegion> &cpgs,
		 vector<T> &meth, vector<U> &reads,
		 vector<size_t> &reset_points) {
  if (VERBOSE)
    cerr << "[SEPARATING BY CPG DESERT]" << endl;
  // eliminate the zero-read cpgs
  size_t j = 0;
  for (size_t i = 0; i < cpgs.size(); ++i)
    if (reads[i] > 0) {
      cpgs[j] = cpgs[i];
      meth[j] = meth[i];
      reads[j] = reads[i];
      ++j;
    }
  cpgs.erase(cpgs.begin() + j, cpgs.end());
  meth.erase(meth.begin() + j, meth.end());
  reads.erase(reads.begin() + j, reads.end());
  
  // segregate cpgs
  size_t prev_cpg = 0;
  for (size_t i = 0; i < cpgs.size(); ++i) {
    const size_t dist = (i > 0 && cpgs[i].same_chrom(cpgs[i - 1])) ? 
      cpgs[i].get_start() - prev_cpg : numeric_limits<size_t>::max();
    if (dist > desert_size)
      reset_points.push_back(i);
    prev_cpg = cpgs[i].get_start();
  }
  reset_points.push_back(cpgs.size());
  if (VERBOSE)
    cerr << "CPGS RETAINED: " << cpgs.size() << endl
	 << "DESERTS REMOVED: " << reset_points.size() - 2 << endl << endl;
}


static void
write_estimated_cpg_meth(const double high_meth, const double low_meth,
			 const vector<SimpleGenomicRegion> &cpgs,
			 const vector<pair<double, double> > &meth, 
			 const vector<size_t> &reads,
			 const vector<double> &scores,
			 std::ostream &out) {
  
  size_t j = 0;
  for (size_t i = 0; i < cpgs.size(); ++i) {
    out << cpgs[i].get_chrom() << '\t'
	<< cpgs[i].get_start() << '\t'
	<< "+\tCpG" << '\t';
    if (reads[i] == 0)
      out << meth[i].first/std::max(1ul, reads[i]) << '\t';
    else {
      out << high_meth*(1.0 - scores[j]) + low_meth*scores[j] << '\t';
      ++j;
    }
    out << reads[i] << endl;
  }
}


static void
load_cpgs(const bool VERBOSE, 
	  const string& cpgs_file, vector<SimpleGenomicRegion> &cpgs,
	  vector<pair<double, double> > &meth, vector<size_t> &reads) {
  if (VERBOSE)
    cerr << "[READING CPGS AND METH PROPS]" << endl;
  vector<GenomicRegion> cpgs_in;
  ReadBEDFile(cpgs_file, cpgs_in);
  if (!check_sorted(cpgs_in))
    throw SMITHLABException("CpGs not sorted in file \"" + cpgs_file + "\"");

  for (size_t i = 0; i < cpgs_in.size(); ++i) {
    cpgs.push_back(SimpleGenomicRegion(cpgs_in[i]));
    meth.push_back(std::make_pair(cpgs_in[i].get_score(), 0.0));
    const string r(cpgs_in[i].get_name());
    reads.push_back(atoi(r.substr(r.find_first_of(":") + 1).c_str()));
    meth.back().first = int(meth.back().first * reads.back());
    meth.back().second = int(reads.back() - meth.back().first);
  }
  if (VERBOSE)
    cerr << "TOTAL CPGS: " << cpgs.size() << endl
	 << "MEAN COVERAGE: " 
	 << accumulate(reads.begin(), reads.end(), 0.0)/reads.size() 
	 << endl << endl;
}



int
main(int argc, const char **argv) {

  try {

    string outfile;
    
    size_t desert_size = 1000;
    size_t max_iterations = 10;
    
    // run mode flags
    bool VERBOSE = false;
    
    // corrections for small values (not parameters):
    double tolerance = 1e-10;
    double min_prob  = 1e-10;

    string params_in_file;
    string params_out_file;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "use an HMM to smooth methylation"
			   " levels at CpG sites", "<cpg-meth-file>");
    opt_parse.add_opt("out", 'o', "output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("desert", 'd', "max dist btwn cpgs with reads in HMR", 
		      false, desert_size);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_iterations); 
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string cpgs_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/
    
    // separate the regions by chrom and by desert
    vector<SimpleGenomicRegion> cpgs;
    // vector<double> meth;
    vector<pair<double, double> > meth;
    vector<size_t> reads;
    if (methpipe::is_methpipe_file_single(cpgs_file)) {
      if (VERBOSE)
        cerr << "[READING CPGS AND METH PROPS]" << endl;
      methpipe::load_cpgs(cpgs_file, cpgs, meth, reads);
      if (VERBOSE)
        cerr << "TOTAL CPGS: " << cpgs.size() << endl
             << "MEAN COVERAGE: " 
             << accumulate(reads.begin(), reads.end(), 0.0)/reads.size() 
             << endl << endl;
    }
    else
      load_cpgs(VERBOSE, cpgs_file, cpgs, meth, reads);

    // separate the regions by chrom and by desert, and eliminate
    // those isolated CpGs
    vector<size_t> reset_points;
    separate_regions(VERBOSE, desert_size, cpgs, meth, reads, reset_points);
    
    vector<double> start_trans(2, 0.5), end_trans(2, 1e-10);
    vector<vector<double> > trans(2, vector<double>(2, 0.25));
    trans[0][0] = trans[1][1] = 0.75;
    
    const TwoStateHMMB hmm(min_prob, tolerance, max_iterations, VERBOSE);
    
    double fg_alpha = 0;
    double fg_beta = 0;
    double bg_alpha = 0;
    double bg_beta = 0;
    
    const double n_reads = 
      accumulate(reads.begin(), reads.end(), 0.0)/reads.size();
    fg_alpha = 0.33*n_reads;
    fg_beta = 0.67*n_reads;
    bg_alpha = 0.67*n_reads;
    bg_beta = 0.33*n_reads;
    
    assert(max_iterations > 0);
    
    hmm.BaumWelchTraining(meth, reset_points, start_trans, trans, 
			  end_trans, fg_alpha, fg_beta, bg_alpha, bg_beta);
    
    /******************************************************************
     * DECODE THE DOMAINS AND IN THE PROCESS GET THE POSTERIORS
     */
    vector<bool> classes;
    vector<double> scores;
    hmm.PosteriorDecoding(meth, reset_points, start_trans, trans, 
			  end_trans, fg_alpha, fg_beta, bg_alpha, 
			  bg_beta, classes, scores);
    
    cpgs.clear();
    meth.clear();
    reads.clear();
    if (VERBOSE)
      cerr << "[RE-READING CPGS AND METH PROPS]" << endl;
    if (methpipe::is_methpipe_file_single(cpgs_file))
      methpipe::load_cpgs(cpgs_file, cpgs, meth, reads);
    else load_cpgs(VERBOSE, cpgs_file, cpgs, meth, reads);
    
    const double high_meth = bg_alpha/(bg_alpha + bg_beta);
    const double low_meth = fg_alpha/(fg_alpha + fg_beta);
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    
    write_estimated_cpg_meth(high_meth, low_meth, cpgs, meth, reads, scores, out);
    
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
