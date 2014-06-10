/*    allelicmeth2:
 *
 *    Copyright (C) 2014 University of Southern California and
 *                       Andrew D. Smith
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

#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <utility>

#include <tr1/cmath>

#include <gsl/gsl_sf_gamma.h>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "MethpipeFiles.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::istringstream;
using std::tr1::unordered_map;
using std::max;
using std::min;


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///////
/////// CODE FOR WORKING WITH EPIREADS BELOW HERE
///////

struct Epiread {
  bool operator<(const Epiread &other) const {
    return (chr < other.chr || (chr == other.chr && pos < other.pos));
  }
  size_t length() const {return seq.length();}
  
  string chr;
  size_t pos;
  string seq;
};


static std::istream&
operator>>(std::istream &in, Epiread &er) {
  string buffer;
  if (getline(in, buffer)) {
    std::istringstream is(buffer);
    if (!(is >> er.chr >> er.pos >> er.seq))
      throw SMITHLABException("malformed epiread line:\n" + buffer);
  }
  return in;
}


std::ostream&
operator<<(std::ostream &out, const Epiread &er) {
  return out << er.chr << '\t' << er.pos << '\t' << er.seq;
}


static inline double
log_sum_log(const double p, const double q) {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}


// p(k) =  C(n1, k) C(n2, t - k) / C(n1 + n2, t)
static double
log_hyper_g(const size_t k, const size_t n1, const size_t n2, const size_t t) {
  return (gsl_sf_lnfact(n1) - gsl_sf_lnfact(k) - gsl_sf_lnfact(n1 - k) +
  	  gsl_sf_lnfact(n2) - gsl_sf_lnfact(t - k) - gsl_sf_lnfact(n2 - (t - k)) -
  	  (gsl_sf_lnfact(n1 + n2) - gsl_sf_lnfact(t) - gsl_sf_lnfact(n1 + n2 - t)));
}


static double
fishers_exact(size_t a, size_t b, size_t c, size_t d) {
  const size_t m = a + c; // sum of first column
  const size_t n = b + d; // sum of second column
  const size_t k = a + b; // sum of first row
  const double observed = log_hyper_g(a, m, n, k);
  double p = 0.0;
  for (size_t i = (n > k ? 0ul : k - n); i <= std::min(k, m); ++i) {
    const double curr = log_hyper_g(i, m, n, k);
    if (curr <= observed)
      p = log_sum_log(p, curr);
  }
  return exp(p);
}


static size_t
state_pair_to_index(const string &s, const size_t idx) {
  assert(idx < s.length() - 1);
  const char a = s[idx];
  if (a == 'C') {
    const char b = s[idx+1];
    if (b == 'C') return 0;
    if (b == 'T') return 1;
    return 4;
  }
  if (a == 'T') {
    const char b = s[idx+1];
    if (b == 'C') return 2;
    if (b == 'T') return 3;
    return 4;
  }
  return 4;
}


template <class T> 
struct PairStateCounter {
  T CC;
  T CT;
  T TC;
  T TT;
  
  double score() const {
    return (CC*TT > CT*TC) ?
      fishers_exact(CC, CT, TC, TT) : fishers_exact(CT, CC, TT, TC);
  }
  double total() const {return CC + CT + TC + TT;}
  
  string tostring() const {
    return toa(CC) + '\t' + toa(CT) + '\t' + toa(TC) + '\t' + toa(TT);
  }

  void increment(const size_t state) {
    if (state == 0) ++CC;
    else if (state == 1) ++CT;
    else if (state == 2) ++TC;
    else if (state == 3) ++TT;
  }
};


template <class T> void
fit_states(const Epiread &er, vector<PairStateCounter<T> > &counts) {
  for (size_t i = 0; i < er.length() - 1; ++i) {
    const size_t pos = er.pos + i;
    assert(pos < counts.size());
    const size_t curr_state = state_pair_to_index(er.seq, i);
    counts[pos].increment(curr_state);
  }
}


static void
get_chrom_sizes(const string &epi_file, 
		unordered_map<string, size_t> &chrom_sizes) {
  
  std::ifstream in(epi_file.c_str());
  if (!in)
    throw SMITHLABException("cannot open input file: " + epi_file);
  
  string chrom;
  Epiread er;
  while (in >> er) {
    if (chrom_sizes.find(er.chr) == chrom_sizes.end())
      chrom_sizes[er.chr] = 0;
    chrom_sizes[er.chr] = std::max(chrom_sizes[er.chr], er.pos + er.length());
  }
}


int 
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;
    string outfile;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "does like methcounts "
                           "except with epireads", "<epireads>");
    opt_parse.add_opt("output", 'o', "output file name (default: stdout)", 
		      false, outfile);
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
    const string epi_file(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/
    
    unordered_map<string, size_t> chrom_sizes;
    get_chrom_sizes(epi_file, chrom_sizes);

    if (VERBOSE)
      cerr << "CHROMS: " << chrom_sizes.size() << endl;
    
    std::ifstream in(epi_file.c_str());
    if (!in)
      throw SMITHLABException("cannot open input file: " + epi_file);
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
    
    vector<PairStateCounter<unsigned short> > counts;
    string chrom;
    Epiread er;
    while (in >> er) {
      if (er.chr != chrom) {
    	if (!chrom.empty())
    	  for (size_t i = 0; i < counts.size(); ++i)
            // this output format assumes it's all CpGs
    	    out << chrom << '\t' << i << "\t+\tCpG\t"
    		<< counts[i].score() << '\t' 
    		<< counts[i].total() << '\t' 
    		<< counts[i].tostring() << endl;
    	counts = vector<PairStateCounter<unsigned short> >(chrom_sizes[er.chr]);
      }
      fit_states(er, counts);
      chrom.swap(er.chr);
    }
    if (!chrom.empty())
      for (size_t i = 0; i < counts.size(); ++i)
    	// this output format assumes it's all CpGs
    	out << chrom << '\t' << i << "\t+\tCpG\t"
    	    << counts[i].score() << '\t' 
    	    << counts[i].total() << '\t' 
    	    << counts[i].tostring() << endl;
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
