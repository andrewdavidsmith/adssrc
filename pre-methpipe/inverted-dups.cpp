/*    hairpinreads: a program for doing quality control relative to
 *    the "hairpin read" problem.
 *
 *    Copyright (C) 2014 University of Southern California and
 *                       Andrew D. Smith, Liz Ji and Jenny Qu
 *
 *    Authors: Andrew D. Smith, Liz Ji and Jenny Qu
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
#include <fstream>
#include <numeric>
#include <sstream>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::min;
using std::accumulate;
using std::tr1::unordered_map;
using std::ptr_fun;


struct FASTQRecord {
  string name;
  string seq;
  string score;

  string tostring() const {
    string tmp;
    std::ostringstream s;
    s << '@' << name << '\n' << seq << '\n' << tmp << '\n' << score;  
    return s.str();
  }
    
  void revcomp();
  
  static bool 
  mates(const size_t to_ignore_at_end, // in case names have #0/1 name ends
	const FASTQRecord &a, const FASTQRecord &b) {
    const string namefield1 = a.name.substr(0, a.name.find_first_of(' '));
    const string namefield2 = b.name.substr(0, b.name.find_first_of(' '));

    const bool same_name = 
      equal(namefield1.begin(), namefield1.begin() + namefield1.length() - to_ignore_at_end,
	    namefield2.begin());
    return same_name && a.seq.length() == b.seq.length();
  }
  
};


void
FASTQRecord::revcomp() {
  string tmp(seq);
  reverse(tmp.begin(), tmp.end());
  revcomp_inplace(seq);
  for (size_t i = 0; i < seq.length(); ++i)
    if (tmp[i] == '_')
      seq[i] = '_';
  reverse(score.begin(), score.end());
}



std::ostream& 
operator<<(std::ostream& s, const FASTQRecord &r) {
  return s << r.tostring();
}



std::istream& 
operator>>(std::istream& s, FASTQRecord &r) {
  if (getline(s, r.name)) {
    
    if (r.name.empty() || r.name[0] != '@') 
      throw SMITHLABException("FASTQ file out of sync at '@'");
    
    size_t j = 0;
    for (size_t i = 1; i < r.name.length() && r.name[i] != ' '
	   && r.name[i] != '\t'; ++i)
      r.name[j++] = r.name[i];
    r.name.erase(r.name.begin() + j, r.name.end());
    
    if (!getline(s, r.seq))
      throw SMITHLABException("FASTQ file truncated expecting seq");
    
    if (!getline(s, r.score))
      throw SMITHLABException("FASTQ file truncated expecting '+' line");
    
    if (r.score.empty() || r.score[0] != '+') 
      throw SMITHLABException("FASTQ file out of sync [missing '+']");
    
    if (!getline(s, r.score))
      throw SMITHLABException("FASTQ file truncated expecting score");
  }
  else s.setstate(std::ios::badbit);

  return s;
}

////////////////////////////////////////////////////////////////////////////////

static bool
similar_letters_bisulfite(const char a, const char b) {
  return (a == b) || (a == 'T' && b == 'C') || (a == 'G' && b == 'A');
}


size_t
hairpin_similarity(FASTQRecord &r1, FASTQRecord &r2) {
  size_t sim = 0;
  string::const_iterator it1(r1.seq.begin());
  string::const_iterator it2(r2.seq.begin());
  while (it1 < r1.seq.end() && it2 < r2.seq.end())
    sim += similar_letters_bisulfite(*it1++, *it2++);
  return sim;
}


int 
main(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    size_t to_ignore_at_end_of_name = 0;
    string outfile;
    double cutoff = 0;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "count the hairpin reads "
			   "in the input files", "<end1-fastq> <end2-fastq>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("cutoff", 'c', "The cutoff for hairpin reads",
		      false, cutoff);
    opt_parse.add_opt("ignore", 'i', "ignore this number of letters "
		      "at end of name", false, to_ignore_at_end_of_name);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
	   << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested() || leftover_args.size() != 2) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    const string reads_file_one = leftover_args.front();
    const string reads_file_two = leftover_args.back();
    /****************** END COMMAND LINE OPTIONS *****************/

    std::ifstream f_read1(reads_file_one.c_str());
    if(!f_read1)
      throw SMITHLABException("cannot open input file " + reads_file_one);
    
    std::ifstream f_read2(reads_file_two.c_str());
    if(!f_read2)
      throw SMITHLABException("cannot open input file " + reads_file_two);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());   

    FASTQRecord end_one, end_two;
    while (f_read1 >> end_one && f_read2 >> end_two) {
      
      if (!FASTQRecord::mates(to_ignore_at_end_of_name, end_one, end_two)) 
	throw SMITHLABException("expected mates, got:" + 
				end_one.tostring() + "\n" +
				end_two.tostring());
      
      const size_t sim = hairpin_similarity(end_one, end_two);
      const double percent_overlap = 
	static_cast<double>(sim)/end_one.seq.length();
      
      out << sim << '\t' << percent_overlap << endl;
    }
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
