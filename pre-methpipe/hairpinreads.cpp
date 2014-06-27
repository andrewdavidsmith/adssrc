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

// #include "bsutils.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::accumulate;
using std::tr1::unordered_map;



struct FASTQRecord {
  string name;
  string seq;
  string score;

  string tostring() const {
    std::ostringstream s;
    s << '@' << name << '\n' << seq << "\n#\n" << score;
    return s.str();
  }
  
  void revcomp();
  
  static bool 
  mates(const size_t to_ignore_at_end, // in case names have #0/1 name ends
	const FASTQRecord &a, const FASTQRecord &b) {
    const bool same_name = 
      equal(a.name.begin(), a.name.begin() + a.name.length() - to_ignore_at_end,
	    b.name.begin());
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




int 
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;
    size_t to_ignore_at_end_of_name = 0;
    string outfile;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "count the hairpin reads "
			   "in the input files", "<end1-fastq> <end2-fastq>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false, outfile);
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
    const string end_one_file = leftover_args.front();
    const string end_two_file = leftover_args.back();
    /****************** END COMMAND LINE OPTIONS *****************/

    std::ifstream in1(end_one_file.c_str());
    if (!in1)
      throw SMITHLABException("bad file: " + end_one_file);
    
    std::ifstream in2(end_two_file.c_str());
    if (!in2)
      throw SMITHLABException("bad file: " + end_two_file);
    

    FASTQRecord end_one, end_two;
    while (in1 >> end_one && in2 >> end_two) {
      
      if (!FASTQRecord::mates(to_ignore_at_end_of_name, end_one, end_two)) 
	throw SMITHLABException("expected mates, got:" + 
				end_one.tostring() + "\n" +
				end_two.tostring());
      
    }

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
    
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
