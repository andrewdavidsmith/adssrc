/*    inverted-dups: a program for doing quality control relative to
 *    the inverted duplicates reads problem and masking "bad" reads if needed.
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
using std::unordered_map;
using std::ptr_fun;


// store each read from one end
struct FASTQRecord {
  string name;
  string seq;
  string seqtag;
  string score;

  string tostring() const {
    string tmp;
    std::ostringstream s;
    s << '@' << name << '\n' << seq << '\n' << tmp << '\n' << score;  
    return s.str();
  }
    
  void revcomp();
  
  // see if two reads from two ends match to each other
  // (they should have the same name)
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


// Reverse the sequence
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


// Read 4 lines one time from fastq and fill in the FASTQRecord structure
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
    
    if (!getline(s, r.seqtag))
      throw SMITHLABException("FASTQ file truncated expecting '+' line");
    
    if (r.seqtag.empty() || r.seqtag[0] != '+') 
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


// Compare two reads to detect the overlapped region
size_t
invdup_similarity(FASTQRecord &r1, FASTQRecord &r2) {
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

    string outfile;
    string statoutfile;
    string new_fqfile;
    double cutoff = 0.9;
    size_t to_ignore_at_end_of_name = 0;
    bool VERBOSE = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "count the invdup reads "
			   "in the input files", "<end1-fastq> <end2-fastq>");
    opt_parse.add_opt("output", 'o', "Name of the scanning results (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("statoutfile", 's', "Name of the output stats file (default: stdout)", 
		      false, statoutfile);
    opt_parse.add_opt("masking", 'm', "Name of the new second-end fastq file if you want to mask the invdup reads", 
		      false, new_fqfile);
    opt_parse.add_opt("cutoff", 'c', "The cutoff for invdup reads (default: 0.9)",
		      false, cutoff);
    opt_parse.add_opt("ignore", 'i', "Ignore this number of letters "
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

	 // Input: paired-end reads with end1 and end2
    std::ifstream f_read1(reads_file_one.c_str());
    if(!f_read1)
      throw SMITHLABException("cannot open input file " + reads_file_one);
    
    std::ifstream f_read2(reads_file_two.c_str());
    if(!f_read2)
      throw SMITHLABException("cannot open input file " + reads_file_two);

	 // Output scanning results: 
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());   

	 // Output the refined fastq files:
	 std::ofstream out_end_two(new_fqfile.c_str());
	 if(!new_fqfile.empty() && !out_end_two)
				throw SMITHLABException("cannot write new fastq file " + new_fqfile);

	 //------------SCAN THE READS------------//
	 size_t num_read = 0;
	 size_t num_bad_read = 0;
	 double sum_percent_overlap = 0;
	 double sum_bad_percent_overlap = 0;
    FASTQRecord end_one, end_two;

    while (f_read1 >> end_one && f_read2 >> end_two) {
      
		// Two reads should be in paired-ends
      if (!FASTQRecord::mates(to_ignore_at_end_of_name, end_one, end_two)) 
			throw SMITHLABException("expected mates, got:" + 
				end_one.tostring() + "\n" + end_two.tostring());
      
		// See if inverted duplicates emerge
      const size_t sim = invdup_similarity(end_one, end_two);
      const double percent_overlap = 
				static_cast<double>(sim)/end_one.seq.length();
      
		if (percent_overlap > cutoff) {
			num_bad_read++;
			sum_bad_percent_overlap += percent_overlap;
		}

		// Write new fastq file if -n 
		if (out_end_two) {
			const string masked_seq(end_two.seq.length(), 'N');
			out_end_two << '@'+end_two.seqtag.substr(1) << "\n";
			if (percent_overlap > cutoff) 
				out_end_two << masked_seq << "\n";
			else
				out_end_two << end_two.seq << "\n";
			out_end_two << end_two.seqtag << "\n"
							<< end_two.score << "\n";
		}

		num_read++;
		sum_percent_overlap += percent_overlap;

      out << sim << '\t' << percent_overlap << endl;
    }

    //------------WRITE STAT INFORMATION------------//
    std::ofstream of_stat;
    if (!statoutfile.empty()) of_stat.open(statoutfile.c_str());
    std::ostream statout(statoutfile.empty() ? cout.rdbuf() : of_stat.rdbuf());   
	 statout << "CUTOFF:\t" << cutoff << "\n"	
			  	<< "TOTAL READ PAIRS:\t" << num_read << "\n"	
			  	<< "SUSPECT INVERTED-DUPLICATED READ PAIRS:\t" << num_bad_read << "\n"	
			  	<< "PERCENTAGE OF GOOD READS:\t" << 1 - static_cast<double>(num_bad_read)/num_read << "\n"	
			  	<< "MEAN OVERLAP PERCENTAGE:\t" << sum_percent_overlap/num_read << "\n"	
			  	<< "MEAN OVERLAP PERCENTAGE OF INVERTED-DUPLICATES:\t" << sum_bad_percent_overlap/num_bad_read << "\n";	
  
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
