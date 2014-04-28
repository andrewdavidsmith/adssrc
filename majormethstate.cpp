/*    majormethstate: alter epireads such that if more states within
 *    the read showed unmethylated, then the states are inverted. 
 *
 *    The purpose of this program is to help understand the complexity
 *    of states independent of whether the methylation is high or low.
 *
 *    Copyright (C) 2014 Andrew D. Smith
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
#include <fstream>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "bsutils.hpp"

using std::string;
using std::cout;
using std::cerr;
using std::endl;

static void
major_state_to_meth(string &seq) {
  
  size_t total = 0, unmeth = 0;
  
  for (std::string::const_iterator i(seq.begin()); i != seq.end(); ++i)
    if (is_cytosine(*i))
      ++total;
    else if (is_thymine(*i)) {
      ++unmeth;
      ++total;
    }
  
  if (static_cast<double>(unmeth)/total > 0.5)
    for (std::string::iterator i(seq.begin()); i != seq.end(); ++i) {
      if (is_cytosine(*i))
	*i = 'T';
      else if (is_thymine(*i))
	*i = 'C';
    }
}


int 
main(int argc, const char **argv) {
  
  try {
    
    // bool VERBOSE = false;
    string outfile;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), 
			   "inverts states in epireads so meth is major",
			   "<epireads>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false, outfile);
    // opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    std::vector<string> leftover_args;
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
    const string epireads_file(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/
    
    std::ifstream in(epireads_file.c_str());
    if (!in) 
      throw SMITHLABException("cannot open input file " + epireads_file);
    
    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
    
    string line;
    while (getline(in, line)) {
      string chrom, seq;
      size_t pos = 0;
      std::istringstream iss(line);
      iss >> chrom >> pos >> seq;
      if (!iss)
	throw SMITHLABException("malformed line: \"" + line + "\"");
      
      major_state_to_meth(seq);
      out << chrom << '\t' << pos << '\t' << seq << endl;
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
