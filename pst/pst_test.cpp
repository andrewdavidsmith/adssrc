/*    pst_test: program for testing my probabilistic suffix tree
 *    implementation, which is in the files ProbSuffTree.hpp and
 *    ProbSuffTree.cpp
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

#include "ProbSuffTree.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;


int 
main(int argc, const char **argv) {

  try {
    
    string outfile;
    bool VERBOSE = false;

    size_t max_depth = 10;

    double pseudocount = 1.0;
    double critical_value = 0.01;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "test PST implementation",
			   "<fasta-file>");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("depth", 'd', "maximum depth", false, max_depth);
    opt_parse.add_opt("crit", 'c', "critical value", false, critical_value);
    opt_parse.add_opt("pseudo", 'p', "pseudocount", false, pseudocount);
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
    if (opt_parse.about_requested() || leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string sequences_file(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/
  
    vector<string> names, sequences;
    read_fasta_file(sequences_file, names, sequences);
    
    PST prob_suff_tree;
    
    for (size_t i = 0; i < sequences.size(); ++i)
      prob_suff_tree.insert(sequences[i], max_depth);
    
    prob_suff_tree.prune(pseudocount, critical_value);
    prob_suff_tree.convert_to_probabilities(pseudocount);
    
    cout << prob_suff_tree.tostring() << endl;
    
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
