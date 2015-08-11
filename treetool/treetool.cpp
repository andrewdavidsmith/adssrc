/*    treetool: a program to manipulate phylogenetic trees in Newick
 *    format and to test parsing code
 *
 *    Copyright (C) 2014 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith and Jenny Qu
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
#include <algorithm>
#include <numeric>
#include <sstream>
#include <tr1/unordered_set>

#include <cctype> // for isspace

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "PhyloTree.hpp"


using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::tr1::unordered_set;
using std::istream_iterator;



int 
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;
    string outfile;
    string species_file;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), 
			   "manipulate Newick format phylogenetic trees",
			   "<newick-input>");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("species", 's', "species defining subtree to extract",
		      false, species_file);
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
    const string newick_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/
    
    std::ifstream in(newick_file.c_str());
    if (!in)
      throw SMITHLABException("problem reading file: " + newick_file);
    
    string nw;
    getline(in, nw);

    PhyloTree t(nw);
    cout << "tree format valid" << endl;
    cout << "number of nodes: " << t.get_size() << endl;
    
    if (!species_file.empty()) {
      std::ifstream spec_in(species_file.c_str());
      unordered_set<string> species((istream_iterator<string>(spec_in)), 
				    istream_iterator<string>());

      PhyloTree u;
      PhyloTree::copy_subtree_with_species(t, species, u);
      
      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
      out << u << endl;
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
