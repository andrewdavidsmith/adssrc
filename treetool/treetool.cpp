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


/*

  ABOUT NEWICK FORMAT
  -------------------

  The following should all be legal, according to wikipedia.

  (,,(,));                               no nodes are named
  (A,B,(C,D));                           leaf nodes are named
  (A,B,(C,D)E)F;                         all nodes are named
  (:0.1,:0.2,(:0.3,:0.4):0.5);           all but root node have a distance to parent
  (:0.1,:0.2,(:0.3,:0.4):0.5):0.0;       all have a distance to parent
  (A:0.1,B:0.2,(C:0.3,D:0.4):0.5);       distances and leaf names (popular)
  (A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;     distances and all names
  ((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;    a tree rooted on a leaf node (rare)
  
  ========================================================================

  (1) It seems like a tree always ends in a semicolon
  (2) Trees need not be binary, but usually are ***
  (3) Since the code needs to have uniform nodes, we will always have
      a distance, even for the root, but this will be 0.0
  (4) The parsing will look for commas, but not within parentheses
  (5) No whitespace will be allowed. It will be removed initially
  (6) There is a question of whether to allow commas or colons inside
      names, for example by using quotes around the names or trying to
      parse intelligently
  (7) Leaves should have proper names as unique identifiers, if not given in the 
      newick string, we should name them properly. 
  (8) Search for nearest common ancestor will return the ancestor's name, for now. 

*/

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <bitset>
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


int 
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;
    string outfile;
    string label_to_check;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "manipulate Newick format "
			   "phylogenetic trees",
			   "<newick-input>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("label", 'l', "check if this label exists", 
		      false, label_to_check);
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
      throw SMITHLABException("bad file: " + newick_file);
    
    PhyloTree t;
    in >> t;
    
    cout << t.tostring() << endl;
    cout << t << endl;

    //Name unnamed nodes
    size_t count = 0;
    t.fill_leaf_names("Leaf", count);
    count = 0;
    t.fill_names("Internal", count);  
 
    //get common ancestor of 2 random selected leaves. 
    string ancestor;
    vector<string> leaf_names;
    t.get_leaf_names(leaf_names); 

    srand(time(0));
    std::random_shuffle(leaf_names.begin(), leaf_names.end());
    vector<string> tmp(leaf_names.begin(), leaf_names.begin()+2);

    cerr << "the common ancestor of "  ;
    copy(tmp.begin(), tmp.end(), std::ostream_iterator<string>(cerr, ","));
    
    t.find_common_ancestor(tmp, ancestor);
    cerr << (ancestor.empty() ? "not found" : ancestor) << endl; 
    
    // print subtree rooted at ancestor
    cout << t.tostring(ancestor) << endl;
    cout << t.Newick_format(ancestor) << endl;
    
    // //Get neighbor relations
    // t.build_neighbor();
    // vector<vector<size_t> > neighbor(t.get_neighbor());
    // for(size_t i=0; i < neighbor.size(); ++i){
    //   cerr << "pattern" << i << " has neighbors" ;
    //   for(size_t j = 0; j < neighbor[i].size(); ++j)
    // 	cerr << neighbor[i][j] <<",";
    //   cerr << endl;
    // }

    if (!label_to_check.empty())
      cout << label_to_check << " "
	   << (t.label_exists(label_to_check) ? 
	       "exists" : "does not exist") << endl;
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
