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

static string 
combine_newick(const string s1, const string s2, 
	       const string rootname,
	       const double branch1, const double branch2, 
	       const double branch0){
  string s1_copy, s2_copy;

  string s = s1; 
  size_t found2 = s.find_last_of(":");
  size_t found3 = s.find_last_of(";");  
  if(found2 > found3){
    s1_copy.assign(s.substr(0, found3));
  }else{
    s1_copy.assign(s.substr(0, found2));
  }
  s.assign(s2); 
  found2 = s.find_last_of(":");
  found3 = s.find_last_of(";");  
  if(found2 > found3){
    s2_copy.assign(s.substr(0, found3));
  }else{
    s2_copy.assign(s.substr(0, found2));
  }

  std::ostringstream oss;
  oss << "(" <<  s1_copy << ":" << branch1 << "," << 
    s2_copy <<  ":" << branch2  << ")" <<
    rootname << ":" << branch0  << ";" ;
  string newick = oss.str();

  return newick; 
}

int 
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;
    string outfile;
    string label_to_check;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), 
			   "manipulate Newick format "
			   "phylogenetic trees",
			   "<newick-input>");
    opt_parse.add_opt("output", 'o', 
		      "Name of output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("label", 'l', "check if this label exists", 
		      false, label_to_check);
    opt_parse.add_opt("verbose", 'v', "print more run info", 
		      false, VERBOSE);
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
    
    string tree_rep;
    in >> tree_rep;
    PhyloTree t(tree_rep);
    
    cout << t.tostring() << endl;
    cout << t << endl;
   
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
    cout << t.tostring() << endl;
    cout << t.Newick_format(ancestor) << endl;
    
    //trim tree
    t.trim_to_keep(tmp);
    cout << "After trimming:" << endl 
	 << t.Newick_format() << endl;

    //combine tree
    string s1 = t.Newick_format(tmp[0]);
    string s2 = t.Newick_format(tmp[1]);
    string sc = combine_newick(s1, s2, "root", 0.1, 0.2, 0.3);
    cerr << "combining "<< s1 << " and " << s2 << endl << sc << endl;
    PhyloTree t2(sc);
    cerr << t2.get_root_branch() << endl;


    //get branch lengths
    vector<double> branches;
    t.get_branches(branches);
    cerr << "branches are:" << endl;
    for(size_t i =0; i < branches.size(); ++i){
      cerr << branches[i] << "\t";
    }
    cerr << endl;

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
