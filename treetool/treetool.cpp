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
    
    string nw;
    in >> nw;
    PhyloTree t(nw);
    size_t namesuf = 1;
    t.fill_names("INT",namesuf);
    
    /***** get common ancestor of 2 random selected leaves. ****/
    // string ancestor;
    // vector<string> leaf_names;
    // t.get_leaf_names(leaf_names); 

    // srand(time(0));
    // std::random_shuffle(leaf_names.begin(), leaf_names.end());
    // vector<string> tmp;
    // for(size_t i =0; i < 2; ++i)
    //   tmp.push_back(leaf_names[i]);

    //  cerr << "the common ancestor of "  ;
    //  copy(tmp.begin(), tmp.end(), 
    //      std::ostream_iterator<string>(cerr, ","));
    
    // t.find_common_ancestor(tmp, ancestor);
    // cerr << "is " <<  (ancestor.empty() ? "not found" : ancestor) << endl; 
    
    // // print subtree rooted at ancestor
    // cout << t.Newick_format(ancestor) << endl;
    
    /**** trim tree *******/
    // t.trim_to_keep(tmp);
    // cout << "After trimming:" << endl 
    // 	 << t.Newick_format() << endl;

    /**** combine tree *****/
    // string s1 = t.Newick_format(tmp[0]);
    // string s2 = t.Newick_format(tmp[1]);
    // string sc = combine_newick(s1, s2, "root", 0.1, 0.2, 0.3);
    // cerr << "combining "<< s1 << " and " << s2 << endl << sc << endl;
    // PhyloTree t2(sc);
    // cerr << t2.get_root_branch() << endl;

    // t.set_branch(leaf_names[2], 3000);
    // cerr << t.Newick_format()<< endl;

    vector<size_t> pa_idx;
    t.get_node_parent_idx(pa_idx);

    vector<string> node_names;
    t.get_node_names(node_names);
    vector<vector<size_t> > child_idx;
    t.get_node_child_idx(child_idx);
    vector<size_t> heights;
    t.get_all_heights(heights);

    for(size_t i = 0; i < node_names.size(); ++i){
      cerr << "node " << node_names[i] << " has height " 
	   << heights[i] << "\tparent\t";
      if( pa_idx[i]==node_names.size())
	cerr << "NA;\t";
      else 
	cerr << node_names[pa_idx[i]] << ";\t";

      cerr << " and children \t";
      if(child_idx[i].size()==0)
	cerr << "NA" << endl;
      else{
	for(size_t j = 0; j < child_idx[i].size(); ++j )
	  cerr << node_names[child_idx[i][j]] << "\t";
	cerr << endl;
      }
    }


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


