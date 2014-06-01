/*    treetool: a program to manipulate phylogenetic trees in Newick
 *    format and to test parsing code
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
#include <fstream>
#include <algorithm>
#include <numeric>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;

class PhyloTreeNode {
public:
  PhyloTreeNode(const double bl) : branch_length(bl) {}
  PhyloTreeNode(const string &subtree_string);
  
  bool has_children() const {return !child.empty();}
  bool is_leaf() const {return child.empty();}
  
  string tostring() const;

private:
  vector<PhyloTreeNode> child;
  string name;
  double branch_length; // distance to parent
  
};

static bool
represents_leaf(const string &subtree_rep) {
  return true;
}

static string
extract_name(const string &subtree_rep) {
  return string();
}

static double
extract_branch_length(const string &subtree_rep) {
  return 0.0;
}

static bool
root_has_name(const string &subtree_rep) {
  return true;
}

static string
extract_left_subtree(const string &subtree_rep) {
  return string();
}

static string
extract_right_subtree(const string &subtree_rep) {
  return string();
}

PhyloTreeNode::PhyloTreeNode(const string &subtree_rep) {
  // This function needs to test the various ways that a string can be
  // passed in to represent a subtree at the current node
  
  branch_length = extract_branch_length(subtree_rep);
  if (root_has_name(subtree_rep))
    name = extract_name(subtree_rep);
  
  if (!represents_leaf(subtree_rep)) {
    const string left_subtree_rep = extract_left_subtree(subtree_rep);
    child.push_back(PhyloTreeNode(left_subtree_rep));
    const string right_subtree_rep = extract_right_subtree(subtree_rep);
    child.push_back(PhyloTreeNode(right_subtree_rep));
  }
}

int 
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;
    string outfile;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "manipulate Newick format "
			   "phylogenetic trees" "<newick-input>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
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
    const string newick_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/
    
    
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
