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

*/

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>

#include <cctype> // for isspace


#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"


using std::string;
using std::vector;
using std::endl;
using std::cerr;


static bool
check_balanced_parentheses(const string &s) {
  int count = 0;
  for (size_t i = 0; i < s.length() && count >= 0; ++i) {
    if (s[i] == '(') ++count;
    if (s[i] == ')') --count;
  }
  return count == 0;
}


class PhyloTreeNode {
public:
  PhyloTreeNode() {}
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
represents_leaf(const string &tree_rep) {
  return true;
}


static string
extract_name(const string &tree_rep) {
  return string();
}


static double
extract_branch_length(const string &tree_rep) {
  return 0.0;
}


static bool
root_has_name(const string &tree_rep) {
  return true;
}


static void
extract_subtrees(const string &tree_rep,
		 vector<string> &sub_subtree_reps) {
  return;
}


PhyloTreeNode::PhyloTreeNode(const string &tree_rep) {
  // This function needs to test the various ways that a string can be
  // passed in to represent a subtree at the current node
  
  branch_length = extract_branch_length(tree_rep);
  if (root_has_name(tree_rep))
    name = extract_name(tree_rep);

  if (!represents_leaf(tree_rep)) {
    vector<string> subtree_reps;
    extract_subtrees(tree_rep, subtree_reps);
    for (size_t i = 0; i < subtree_reps.size(); ++i)
      child.push_back(PhyloTreeNode(subtree_reps[i]));
  }
}


class PhyloTree {
public:
  PhyloTree(string tree_rep) {
    check_balanced_parentheses(tree_rep);
    // remove whitespace
    string::iterator w = std::remove_copy_if(tree_rep.begin(), tree_rep.end(),
					     tree_rep.begin(), &isspace);
    assert(w != tree_rep.begin());
    tree_rep.erase(--w, tree_rep.end()); // The "--" is for the ";"
    root = PhyloTreeNode(tree_rep);
  }
private:
  PhyloTreeNode root;
};


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
