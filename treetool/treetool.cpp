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

#include <cctype> // for isspace

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::cout;


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
  PhyloTreeNode(const string &subtree_string);
  
  bool has_children() const {return !child.empty();}
  bool is_leaf() const {return child.empty();}
  
  bool label_exists(const string &label) const;
  string tostring(const size_t depth = 0) const;
  string tostring(const string &label) const;
  string Newick_format() const;
  string Newick_format(const string & label) const;

  void fill_leaf_names(const string prefix, size_t &count);
  void fill_names(const string prefix, size_t &count);
  
  void get_leaf_names(vector<string> &leaf_names);
  bool find_common_ancestor(const vector<string> &leaf_names, 
			    size_t &count, string &ancestor);

private:
  vector<PhyloTreeNode> child;
  string name;
  double branch_length; // distance to parent
};



bool
PhyloTreeNode::label_exists(const string &label) const {
  if (name == label) return true;
  else {
    for (size_t i = 0; i < child.size(); ++i)
      if (child[i].label_exists(label))
	return true;
    return false;
  }
}


string
PhyloTreeNode::tostring(const size_t depth) const {
  std::ostringstream oss;
  oss << string(depth, '\t') << branch_length << ':' << name;
  for (size_t i = 0; i < child.size(); ++i)
    oss << endl << child[i].tostring(depth + 1);
  return oss.str();
}


string
PhyloTreeNode::tostring(const string &label)const{
  assert(label_exists(label));
  if(label == name)
    return( tostring());
  
  string stringrep;
  for(size_t i = 0; i < child.size(); ++i){
    if(child[i].label_exists(label))
      stringrep = child[i].tostring(label);
  }
  return stringrep;
}


string
PhyloTreeNode::Newick_format() const {
  std::ostringstream oss;
  if (!child.empty()) {
    oss << '(';
    oss << child.front().Newick_format();
    for (size_t i = 1; i < child.size(); ++i)
      oss << ',' << child[i].Newick_format();
    oss << ')';
  }
  oss << name << ':' << branch_length;
  return oss.str();
}

string
PhyloTreeNode::Newick_format(const string &label) const {
  string nf;
  if (name == label){
    nf = Newick_format();
  }  else {
    for (size_t i=0; i < child.size(); ++i)
      if(child[i].label_exists(label))
	nf = child[i].Newick_format(label);
  } 
  return nf;
}

void 
PhyloTreeNode::fill_leaf_names(const string prefix, size_t &count)  {
  if(is_leaf()){
    if(name.length()==0){
      std::stringstream ss; ss << prefix << count;
      name.assign( ss.str());
      count ++; 
    }
  }else {
    for(size_t i = 0; i < child.size(); ++i)
      child[i].fill_leaf_names(prefix, count);
  }
}


void 
PhyloTreeNode::fill_names(const string prefix, size_t &count)  {
  if(name.empty()){
    std::stringstream ss; ss << prefix << count;
    name.assign( ss.str());
    count ++; 
  } 
  if(!is_leaf()) {
    for(size_t i = 0; i < child.size(); ++i)
      child[i].fill_names(prefix, count);
  }
}

void 
PhyloTreeNode::get_leaf_names(vector<string> &leaf_names){
  if(is_leaf()){
    assert(name.length()>0);
    leaf_names.push_back(name);
  }
  else{
    for(size_t i = 0; i < child.size(); ++i)
      child[i].get_leaf_names(leaf_names);
  }
}

bool 
PhyloTreeNode::find_common_ancestor(const vector<string> &leaf_names, 
				    size_t &count, string &ancestor){
  bool found=false;
  if(count == leaf_names.size()){
    assert(!ancestor.empty());
    found = true; 
  }else{
    if(is_leaf()){
      vector<string>::const_iterator it;
      it = find(leaf_names.begin(),leaf_names.end(), name);
      if( it != leaf_names.end()){
	count = 1; 
      }
    } else{
      vector<size_t> counts(child.size(), 0);
      for(size_t i =0; i< child.size(); ++i)
	found = child[i].find_common_ancestor(leaf_names, counts[i], ancestor);
      count = std::accumulate(counts.begin(), counts.end(), 0); 
    }
    if(count == leaf_names.size()) found = true;
    if(found && ancestor.empty())
      ancestor = name; 
  }
  return found;
}



static bool
represents_leaf(const string &tree_rep) {
  return (tree_rep.find_first_of(',') == string::npos);
}


static string
extract_name(const string &tree_rep) {
  const size_t final_parenthesis = tree_rep.find_last_of(")");
  const size_t possible_name_start = 
    (final_parenthesis == string::npos) ? 0 : final_parenthesis + 1;
  const size_t branch_len_start = 
    std::min(tree_rep.find_first_of(":", possible_name_start), 
	     tree_rep.length());
  return string(tree_rep.begin() + possible_name_start,
		tree_rep.begin() + branch_len_start);
}


static double
extract_branch_length(const string &tree_rep) {
  const size_t final_parenthesis = tree_rep.find_last_of(")");
  const size_t colon_pos = tree_rep.find_first_of(":", final_parenthesis + 1);
  if (colon_pos == string::npos)
    return 0.0;
  const size_t non_numeric = 
    tree_rep.find_first_not_of("0123456789.", colon_pos + 1);
  if (non_numeric == 0)
    return 0.0;
  return atof(tree_rep.substr(colon_pos + 1).c_str());
}


static bool
root_has_name(const string &tree_rep) {
  const size_t final_parenthesis = tree_rep.find_last_of(")");
  const size_t possible_name_start = 
    (final_parenthesis == string::npos) ? 0 : final_parenthesis;
  return (tree_rep[possible_name_start] != ':');
}


static void
extract_subtrees(const string &tree_rep,
		 vector<string> &subtree_reps) {
  const size_t offset = (tree_rep[0] == '(') ? 1 : 0;
  const string tmp(tree_rep.begin() + offset, 
		   tree_rep.begin() + tree_rep.find_last_of(")"));
  vector<size_t> split_points;
  size_t p_count = 0;
  for (size_t i = 0; i < tmp.length(); ++i) {
    if (p_count == 0 && tmp[i] == ',')
      split_points.push_back(i);
    if (tmp[i] == '(') ++p_count;
    if (tmp[i] == ')') --p_count;
  }
  split_points.push_back(tmp.length());
  subtree_reps.push_back(tmp.substr(0, split_points[0]));
  for (size_t i = 0; i < split_points.size() - 1; ++i)
    subtree_reps.push_back(string(tmp.begin() + split_points[i] + 1,
				  tmp.begin() + split_points[i + 1]));
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
  PhyloTree() {}
  PhyloTree(string tree_rep) {
    check_balanced_parentheses(tree_rep);
    // remove whitespace
    string::iterator w = 
      std::remove_copy_if(tree_rep.begin(), tree_rep.end(),
			  tree_rep.begin(), &isspace);
    assert(w != tree_rep.begin());
    tree_rep.erase(--w, tree_rep.end()); // The "--w" is for the ";"
    root = PhyloTreeNode(tree_rep);
  }
  string tostring() const {return root.tostring();}
  string tostring(const string &label) const;
  string Newick_format() const {return root.Newick_format() + ";";}
  string Newick_format(const string &label) const;
  bool label_exists(const string &label) const {
    return root.label_exists(label);
  }
  void fill_leaf_names(const string prefix, size_t &count); 
  void fill_names(const string prefix, size_t &count);
  void get_leaf_names(vector<string> &leaf_names );
  string find_common_ancestor(const vector<string> &leaf_names); 

private:
  PhyloTreeNode root;
};

string
PhyloTree::tostring(const string &label) const{
  return root.tostring(label);
}

string 
PhyloTree::Newick_format(const string &label) const{
  assert(label_exists(label));
  return root.Newick_format(label)+ ";";
}

void
PhyloTree::fill_leaf_names(const string prefix, size_t &count) {
  root.fill_leaf_names(prefix, count);
}

void
PhyloTree::fill_names(const string prefix, size_t &count) {
  root.fill_names(prefix, count);
}

void
PhyloTree::get_leaf_names(vector<string> & leaf_names){
  root.get_leaf_names(leaf_names);
}

string
PhyloTree::find_common_ancestor(const vector<string> &leaf_names){
  vector<string> copy_names = leaf_names;
  vector<string>::iterator it;
  it = std::unique(copy_names.begin(), copy_names.end());
  if(it != copy_names.end()){
    throw SMITHLABException("Names are not unique in query for common ancestor");
  }
  size_t count = 0; //# of leaves found
  string ancestor;
  bool found = root.find_common_ancestor(leaf_names, count, ancestor);
  if (!found)
    throw SMITHLABException("Ancestor not found");
  return ancestor;
}

std::istream&
operator>>(std::istream &in, PhyloTree &t) {
  /* doing it this way to just get one tree at a time */
  string r;
  char c;
  bool found_end = false;
  while (in >> c && !found_end) {
    r += c;
    if (c == ';')
      found_end = true;
  }
  if (!found_end)
    throw SMITHLABException("bad tree format");
  t = PhyloTree(r);
  return in;
}


std::ostream&
operator<<(std::ostream &out, const PhyloTree &t) {
  return out << t.Newick_format();
}


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
   

    
 
   //get common ancestor of two random selected leaves. 
    string ancestor;
    vector<string> tmp;
    vector<string> leaf_names;
    t.get_leaf_names(leaf_names); 
    srand(time(NULL));
    std::random_shuffle(leaf_names.begin(), leaf_names.end());
    tmp.assign(leaf_names.begin(), leaf_names.begin()+3);
    cerr << "the common ancestor of "  ;
    for(size_t i =0; i < tmp.size()-1; ++i) 
      cerr << tmp[i] << ","; 
    cerr <<  tmp[tmp.size()-1] << " is " ;
    ancestor = t.find_common_ancestor(tmp);
    cerr << ancestor << endl; 
    // print subtree rooted at ancestor
    cout << t.tostring(ancestor) << endl;
    cout << t.Newick_format(ancestor) << endl;


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
