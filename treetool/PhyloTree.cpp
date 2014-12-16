/*    Copyright (C) 2014 University of Southern California and
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


#include "PhyloTree.hpp"

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

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::tr1::unordered_set;


static bool
check_balanced_parentheses(const string &s) {
  int count = 0;
  for (size_t i = 0; i < s.length() && count >= 0; ++i) {
    if (s[i] == '(') ++count;
    if (s[i] == ')') --count;
  }
  return count == 0;
}

/// ADS: this function is really poorly coded
//check if the names in newick string are unique
static bool
check_unique_names(const string &s) {
  bool is_unique = true;
  string n, s_copy = s;
  unordered_set<string> names;
  size_t pos;
  while(!s_copy.empty() && is_unique){
    pos = s_copy.find_first_not_of("()0123456789:.,;");
    s_copy.erase(0,pos);
    if(!s_copy.empty()){
      pos = s_copy.find_first_of("()0123456789:.,;");
      n = s_copy.substr(0,pos);
      s_copy.erase(0,pos);
      if (!names.insert(n).second)
	is_unique = false; 
    }
  }
  return is_unique;
}





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
  if (label == name){
    return tostring();
  }
  string stringrep;
  for (size_t i = 0; i < child.size(); ++i){
    if (child[i].label_exists(label))
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

/*Name unnamed leaf nodes in the descendants
 *in the format of [prefix][index], keeping track
 *of next available index with variable count
 */
void 
PhyloTreeNode::fill_leaf_names(const string prefix, size_t &count)  {
  if (is_leaf()) {
    if (name.empty()) {
      name = prefix + toa(count);
      count++; 
    }
  }
  else
    for (size_t i = 0; i < child.size(); ++i)
      child[i].fill_leaf_names(prefix, count);
}

/*Name all unnamed nodes in the descendants
 *in the format of [prefix][index], keeping track
 *of next available index with variable count */
void 
PhyloTreeNode::fill_names(const string prefix, size_t &count)  {
  if (name.empty()) {
    name = prefix + toa(count);
    count++; 
  } 
  if (!is_leaf()) {
    for (size_t i = 0; i < child.size(); ++i)
      child[i].fill_names(prefix, count);
  }
}


/*Put names of all descendant leaf nodes into vector leaf_names*/
void 
PhyloTreeNode::get_leaf_names(vector<string> &leaf_names){
  if (is_leaf()) {
    assert(!name.empty());
    leaf_names.push_back(name);
  }
  else {
    for (size_t i = 0; i < child.size(); ++i)
      child[i].get_leaf_names(leaf_names);
  }
}


//Return number of hits of target ``names'' in the clade rooted here
//If nearest common ancestor of ``names'' is found in the subtree,
//``ancestor'' would be assigned with the ancestor's name 
size_t
PhyloTreeNode::find_common_ancestor(const vector<string> &names, 
				    string &ancestor) {
  size_t count = 0;
  for (size_t i = 0; i < child.size() && ancestor.empty(); ++i)
    count += child[i].find_common_ancestor(names, ancestor);
  
  count += find(names.begin(), names.end(), name) != names.end();
  
  if (ancestor.empty() && count == names.size())
    ancestor = name;
  
  return count;
}



/*Return total number of leaf nodes in the descendants */
size_t 
PhyloTreeNode::get_leaf_num() const{
  size_t num = 0;
  if (is_leaf()){
    num = 1;
  } else{
    for(size_t i =0; i < child.size(); ++i)
      num += child[i].get_leaf_num();
  }
  return num;

}



/*Put names of all children names into vector child_names*/
void 
PhyloTreeNode::get_child_names(vector<string> &child_names){
  if (!is_leaf()) 
    for (size_t i = 0; i < child.size(); ++i)
      child_names.push_back(child[i].get_name());
}



/*Return true if all names in the subtree are unique, 
 *and none of which have appeared in ``existing_names''.
 *Add all names in the subtree to ``existing_names''.
 */
bool 
PhyloTreeNode::unique_names(unordered_set<string> &existing_names) const{
  bool is_unique =  existing_names.count(name) ==0 ;
  existing_names.insert(name);

  if (is_unique && !is_leaf()){
      for(size_t i =0; i < child.size(); ++i)
	is_unique = child[i].unique_names(existing_names);
  }

  return is_unique;
}


/*Get all leaves in the subtree, and add the set to ``clade_leaves'' */
void 
PhyloTreeNode::get_clade_leaves(vector<unordered_set<string> > &clade_leaves){
  unordered_set<string> clade;
  unordered_set<string> tmp;
  if(is_leaf()){
    clade.insert(name);
  } else{
    for(size_t i =0; i < child.size(); ++i){
      child[i].get_clade_leaves(clade_leaves);
      tmp = clade_leaves.back();
      clade.insert(tmp.begin(), tmp.end());
    }
  }
  clade_leaves.push_back(clade);
}

/*Get all node names in the subtree*/
void 
PhyloTreeNode::get_node_names(std::vector<std::string> &node_names){
  node_names.push_back(get_name());
  if(!is_leaf()){
    for(size_t i=0; i < child.size(); ++i){
      child[i].get_node_names(node_names);
    }
  }
}

/*Get all node names in the subtree rooted at the node with name label*/
void 
PhyloTreeNode::get_node_names(const std::string label, std::vector<std::string> &node_names){
  if(name == label)
    get_node_names(node_names);
  else if (label_exists(label)){
    for(size_t i =0; i < child.size(); ++i )
      child[i].get_node_names(label, node_names);
  }
}


bool 
PhyloTreeNode::trim_to_keep(const std::vector<std::string>& leaves){
  bool keep = true;
  if(is_leaf() ){
    if(std::find(leaves.begin(), leaves.end(), name)!=leaves.end())
      keep = true;
  } else{
    size_t j = 0; 
    for(size_t i =0; i < child.size(); ++i){
      if(child[i].trim_to_keep(leaves) ){
	child[j] = child[i];
	++j;
      }
    }
    child.resize(j);
    keep = (j>0); 
    if (j ==1) { //only one leaf stays
      set_branch_length(branch_length + child[0].get_branch_length());
      set_name(child[0].get_name());
      std::vector<PhyloTreeNode> newchild;
      child[0].get_child(newchild);
      set_child(newchild);
    }
  }
  return keep;
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

// PhyloTreeNode::PhyloTreeNode(const PhyloTreeNode &node1,
// 			     const PhyloTreeNode &node2,
// 			     const string label,
// 			     const double bl){
//   name = label;
//   branch_length = bl; 
//   child.push_back(node1);
//   child.push_back(node2);
//}




////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////  PHYLO TREE CLASS BELOW HERE                                   ////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


PhyloTree::PhyloTree(string newick) {
  if(!check_balanced_parentheses(newick))
    throw SMITHLABException("Unbalanced parentheses in the string representation.");
  if(!check_unique_names(newick))
    throw SMITHLABException("Repeated names in string representation.");
  // remove whitespace
  string::iterator w = 
    std::remove_copy_if(newick.begin(), newick.end(),
			newick.begin(), &isspace);
  assert(w != newick.begin());
  newick.erase(--w, newick.end()); // The "--w" is for the ";"
  
  if(! check_unique_names(newick))
    throw SMITHLABException("Duplicated names in the string representation.");
  
  root = PhyloTreeNode(newick);
}

// PhyloTree::PhyloTree(PhyloTree &t1, PhyloTree &t2){
  
  
// }

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
PhyloTree::get_leaf_names(vector<string> &leaf_names){
  root.get_leaf_names(leaf_names);
}

void
PhyloTree::get_child_names(vector<string> &child_names){
  root.get_child_names(child_names);
}

void 
PhyloTree::get_node_names(std::vector<std::string> &node_names){
  root.get_node_names(node_names);
}

void 
PhyloTree::get_node_names(const std::string label, std::vector<std::string> &node_names){
  root.get_node_names(label, node_names);
}



bool 
PhyloTree::unique_names(){
  unordered_set<string> names;
  names.insert("");
  bool is_unique = root.unique_names(names);
  return is_unique;
}


/*Find the nearest common ancestor for a set of nodes in the tree
 *Return true if found;
 */
void
PhyloTree::find_common_ancestor(const vector<string> &names, string &ancestor) {
  root.find_common_ancestor(names, ancestor);
}

/*Get a collection of complete set of leaves 
 *for all clades in the tree
 */
void
PhyloTree::get_clade_leaves(vector<unordered_set<string> > &clade_leaves) {
  assert(unique_names());
  root.get_clade_leaves(clade_leaves);
}

void 
PhyloTree::trim_to_keep(const std::vector<std::string>& leaves){
  for(size_t i = 0; i < leaves.size(); ++i)
    assert(label_exists(leaves[i]));
  root.trim_to_keep(leaves);
}

std::istream&
operator>>(std::istream &in, PhyloTree &t) {
  /* doing it this way to just get one tree at a time */
  std::string r;
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
