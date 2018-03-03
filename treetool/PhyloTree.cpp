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
*/


#include "PhyloTree.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <unordered_set>

#include <cctype> // for isspace

#include "smithlab_utils.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::unordered_set;


static bool
check_balanced_parentheses(const string &s) {
  int count = 0;
  for (size_t i = 0; i < s.length(); ++i) {
    if (s[i] == '(') ++count;
    if (s[i] == ')') --count;
  }
  return count == 0;
}


// static bool
// check_unique_names(const string &s) {
//   /// ADS: this function is really poorly coded
//   bool is_unique=true;
//   string s_copy =s;
//   unordered_set<string> names;
//   size_t f1, f2;
//   while (!s_copy.empty()) {
//     f1 = s_copy.find_first_not_of("()0123456789:,;");
//     s_copy.erase(0,f1);
//     if (s_copy.empty()) break;
//     f2 = s_copy.find_first_of("()0123456789:,;");
//     string name=s_copy.substr(0,f2);
//     if (names.find(name)== names.end()){
//       names.insert(name);
//     }
//     else {
//       is_unique=false;
//       break;
//     }
//   }
//   return is_unique;
// }


void
PhyloTree::PTNode::swap(PTNode &other) {
  std::swap(child, other.child);
  std::swap(name, other.name);
  std::swap(branch_length, other.branch_length);
}

size_t
PhyloTree::PTNode::get_size() const {
  size_t sz = 1;
  for (size_t i = 0; i < child.size(); ++i)
    sz += child[i].get_size();
  return sz;
}


string
PhyloTree::PTNode::tostring(const size_t depth) const {
  std::ostringstream oss;
  oss << string(depth, '\t') << branch_length << ':' << name;
  for (size_t i = 0; i < child.size(); ++i)
    oss << endl << child[i].tostring(depth + 1);
  return oss.str();
}


string
PhyloTree::PTNode::Newick_format() const {
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


PhyloTree::PTNode::PTNode(const string &tree_rep) {
  // This function needs to test the various ways that a string can be
  // passed in to represent a subtree at the current node

  branch_length = extract_branch_length(tree_rep);
  if (root_has_name(tree_rep))
    name = extract_name(tree_rep);

  if (!represents_leaf(tree_rep)) {
    vector<string> subtree_reps;
    extract_subtrees(tree_rep, subtree_reps);
    for (size_t i = 0; i < subtree_reps.size(); ++i)
      child.push_back(PTNode(subtree_reps[i]));
  }
}

void
PhyloTree::PTNode::copy_subtree_with_species(const PTNode &t,
                                             const
                                             unordered_set<string> &species,
                                             PTNode &u) {
  if (t.is_leaf()) {
    if (species.find(t.name) != species.end()) {
      u.name = t.name;
      u.branch_length = t.branch_length;
    }
  }
  else {

    vector<PTNode> children;
    for (size_t i = 0; i < t.child.size(); ++i) {
      PTNode v;
      copy_subtree_with_species(t.child[i], species, v);
      if (v.has_children() || !v.name.empty()) {
        children.push_back(PTNode());
        children.back().swap(v);
      }
    }

    if (children.size() > 0) {
      if (children.size() == 1) {
        u.swap(children.front());
        u.branch_length += t.branch_length;
      }
      else {
        u.name = t.name;
        u.branch_length = t.branch_length;
        u.child.swap(children);
      }
    }
  }
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////  PHYLO TREE CLASS BELOW HERE                                   ////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


PhyloTree::PhyloTree(string tree_rep) {
  if (!check_balanced_parentheses(tree_rep))
    throw SMITHLABException("Unbalanced parentheses in Newick format: " +
                            tree_rep);
  // remove whitespace
  string::iterator w =
    std::remove_copy_if(tree_rep.begin(), tree_rep.end(),
                        tree_rep.begin(), &isspace);
  assert(w != tree_rep.begin());
  if (*(w - 1) == ';')
    tree_rep.erase(--w, tree_rep.end()); // The "--w" is for the ";"

  /* ADS: requirement for unique names should only apply to leafs
   */
  // if (!check_unique_names(tree_rep))
  //   throw SMITHLABException("duplicate names in: " + tree_rep);

  root = PTNode(tree_rep);
}


void
PhyloTree::copy_subtree_with_species(const PhyloTree &t,
                                     const unordered_set<string> &species,
                                     PhyloTree &u) {
  PTNode::copy_subtree_with_species(t.root, species, u.root);
}


std::istream&
operator>>(std::istream &in, PhyloTree &t) {
  /* doing it this way to just get one tree at a time */
  string r;
  char c;
  bool found_end = false;
  while (!found_end && in >> c) {
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
