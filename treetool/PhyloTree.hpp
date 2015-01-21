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

#ifndef PHYLOTREE_HPP
#define PHYLOTREE_HPP

#include <string>
#include <vector>
#include <tr1/unordered_set>


class PhyloTreeNode {
public:
  PhyloTreeNode(){}
  PhyloTreeNode(const std::string &subtree_string);
  
  //utilities
  bool has_children() const {return !child.empty();}
  bool is_leaf() const {return child.empty();}
  bool unique_names(std::tr1::unordered_set<std::string> &names) const;
  bool label_exists(const std::string &label) const;
  std::string tostring(const size_t depth = 0) const;
  std::string tostring(const std::string &label) const;
  std::string treerep() const;
  std::string treerep(const std::string &label) const;
  std::string Newick_format() const{ return treerep() + ";";}
  std::string Newick_format(const std::string &label) const {
    return treerep(label) + ";";
  }
  size_t find_common_ancestor(const std::vector<std::string> &names, 
			      std::string &ancestor, bool &found)const;
  //mutator
  void fill_leaf_names(const std::string prefix, size_t &count);
  void fill_names(const std::string prefix, size_t &count);
  void set_branch(const double newlength) { branch = newlength;}
  bool set_branch(const std::string label, const double newlength);
  void set_name(const std::string newname) { name = newname;}
  void set_child(std::vector<PhyloTreeNode> &newchild){
    child.assign(newchild.begin(), newchild.end());
  }
  bool trim_to_keep(const std::vector<std::string>& leaves); 

  //accessor
  std::string get_name() const {return name;}
  double get_branch() const{return branch;}
  void get_child(std::vector<PhyloTreeNode> &newchild) const{
    newchild.assign(child.begin(), child.end());
  }
  size_t get_leaf_num() const;
  size_t get_child_size() const {return child.size();}
  void get_child_names(std::vector<std::string> &child_names)const;
  void get_leaf_names(std::vector<std::string> &leaf_names)const;
  void get_clade_leaves(std::vector<std::tr1::unordered_set<std::string> > 
			&clade_leaves)const;
  void get_node_names(std::vector<std::string> &node_names)const;
  void get_node_names(const std::string label, 
		      std::vector<std::string> &node_names)const;
  void get_branches(std::vector<double> &branches)const;

private:
  std::vector<PhyloTreeNode> child;
  std::string name;
  double branch; // distance to parent
};

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////  PHYLO TREE CLASS BELOW HERE                                   ////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

class PhyloTree {
public:
  PhyloTree() {}
  PhyloTree(std::string newick);
  PhyloTree(const PhyloTree &rhs){
    root = PhyloTreeNode(rhs.Newick_format());
  }

  std::string tostring() const {return root.tostring();}
  std::string Newick_format() const {return root.Newick_format() ;}
  std::string Newick_format(const std::string &label) const;
  bool label_exists(const std::string &label) const {
    return root.label_exists(label);
  }
  bool unique_names() const;

  void fill_leaf_names(const std::string prefix, size_t &count); 
  void fill_names(const std::string prefix, size_t &count);
  bool set_branch(const std::string label, const double newlength);


  void get_leaf_names(std::vector<std::string> &leaf_names)const;
  void get_child_names(std::vector<std::string> &child_names)const;
  void get_node_names(std::vector<std::string> &node_names)const;
  void get_node_parent_idx(std::vector<size_t> &pa_idx)const;
  void get_node_child_idx(std::vector<std::vector<size_t> > &child_idx )const;
  void get_node_names(const std::string label, 
		      std::vector<std::string> &node_names)const;
  void get_clade_leaves(std::vector<std::tr1::unordered_set<std::string> > 
			&clade_leaves)const;
  std::string get_root_name() const{ return root.get_name();}
  double get_root_branch() const{ return root.get_branch();}
  void get_branches(std::vector<double> &branches) const;
  size_t get_child_size() const{return root.get_child_size();}
  
  bool find_common_ancestor(const std::vector<std::string> &names, 
			    std::string &ancestor) const;
  void trim_to_keep(const std::vector<std::string>& leaves);
  
private:
  PhyloTreeNode root;
};

std::istream&
operator>>(std::istream &in, PhyloTree &t);

std::ostream&
operator<<(std::ostream &out, const PhyloTree &t);


std::string 
combine_newick(const std::string s1, const std::string s2, 
	       const std::string rootname,
	       const double branch1, const double branch2, 
	       const double branch0);


#endif
