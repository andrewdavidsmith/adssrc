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


class PhyloTree {
public:
  PhyloTree() {}
  PhyloTree(std::string tree_rep);

  size_t get_size() const {return root.get_size();}
  std::string tostring() const {return root.tostring();}
  std::string Newick_format() const {return root.Newick_format() + ";";}

  static void
  copy_subtree_with_species(const PhyloTree &t, 
			    const std::tr1::unordered_set<std::string> &species,
			    PhyloTree &u);
  
protected:
  
  struct PTNode {
    PTNode() {}
    PTNode(const std::string &subtree_string);
    void swap(PTNode &other);
    
    bool has_children() const {return !child.empty();}
    bool is_leaf() const {return child.empty();}
    size_t get_size() const;

    std::string tostring(const size_t depth = 0) const;
    std::string Newick_format() const;
    
    static void
    copy_subtree_with_species(const PTNode &t, 
			      const std::tr1::unordered_set<std::string> 
			      &species,
			      PTNode &u);
    
    std::vector<PTNode> child;
    std::string name;
    double branch_length; // distance to parent
  };
  
  PTNode root;
};

std::istream&
operator>>(std::istream &in, PhyloTree &t);

std::ostream&
operator<<(std::ostream &out, const PhyloTree &t);

#endif
