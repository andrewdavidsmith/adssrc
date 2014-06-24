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

#ifndef SINGLE_CHANGE_NEIGHBORS_CPP
#define SINGLE_CHANGE_NEIGHBORS_CPP

#include <string>
#include <vector>
#include <tr1/unordered_set>

#include "PhyloTree.hpp"

struct MethPhyloTree {
  MethPhyloTree(){}
  MethPhyloTree(std::string tree_rep){
    tree = PhyloTree(tree_rep);
    //Name unnamed nodes
    size_t count = 0;
    tree.fill_leaf_names("Leaf", count);
    count = 0;
    tree.fill_names("Internal", count);  
  };
  
  void build_states();
  void build_neighbor();
  void find_common_ancestor(const std::vector<std::string> &names, 
			    std::string &ancestor);

  std::vector<std::tr1::unordered_set<std::string> > states;
  std::vector<std::vector<size_t> > neighbor;
  PhyloTree tree;

};


std::ostream&
operator<<(std::ostream &out, const MethPhyloTree &t);


#endif
