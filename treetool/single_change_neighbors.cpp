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

#include "single_change_neighbors.hpp"

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

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::tr1::unordered_set;


bool operator==( const unordered_set<string> &lhs, const unordered_set<string> &rhs ){
  bool equal = (lhs.size()==rhs.size()) ;
  for(unordered_set<string>::const_iterator it = lhs.begin(); equal && it != lhs.end() ; ++it){
    const string tmpkey= *it; 
    equal = (rhs.count(tmpkey)>0);
  }
  return equal;
}

/* Binary combination of states among all leaf nodes,
 * leaf name is present in the set if its state is 'on' 
 * in the combination
*/
void
MethPhyloTree::build_states() {
  assert(tree.unique_names());
  vector<string> leaf_names;
  tree.get_leaf_names(leaf_names);
  states.clear();
  unordered_set<string> EmptySet;
  unordered_set<string> tmp;
  states.push_back(EmptySet);
  for (size_t i =0; i< leaf_names.size(); ++i) {
    tmp.clear();
    size_t N = states.size();
    for (size_t j = 0; j < N; ++j) {
      tmp = states[j];
      tmp.insert(leaf_names[i]);
      states.push_back(tmp);
    }
  }
}

/*Assign Neighbor relations*/
void
MethPhyloTree::build_neighbor() {
  assert(!states.empty());
  vector<unordered_set<string> > clade_leaves;
  tree.get_clade_leaves(clade_leaves);
  neighbor.clear();
  for (size_t i = 0; i < states.size(); ++i)
    neighbor.push_back(vector<size_t>(1, i)); //neighbor with it self
  
  for (size_t i = 0; i < states.size(); ++i) {
    for (size_t j = 0; j < clade_leaves.size(); ++j) {
      unordered_set<string> tmp  = states[i];
      size_t nb;
      // First exclude these leaves
      
      for (const string &x: clade_leaves[j])
	tmp.erase(x);

      vector<unordered_set<string> >::iterator it;
      bool found = false; 
      for( it = states.begin(); !found && it !=states.end(); ++ it){ 
	found = (*it == tmp);
      }
      assert(found);
      it --;	    
      nb = std::distance(states.begin(), it);
      neighbor[i].push_back(nb);
      neighbor[nb].push_back(i);
          
      //Then include these leaves
      tmp.insert(clade_leaves[j].begin(), clade_leaves[j].end());
      found = false;
      for( it = states.begin(); !found && it != states.end(); ++ it){ 
	found = (*it == tmp);
      }
      assert(found);
      it --;
      nb = std::distance(states.begin(), it);
      neighbor[i].push_back(nb);
      neighbor[nb].push_back(i);
    }
  }
  
  //remove redundancy in neighbor
  vector<size_t>::iterator it;
  for(size_t i = 0; i < neighbor.size(); ++i){
    sort(neighbor[i].begin(), neighbor[i].end());
    it = unique(neighbor[i].begin(), neighbor[i].end()); //                
    neighbor[i].resize( std::distance(neighbor[i].begin(),it) );
  }
}


