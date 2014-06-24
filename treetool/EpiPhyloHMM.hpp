
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


/****************************************
 A set of aligned methylomes
 MethPhyloTree containing the PhyloTree for the methylomes,
 and the neighbor relations among states of methylation patterns
 ***************************************/





#ifndef PHYLO_HMM_HPP
#define PHYLO_HMM_HPP

#include "smithlab_utils.hpp"
#include "PhyloTree.hpp"
#include "single_change_neighbors.hpp"
//#include<memory>

struct betabin;

class EpiPhyloHMM{
public:
  
  EpiPhyloHMM(const double mp, const double tol, 
	      const size_t max_itr, const bool v, bool d = false) :
    MIN_PROB(mp), tolerance(tol), max_iterations(max_itr),
    VERBOSE(v), DEBUG(d) {}

private:
  double MIN_PROB;
  double tolerance;
  size_t max_iterations;
  bool VERBOSE;
  bool DEBUG;
  MethPhyloTree mt;

};


double 
tree_prob(std::pair<double, double> init_distr,
	const MethPhyloTree mt,
	const double U2Mrate, const double M2Urate, 
	  const std::tr1::unordered_set<std::string> &state);


#endif
