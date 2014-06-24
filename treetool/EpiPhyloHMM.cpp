
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

#include "EpiPhyloHMM.hpp"

#include <iomanip>
#include <numeric>
#include <limits>
#include <cmath>

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>

// #pragma omp <rest of pragma>

using std::vector;
using std::pair;
using std::setw;
using std::max;
using std::cerr;
using std::endl;
using std::string;
using std::setprecision;
using std::isfinite;
using std::tr1::unordered_set;

double TOL = 1.0e-10;

struct betabin {
  betabin(const double a, const double b) :
    alpha(a), beta(b), lnbeta_helper(gsl_sf_lnbeta(a, b)) {}
  double operator()(const pair<double, double> &val) const;
  void fit(const vector<double> &vals_a,
const vector<double> &vals_b,
const vector<double> &p);
  string tostring() const;
  double alpha;
  double beta;
  double lnbeta_helper;
};

string
betabin::tostring() const {
  std::ostringstream os;
  os << setprecision(4) << alpha << " " << setprecision(4) << beta;
  return os.str();
}

double
betabin::operator()(const pair<double, double> &val) const {
  const size_t x = static_cast<size_t>(val.first);
  const size_t n = static_cast<size_t>(x + val.second);
  return gsl_sf_lnchoose(n, x) +
    gsl_sf_lnbeta(alpha + x, beta + val.second) - lnbeta_helper;
}

inline static double
sign(double x) {
  return (x >= 0) ? 1.0 : -1.0;
}
static const double tolerance = 1e-10;

inline static double
invpsi(const double tolerance, const double x) {
  double L = 1.0, Y = std::exp(x);
  while (L > tolerance) {
    Y += L*sign(x - gsl_sf_psi(Y));
    L /= 2.0;
  }
  return Y;
}

static double
movement(const double curr, const double prev) {
  return std::abs(curr - prev)/std::max(std::fabs(curr), std::fabs(prev));
}

void
betabin::fit(const vector<double> &vals_a, const vector<double> &vals_b,
const vector<double> &p) {
  const double p_total = std::accumulate(p.begin(), p.end(), 0.0);
  const double alpha_rhs = inner_product(vals_a.begin(), vals_a.end(),
p.begin(), 0.0)/p_total;
  const double beta_rhs = inner_product(vals_b.begin(), vals_b.end(),
p.begin(), 0.0)/p_total;
  double prev_alpha = 0.0, prev_beta = 0.0;
  alpha = beta = 0.01;
  while (movement(alpha, prev_alpha) > tolerance &&
movement(beta, prev_beta) > tolerance) {
    prev_alpha = alpha;
    prev_beta = beta;
    alpha = invpsi(tolerance, gsl_sf_psi(prev_alpha + prev_beta) + alpha_rhs);
    beta = invpsi(tolerance, gsl_sf_psi(prev_alpha + prev_beta) + beta_rhs);
  }
  lnbeta_helper = gsl_sf_lnbeta(alpha, beta);
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
void
TransProb(const double tol, 
	  const double U2Mrate, const double M2Urate,
	  const double t, double &U2Mlp, double &M2Ulp){
  assert(U2Mrate > 0 && M2Urate > 0 && t > 0); 
  double h = std::min(1-tol, std::max(tol, exp(-t*(U2Mrate+M2Urate))));
  U2Mlp = log(U2Mrate)+ log(1-h) - log(U2Mrate+M2Urate);
  M2Ulp = log(M2Urate)+ log(1-h) - log(U2Mrate+M2Urate);
}


//need an unordered_map to hold species label, and their index in the read-count vectors

std::pair<double, double>
tree_prob_aux(const MethPhyloTree mt, 
	      const double U2Mrate, const double M2Urate, 
	      const unordered_set<string> &state,
	      const string label){
  // this will return the tree likelihood at position i 
  // for the subtree rooted at the node with name ``label'' 
  // given that the subtree root has methylation state
  // first: M, second: U
  assert(mt.tree.label_exists(label));
  double first = 1, second = 1;
  string subtreenk = mt.tree.Newick_format(label);
  PhyloTree tree= PhyloTree(subtreenk);
  if (tree.get_child_size() ==0){ //base cases
    first = ( state.count(label) >0 )? 1 : 0;
    second = 1 - first;
  } else{
    vector<string> child_names;
    tree.get_child_names(child_names);
    for(size_t i = 0; i < child_names.size(); ++i){
      std::pair<double, double> aux = tree_prob_aux( mt, U2Mrate, M2Urate, state, child_names[i]);
      PhyloTreeNode tmp = PhyloTreeNode(mt.tree.Newick_format(child_names[i]));
      double branch = tmp.get_branch_length();
      double U2Mlp, M2Ulp;
      TransProb(TOL, U2Mrate, M2Urate,branch, U2Mlp, M2Ulp);
      first = first*((1.0-exp(M2Ulp))*aux.first + exp(M2Ulp)*aux.second);
      second = second*(exp(U2Mlp)*aux.first + (1.0-exp(U2Mlp))*aux.second);
    }  
  }
  return std::make_pair(first, second);  
}

double 
tree_prob(std::pair<double, double> init_distr,
	const MethPhyloTree mt,
	const double U2Mrate, const double M2Urate, 
	const unordered_set<string> &state){
  string rootname = mt.tree.get_root_name();
  std::pair<double, double> root_aux = tree_prob_aux(mt, U2Mrate, M2Urate, state, rootname); 
  return root_aux.first*init_distr.first + root_aux.second*init_distr.second; 
}
