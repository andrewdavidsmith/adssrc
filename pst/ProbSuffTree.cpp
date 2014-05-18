/*
 *    Copyright (C) 2014 Andrew D. Smith
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

#include "ProbSuffTree.hpp"

#include <algorithm>
#include <numeric>
#include <iomanip>

#include <gsl/gsl_cdf.h>

#include "smithlab_utils.hpp"
#include "RNG.hpp"

using std::vector;
using std::endl;
using std::divides;
using std::string;
using std::ostringstream;
using std::bind2nd;
using std::plus;
using std::cerr;

using std::ios_base;

using smithlab::alphabet_size;

Runif PST::rng;


static string
format_probs(const vector<double> &p) {
  ostringstream oss;
  oss << '(' << std::setprecision(2) << std::fixed << p.front();
  for (size_t i = 1; i < p.size(); ++i)
    oss << ',' << std::setprecision(2) << std::fixed << p[i];
  oss << ')';
  return oss.str();
}


static char
sample_base(const vector<double> &probs, const Runif &rng) {
  assert(!probs.empty());

  const double random = rng.runif(0.0, 1.0);
  size_t i = 1;
  double curr_sum = probs[0];
  while (i < alphabet_size && curr_sum < random)
    curr_sum += probs[i++];
  return int2base(i);
}


static bool
sig_diff_probs(const vector<double> &a, const vector<double> &b, 
	       const double pseudo, const double crit) {
  assert(a.size() == b.size() && a.size() == alphabet_size);
  
  const double a_total = accumulate(a.begin(), a.end(), 0.0) + a.size()*pseudo;
  const double b_total = accumulate(b.begin(), b.end(), 0.0) + a.size()*pseudo;
  
  vector<double> expected(a);
  transform(expected.begin(), expected.end(), expected.begin(), 
	    bind2nd(plus<double>(), pseudo));
  transform(expected.begin(), expected.end(), expected.begin(),
	    bind2nd(divides<double>(), a_total/b_total));

  vector<double> observed(b);
  transform(observed.begin(), observed.end(), observed.begin(), 
	    bind2nd(plus<double>(), pseudo));
  
  double test_stat = 0.0;
  for (size_t i = 0; i < a.size(); ++i)
    test_stat += 
      (observed[i] - expected[i])*(observed[i] - expected[i])/expected[i];
  
  return gsl_cdf_chisq_Q(test_stat, alphabet_size - 1) < crit;
}


/******************************* "PSTNode" CLASS *********************************/


bool
PSTNode::has_child(const size_t i) const {
  assert(!child.empty() && child.size() > i); // sanity check
  return !child[i].probs.empty();
}


string
PSTNode::tostring(size_t depth) const {
  ostringstream oss;
  if (!child.empty())
    for (size_t i = 0; i < alphabet_size; i++)
      if (has_child(i))
	oss << string(depth, ' ') << int2base(i) << ' ' 
	    << format_probs(child[i].probs) << endl
	    << child[i].tostring(depth + 1);
  return oss.str();
}


char
PSTNode::generate(const string &s, const size_t pos, const Runif &rng) const {
  if (pos == s.length() || !isvalid(s[pos]) ||
      child.empty() || !has_child(base2int(s[pos])))
    return sample_base(probs, rng);
  else return child[base2int(s[pos])].generate(s, pos + 1, rng);
}


double
PSTNode::score_emission(const string &s, const size_t pos, 
			const char emission) const {
  assert(isvalid(emission));
  if (pos == s.length() || !isvalid(s[pos]) ||
      child.empty() || !has_child(base2int(s[pos])))
    return log(probs[base2int(emission)]);
  else return child[base2int(s[pos])].score_emission(s, pos + 1, emission);
}


void
PSTNode::insert(const string &s, const size_t pos,
		const char next_base, const size_t current_depth, 
		const size_t max_depth) {
  // at this depth is the creating and updating of probabilities
  if (probs.empty())
    probs.resize(alphabet_size, 0.0);
  probs[base2int(next_base)]++;
  
  // recursion adds or updates children
  if (current_depth < max_depth && pos < s.length() && isvalid(s[pos])) {
    if (child.empty()) 
      child.resize(alphabet_size);
    child[base2int(s[pos])].insert(s, pos + 1, next_base, 
				   current_depth + 1, max_depth);
  }
}


void
PSTNode::convert_to_probabilities(const double pseudocount) {
  assert(!probs.empty());
  transform(probs.begin(), probs.end(), probs.begin(),
	    bind2nd(plus<double>(), pseudocount));
  const double total = accumulate(probs.begin(), probs.end(), 0.0);
  transform(probs.begin(), probs.end(), probs.begin(), 
	    bind2nd(divides<double>(), total));
  
  if (!child.empty())
    for (size_t i = 0; i < child.size(); ++i)
      if (has_child(i))
	child[i].convert_to_probabilities(pseudocount);
}


void
PSTNode::prune(const double pseudo, const double crit) {
  if (child.empty()) return;
  
  size_t child_count = 0ul;
  for (size_t i = 0; i < alphabet_size; i++)
    if (has_child(i)) {
      child[i].prune(pseudo, crit);
      ++child_count;
    }
  for (size_t i = 0; i < alphabet_size; i++)
    if (has_child(i) && 
	!sig_diff_probs(probs, child[i].probs, pseudo, crit)) {
      child[i].probs.clear();
      --child_count;
    }
  if (child_count == 0)
    child.clear();
}


/****************************** 'PST' CLASS ************************************/


string
PST::tostring() const {
  ostringstream s;
  s << root.tostring(0);
  return s.str();
}


void
PST::insert(const string &seq, const size_t max_depth) {
  if (ALREADY_PRUNED || HAS_PROBABILITIES)
    throw SMITHLABException("PST no longer insertable");
  
  string rev(seq);
  reverse(rev.begin(), rev.end());
  for (size_t j = 0; j < rev.length(); ++j)
    if (isvalid(rev[j]))
      root.insert(rev, j + 1, rev[j], 0, max_depth);
}


void
PST::prune(const double pseudocount, const double critical_value) {
  if (ALREADY_PRUNED || HAS_PROBABILITIES)
    throw SMITHLABException("trying to prune PST twice");
  
  ALREADY_PRUNED = true;
  root.prune(pseudocount, critical_value);
}
  

void
PST::convert_to_probabilities(const double pseudocount) {
  HAS_PROBABILITIES = true;
  root.convert_to_probabilities(pseudocount);
}


void
PST::generate(const size_t len, string &to_generate) const {
  if (!HAS_PROBABILITIES)
    throw SMITHLABException("in generate but PST probabilities not set");
  
  to_generate = string(' ', len);
  for (size_t i = len; i > 0; --i) 
    to_generate[i - 1] = root.generate(to_generate, i, rng);
  std::reverse(to_generate.begin(), to_generate.end());
}


double
PST::score_sequence(string s) const {
  if (!HAS_PROBABILITIES)
    throw SMITHLABException("in score_sequence but PST probabilities not set");
  
  std::reverse(s.begin(), s.end());
  double score = 0.0;
  for (size_t j = 0; j < s.length(); j++)
    score += root.score_emission(s, j + 1, s[j]);
  return score;
}
