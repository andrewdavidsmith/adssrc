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

#ifndef PROBSUFFTREE_H
#define PROBSUFFTREE_H

#include <string>
#include <vector>

class Runif;

class PSTNode {

public:
  
  std::string tostring(size_t) const;
  void insert(const std::string &s, const size_t pos,
	      const char next_base, const size_t current_depth,
	      const size_t max_depth);
  double score_emission(const std::string &context_seq, 
			const size_t current_position,
			const char emission) const;
  char generate(const std::string &s, const size_t pos,
		const Runif &rng) const;
  
  void convert_to_probabilities(const double pseudocount);
  void prune(const double pseudocount, const double critical_value);
  
private:
  std::vector<PSTNode> child;
  std::vector<double> probs;
  
  bool has_child(const size_t base) const;
};


class PST {
public:

  void prune(const double pseudocount, const double critical_value);
  void insert(const std::string &training_sequence,
	      const size_t max_depth);
  double score_sequence(std::string to_score) const;
  void generate(const size_t len, std::string &to_generate) const;
  void convert_to_probabilities(const double pseudocount);
  std::string tostring() const;
  
private:
  PSTNode root;
  
  bool HAS_PROBABILITIES;
  bool ALREADY_PRUNED;
  
  static Runif rng;
  
};

#endif
