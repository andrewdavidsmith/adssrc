/*    treetool: a program to manipulate phylogenetic trees in Newick
 *    format and to test parsing code
 *
 *    Copyright (C) 2014 University of Southern California and
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

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

using std::string;
using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using std::tr1::unordered_set;

static const size_t MAX_LEAF_NUM = 10;

static bool
check_balanced_parentheses(const string &s) {
  int count = 0;
  for (size_t i = 0; i < s.length() && count >= 0; ++i) {
    if (s[i] == '(') ++count;
    if (s[i] == ')') --count;
  }
  return count == 0;
}

static bool
check_unique_names(const string &s) {
  bool is_unique=true;
  string s_copy =s;
  unordered_set<std::string> names;
  size_t f1, f2;
  while(!s_copy.empty()){
    f1 = s_copy.find_first_not_of("()0123456789:,;");
    s_copy.erase(0,f1);
    if(s_copy.empty()) break;
    f2 =s_copy.find_first_of("()0123456789:,;");
    string name=s_copy.substr(0,f2);
    if (names.find(name)== names.end()){
      names.insert(name);
    }else{
      is_unique=false; 
      break;
    }
  }
  return is_unique;
}

class PhyloTreeNode {
public:
  PhyloTreeNode() {}
  PhyloTreeNode(const string &subtree_string);
  
  bool has_children() const {return !child.empty();}
  bool is_leaf() const {return child.empty();}
  
  bool label_exists(const string &label) const;
  string tostring(const size_t depth = 0) const;
  string tostring(const string &label) const;
  string Newick_format() const;
  string Newick_format(const string & label) const;
  string get_bitlabel() const {return bitlabel;}
  size_t get_leaf_num() const;

  void fill_leaf_names(const string prefix, size_t &count);
  void fill_names(const string prefix, size_t &count);
  
  void get_leaf_names(vector<string> &leaf_names);
  size_t find_common_ancestor(const vector<string> &names, string &ancestor);
  void set_bitlabel(size_t &leaftracker, vector<string> &bitlabels);

private:
  vector<PhyloTreeNode> child;
  string name;
  string bitlabel;
  double branch_length; // distance to parent
};



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
  if(label == name)
    return( tostring());
  
  string stringrep;
  for(size_t i = 0; i < child.size(); ++i){
    if(child[i].label_exists(label))
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

/*Make bitlabel:
 *Label leaf nodes with bitstrings, each having exactly
 *one bit set to 1, depending on a DF traverse oreder;
 *A parent's bitlabel is the OR operation result of all the children;
 *The first argument keeps track of total number of labeled leaves.  
 *The second argument is a collection of all bitlabels assigned so far.
 */
void 
PhyloTreeNode::set_bitlabel(size_t &leaftracker, vector<string> &bitlabels){
  std::bitset<MAX_LEAF_NUM> bits;
  if(is_leaf()){
    bits.set(leaftracker);
    leaftracker ++;
  }else{
    for(size_t i =0; i< child.size();++i){
      child[i].set_bitlabel(leaftracker, bitlabels);
      std::bitset<MAX_LEAF_NUM> tmp(child[i].get_bitlabel());
      bits |= tmp;
    }
  }
  bitlabel = bits.to_string<char,std::string::traits_type,std::string::allocator_type>();
  bitlabels.push_back(bitlabel);
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


class PhyloTree {
public:
  PhyloTree() {}
  PhyloTree(string tree_rep) {
    if(!check_balanced_parentheses(tree_rep))
      throw SMITHLABException("Unbalanced parentheses in the string representation.");
    // remove whitespace
    string::iterator w = 
      std::remove_copy_if(tree_rep.begin(), tree_rep.end(),
			  tree_rep.begin(), &isspace);
    assert(w != tree_rep.begin());
    tree_rep.erase(--w, tree_rep.end()); // The "--w" is for the ";"
    
    if(! check_unique_names(tree_rep))
      throw SMITHLABException("Duplicated names in the string representation.");
    
    root = PhyloTreeNode(tree_rep);
  }
  string tostring() const {return root.tostring();}
  string tostring(const string &label) const;
  string Newick_format() const {return root.Newick_format() + ";";}
  string Newick_format(const string &label) const;
  bool label_exists(const string &label) const {
    return root.label_exists(label);
  }
  void fill_leaf_names(const string prefix, size_t &count); 
  void fill_names(const string prefix, size_t &count);
  void get_leaf_names(vector<string> &leaf_names );

  void find_common_ancestor(const vector<string> &names, string &ancestor);
  
  void build_neighbor();
  vector<vector<size_t> > get_neighbor()const{return neighbor;}
private:
  PhyloTreeNode root;
  size_t leaf_num;
  vector<vector<size_t> > neighbor;
  void set_bitlabel(vector<string> &bitlabels);
  void set_leaf_num();

};

string
PhyloTree::tostring(const string &label) const{
  return root.tostring(label);
}

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
PhyloTree::get_leaf_names(vector<string> & leaf_names){
  root.get_leaf_names(leaf_names);
}

/*Find the nearest common ancestor for a set of nodes in the tree
 *Return true if found;
 */
void
PhyloTree::find_common_ancestor(const vector<string> &names, string &ancestor){
  root.find_common_ancestor(names, ancestor);
}

/*Assign value to the private member leaf_num*/
void
PhyloTree::set_leaf_num(){
  leaf_num = root.get_leaf_num();
  if(leaf_num > MAX_LEAF_NUM)
    throw SMITHLABException("Number of leaves over MAX_LEAF_NUM");
}

/*Assign bitlabels to nodes*/
void 
PhyloTree::set_bitlabel(vector<string> &bitlabels){
  size_t tracker = 0;
  root.set_bitlabel(tracker, bitlabels);
}

/*Assign Neighbor relations basing on bitlabels*/
void 
PhyloTree::build_neighbor(){
  set_leaf_num();
  size_t N = pow(2, leaf_num);
  
  std::bitset<MAX_LEAF_NUM> modifier;
  for(size_t i=leaf_num; i<MAX_LEAF_NUM ; ++i)
    modifier.set(i);
  modifier.flip(); //0...0111...1  
  
  vector<string> bitlabels;
  set_bitlabel(bitlabels);

  neighbor.clear();
  for(size_t i =0; i < N; ++i)
    neighbor.push_back(vector<size_t>(1, i)); //each pattern is a neighbor of itself

   for(size_t i = 0; i < N; ++i){
     std::bitset<MAX_LEAF_NUM> pattern(i); //binary representation of the i-th pattern
     for(size_t j =0; j < bitlabels.size(); ++j){
       std::bitset<MAX_LEAF_NUM> subtree_bit_rep(bitlabels[j]);
       size_t n1 = (pattern | subtree_bit_rep).to_ulong();
       size_t n2 = ((~pattern | subtree_bit_rep)& modifier).to_ulong();
       neighbor[i].push_back(n1);
       neighbor[i].push_back(n2);
       neighbor[n1].push_back(i);
       neighbor[n2].push_back(i);
     }
   }
   //remover redundancy
   for(size_t i = 0; i < N; ++i){
     std::sort(neighbor[i].begin(), neighbor[i].end());
     vector<size_t>::iterator it;
     it = std::unique(neighbor[i].begin(), neighbor[i].end());
     neighbor[i].resize( std::distance(neighbor[i].begin(),it));
   }

}



std::istream&
operator>>(std::istream &in, PhyloTree &t) {
  /* doing it this way to just get one tree at a time */
  string r;
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


int 
main(int argc, const char **argv) {
  
  try {
    
    bool VERBOSE = false;
    string outfile;
    string label_to_check;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "manipulate Newick format "
			   "phylogenetic trees",
			   "<newick-input>");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false, outfile);
    opt_parse.add_opt("label", 'l', "check if this label exists", 
		      false, label_to_check);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
	   << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string newick_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/
    
    std::ifstream in(newick_file.c_str());
    if (!in)
      throw SMITHLABException("bad file: " + newick_file);
    
    PhyloTree t;
    in >> t;
    
    cout << t.tostring() << endl;
    cout << t << endl;

    //Name unnamed nodes
    size_t count = 0;
    t.fill_leaf_names("Leaf", count);
    count = 0;
    t.fill_names("Internal", count);  
 
    //get common ancestor of 2 random selected leaves. 
    string ancestor;
    vector<string> leaf_names;
    t.get_leaf_names(leaf_names); 

    srand(time(0));
    std::random_shuffle(leaf_names.begin(), leaf_names.end());
    vector<string> tmp(leaf_names.begin(), leaf_names.begin()+2);

    cerr << "the common ancestor of "  ;
    copy(tmp.begin(), tmp.end(), std::ostream_iterator<string>(cerr, ","));
    
    t.find_common_ancestor(tmp, ancestor);
    cerr << (ancestor.empty() ? "not found" : ancestor) << endl; 
    
    // print subtree rooted at ancestor
    cout << t.tostring(ancestor) << endl;
    cout << t.Newick_format(ancestor) << endl;
    
    // //Get neighbor relations
    // t.build_neighbor();
    // vector<vector<size_t> > neighbor(t.get_neighbor());
    // for(size_t i=0; i < neighbor.size(); ++i){
    //   cerr << "pattern" << i << " has neighbors" ;
    //   for(size_t j = 0; j < neighbor[i].size(); ++j)
    // 	cerr << neighbor[i][j] <<",";
    //   cerr << endl;
    // }

    if (!label_to_check.empty())
      cout << label_to_check << " "
	   << (t.label_exists(label_to_check) ? 
	       "exists" : "does not exist") << endl;
  }
  catch (const SMITHLABException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
