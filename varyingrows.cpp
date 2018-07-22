/* varyingrows: get (or remove) the most varying rows
 *
 * Copyright (C) 2018 Andrew D. Smith
 *
 * Authors: Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#include <string>
#include <vector>
#include <iostream>
#include <unordered_set>
#include <queue>
#include <fstream>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include <gsl/gsl_statistics_double.h>

using std::string;
using std::vector;
using std::ifstream;
using std::cerr;
using std::cout;
using std::endl;
using std::unordered_set;
using std::istringstream;
using std::priority_queue;
using std::pair;
using std::make_pair;

using std::greater;
using std::runtime_error;

static void
parse_table_row(const string &row, vector<double> &values) {
  std::istringstream is;
  is.rdbuf()->pubsetbuf(const_cast<char*>(row.c_str()), row.size());

  string dummy;
  is >> dummy; //eliminate the row name

  values.clear();
  double val = 0.0;
  while (is >> val)
    values.push_back(val);
}


static void
get_row_name(const string &line, string &rowname) {
  rowname = line.substr(0, line.find_first_of(" \t"));
}


static double
get_row_variance(const string &line) {
  vector<double> vals;
  parse_table_row(line, vals);
  return gsl_stats_variance(&vals[0], 1, vals.size());
}


int
main(int argc, const char **argv) {

  try {

    string outfile;
    bool VERBOSE = false;
    bool INVERT = false;

    size_t n_top = 10000;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "select to keep or remove the top "
                           "most varying rows in the matrix", "<matrix>");
    opt_parse.add_opt("outfile", 'o', "output file", false, outfile);
    opt_parse.add_opt("top", 't', "number of top varying to get", false, n_top);
    opt_parse.add_opt("invert", 'I', "invert the selection", false, INVERT);
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
    const string table_file(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/


    if (VERBOSE)
      cerr << "[extracting row variances]" << endl;
    ifstream in(table_file);
    if (!in)
      throw runtime_error("could not open file: " + table_file);

    typedef pair<double, string> var_row;
    priority_queue<var_row, vector<var_row>, greater<var_row>  > top_row_sorter;

    string header;
    getline(in, header); // first remove header

    string line;
    size_t lines_read = 0;
    while (getline(in, line)) {
      lines_read++;
      string curr_row;
      get_row_name(line, curr_row);
      const double var = get_row_variance(line);
      top_row_sorter.push(make_pair(var, curr_row));
      if (top_row_sorter.size() > n_top)
        top_row_sorter.pop();
      if (VERBOSE && (lines_read % 10000 == 0))
        cerr << "lines read: " << lines_read << '\r';
    }
    if (VERBOSE)
      cerr << "lines read: " << ++lines_read << endl;
    in.close();

    if (VERBOSE)
      cerr << "[identifying top varying rows]" << endl;
    unordered_set<string> good_rows;
    while (!top_row_sorter.empty()) {
      good_rows.insert(top_row_sorter.top().second);
      top_row_sorter.pop();
    }

    in.open(table_file);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    getline(in, line);
    out << line << endl;

    if (VERBOSE)
      cerr << "[selecting top varying rows]" << endl;
    while (getline(in, line)) {
      string curr_row;
      get_row_name(line, curr_row);
      const bool found = (good_rows.find(curr_row) != good_rows.end());
      if ((!INVERT && found) || (INVERT && !found))
        out << line << endl;
    }

  }
  catch (const std::runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
