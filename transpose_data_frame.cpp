/* transpose_data_frame: transpose a data frame (hopefully fast and
 * light)
 *
 * Copyright (C) 2020 Andrew D. Smith
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
#include <fstream>
#include <iterator>
#include <sstream>

#include "OptionParser.hpp"
#include "smithlab_os.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::runtime_error;

static void
parse_table_row(const string &row, vector<double> &values) {
  // ADS: not safe c++! This function is an attempt to be fast,
  // knowing "row" is entirely floating point numbers
  if (row.empty() || !isalnum(row[0]))
    throw runtime_error("rownames must begin with alphanumeric");

  values.clear(); // keep the reserve, but use push_back

  // skip the row's "name"
  const char* a = &row[row.find_first_of(" \t")];
  char* b = 0;
  double val = strtod(a, &b);

  while (a != b) {
    values.push_back(val);
    a = b;
    val = strtod(a, &b);
  }
}


int
main(int argc, const char **argv) {

  try {

    static const string description =
      "program to transpose a matrix";

    bool verbose = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<input-file> <output-file>");
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
    opt_parse.set_show_defaults();
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
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string infile(leftover_args.front());
    const string outfile(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    std::ifstream in(infile);
    if (!in)
      throw runtime_error("could not open file: " + infile);

    if (verbose)
      cerr << "[reading input]" << endl;

    string header;
    getline(in, header);
    std::istringstream h(header);
    vector<string> colnames = {std::istream_iterator<string>{h}, {}};

    if (verbose)
      cerr << "[n columns in data matrix: " << colnames.size() << "]" << endl;

    vector<vector<double> > the_matrix;

    string line;
    vector<string> rownames; // name for each row
    vector<double> curr_row(colnames.size());
    while (getline(in, line)) {

      rownames.push_back(line.substr(0, line.find_first_of(" \t")));
      parse_table_row(line, curr_row);

      if (curr_row.size() != colnames.size())
        throw runtime_error("wrong number of values: " + line);

      the_matrix.push_back(curr_row);
    }

    if (verbose)
      cerr << "[n rows in data matrix: " << rownames.size() << "]" << endl;

    std::ofstream out(outfile);
    if (!out)
      throw runtime_error("could not open output file: " + outfile);
    out.precision(std::numeric_limits<double>::max_digits10);

    if (verbose)
      cerr << "[writing output]" << endl;

    copy(std::begin(rownames), std::end(rownames),
         std::ostream_iterator<string>(out, "\t"));
    out << endl;

    for (size_t i = 0; i < colnames.size(); ++i) {
      out << colnames[i];
      for (size_t j = 0; j < rownames.size(); ++j)
        out << '\t' << the_matrix[j][i];
      out << '\n';
    }
  }
  catch (const std::runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
