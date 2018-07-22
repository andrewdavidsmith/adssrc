/* dataframe-euclidean-dist: produce a distance matrix for all pairs
 * of columns in a matrix with very many rows
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
#include <fstream>
#include <exception>
#include <sstream>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

using std::string;
using std::to_string;
using std::vector;
using std::ifstream;
using std::cerr;
using std::cout;
using std::endl;
using std::istringstream;
using std::runtime_error;

#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>

static size_t
count_lines_fast(const string &filename) {

  struct stat st;
  stat(filename.c_str(), &st);

  int fd = open(filename.c_str(), O_RDONLY, 0);
  if (fd < 0)
    throw runtime_error("bad file: " + filename);

  char *mmap_data = static_cast<char *>(mmap(NULL, st.st_size, PROT_READ,
                                             MAP_PRIVATE | MAP_NORESERVE, fd, 0));
  if (mmap_data == MAP_FAILED)
    throw runtime_error("failed mmap for: " + filename);

  const size_t n_lines = std::count(mmap_data, mmap_data + st.st_size, '\n');

  if (munmap(static_cast<void *>(mmap_data), st.st_size) != 0)
    throw runtime_error("failed to release mmap for: " + filename);

  close(fd);
  return n_lines;
}

template <typename T> T
euclidean_dist(const vector<T> &a, const vector<T> &b, vector<T> &aux) {
  transform(begin(a), end(a), begin(b), begin(aux),
            [](T elmnt1, T elmnt2) {
              const T x = elmnt1-elmnt2; return x*x;});
  return std::sqrt(accumulate(begin(aux), end(aux), 0.0));
}

template <typename T> static void
parse_table_row(const string &row, const size_t id, vector<vector<T> > &values) {
  std::istringstream is;
  is.rdbuf()->pubsetbuf(const_cast<char*>(row.c_str()), row.size());

  string dummy;
  is >> dummy; //eliminate first column

  for (size_t i = 0; i < values.size(); ++i)
    if (!(is >> values[i][id]))
      throw runtime_error("bad line: " + row + "[" + to_string(id + 1) + "]");
}

template <typename T> static void
parse_strings_whitespace(const string &line, vector<T> &parts) {
  istringstream parser;
  parser.rdbuf()->pubsetbuf(const_cast<char*>(line.c_str()), line.size());

  parts.clear();
  string buffer;
  while (parser >> buffer)
    parts.push_back(buffer);
}


int
main(int argc, const char **argv) {

  try {

    string outfile;
    bool VERBOSE = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "make a dist matrix "
                           "from a data frame", "<data-frame>");
    opt_parse.add_opt("outfile", 'o', "output file", false, outfile);
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

    ifstream in(table_file);
    if (!in)
      throw runtime_error("cannot open: " + table_file);

    string header;
    if (!getline(in, header))
      throw runtime_error("could not extract header from: " + table_file);

    const size_t n_rows = count_lines_fast(table_file) - 1;
    if (VERBOSE)
      cerr << "n_rows=" << n_rows << endl;

    vector<string> column_names;
    parse_strings_whitespace(header, column_names);

    const size_t n_columns = column_names.size();
    if (VERBOSE)
      cerr << "n_columns=" << n_columns << endl;

    vector<vector<double> > the_table(n_columns, vector<double>(n_rows));

    string line;
    size_t row_idx = 0;
    while (getline(in, line))
      parse_table_row(line, row_idx++, the_table);

    vector<double> auxiliary(n_rows, 0.0);
    vector<vector<double> > dist(n_columns, vector<double>(n_columns));
    for (size_t i = 0; i < the_table.size(); ++i)
      for (size_t j = 0; j < i; ++j)
        dist[i][j] = euclidean_dist(the_table[i], the_table[j], auxiliary);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    out << header;
    for (size_t i = 0; i < dist.size(); ++i) {
      out << column_names[i];
      for (size_t j = 0; j < i; ++j)
        out << '\t' << dist[i][j];
      out << endl;
    }

  }
  catch (const std::runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
