/* Copyright (C) 2017 University of Southern California and
 *                    Andrew D. Smith
 *
 * Authors: Andrew D. Smith
 *
 * expand_chains: expand chains to intervals. The chains should be
 * in available from the UCSC Genome Browser downloads and in a
 * directory called liftOver for a particular species. The names
 * should look like mm10ToHg19.over.chain.gz. An example can be
 * found here:
 *
 * http://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToHg19.over.chain.gz
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <exception>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "GenomicRegion.hpp"

#include "chain_file_utils.hpp"

using std::string;
using std::ios_base;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::runtime_error;

int
main(int argc, const char **argv) {
  try{
    bool VERBOSE = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "expand a chain format file",
                           "<input> <output>");
    opt_parse.add_opt("verbose", 'v', "print more information",
                      false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
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
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    const string chain_file(leftover_args.front());
    const string outfile(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile.c_str());
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
    if (!out)
      throw runtime_error("cannot open output file: " + outfile);

    chain tmp_chain;
    std::ifstream in_chains(chain_file.c_str());
    if (!in_chains)
      throw runtime_error("cannot open input file: " + chain_file);

    while (in_chains >> tmp_chain)
      out << tmp_chain;
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
