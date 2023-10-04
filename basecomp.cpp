/* basecomp: get the base comp
 *
 * Copyright (C) 2023 Andrew D. Smith
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
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstdint>

using std::cout;
using std::cerr;
using std::tolower;

constexpr std::uint8_t nuc_to_idx[] = {
 /*  0*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 16*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 32*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 48*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 64*/  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 80*/  4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 96*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /*112*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /*128*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
};

int
main(int argc, const char **argv) {
  static constexpr auto buf_size = 64u * 1024u;

  if (argc != 2) {
    std::perror("basecomp <fasta-file>");
    return EXIT_FAILURE;
  }

  const auto genome_file = argv[1];

  std::array<std::uint_fast64_t, 256> counter;
  std::fill_n(begin(counter), 256, 0);

  int is_ok = EXIT_FAILURE;
  FILE *fp = std::fopen(genome_file, "r");
  if (!fp) {
    std::perror(genome_file);
    return is_ok;
  }

  std::array<uint8_t, 64 * 1024> buf;

  std::int_fast32_t n = 0;
  while ((n = std::fread(buf.data(), 1, buf_size, fp)) > 0) {
    auto c = buf.data();
    // not bother to exclude name lines that start with '>'
    while (n-- > 0) counter[*c]++;
  }
  if (std::ferror(fp)) {
    std::perror(genome_file);
    return is_ok;
  }

  std::fclose(fp);

  std::uint_fast32_t tot = 0;
  for (auto i : {'A', 'C', 'G', 'T'}) tot += counter[i] + counter[tolower(i)];

  if (tot == 0u) tot = 1u;

  for (auto i : {'A', 'C', 'G', 'T'}) {
    const double denom = counter[i] + counter[tolower(i)];
    cout << static_cast<char>(i) << '\t' << denom / tot << '\n';
  }
  return EXIT_SUCCESS;
}
