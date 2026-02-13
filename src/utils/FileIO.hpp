/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 // FileIO.hpp

#pragma once

#include <string>
#include <vector>
#include <functional>
#include <iosfwd>   // forward declare std::ostream

namespace utils::io {

// ---------- Binary reader ----------

// Reads a rows x cols float array from a binary file into `array`.
// (Assumes row-major layout in the file)
void readFloatBinaryfileData(const std::string& filename,
                             float* array,
                             int rows,
                             int cols);

// ---------- ASCII .dat table writer ----------

struct MetaLine {
    std::string key;
    std::string value;
};

// ISO-8601 UTC timestamp
std::string now_utc_iso8601();

// Fixed-width columns (all columns same width)
bool write_dat_table_fixed_width(
    const std::string& filename,
    const std::string& title,
    const std::vector<MetaLine>& meta,
    const std::vector<std::string>& column_names,
    int nrows,
    const std::function<void(int row, std::ostream& os, int colw)>& write_row,
    int col_width = 16,
    int precision = 10
);

// Per-column widths (widths.size() should match column_names.size())
bool write_dat_table_column_widths(
    const std::string& filename,
    const std::string& title,
    const std::vector<MetaLine>& meta,
    const std::vector<std::string>& column_names,
    const std::vector<int>& column_widths,
    int nrows,
    const std::function<void(int row, std::ostream& os, const std::vector<int>& w)>& write_row,
    int precision = 10
);

// Writes a header (comment lines) then lets the caller write the body.
// Does NOT set scientific/fixed/precision/widths — it leaves stream defaults alone.
bool write_dat_with_header(
    const std::string& filename,
    const std::string& title,
    const std::vector<MetaLine>& meta,
    const std::vector<std::string>& column_names,
    const std::function<void(std::ostream& os)>& write_body
);


} // namespace utils::io
