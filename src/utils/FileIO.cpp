/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

 // FileIO.cpp


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <chrono>
#include <ctime>

using namespace std;

#include "FileIO.hpp"

namespace utils::io {


    // ---------- Binary reader ----------

    void readFloatBinaryfileData(const string& filename, 
                                        float* array, 
                                        const int rows, 
                                        const int cols) {
        ifstream file(filename, ios::binary);
        if (!file.is_open()) {
        cerr << "Unable to open file: " << filename << endl;
        exit(1);
        }
        // cout << "Reading binary floating point data from file: " << filename << endl;
        file.read(reinterpret_cast<char*>(array), rows * cols * sizeof(float));
        
        file.close();
    }

    // ---------- ASCII .dat table writer ----------

    std::string now_utc_iso8601() {
        using namespace std::chrono;
        auto t = system_clock::to_time_t(system_clock::now());
        std::tm tm{};
    #if defined(_WIN32)
        gmtime_s(&tm, &t);
    #else
        gmtime_r(&t, &tm);
    #endif
        std::ostringstream ss;
        ss << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
        return ss.str();
    }

    static std::string join_cols(const std::vector<std::string>& cols, const char* sep = " ") {
        std::ostringstream ss;
        for (size_t i = 0; i < cols.size(); ++i) {
            if (i) ss << sep;
            ss << cols[i];
        }
        return ss.str();
    }

    static void write_columns_line_fixed(std::ostream& os,
                                        const std::vector<std::string>& cols,
                                        int colw) {
        os << "# ";
        for (const auto& c : cols) {
            os << " " << std::left << std::setw(colw) << c;
        }
        os << "\n";
    }

    static void write_columns_line_per(std::ostream& os,
                                    const std::vector<std::string>& cols,
                                    const std::vector<int>& w) {
        os << "# ";
        for (size_t i = 0; i < cols.size(); ++i) {
            const int wi = (i < w.size() ? w[i] : 16);
            os << " " << std::left << std::setw(wi) << cols[i];
        }
        os << "\n";
    }

    bool write_dat_table_fixed_width(
        const std::string& filename,
        const std::string& title,
        const std::vector<MetaLine>& meta,
        const std::vector<std::string>& column_names,
        int nrows,
        const std::function<void(int row, std::ostream& os, int colw)>& write_row,
        int col_width,
        int precision
    ) {
        try {
            std::filesystem::path p(filename);
            if (p.has_parent_path()) std::filesystem::create_directories(p.parent_path());

            std::ofstream f(filename, std::ios::trunc);
            if (!f) {
                std::cerr << "ERROR: unable to open file for writing: " << filename << "\n";
                return false;
            }

            f << "# " << title << "\n";
            f << "# generated_utc: " << now_utc_iso8601() << "\n";
            for (const auto& kv : meta) f << "# " << kv.key << ": " << kv.value << "\n";
            if (!column_names.empty()) write_columns_line_fixed(f, column_names, col_width);

            // Numeric format
            f.setf(std::ios::scientific);
            f << std::setprecision(precision);

            for (int r = 0; r < nrows; ++r) {
                write_row(r, f, col_width);
                f << "\n";
            }
            return true;

        } catch (const std::exception& e) {
            std::cerr << "ERROR: exception while writing " << filename << ": " << e.what() << "\n";
            return false;
        }
    }

    bool write_dat_table_column_widths(
        const std::string& filename,
        const std::string& title,
        const std::vector<MetaLine>& meta,
        const std::vector<std::string>& column_names,
        const std::vector<int>& column_widths,
        int nrows,
        const std::function<void(int row, std::ostream& os, const std::vector<int>& w)>& write_row,
        int precision
    ) {
        try {
            std::filesystem::path p(filename);
            if (p.has_parent_path()) std::filesystem::create_directories(p.parent_path());

            std::ofstream f(filename, std::ios::trunc);
            if (!f) {
                std::cerr << "ERROR: unable to open file for writing: " << filename << "\n";
                return false;
            }

            f << "# " << title << "\n";
            f << "# generated_utc: " << now_utc_iso8601() << "\n";
            for (const auto& kv : meta) f << "# " << kv.key << ": " << kv.value << "\n";
            if (!column_names.empty()) write_columns_line_per(f, column_names, column_widths);

            f.setf(std::ios::scientific);
            f << std::setprecision(precision);

            for (int r = 0; r < nrows; ++r) {
                write_row(r, f, column_widths);
                f << "\n";
            }
            return true;

        } catch (const std::exception& e) {
            std::cerr << "ERROR: exception while writing " << filename << ": " << e.what() << "\n";
            return false;
        }
    }

    bool write_dat_with_header(
    const std::string& filename,
    const std::string& title,
    const std::vector<MetaLine>& meta,
    const std::vector<std::string>& column_names,
    const std::function<void(std::ostream& os)>& write_body
) {
    try {
        std::filesystem::path p(filename);
        if (p.has_parent_path()) std::filesystem::create_directories(p.parent_path());

        std::ofstream f(filename, std::ios::trunc);
        if (!f) {
            std::cerr << "ERROR: unable to open file for writing: " << filename << "\n";
            return false;
        }

        // Header
        f << "# " << title << "\n";
        f << "# generated_utc: " << now_utc_iso8601() << "\n";
        for (const auto& kv : meta) {
            f << "# " << kv.key << ": " << kv.value << "\n";
        }
        if (!column_names.empty()) {
            f << "# cols:" << join_cols(column_names) << "\n";
        }

        // Body (caller controls everything; defaults preserved)
        write_body(f);

        return true;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: exception while writing " << filename << ": " << e.what() << "\n";
        return false;
    }
}


} // namespace utils::io

