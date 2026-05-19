/*
 * Copyright (c) 2026 J. D. Nichols, University of Leicester, UK.
 * SPDX-License-Identifier: MIT
 *
 * This file is part of the Jupiter Auroral Ionosphere Code (JAIC) project.
 * Licensed under the MIT License. See the LICENSE file in the project root for full license information.
 */

/*
This file defines utility functions for working with std::vector arrays: copy,
copy_as, arange, linspace, logspace, and interp.
*/

// ArrayOps.hpp

#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace utils::array {

namespace detail {

    template <typename TOut, typename... Ts>
    struct output_type {
        using type = TOut;
    };

    template <typename... Ts>
    struct output_type<void, Ts...> {
        using type = std::common_type_t<Ts...>;
    };

    template <typename TOut, typename... Ts>
    using output_type_t = typename output_type<TOut, Ts...>::type;

    template <typename T>
    constexpr bool is_numeric_v = std::is_arithmetic_v<T>;

    template <typename T>
    constexpr bool is_real_v = std::is_floating_point_v<T>;

    template <typename T>
    bool is_strictly_increasing(const std::vector<T>& x) {
        for (std::size_t i = 1; i < x.size(); ++i) {
            if (!(x[i] > x[i - 1])) {
                return false;
            }
        }
        return true;
    }

    template <typename T>
    bool nearly_equal(T a, T b) {
        const T scale = std::max(static_cast<T>(1), std::max(std::abs(a), std::abs(b)));
        const T tol = static_cast<T>(100) * std::numeric_limits<T>::epsilon() * scale;
        return std::abs(a - b) <= tol;
    }

} // namespace detail

// -----------------------------------------------------------------------------
// Copy helpers
// -----------------------------------------------------------------------------

// By default, copy preserves the input value type:
//     std::vector<float>  b = copy(a_float);
//     std::vector<double> d = copy(a_double);
//
// To convert, specify the output type explicitly:
//     std::vector<double> d = copy<double>(a_float);
//     std::vector<float>  f = copy<float>(a_double);
template <typename TOut = void, typename TIn>
auto copy(const std::vector<TIn>& src) {
    static_assert(detail::is_numeric_v<TIn>, "copy: input type must be arithmetic");

    using ReturnT = detail::output_type_t<TOut, TIn>;
    static_assert(detail::is_numeric_v<ReturnT>, "copy: output type must be arithmetic");

    std::vector<ReturnT> dst;
    dst.reserve(src.size());

    for (const auto& value : src) {
        dst.push_back(static_cast<ReturnT>(value));
    }

    return dst;
}

template <typename TOut = void, typename TIn, std::size_t N>
auto copy(const TIn (&src)[N]) {
    static_assert(detail::is_numeric_v<TIn>, "copy: input type must be arithmetic");

    using ReturnT = detail::output_type_t<TOut, TIn>;
    static_assert(detail::is_numeric_v<ReturnT>, "copy: output type must be arithmetic");

    std::vector<ReturnT> dst;
    dst.reserve(N);

    for (std::size_t i = 0; i < N; ++i) {
        dst.push_back(static_cast<ReturnT>(src[i]));
    }

    return dst;
}

// Backwards-readable alias
template <typename TOut, typename TIn>
auto copy_as(const std::vector<TIn>& src) {
    return copy<TOut>(src);
}

template <typename TOut, typename TIn, std::size_t N>
auto copy_as(const TIn (&src)[N]) {
    return copy<TOut>(src);
}

// -----------------------------------------------------------------------------
// Grid generation
// -----------------------------------------------------------------------------

template <typename TOut = void, typename TStart, typename TStop, typename TStep>
auto arange(TStart start, TStop stop, TStep step) {
    static_assert(detail::is_numeric_v<TStart>, "arange: start must be arithmetic");
    static_assert(detail::is_numeric_v<TStop>,  "arange: stop must be arithmetic");
    static_assert(detail::is_numeric_v<TStep>,  "arange: step must be arithmetic");

    using WorkT   = std::common_type_t<TStart, TStop, TStep>;
    using ReturnT = detail::output_type_t<TOut, TStart, TStop, TStep>;

    static_assert(detail::is_numeric_v<ReturnT>, "arange: output type must be arithmetic");

    const WorkT start_w = static_cast<WorkT>(start);
    const WorkT stop_w  = static_cast<WorkT>(stop);
    const WorkT step_w  = static_cast<WorkT>(step);

    if (step_w == static_cast<WorkT>(0)) {
        throw std::invalid_argument("arange: step must be non-zero");
    }

    std::vector<ReturnT> out;

    if (step_w > static_cast<WorkT>(0)) {
        for (WorkT value = start_w; value < stop_w; value += step_w) {
            out.push_back(static_cast<ReturnT>(value));
        }
    } else {
        for (WorkT value = start_w; value > stop_w; value += step_w) {
            out.push_back(static_cast<ReturnT>(value));
        }
    }

    return out;
}

template <typename TOut = void, typename TStart, typename TEnd>
auto linspace(TStart start, TEnd end, int num) {
    static_assert(detail::is_real_v<TStart>, "linspace: start must be floating point");
    static_assert(detail::is_real_v<TEnd>,   "linspace: end must be floating point");

    using WorkT   = std::common_type_t<TStart, TEnd>;
    using ReturnT = detail::output_type_t<TOut, TStart, TEnd>;

    static_assert(detail::is_real_v<ReturnT>, "linspace: output type must be floating point");

    std::vector<ReturnT> out;

    if (num <= 0) {
        return out;
    }

    out.reserve(static_cast<std::size_t>(num));

    if (num == 1) {
        out.push_back(static_cast<ReturnT>(start));
        return out;
    }

    const WorkT start_w = static_cast<WorkT>(start);
    const WorkT end_w   = static_cast<WorkT>(end);
    const WorkT step_w  = (end_w - start_w) / static_cast<WorkT>(num - 1);

    for (int i = 0; i < num; ++i) {
        out.push_back(static_cast<ReturnT>(start_w + static_cast<WorkT>(i) * step_w));
    }

    return out;
}

template <typename TOut = void, typename TStart, typename TEnd>
auto logspace(TStart start, TEnd end, int num) {
    static_assert(detail::is_real_v<TStart>, "logspace: start must be floating point");
    static_assert(detail::is_real_v<TEnd>,   "logspace: end must be floating point");

    using WorkT   = std::common_type_t<TStart, TEnd>;
    using ReturnT = detail::output_type_t<TOut, TStart, TEnd>;

    static_assert(detail::is_real_v<ReturnT>, "logspace: output type must be floating point");

    const WorkT start_w = static_cast<WorkT>(start);
    const WorkT end_w   = static_cast<WorkT>(end);

    if (start_w <= static_cast<WorkT>(0) || end_w <= static_cast<WorkT>(0)) {
        throw std::invalid_argument("logspace: start and end must be positive");
    }

    std::vector<ReturnT> out;

    if (num <= 0) {
        return out;
    }

    out.reserve(static_cast<std::size_t>(num));

    if (num == 1) {
        out.push_back(static_cast<ReturnT>(start_w));
        return out;
    }

    const WorkT log_start = std::log10(start_w);
    const WorkT log_end   = std::log10(end_w);
    const WorkT step_w    = (log_end - log_start) / static_cast<WorkT>(num - 1);

    for (int i = 0; i < num; ++i) {
        const WorkT exponent = log_start + static_cast<WorkT>(i) * step_w;
        out.push_back(static_cast<ReturnT>(std::pow(static_cast<WorkT>(10), exponent)));
    }

    return out;
}

// -----------------------------------------------------------------------------
// Linear interpolation
// -----------------------------------------------------------------------------

namespace detail {

    template <typename TOut = void, typename TX, typename TY, typename TXNew>
    auto interp_vector_impl(
        const std::vector<TX>& x,
        const std::vector<TY>& y,
        const std::vector<TXNew>& x_new
    ) {
        static_assert(is_real_v<TX>,    "interp: x values must be floating point");
        static_assert(is_real_v<TY>,    "interp: y values must be floating point");
        static_assert(is_real_v<TXNew>, "interp: x_new values must be floating point");

        using WorkT   = std::common_type_t<TX, TY, TXNew>;
        using ReturnT = output_type_t<TOut, TX, TY, TXNew>;

        static_assert(is_real_v<ReturnT>, "interp: output type must be floating point");

        if (x.size() != y.size()) {
            throw std::invalid_argument("interp: x and y must have the same size");
        }

        if (x.empty()) {
            throw std::invalid_argument("interp: input vectors cannot be empty");
        }

        if (!is_strictly_increasing(x)) {
            throw std::invalid_argument("interp: x must be strictly increasing");
        }

        std::vector<ReturnT> y_new;
        y_new.reserve(x_new.size());

        for (const TXNew& xq_raw : x_new) {
            const WorkT xq = static_cast<WorkT>(xq_raw);

            if (xq <= static_cast<WorkT>(x.front())) {
                y_new.push_back(static_cast<ReturnT>(y.front()));
                continue;
            }

            if (xq >= static_cast<WorkT>(x.back())) {
                y_new.push_back(static_cast<ReturnT>(y.back()));
                continue;
            }

            for (std::size_t i = 1; i < x.size(); ++i) {
                if (xq < static_cast<WorkT>(x[i])) {
                    const WorkT x0 = static_cast<WorkT>(x[i - 1]);
                    const WorkT x1 = static_cast<WorkT>(x[i]);
                    const WorkT y0 = static_cast<WorkT>(y[i - 1]);
                    const WorkT y1 = static_cast<WorkT>(y[i]);

                    const WorkT t  = (xq - x0) / (x1 - x0);
                    const WorkT yq = (static_cast<WorkT>(1) - t) * y0 + t * y1;

                    y_new.push_back(static_cast<ReturnT>(yq));
                    break;
                }
            }
        }

        return y_new;
    }

} // namespace detail

// Simple linear interpolation. x must be sorted in strictly increasing order.
// Values outside the input domain are clamped to the first/last y value
//
// By default the output type is std::common_type_t<TX, TY, TXNew>.
// Use interp<double>(x, y, x_new) or interp<float>(x, y, x_new) to force it.
template <typename TOut = void, typename TX, typename TY, typename TXNew>
auto interp(
    const std::vector<TX>& x,
    const std::vector<TY>& y,
    const std::vector<TXNew>& x_new
) {
    return detail::interp_vector_impl<TOut, TX, TY, TXNew>(x, y, x_new);
}

// Convenience overload for a single interpolation point.
template <
    typename TOut = void,
    typename TX,
    typename TY,
    typename TXNew,
    std::enable_if_t<detail::is_real_v<TXNew>, int> = 0
>
auto interp(
    const std::vector<TX>& x,
    const std::vector<TY>& y,
    TXNew x_new
) {
    std::vector<TXNew> x_query{ x_new };
    auto y_query = detail::interp_vector_impl<TOut, TX, TY, TXNew>(x, y, x_query);
    return y_query.front();
}

// -----------------------------------------------------------------------------
// Conservative rebinning (preserves total)
// -----------------------------------------------------------------------------

// Rebins piecewise-constant bin totals y_old defined on old_edges onto new_edges.
// y_old.size() must be old_edges.size() - 1.
// The input/output outer edges must match (within floating-point tolerance).
// Each output bin receives overlap-weighted contributions so sum(y_new) == sum(y_old).
template <typename TOut = void, typename TEdgeOld, typename TYOld, typename TEdgeNew>
auto rebin_conservative(
    const std::vector<TEdgeOld>& old_edges,
    const std::vector<TYOld>& y_old,
    const std::vector<TEdgeNew>& new_edges
) {
    static_assert(detail::is_real_v<TEdgeOld>, "rebin_conservative: old_edges must be floating point");
    static_assert(detail::is_real_v<TEdgeNew>, "rebin_conservative: new_edges must be floating point");
    static_assert(detail::is_numeric_v<TYOld>, "rebin_conservative: y_old must be arithmetic");

    using WorkT   = std::common_type_t<TEdgeOld, TYOld, TEdgeNew>;
    using ReturnT = detail::output_type_t<TOut, TEdgeOld, TYOld, TEdgeNew>;

    static_assert(detail::is_real_v<ReturnT>, "rebin_conservative: output type must be floating point");

    if (old_edges.size() < 2 || new_edges.size() < 2) {
        throw std::invalid_argument("rebin_conservative: edge vectors must have at least 2 elements");
    }

    if (y_old.size() + 1 != old_edges.size()) {
        throw std::invalid_argument("rebin_conservative: y_old.size() must be old_edges.size() - 1");
    }

    if (!detail::is_strictly_increasing(old_edges)) {
        throw std::invalid_argument("rebin_conservative: old_edges must be strictly increasing");
    }

    if (!detail::is_strictly_increasing(new_edges)) {
        throw std::invalid_argument("rebin_conservative: new_edges must be strictly increasing");
    }

    const WorkT old_lo = static_cast<WorkT>(old_edges.front());
    const WorkT old_hi = static_cast<WorkT>(old_edges.back());
    const WorkT new_lo = static_cast<WorkT>(new_edges.front());
    const WorkT new_hi = static_cast<WorkT>(new_edges.back());

    if (!detail::nearly_equal(old_lo, new_lo) || !detail::nearly_equal(old_hi, new_hi)) {
        throw std::invalid_argument("rebin_conservative: old and new outer edges must match");
    }

    std::vector<ReturnT> y_new(new_edges.size() - 1, static_cast<ReturnT>(0));

    std::size_t i_old = 0;

    for (std::size_t i_new = 0; i_new + 1 < new_edges.size(); ++i_new) {
        const WorkT nl = static_cast<WorkT>(new_edges[i_new]);
        const WorkT nr = static_cast<WorkT>(new_edges[i_new + 1]);

        while (i_old + 1 < old_edges.size() && static_cast<WorkT>(old_edges[i_old + 1]) <= nl) {
            ++i_old;
        }

        std::size_t k = i_old;
        WorkT accum = static_cast<WorkT>(0);

        while (k + 1 < old_edges.size()) {
            const WorkT ol = static_cast<WorkT>(old_edges[k]);
            const WorkT orr = static_cast<WorkT>(old_edges[k + 1]);

            if (ol >= nr) {
                break;
            }

            const WorkT overlap_l = std::max(ol, nl);
            const WorkT overlap_r = std::min(orr, nr);
            const WorkT overlap = overlap_r - overlap_l;

            if (overlap > static_cast<WorkT>(0)) {
                const WorkT old_width = orr - ol;
                const WorkT fraction = overlap / old_width;
                accum += static_cast<WorkT>(y_old[k]) * fraction;
            }

            if (orr <= nr) {
                ++k;
            } else {
                break;
            }
        }

        y_new[i_new] = static_cast<ReturnT>(accum);
    }

    return y_new;
}









} // namespace utils::array

