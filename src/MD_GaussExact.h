/// MD DUST: exact maximum for the Gaussian model, k >= 3 (method "exact")

#ifndef MD_GAUSS_EXACT_H
#define MD_GAUSS_EXACT_H

#include <functional>
#include <numeric>

#include "MD_Decision.h"

namespace dust_md {

////////////////////////////////////////////////////////////////////////////////
/// GAUSS: D(x) = D(0) + h'x - x'Gx/2, h = -M'S - u, G = M'M
/// on each face: G_II x_I = h_I + KKT conditions
enum class GaussMaximum { finite, unbounded, unresolved };

struct GaussMaximumResult
{
  GaussMaximum kind = GaussMaximum::unresolved;
  std::vector<double> point;
  double value = -std::numeric_limits<double>::infinity();
};

inline GaussMaximumResult gauss_joint_maximum(const Decision<GaussPolicy>& test)
{
  const size_t p = test.constraints;
  const size_t d = test.dimension;
  const double eps = std::numeric_limits<double>::epsilon();
  std::vector<double> h(p), gram(p * p);
  for (size_t j = 0; j < p; ++j)
  {
    h[j] = -test.u[j];
    for (size_t row = 0; row < d; ++row)
      h[j] -= test.a[row] * test.matrix[j * d + row];
    for (size_t k = 0; k <= j; ++k)
    {
      double dot = 0.0;
      for (size_t row = 0; row < d; ++row)
        dot += test.matrix[j * d + row] * test.matrix[k * d + row];
      gram[j * p + k] = gram[k * p + j] = dot;
    }
  }

  GaussMaximumResult result;
  std::vector<size_t> face;
  std::vector<double> unconstrained;
  std::vector<double> face_point, face_gradient;
  bool face_solved = false;
  std::function<bool(size_t)> visit = [&](size_t next) -> bool {
    face_solved = false;
    const size_t k = face.size();
    std::vector<double> q(d * k, 0.0), r(k * k, 0.0), x(p, 0.0);
    for (size_t col = 0; col < k; ++col)
    {
      const size_t j = face[col];
      double column_norm2 = 0.0;
      for (size_t row = 0; row < d; ++row)
      {
        q[col * d + row] = test.matrix[j * d + row];
        column_norm2 += q[col * d + row] * q[col * d + row];
      }
      // reorthogonalization
      for (unsigned int pass = 0; pass < 2; ++pass)
        for (size_t previous = 0; previous < col; ++previous)
        {
          double projection = 0.0;
          for (size_t row = 0; row < d; ++row)
            projection += q[previous * d + row] * q[col * d + row];
          r[previous * k + col] += projection;
          for (size_t row = 0; row < d; ++row)
            q[col * d + row] -= projection * q[previous * d + row];
        }
      double remainder2 = 0.0;
      for (size_t row = 0; row < d; ++row)
        remainder2 += q[col * d + row] * q[col * d + row];
      const double diagonal = std::sqrt(remainder2);
      const double rank_guard = 64.0 * eps *
        static_cast<double>(std::max(d, p)) * std::sqrt(column_norm2);
      if (!std::isfinite(diagonal) || diagonal <= rank_guard)
      {
        // dependent face: nonnegative direction where D increases?
        if (col + 1 == k)
        {
          std::vector<double> ray(p, 0.0);
          ray[j] = 1.0;
          for (size_t back = col; back-- > 0; )
          {
            double coefficient = r[back * k + col];
            for (size_t later = back + 1; later < col; ++later)
              coefficient += r[back * k + later] * ray[face[later]];
            ray[face[back]] = -coefficient / r[back * k + back];
          }
          bool nonnegative = true;
          for (size_t index : face)
          {
            if (ray[index] < -64.0 * eps) nonnegative = false;
            if (ray[index] < 0.0) ray[index] = 0.0;
          }
          double slope = 0.0, slope_scale = 0.0;
          for (size_t index : face)
          {
            slope += h[index] * ray[index];
            slope_scale += std::abs(h[index] * ray[index]);
          }
          if (nonnegative && dustib::positive(slope, slope_scale))
          {
            double scale = 0.0;
            const double at_zero = test.value(x, scale);
            double multiplier = (std::abs(at_zero) + scale + 1.0) / slope;
            if (!std::isfinite(multiplier) || multiplier <= 0.0)
              multiplier = 1.0;
            for (unsigned int attempt = 0; attempt < 64; ++attempt)
            {
              for (size_t index : face) x[index] = multiplier * ray[index];
              double witness_scale = 0.0;
              const double witness = test.value(x, witness_scale);
              if (dustib::positive(witness, witness_scale))
              {
                result = {GaussMaximum::unbounded, x, witness};
                return true;
              }
              if (!std::isfinite(multiplier * 2.0)) break;
              multiplier *= 2.0;
            }
          }
        }
        return false;
      }
      r[col * k + col] = diagonal;
      for (size_t row = 0; row < d; ++row)
        q[col * d + row] /= diagonal;
    }

    // R'R x_I = h_I
    std::vector<double> intermediate(k, 0.0);
    for (size_t col = 0; col < k; ++col)
    {
      double rhs = h[face[col]];
      for (size_t previous = 0; previous < col; ++previous)
        rhs -= r[previous * k + col] * intermediate[previous];
      intermediate[col] = rhs / r[col * k + col];
    }
    bool feasible = true;
    for (size_t col = k; col-- > 0; )
    {
      double rhs = intermediate[col];
      for (size_t later = col + 1; later < k; ++later)
        rhs -= r[col * k + later] * x[face[later]];
      x[face[col]] = rhs / r[col * k + col];
      if (!std::isfinite(x[face[col]]) || x[face[col]] < 0.0)
        feasible = false;
    }
    if (p > 3)
    {
      if (k == p) unconstrained = x;
      face_point = x;
      face_gradient = h;
      for (size_t j = 0; j < p; ++j)
        for (size_t index : face)
          face_gradient[j] -= gram[j * p + index] * x[index];
      face_solved = true;
    }
    if (feasible)
    {
      for (size_t j = 0; j < p; ++j)
      {
        double gradient = p > 3 ? face_gradient[j] : h[j];
        double gradient_scale = std::abs(h[j]);
        for (size_t index : face)
        {
          const double term = gram[j * p + index] * x[index];
          if (p <= 3) gradient -= term;
          gradient_scale += std::abs(term);
        }
        const double tolerance = 128.0 * eps * (1.0 + gradient_scale);
        const bool active = std::find(face.begin(), face.end(), j) != face.end();
        if (!std::isfinite(gradient) ||
            (active ? std::abs(gradient) > tolerance :
                      gradient > tolerance))
        {
          feasible = false;
          break;
        }
      }
    }
    if (feasible)
    {
      double scale = 0.0;
      const double maximum = test.value(x, scale);
      if (std::isfinite(maximum))
      {
        result = {GaussMaximum::finite, x, maximum};
        return true;
      }
    }
    for (size_t j = next; j < p; ++j)
    {
      face.push_back(j);
      const bool found = visit(j + 1);
      face.pop_back();
      if (found) return true;
    }
    return false;
  };
  // full face first
  if (p > 1)
  {
    face.resize(p);
    std::iota(face.begin(), face.end(), 0);
    if (visit(p)) return result;
    face.clear();
    // active set from the unconstrained solution, then all the faces
    if (p > 3 && unconstrained.size() == p)
    {
      for (size_t j = 0; j < p; ++j)
        if (std::isfinite(unconstrained[j]) && unconstrained[j] > 0.0)
          face.push_back(j);
      std::vector<std::vector<size_t>> seen;
      for (size_t step = 0; step < 3 * p + 3; ++step)
      {
        if (std::find(seen.begin(), seen.end(), face) != seen.end()) break;
        seen.push_back(face);
        if (visit(p)) return result;
        if (!face_solved) break;
        size_t remove = p;
        double most_negative = 0.0;
        for (size_t j : face)
          if (face_point[j] < most_negative)
          {
            most_negative = face_point[j];
            remove = j;
          }
        if (remove < p)
        {
          face.erase(std::find(face.begin(), face.end(), remove));
          continue;
        }
        size_t add = p;
        double most_positive = 0.0;
        for (size_t j = 0; j < p; ++j)
          if (std::find(face.begin(), face.end(), j) == face.end() &&
              face_gradient[j] > most_positive)
          {
            most_positive = face_gradient[j];
            add = j;
          }
        if (add == p) break;
        face.insert(std::lower_bound(face.begin(), face.end(), add), add);
      }
      face.clear();
    }
  }
  visit(0);
  return result;
}

inline bool gauss_exact_search(const Decision<GaussPolicy>& test)
{
  const auto maximum = gauss_joint_maximum(test);
  return maximum.kind != GaussMaximum::unresolved &&
    test.positive(maximum.point);
}

} // namespace dust_md

#endif
