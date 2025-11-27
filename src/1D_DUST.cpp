#include <Rcpp.h>
#include <cmath>
#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL
#include "1D_DUST.h"
#include "preProcessing.h"

#include <fstream>
#include <iostream>

using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// --- // Constructor // --- //
DUST_1D::DUST_1D(std::string dualmax_algo,
                 std::string constr_index,
                 Nullable<int> nbLoops)
  : dualmax_algo(dualmax_algo),
    constr_index(constr_index),
    indices(nullptr),
    n(0)
{
  if(nbLoops.isNull()){nb_Loops = 10;}else{nb_Loops = as<int>(nbLoops);}
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

DUST_1D::~DUST_1D()
{
  delete indices;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void DUST_1D::pruning_method()
{
  delete indices;

  /// /// ///
  /// /// /// index METHOD
  /// /// ///
  if(constr_index == "rand"){indices = new Indices_1D_Rand();}
  else{indices = new Indices_1D_Det;}

  /// /// ///
  /// /// /// dualmax_algo METHOD
  /// /// ///
  if(dualmax_algo == "DUSTr"){current_test = &DUST_1D::dualMaxAlgo0;}
  if(dualmax_algo == "DUST"){current_test = &DUST_1D::dualMaxAlgo1;}
  if(dualmax_algo == "DUSTgs"){current_test = &DUST_1D::dualMaxAlgo2;}
  if(dualmax_algo == "DUSTbs"){current_test = &DUST_1D::dualMaxAlgo3;}
  if(dualmax_algo == "DUSTqn"){current_test = &DUST_1D::dualMaxAlgo4;}
  if(dualmax_algo == "PELT"){current_test = &DUST_1D::dualMaxAlgo5;}
  if(dualmax_algo == "OP"){current_test = &DUST_1D::dualMaxAlgo6;}

  /// /// ///
  /// /// /// INIT RANDOM GENERATOR
  /// /// ///
  engine.seed(std::random_device{}());
  dist = std::uniform_real_distribution<double>(0.0, 1.0);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// --- // append_data, i. e. initializes all data-dependent vectors // --- //
// --- // append_data, i. e. initializes all data-dependent vectors // --- //
// --- // append_data, i. e. initializes all data-dependent vectors // --- //

void DUST_1D::append_data(std::vector<double>& inData,
                       Nullable<double> inPenalty)
{
  bool first_execution = (n == 0);
  n += inData.size();   // update total size

  // Reserve memory for all vectors
  cumsum.reserve(n + 1);
  chptRecord.reserve(n + 1);
  costRecord.reserve(n + 1);
  nb_indices.reserve(n + 1);

  // On first execution: initialize all vectors
  if (first_execution)
  {
    if (inPenalty.isNull()) { penalty = 2 * std::log(n); } else { penalty = as<double>(inPenalty); }

    cumsum.push_back(0);
    costRecord.push_back(-penalty);
    chptRecord.push_back(0);
    nb_indices.push_back(1);

    pruning_method();
    indices->add_first(0);
  }

  // store data by updating the cumsum
  for (auto current_datum = inData.begin(); current_datum != inData.end(); ++current_datum)
    cumsum.push_back(cumsum.back() + statistic(*current_datum));
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void DUST_1D::update_partition()
{
  // Initialize OP step value
  double lastCost;
  unsigned int nbt = nb_indices.back();
  for (unsigned t = indices->get_first() + 1; t <= n; t++)
  {
    ///////////// OP step /////////////
    ///////////// OP step /////////////
    indices->reset();
    double minCost = std::numeric_limits<double>::infinity();
    unsigned int argMin = 0;
    do
    {
      unsigned int s = indices->get_current();
      lastCost = costRecord[s] + Cost(t, s); // without the penalty beta
      if (lastCost < minCost)
      {
        minCost = lastCost;
        argMin = s;
      }
      indices->next();
    }
    while(indices->is_not_the_last());
    //////// END (OP step) ////////
    //////// END (OP step) ////////

    // OP update minCost and save values
    // minCost = Q_t. Here + beta to get Q_t = Q_i + C(y_it) + beta
    minCost += penalty;
    costRecord.push_back(minCost);
    chptRecord.push_back(argMin);

    ///////////// DUST step /////////////
    ///////////// DUST step /////////////

    indices->reset_pruning();

    //////// DUST loop
    //////// DUST loop
    // we use :
    // is_not_the_last_pruning
    // current_test
    // prune_current
    // next_pruning
    // prune_last
    while (indices->is_not_the_last_pruning()) // is true, while we are not on the smallest index
    {
      // valid if we prune the index get_current using index in get_constraint
      if ((this->*current_test)(minCost, t, indices->get_current(), indices->get_constraint()))
      {
        // remove the pruned index and its pointer
        // removing the elements increments the cursors i and pointersCurrent, while before stands still
        indices->prune_current();
        nbt--;
      }
      else
      {
        // increment all cursors
        indices->next_pruning();
      }
    }
    //////// END (DUST loop)
    //////// END (DUST loop)

    // Prune the last index (analogous with a "mu* = 0" duality simple test)
    // this is the smallest available index = PELT RULE
    if (lastCost > minCost)
    {
      indices->prune_last();
      nbt--;
    }

    ///////////// Update to next index /////////////
    ///////////// Update to next index /////////////
    indices->add_first(t);
    nb_indices.push_back(nbt);
    nbt++;
  }
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// --- // Wrapper method : all included  // --- //
List DUST_1D::dust(std::vector<double>& inData, Nullable<double> inPenalty)
{
  append_data(inData, inPenalty);
  update_partition();
  return get_partition();
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// --- // Builds changepoints // --- //
// --- // Builds changepoints // --- //
// --- // Builds changepoints // --- //
std::forward_list<unsigned int> DUST_1D::backtrack_changepoints()
{
  std::forward_list<unsigned int> changepoints {n};
  for (int tau = chptRecord[n]; tau != 0; tau = chptRecord[tau])
  {
    changepoints.push_front(tau);
  }
  return changepoints;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// --- // optimal partition info // --- //
// --- // optimal partition info // --- //
// --- // optimal partition info // --- //
List DUST_1D::get_partition()
{
  std::forward_list<unsigned int> chpts = backtrack_changepoints();

  return List::create(
    _["changepoints"] = chpts,
    _["lastIndexSet"] = indices->get_list(),
    _["nb"] = std::vector<unsigned>(nb_indices.begin() + 1, nb_indices.end()),
    _["costQ"] = std::vector<double>(costRecord.begin() + 1, costRecord.end())
  );
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// --- // get_info object // --- //
// --- // get_info object // --- //
// --- // get_info object // --- //
List DUST_1D::get_info()
{
  return List::create(
    _["data_statistic"] = cumsum,
    _["data_length"] = n,
    _["current_penalty"] = penalty,
    _["model"] = get_model(),
    _["pruning_algo"] = dualmax_algo,
    _["pruning_constr_index"] = constr_index,
    _["pruning_nb_loops"] = nb_Loops
  );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// SPECIAL CASE BEFORE DUAL EVALUATION / OPTIMIZATION

bool DUST_1D::isOnePointOrLinear(double a, double b)
{
  if(isBoundary(a) == true){return true;}
  if(a == b){return true;}
  return false;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// if special_cases == true => return true;

bool DUST_1D::specialCasePruning(double a,
                                 double b,
                                 double c,
                                 double d,
                                 double mu_max)
{
  ///
  /// case "objectiveMean = a = left bound"
  ///
  if(isBoundary(a) == true)
  {
    if (-c > 0) {return true;} /// test in mu = left bound
    if(isBoundary(b) == false){return false;} /// case mu_max = 0
    else /// case linear
    {
      if(-(c - mu_max * d)  > 0){return true;} /// right boun.  case bern / binom : mu_max != 1
      return false;
    }
  }

  ///
  /// case "a = b", here not on the boundary
  ///
  if(a == b)
  {
    if ( -c - Dstar(a) > 0){return true;} /// left bound
    if (-(c - mu_max * d) - Dstar(a) > 0){return true;} /// right bound
  }
  return false;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// DUAL = -(c - mu * d)  - (1- mu) * Dstar((a - mu * b) / (1 - mu));

// a = (cumsum[t] - cumsum[s]) / (t - s)
// b = (cumsum[s] - cumsum[r]) / (s - r)

// c = (Qt - Qs) / (t - s) =  (minCost - costRecord[s]) / (t - s);
// d = (Qs - Qr) / (s - r)  = (costRecord[s] - costRecord[r]) / (s - r);

/// RANDOM EVAL


bool DUST_1D::dualMaxAlgo0(double minCost,
                           unsigned int t,
                           unsigned int s,
                           unsigned int r)
{
  return (dualEval(dist(engine), minCost, t, s, r) > 0);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


/// EXACT EVAL

bool DUST_1D::dualMaxAlgo1(double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  double a = (cumsum[t] - cumsum[s]) / (t - s);
  double b = (cumsum[s] - cumsum[r]) / (s - r);
  double c = (minCost - costRecord[s]) / (t - s);
  double d = (costRecord[s] - costRecord[r]) / (s - r);
  double mu_max = muMax(a, b);

  ///
  /// case dual domain = one point OR dual = linear function
  ///
  if(isOnePointOrLinear(a,b) == true){return specialCasePruning(a,b,c,d,mu_max);}
  return (dualMax(minCost, t, s, r) > 0);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// Golden-section search

bool DUST_1D::dualMaxAlgo2(double minCost,
                           unsigned int t,
                           unsigned int s,
                           unsigned int r)
{
  double a = (cumsum[t] - cumsum[s]) / (t - s);
  double b = (cumsum[s] - cumsum[r]) / (s - r);

  double lt = 0.0;
  double rt = muMax(a, b);

  double cs1 = (1 - 1/phi) * rt;
  double cs2 = rt/phi;

  double d = (costRecord[s] - costRecord[r]) / (s - r);
  double c = (minCost - costRecord[s]) / (t - s);

  double fl =  - (c - cs1*d) - (1.0 - cs1) * Dstar((a - cs1*b)/(1 - cs1));
  double fr =  - (c - cs2*d) - (1.0 - cs2) * Dstar((a - cs2*b)/(1 - cs2));
  if(fl > 0 || fr > 0){return true;}
  double max_val = std::max(fl, fr);

  for (int i = 0; i < nb_Loops; i++)
  {
    if (fl > fr)
    {
      rt = cs2;
      cs2 = cs1;
      fr = fl;
      cs1 = rt - (rt - lt) / phi;
      fl =  - (c - cs1*d) - (1.0 - cs1) * Dstar((a - cs1*b)/(1 - cs1));
    }
    else
    {
      lt = cs1;
      cs1 = cs2;
      fl = fr;
      cs2 = lt + (rt - lt) / phi;
      fr = - (c - cs2*d) - (1.0 - cs2) * Dstar((a - cs2*b)/(1 - cs2));
    }
    max_val = std::max(max_val, std::max(fl, fr));
    if(max_val > 0){return true;}
  }
  return false;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// binary search. At each step, we evaluate the tangent line to the current point
// at its max to stop the search at early step (when possible)

bool DUST_1D::dualMaxAlgo3(double minCost,
                           unsigned int t,
                           unsigned int s,
                           unsigned int r)
{
  double a = (cumsum[t] - cumsum[s]) / (t - s);
  double b = (cumsum[s] - cumsum[r]) / (s - r);
  double c = (minCost - costRecord[s]) / (t - s);
  double d = (costRecord[s] - costRecord[r]) / (s - r);

  const double mu_max = muMax(a, b);

  ///
  /// case dual domain = one point OR dual = linear function
  ///
  if(isOnePointOrLinear(a,b) == true){return specialCasePruning(a,b,c,d,mu_max);}

  ///
  /// PELT test (dual at mu = 0)
  ///
  const double Dstar_a = Dstar(a);
  double test_value = - Dstar_a - c;  // dual(0)
  if (test_value > 0.0) {return true;}  // We already found a positive dual value at mu = 0

  // if non pruning at mu = 0
  // AND
  // If the tangent at 0 is already ≤ 0 on [0, mu_max], we can prune everything
  const double meanGap = a - b;
  double grad = Dstar_a - meanGap * DstarPrime(a) + d;
  if (test_value + mu_max * grad <= 0) {return false;} //check value in mu = mu_max

  /// iterative algo
  /// iterative algo
  /// iterative algo
  ///
  /// Iterative derivative-aided bisection on [0, mu_max]
  ///
  double lt = 0.0;
  double rt = mu_max;
  double mu = 0.5 * (lt + rt);

  const double grad_eps = 1e-12;  // tolerance for "zero" gradient

  for (int i = 0; i < nb_Loops; ++i)
  {
    const double one_minus_mu = 1.0 - mu;
    const double R = (a - mu * b) / one_minus_mu;
    const double Dstar_R      = Dstar(R);
    const double DstarPrime_R = DstarPrime(R);

    // value of the dual at mu
    test_value = - (c - mu * d) - one_minus_mu * Dstar_R;
    if (test_value > 0.0) {return true;}

    // derivative of f at mu
    grad = Dstar_R - (a - b) / one_minus_mu * DstarPrime_R + d;

    // If gradient is ~0 and value ≤ 0, for a concave f this implies f ≤ 0 everywhere
    if (std::abs(grad) < grad_eps) {
      return false;
    }

    if (grad > 0.0)
    {
      // Candidate max can only be on the right side for a concave f
      if (test_value + (rt - mu) * grad <= 0.0) {
        // Upper bound on [mu, rt] is ≤ 0 -> nothing positive to the right
        return false;
      }
      lt = mu;
    }
    else // grad < 0
    {
      // Candidate max can only be on the left side for a concave f
      if (test_value + (lt - mu) * grad <= 0.0) {
        // Upper bound on [lt, mu] is ≤ 0 -> nothing positive to the left
        return false;
      }
      rt = mu;
    }

    // New midpoint of the remaining interval
    mu = 0.5 * (lt + rt);
  }
  return false;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

bool DUST_1D::dualMaxAlgo4(double minCost,
                           unsigned int t,
                           unsigned int s,
                           unsigned int r)
{
  // Basic sanity: if these don't hold, something is structurally wrong.
  if (!(r < s && s < t)) {
    return false; // or assert(false);
  }

  // --- Small numerical constants (local to this function) ---
  constexpr double kEps            = 1e-12;
  constexpr double kMuMargin       = 1e-9;   // keep μ away from {0, μ_max}
  constexpr double kGradTol        = 1e-10;  // early stop if |grad| is tiny
  constexpr int    kMaxBacktracking = 10;

  const double inv_ts = 1.0 / static_cast<double>(t - s);
  const double inv_sr = 1.0 / static_cast<double>(s - r);

  // --- Notation ---
  // a = mean on [s+1, ..., t]
  // b = mean on [r+1, ..., s]
  // c = (minCost - Q_s) / (t - s)
  // d = (Q_s - Q_r) / (s - r)
  const double a = (cumsum[t] - cumsum[s]) * inv_ts;
  const double b = (cumsum[s] - cumsum[r]) * inv_sr;
  const double c = (minCost - costRecord[s]) * inv_ts;
  const double d = (costRecord[s] - costRecord[r]) * inv_sr;

  // 1) PELT test at μ = 0
  double nonLinear = Dstar(a);
  double test_value = -nonLinear - c;   // = -D*(a) - c

  if (test_value > 0.0) {
    // PELT already prunes at μ = 0
    return true;
  }

  // 2) μ_max: if zero (or negative), DUST does nothing beyond PELT
  const double raw_mu_max = muMax(a, b);
  if (raw_mu_max <= 0.0) {
    return false;
  }

  // Keep μ in [μ_min, μ_max_eff] with a small margin
  const double mu_min     = kMuMargin;
  const double mu_max_eff = std::max(raw_mu_max - kMuMargin, mu_min + kMuMargin);
  if (mu_max_eff <= mu_min) {
    // numerical degeneracy: effectively no room for μ
    return false;
  }

  // 3) Gradient at μ = 0
  const double meanGap = a - b;
  double grad = nonLinear - meanGap * DstarPrime(a) + d;

  // If grad < 0, maximum is at μ = 0 → PELT was already tested
  if (grad < 0.0) {
    return false;
  }

  // Quick concavity check using tangent from μ=0 to μ_max
  if (test_value + mu_max_eff * grad <= 0.0) {
    return false;
  }

  // --- Quasi-Newton state ---
  double mu = 0.0;            // current μ
  double direction = 0.0;     // search direction
  double mu_diff = 0.0;       // step in μ
  double m_value = a;         // m(μ)
  double grad_diff = -grad;   // will store g_{k+1} - g_k
  double inverseHessian = -1.0; // BFGS inverse Hessian in 1D

  // --- Helpers ---

  // m(μ) = (a - μ b) / (1 - μ)
  auto getM = [&](double mu_val) -> double {
    double denom = 1.0 - mu_val;
    if (std::fabs(denom) < kEps) {
      denom = (denom >= 0.0 ? kEps : -kEps);
    }
    return (a - mu_val * b) / denom;
  };

  // Ensure μ stays in [μ_min, μ_max_eff]
  auto clampMu = [&](double & mu_val) {
    if (mu_val < mu_min)       mu_val = mu_min;
    else if (mu_val > mu_max_eff) mu_val = mu_max_eff;
  };

  // Update direction and immediately clamp where we land
  auto updateDirection = [&]() {
    direction = -inverseHessian * grad;
    double proposed = mu + direction;

    clampMu(proposed);
    direction = proposed - mu;
  };

  // Recompute test_value at new μ with backtracking (Armijo-like condition)
  auto updateTestValue = [&]() {
    const double gradCondition = m1 * grad;
    mu_diff = direction;
    mu += mu_diff;
    clampMu(mu);

    m_value = getM(mu);
    nonLinear = Dstar(m_value);

    double new_test = -(1.0 - mu) * nonLinear + mu * d - c;

    int ls_it = 0;
    while (new_test < test_value + mu_diff * gradCondition && ls_it < kMaxBacktracking) {
      // Backtrack
      mu -= 0.5 * mu_diff;
      mu_diff *= 0.5;
      clampMu(mu);

      m_value = getM(mu);
      nonLinear = Dstar(m_value);
      new_test = -(1.0 - mu) * nonLinear + mu * d - c;

      ++ls_it;
    }

    test_value = new_test;
  };

  // Gradient at current μ
  auto updateGrad = [&]() {
    double denom = 1.0 - mu;
    if (std::fabs(denom) < kEps) {
      denom = (denom >= 0.0 ? kEps : -kEps);
    }
    grad = nonLinear - meanGap * DstarPrime(m_value) / denom + d;
  };

  // 1D BFGS: H_{k+1} = s_k / y_k with a bit of regularization
  auto updateHessian = [&]() {
    grad_diff += grad; // now grad_diff = g_{k+1} - g_k

    if (std::fabs(grad_diff) < kEps) {
      grad_diff = (grad >= 0.0 ? kEps : -kEps);
    }

    double newH = mu_diff / grad_diff;

    // keep |H| within reasonable bounds to avoid exploding steps
    if (newH > 0.0) {
      // sign should be negative for concave maximization; reflect if needed
      newH = -newH;
    }
    const double maxAbsH = 1e6;
    if (std::fabs(newH) > maxAbsH) {
      newH = (newH > 0.0 ? -maxAbsH : maxAbsH * -1.0);
    }

    inverseHessian = newH;
    grad_diff = -grad; // prepare for next iteration
  };

  // 4) Quasi-Newton iterations
  for (int iter = 0; iter < nb_Loops; ++iter) {
    // If gradient almost zero, we are at (or near) a stationary point.
    if (std::fabs(grad) < kGradTol) {
      // If test_value is still ≤ 0, we can't prune
      return test_value > 0.0;
    }

    // Concavity-based early exits
    if (grad > 0.0) {
      // Tangent from μ to μ_max cannot cross zero
      if (test_value + (mu_max_eff - mu) * grad <= 0.0) {
        return false;
      }
    } else {
      // Tangent from 0 to μ cannot cross zero
      if (test_value - mu * grad <= 0.0) {
        return false;
      }
    }

    // Quasi-Newton step
    updateDirection();
    if (std::fabs(direction) < kEps) {
      // Step too small: we are effectively stuck
      return test_value > 0.0;
    }

    updateTestValue();
    if (test_value > 0.0) {
      // Dual decision function is positive → prune s
      return true;
    }

    updateGrad();
    updateHessian();
  }

  // No pruning found within nb_Loops iterations
  return false;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// PELT = dual eval in mu = 0

// DUAL = -(c - mu * d)  - (1- mu) Dstar((a - mu b) / (1 - mu));

// a = (cumsum[t] - cumsum[s]) / (t - s)
// b = (cumsum[s] - cumsum[r]) / (s - r)

// c = (Qt - Qs) / (t - s) =  (minCost - costRecord[s]) / (t - s);
// d = (Qs - Qr) / (s - r)  = (costRecord[s] - costRecord[r]) / (s - r);

bool DUST_1D::dualMaxAlgo5(double minCost,
                           unsigned int t,
                           unsigned int s,
                           unsigned int r)
{
  double a = (cumsum[t] - cumsum[s]) / (t - s);
  double c = (minCost - costRecord[s]) / (t - s);
  if (- Dstar(a) - c > 0) {return true;} // PELT test
  return false;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// OP

bool DUST_1D::dualMaxAlgo6(double minCost,
                           unsigned int t,
                           unsigned int s,
                           unsigned int r)
{
  return false;
}




