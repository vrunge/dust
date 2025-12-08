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
  if(dualmax_algo == "DUSTib"){current_test = &DUST_1D::dualMaxAlgo7;}

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
    double minCost_t = std::numeric_limits<double>::infinity();
    unsigned int argMin = 0;
    do
    {
      unsigned int s = indices->get_current();
      lastCost = costRecord[s] + costMin(t, s); // without the penalty beta
      if (lastCost < minCost_t)
      {
        minCost_t = lastCost;
        argMin = s;
      }
      indices->next();
    }
    while(indices->is_not_the_last());
    //////// END (OP step) ////////
    //////// END (OP step) ////////

    // OP update minCost_t and save values
    // minCost_t = Q_t. Here + beta to get Q_t = Q_i + C(y_it) + beta
    minCost_t += penalty;
    costRecord.push_back(minCost_t);
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
      if ((this->*current_test)(minCost_t, t, indices->get_current(), indices->get_constraint()))
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
    if (lastCost > minCost_t)
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

  //Rcout << "NB : " << nb00 << " " << nb0 << " " << nb1 << " " << " " << nb2 << std::endl;
  //Rcout << "NBT: " << nb00T << " "<< nb0T << " " << nb1T << " " << " " << nb2T << " S: " << nb00T + nb0T + nb1T + nb2T << std::endl;
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
  if(isLeftBoundary(a) == true){return true;}
  if(a == b){return true;}
  return false;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


/// cases Bern/binom with right bound, not really considered here

bool DUST_1D::specialCasePruning(double a,
                                 double b,
                                 double c,
                                 double d,
                                 double mu_max)
{
  /// Dual in one point
  ///
  // DUAL = -(c)  - (1) * Dstar(a);
  ///
  if(isLeftBoundary(a) == true)
  {
    if (-c  - Dstar_leftboundary() > 0) {return true;} /// test in mu = 0
  }
  ///
  /// case "a = b", here "a" not on the left boundary
  ///
  if(a == b)
  {
    if ( -c - Dstar(a) > 0){return true;} /// mu = 0
    if (-(c - mu_max * d)  - Dstar(a) > 0){return true;} /// mu = mu_max
  }
  return false;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// DUAL = -(c - mu * d)  - (1- mu) * Dstar((a - mu * b) / (1 - mu));

// a = (cumsum[t] - cumsum[s]) / (t - s)
// b = (cumsum[s] - cumsum[r]) / (s - r)

// c = (Qt - Qs) / (t - s) =  (minCost_t - costRecord[s]) / (t - s);
// d = (Qs - Qr) / (s - r)  = (costRecord[s] - costRecord[r]) / (s - r);

/// RANDOM EVAL

bool DUST_1D::dualMaxAlgo0(double minCost_t,
                           unsigned int t,
                           unsigned int s,
                           unsigned int r)
{
  return (dualEval(dist(engine), minCost_t, t, s, r) > 0);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// EXACT EVAL

bool DUST_1D::dualMaxAlgo1(double minCost_t, unsigned int t, unsigned int s, unsigned int r)
{

  double a = (cumsum[t] - cumsum[s]) / (t - s);
  double b = (cumsum[s] - cumsum[r]) / (s - r);
  double c = (minCost_t - costRecord[s]) / (t - s);
  double d = (costRecord[s] - costRecord[r]) / (s - r);
  double mu_max = muMax(a, b);
  ///
  /// case dual domain = one point OR dual = linear function
  ///
  //if(isOnePointOrLinear(a,b) == true)
  //{
  //nb00 = nb00 + 1;
  //  if (specialCasePruning(a,b,c,d,mu_max)) nb00T = nb00T + 1;
  //}
  if(isOnePointOrLinear(a,b) == true){return specialCasePruning(a,b,c,d,mu_max);}


  // value a visible, not hidden by the r function
  // PELT pruning
  //Rcout << DstarPrime(a) << " "<<  -(a-b)*DstarPrime(a) << " " << -(a-b)*DstarPrime(a) - (c-d) << std::endl;
  //
  // if derivative in 0 is negative
  //
  if(-(a-b)*DstarPrime(a) - (c-d) < 0)
  {
    //nb0 = nb0 + 1;
    //if(- Dstar(a) -c > 0) nb0T = nb0T + 1;
    return (- Dstar(a) -c > 0);
  }

  //
  // if derivative in +inf is positive
  //
  if( a > b)
  {
    if(-DstarPrime(std::numeric_limits<double>::infinity()) - (c-d) > 0)
    {
      return true;
    }
  }

  return(costEval(-(c-d)/(a-b), t, s) > c);

  //double x_max = xMax(a, b);
  //if(x_max == std::numeric_limits<double>::infinity())
  //  {
    //Rcout << "MAX "<< -(a-b) << " "<< a+ x_max*(a-b) << " " << DstarPrime(a+ x_max*(a-b)) << " " << -(a-b)*DstarPrime(a+ x_max*(a-b)) << std::endl;
    //  if(-(a-b)*DstarPrime(a+ x_max*(a-b)) - (c-d) > 0)
    //{
    //  return true;
      //nb2 = nb2 + 1;
      //if(-Dstar_superLinearLimit() - (c-d) > 0) nb2T = nb2T + 1;
      //return(-Dstar_superLinearLimit() - (c-d) > 0);
      //}
      //  }
  // Rcout << "COST1 "<<   -(c-d)/(a-b) << " "<<  costEval(-(c-d)/(a-b), t, s) << " "<< cumsum[t] - cumsum[s] << " "<< t << " "<<  s << std::endl;
  //nb1 = nb1 + 1;
  //if(costEval(-(c-d)/(a-b), t, s) + costRecord[s] > costRecord[t]) nb1T = nb1T + 1;

  // (t-s) costEval(-(c-d)/(a-b), t, s) + Qs + beta > Qt + beta
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// Golden-section search

bool DUST_1D::dualMaxAlgo2(double minCost_t,
                           unsigned int t,
                           unsigned int s,
                           unsigned int r)
{
  double a = (cumsum[t] - cumsum[s]) / (t - s);
  double b = (cumsum[s] - cumsum[r]) / (s - r);
  double c = (minCost_t - costRecord[s]) / (t - s);
  double d = (costRecord[s] - costRecord[r]) / (s - r);
  double mu_max = muMax(a, b);

  ///
  /// case dual domain = one point OR dual = linear function
  if(isOnePointOrLinear(a,b) == true){return specialCasePruning(a,b,c,d,mu_max);}



  double lt = 0.0;
  double rt = mu_max;

  double cs1 = (1 - 1/phi) * rt;
  double cs2 = rt/phi;

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

bool DUST_1D::dualMaxAlgo3(double minCost_t,
                           unsigned int t,
                           unsigned int s,
                           unsigned int r)
{
  double a = (cumsum[t] - cumsum[s]) / (t - s);
  double b = (cumsum[s] - cumsum[r]) / (s - r);
  double c = (minCost_t - costRecord[s]) / (t - s);
  double d = (costRecord[s] - costRecord[r]) / (s - r);

  double mu_max = muMax(a, b);

  ///
  /// case dual domain = one point OR dual = linear function
  ///
  if(isOnePointOrLinear(a,b) == true){return specialCasePruning(a,b,c,d,mu_max);}

  ///
  /// PELT test (dual at mu = 0)
  ///
  const double Dstar_a = Dstar(a);
  double test_value = - Dstar_a - c;  // dual(mu = 0)
  if (test_value > 0.0) {return true;}  // We already found a positive dual value at mu = 0

  // if non pruning at mu = 0
  // AND
  // If the value on the tangent at 0 in mu_max is negative we can prune everything
  double grad = Dstar_a - (a - b) * DstarPrime(a) + d;
  if (test_value + mu_max * grad <= 0) {return false;} //check value in mu = mu_max

  /// iterative algo
  /// iterative algo
  /// iterative algo
  ///
  /// Iterative derivative-aided bisection on [0, mu_max]
  ///
  double lt = 0.0;
  double rt = mu_max;
  double mu = 0.5 * rt;

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

    if (grad > 0.0)
    {
      // Candidate max can only be on the right side for a concave f
      if (test_value + (rt - mu) * grad <= 0.0) {
        // Upper bound on [mu, rt] is ≤ 0 -> nothing positive to the right
        return false;
      }
      lt = mu;
    }
    else // grad <= 0
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



bool DUST_1D::dualMaxAlgo4(double minCost_t,
                           unsigned int t,
                           unsigned int s,
                           unsigned int r)
{
  double a = (cumsum[t] - cumsum[s]) / (t - s);
  double b = (cumsum[s] - cumsum[r]) / (s - r);
  double c = (minCost_t - costRecord[s]) / (t - s);
  double d = (costRecord[s] - costRecord[r]) / (s - r);
  double mu_maxx = muMax(a, b);

  ///
  /// case dual domain = one point OR dual = linear function
  ///
  if(isOnePointOrLinear(a,b) == true){return specialCasePruning(a,b,c,d,mu_maxx);}


  ///
  /// TEST PELT in MU = 0
  ///
  double nonLinear = Dstar(a);
  double test_value = - c - nonLinear;

  if (test_value > 0) {return true;} // PELT test

  ///
  /// mu_max (stop ALGO if mu_max == 0) !!!
  ///

  double mu_max = muMax(a, b);
  if (mu_max == 0) {return false;}

  ///
  /// INIT Quasi Newton algo in zero
  ///
  double meanGap = a - b;
  double grad = nonLinear - meanGap * DstarPrime(a) + d;

  /// if grad < 0 the max is in zero, no pelt pruning at this step, return false
  if (grad < 0) {return false;}

  // the duality function is concave, meaning we can check if the tangent at mu ever reaches the desired value.
  if (test_value + mu_max * grad <= 0) {return false;}
  double mu = 0;
  double direction;
  double mu_diff;
  double m_value;
  double grad_diff = -grad; // stores grad difference between two steps, initialized each step as g_k, then once g_{k+1} is computed, y += g_{k+1}
  double inverseHessian = -1;

  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  auto updateDirection = [&] () // update and clip direction
  {
    direction = - inverseHessian * grad;
    if (mu + direction > mu_max) { direction = mu_max - 1e-9 - mu; } // clip ascent direction, as dual may increase beyond the bound
    else if (mu + direction < 0) { direction = -mu + 1e-9; }
  };

  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  auto getM = [&] (double point) // for use in Dstar, DstarPrime
  {
    return (a - point * b) / (1 - point);
  };

  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  auto updateTestValue = [&] ()
  {
    double gradCondition = m1 * grad;

    // Initialize all values
    mu_diff = direction;
    mu += mu_diff;
    m_value = getM(mu);
    nonLinear = Dstar(m_value);
    double new_test = - (1 - mu) * nonLinear + mu * d -c;

    int i = 0;
    while(new_test < test_value + mu_diff * gradCondition)
    {
      mu_diff *= .5; // shrink if unsuitable stepsize
      mu -= mu_diff; // relay shrinking
      m_value = getM(mu); // update values
      nonLinear = Dstar(m_value); // update values
      new_test = - (1 - mu) * nonLinear + mu * d -c; // update values
      i++;
      if (i == 10) { break; }
    }
    test_value = new_test;
  };

  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  auto updateGrad = [&] ()
  {
    grad = nonLinear - meanGap * DstarPrime(m_value) / (1 - mu) + d;
  };

  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  auto updateHessian = [&] () // uses BFGS formula (in 1D, the formula simplifies to s / y)
  {
    grad_diff += grad;
    if(grad_diff == 0) { grad_diff = 1e-9; }

    inverseHessian = mu_diff / grad_diff;
    grad_diff = - grad; // setting up y for next step
  };

  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  int i = 0;
  do
  {
    updateDirection();
    updateTestValue();
    if(test_value > 0) {return true;} // index s is pruned
    updateGrad();

    if (grad > 0) // the duality function is concave, meaning we can check if the tangent at mu ever reaches the desired value.
    {
      if (test_value + (mu_max - mu) * grad <= 0) {return false;}
    }
    else if (test_value - mu * grad <= 0) {return false;}
    updateHessian();
    i++;
  } while (i < nb_Loops);
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

// c = (Qt - Qs) / (t - s) =  (minCost_t - costRecord[s]) / (t - s);
// d = (Qs - Qr) / (s - r)  = (costRecord[s] - costRecord[r]) / (s - r);

bool DUST_1D::dualMaxAlgo5(double minCost_t,
                           unsigned int t,
                           unsigned int s,
                           unsigned int r)
{
  double a = (cumsum[t] - cumsum[s]) / (t - s);
  double c = (minCost_t - costRecord[s]) / (t - s);
  if (- Dstar(a) - c > 0) {return true;} // PELT test
  return false;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// OP

bool DUST_1D::dualMaxAlgo6(double minCost_t,
                           unsigned int t,
                           unsigned int s,
                           unsigned int r)
{
  return false;
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// DUST exact with decision function inequality based test

bool DUST_1D::dualMaxAlgo7(double minCost_t,
                           unsigned int t,
                           unsigned int s,
                           unsigned int r)
{
  double a = (cumsum[t] - cumsum[s]) / (t - s);
  double b = (cumsum[s] - cumsum[r]) / (s - r);
  double c = (minCost_t - costRecord[s]) / (t - s);
  double d = (costRecord[s] - costRecord[r]) / (s - r);
  double mu_max = muMax(a, b);
  ///
  /// case dual domain = one point OR dual = linear function
  ///
  if(isOnePointOrLinear(a,b) == true){return specialCasePruning(a,b,c,d,mu_max);}

  double R = -(c-d)/(a-b);
  double x_max = xMax(a, b);

  double xstar = 1.0/(a - b) * (DstarPrimeInv(R) - a);

  /// case PELT
  if(xstar <= 0 && xstar > -1)
  {
    return (- Dstar(a) -c > 0);
  }

  /// INTERIOR POINT
  if(xstar >= 0 && xstar < x_max)
  {
    return  (- Dstar(a + xstar*(a-b)) - (c + xstar*(c-d))) > 0;
  }

  if((xstar <= -1) && (x_max < std::numeric_limits<double>::infinity()))
  {
    return (-Dstar_leftboundary() - (c + x_max*(c-d))) > 0;
  }

  if(xstar >= x_max) // case finite only
  {
    return  (- Dstar_leftboundary() - (c + x_max*(c-d))) > 0;
  }
  return (-Dstar_superLinearLimit() - (c-d) > 0);

}





