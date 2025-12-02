#include <Rcpp.h>
#include <cmath>
#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL
#include "2DmeanVar_DUST.h"
#include "preProcessing.h"

#include <fstream>
#include <iostream>

using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// --- // Constructor // --- //
DUST_meanVar::DUST_meanVar(std::string dualmax_algo,
                 std::string constr_index)
  : dualmax_algo(dualmax_algo),
    constr_index(constr_index),
    indices(nullptr),
    n(0)
{}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

DUST_meanVar::~DUST_meanVar()
{
  delete indices;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void DUST_meanVar::pruning_method()
{
  delete indices;

  /// /// ///
  /// /// /// index METHOD
  /// /// ///
  if(constr_index == "det1"){indices = new Indices_2D_Det1();}
  else{indices = new Indices_2D_Det2();}

  /// /// ///
  /// /// /// dualmax_algo METHOD
  /// /// ///
  if(dualmax_algo == "DUSTr"){current_test = &DUST_meanVar::dualMaxAlgoRand;}
  if(dualmax_algo == "PELT"){current_test = &DUST_meanVar::dualMaxPelt;}
  if(dualmax_algo == "OP"){current_test = &DUST_meanVar::dualOP;}
  if(dualmax_algo == "DUST1"){current_test = &DUST_meanVar::decisionTest1;}
  if(dualmax_algo == "DUST2"){current_test = &DUST_meanVar::decisionTest2;}

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

void DUST_meanVar::append_data(std::vector<double>& inData, Nullable<double> inPenalty)
{
  bool first_execution = (n == 0);
  n += inData.size();   // update total size

  // Reserve memory for all vectors
  cumsum.reserve(n + 1);
  cumsum2.reserve(n + 1);
  chptRecord.reserve(n + 1);
  costRecord.reserve(n + 1);
  nb_indices.reserve(n + 1);

  // On first execution: initialize all vectors
  if (first_execution)
  {
    if (inPenalty.isNull()) { penalty = 2 * std::log(n); } else { penalty = as<double>(inPenalty); }

    cumsum.push_back(0);
    cumsum2.push_back(0);
    costRecord.push_back(-penalty);
    chptRecord.push_back(0);
    nb_indices.push_back(1);

    pruning_method();
    indices->add_first(0);
  }

  // store data by updating the cumsum
  for (auto current_datum = inData.begin(); current_datum != inData.end(); ++current_datum)
  {
    cumsum.push_back(cumsum.back() + *current_datum);
    cumsum2.push_back(cumsum2.back() + (*current_datum) * (*current_datum));
  }
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void DUST_meanVar::update_partition()
{
  // Initialize OP step value
  double lastCost;
  unsigned int nbt = nb_indices.back();
  double minCost_tm1 = std::numeric_limits<double>::infinity();

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
      if ((this->*current_test)(minCost_tm1, t-1,
           indices->get_current(), indices->get_constraint_r1(), indices->get_constraint_r2()))
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
      //indices->prune_last();
      //nbt--;
    }

    minCost_tm1 = minCost_t;
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
List DUST_meanVar::dust(std::vector<double>& inData, Nullable<double> inPenalty)
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
std::forward_list<unsigned int> DUST_meanVar::backtrack_changepoints()
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
List DUST_meanVar::get_partition()
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
List DUST_meanVar::get_info()
{
  return List::create(
    _["data_statistic"] = cumsum,
    _["data_statistic2"] = cumsum2,
    _["data_length"] = n,
    _["current_penalty"] = penalty,
    _["pruning_algo"] = dualmax_algo,
    _["pruning_constr_index"] = constr_index
  );
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool DUST_meanVar::dualMaxAlgoRand(double minCost_t,
                           unsigned int t,
                           unsigned int s,
                           unsigned int r1,
                           unsigned int r2)
{
  return (dualEval(dist(engine), minCost_t, t, s, r1) > 0);
}


bool DUST_meanVar::dualOP(double minCost_t,
                          unsigned int t,
                          unsigned int s,
                          unsigned int r1,
                          unsigned int r2)
{
  return false;
}


bool DUST_meanVar::dualMaxPelt(double minCost_t,
                              unsigned int t,
                              unsigned int s,
                              unsigned int r1,
                              unsigned int r2)
{
  if(s + 1 >= t){return false;}
  //if(costRecord[s] == std::numeric_limits<double>::infinity()){return false;}
  //if(costMin(t,s) == std::numeric_limits<double>::infinity()){return false;}
  //Rcout << t << " " << s << " Q[s]: " << costRecord[s] << " c(t,s): " << costMin(t,s) << " Q[t]: " << minCost_t << " PELT: "<< (costRecord[s] + costMin(t,s) > minCost_t) << " val: "<< (costRecord[s] + costMin(t,s)) << std::endl;
  return ((costRecord[s] + costMin(t,s)) > minCost_t);
}



bool DUST_meanVar::decisionTest1(double minCost_t,
                                 unsigned int t,
                                 unsigned int s,
                                 unsigned int r1,
                                 unsigned int r2)
{
  if(s + 1 >= t){return(false);}
  double a = (cumsum[t] - cumsum[s]) / (t - s);
  double a2 = (cumsum2[t] - cumsum2[s]) / (t - s);
  double b = (cumsum[s] - cumsum[r1]) / (s - r1);
  double b2 = (cumsum2[s] - cumsum2[r1]) / (s - r1);

  double c = (minCost_t - costRecord[s]) / (t - s);
  double d = (costRecord[s] - costRecord[r1]) / (s - r1);

  // Compute variance terms
  double Va = a2 - std::pow(a, 2);
  double Vb = b2 - std::pow(b, 2);
  double amb2 = ((a-b)*(a-b));
  double x0 = 0.5*((Va - Vb)/amb2 - 1);
  double x1 = x0*x0 + Va/amb2;
  double x2 = c - d;

  double sign_x2 = (x2 > 0.0) ? 1.0 : (x2 < 0.0 ? -1.0 : 0.0);
  double xstar = std::max(0.0, x0 + 1.0 / (2.0 * x2) - sign_x2 * std::sqrt(x1 + 1.0 / (4.0 * x2* x2)));
  double A = a2 + xstar * (a2 - b2);
  double B = a + xstar * (a - b);
  //Rcout << t << " " << s << " " << r << " xstar: " << xstar << " A: " << A << " B: " << B << " res: "<< (0.5 * (1 + std::log(A - B*B)) - (c + xstar * (c - d))) << std::endl;
  return  0.5 * (1 + std::log(A - B*B)) - (c + xstar * (c - d)) > 0;
}

bool DUST_meanVar::decisionTest2(double minCost_t,
                                 unsigned int t,
                                 unsigned int s,
                                 unsigned int r1,
                                 unsigned int r2)
{
  if(s + 1 >= t){return(false);}

  const double inv_ts   = 1.0 / static_cast<double>(t - s);
  const double inv_s_r1 = 1.0 / static_cast<double>(s - r1);
  const double inv_s_r2 = 1.0 / static_cast<double>(s - r2);

  double st  = (cumsum[t] - cumsum[s]) * inv_ts;
  double st2 = (cumsum2[t] - cumsum2[s]) * inv_ts;

  double r1s  = (cumsum[s] - cumsum[r1])  * inv_s_r1;
  double r1s2 = (cumsum2[s] - cumsum2[r1])  * inv_s_r1;

  double r2s  = (cumsum[s]  - cumsum[r2])  * inv_s_r2;
  double r2s2 = (cumsum2[s] - cumsum2[r2]) * inv_s_r2;

  // Q_st and Q_rk,s
  double c0 = (minCost_t     - costRecord[s])  * inv_ts;    // Q_st
  double d1 = (costRecord[s] - costRecord[r1]) * inv_s_r1;  // Q_r1s
  double d2 = (costRecord[s] - costRecord[r2]) * inv_s_r2;  // Q_r2s

  // Compute variance terms
  double Va = st2 - std::pow(st, 2);

  double Vb, amb2, x0, x1, x2, sign_x2, xstar, A, B;

  Vb = r1s2 - std::pow(r1s, 2);
  amb2 = ((st-r1s)*(st-r1s));
  x0 = 0.5*((Va - Vb)/amb2 - 1);
  x1 = x0*x0 + Va/amb2;
  x2 = c0 - d1;
  sign_x2 = (x2 > 0.0) ? 1.0 : (x2 < 0.0 ? -1.0 : 0.0);
  xstar = std::max(0.0, x0 + 1.0 / (2.0 * x2) - sign_x2 * std::sqrt(x1 + 1.0 / (4.0 * x2* x2)));
  A = st2 + xstar * (st2 - r1s2);
  B = st + xstar * (st - r1s);
  if(0.5 * (1 + std::log(A - B*B)) - (c0 + xstar * (c0 - d1)) > 0){return true;}

  Vb = r2s2 - std::pow(r2s, 2);
  amb2 = ((st-r2s)*(st-r2s));
  x0 = 0.5*((Va - Vb)/amb2 - 1);
  x1 = x0*x0 + Va/amb2;
  x2 = c0 - d2;
  sign_x2 = (x2 > 0.0) ? 1.0 : (x2 < 0.0 ? -1.0 : 0.0);
  xstar = std::max(0.0, x0 + 1.0 / (2.0 * x2) - sign_x2 * std::sqrt(x1 + 1.0 / (4.0 * x2* x2)));
  A = st2 + xstar * (st2 - r2s2);
  B = st + xstar * (st - r2s);
  if(0.5 * (1 + std::log(A - B*B)) - (c0 + xstar * (c0 - d2)) > 0){return true;}


  // 2) Interior candidate from the first-order conditions
  const double eps = 1e-12;

  // ΔS and ΔS^2 and ΔQ
  const double b0 = st;
  const double b1 = st - r1s;      // ΔS_r1st
  const double b2 = st - r2s;      // ΔS_r2st

  const double a0 = st2;
  const double a1 = st2 - r1s2;    // ΔS^2_r1st
  const double a2 = st2 - r2s2;    // ΔS^2_r2st

  const double c1 = c0 - d1;       // ΔQ_r1st
  const double c2 = c0 - d2;       // ΔQ_r2st

  // A_* (variance parameter) from first-order conditions
  const double denomA = 2.0 * (b1 * c2 - b2 * c1);
  if (std::fabs(denomA) < eps) {
    // Degenerate system: we already know the boundary does not give a positive
    // value, so we just check the (0,0) point.
    return decisionEval2(0.0, 0.0, minCost_t, t, s, r1, r2) > 0.0;
  }
  const double A_star = (b1 * a2 - b2 * a1) / denomA;

  // y_* : use the more stable expression among the two (index 1 or 2)
  double y_star;
  if (std::fabs(b1) >= std::fabs(b2) && std::fabs(b1) > eps) {
    y_star = (a1 - 2.0 * c1 * A_star) / (2.0 * b1);
  } else if (std::fabs(b2) > eps) {
    y_star = (a2 - 2.0 * c2 * A_star) / (2.0 * b2);
  } else {
    // Both b1 and b2 are ~0: effectively 1D, and boundaries are already checked
    return decisionEval2(0.0, 0.0, minCost_t, t, s, r1, r2) > 0.0;
  }

  // Linear system for x1*, x2*:
  //   a1 x1 + a2 x2 = k1
  //   b1 x1 + b2 x2 = k2
  const double k1 = A_star + y_star * y_star - a0;
  const double k2 = y_star - b0;

  const double Delta = a1 * b2 - a2 * b1;
  if (std::fabs(Delta) < eps) {
    // Another degeneracy: again, we fall back to (0,0)
    return decisionEval2(0.0, 0.0, minCost_t, t, s, r1, r2) > 0.0;
  }

  double x1_star = (k1 * b2 - k2 * a2) / Delta;
  double x2_star = (a1 * k2 - b1 * k1) / Delta;

  // Check domain: x1 >= 0, x2 >= 0, and A(x1*, x2*) > 0 (for the log)
  const double y_val = b0 + b1 * x1_star + b2 * x2_star;
  const double A_val = a0 + a1 * x1_star + a2 * x2_star - y_val * y_val;

  if (A_val <= 0.0 || x1_star <= 0.0 || x2_star <= 0.0) {
    // Interior candidate is not admissible; optimum is on the boundary, and
    // we already know the boundary does not produce a positive value.
    return decisionEval2(0.0, 0.0, minCost_t, t, s, r1, r2) > 0.0;
  }

  // Valid interior maximiser: test its value
  return decisionEval2(x1_star, x2_star, minCost_t, t, s, r1, r2) > 0.0;
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////





double DUST_meanVar::costEval(double point1, double point2, unsigned int t, unsigned int s)
{
  return - point1*point1/(4.0*point2) + 0.5*std::log(-1.0/(2.0*point2))
  - point1 *(cumsum[t] - cumsum[s]) - point2 * (cumsum2[t] - cumsum2[s]);
}

double DUST_meanVar::costMin(unsigned int t, unsigned int s)
{
  if(s + 1 == t){return(std::numeric_limits<double>::infinity());} // infinite cost segment if one data point only
  double m = (cumsum[t] - cumsum[s]) / (t - s);
  double var = (cumsum2[t] - cumsum2[s]) / (t - s) - (m * m);
  if(var <= 0){return(std::numeric_limits<double>::infinity());} /// choice  1e-100 to avoid -Inf /// THIS IS IMPORTANT
  return 0.5 * (t - s) * (1.0 + std::log(var));
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



double DUST_meanVar::decisionEval1(double point,
                                   double minCost_t,
                                   unsigned int t,
                                   unsigned int s,
                                   unsigned int r)
{
  return (costRecord[s] - minCost_t) / (t - s) + point * (costRecord[s] - costRecord[r]) / (s - r)
  - 0.5 * std::pow((cumsum[t] - cumsum[s]) / (t - s) - point * ((cumsum[s] - cumsum[r]) / (s - r)), 2) / (1 - point);
}


double DUST_meanVar::decisionEval2(double point1,
                                   double point2,
                                   double minCost_t,
                                   unsigned int t,
                                   unsigned int s,
                                   unsigned int r1,
                                   unsigned int r2)
{
  // Precompute reciprocals once (avoids repeated divisions and is
  // slightly more numerically consistent).
  double inv_ts   = 1.0 / static_cast<double>(t - s);
  double inv_s_r1 = 1.0 / static_cast<double>(s - r1);
  double inv_s_r2 = 1.0 / static_cast<double>(s - r2);

  // Means of y and y^2 on [s+1, t], [r1+1, s], [r2+1, s]
  double st   = (cumsum[t]  - cumsum[s])  * inv_ts;
  double st2  = (cumsum2[t] - cumsum2[s]) * inv_ts;
  double r1s  = (cumsum[s]  - cumsum[r1])  * inv_s_r1;
  double r1s2 = (cumsum2[s] - cumsum2[r1]) * inv_s_r1;
  double r2s  = (cumsum[s]  - cumsum[r2])  * inv_s_r2;
  double r2s2 = (cumsum2[s] - cumsum2[r2]) * inv_s_r2;

  // Q_st and Q_rk,s
  double c0 = (minCost_t     - costRecord[s])  * inv_ts;    // Q_st
  double d1 = (costRecord[s] - costRecord[r1]) * inv_s_r1;  // Q_r1s
  double d2 = (costRecord[s] - costRecord[r2]) * inv_s_r2;  // Q_r2s

  // ΔS and ΔS^2 and ΔQ, in the notations of the paper
  const double b0 = st;
  const double b1 = st - r1s;      // ΔS_r1st
  const double b2 = st - r2s;      // ΔS_r2st

  const double a0 = st2;
  const double a1 = st2 - r1s2;    // ΔS^2_r1st
  const double a2 = st2 - r2s2;    // ΔS^2_r2st

  const double c1 = c0 - d1;       // ΔQ_r1st
  const double c2 = c0 - d2;       // ΔQ_r2st

  // y = S_st + x1 ΔS_r1st + x2 ΔS_r2st
  double y = b0 + point1 * b1 + point2 * b2;
  // D(x1, x2) = 0.5 * (1 + log(A)) - (Q_st + x1 ΔQ_r1st + x2 ΔQ_r2st)
  return 0.5 * (1.0 + std::log(a0 + point1 * a1 + point2 * a2 - y * y)) - (c0 + point1 * c1 + point2 * c2);
}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


double DUST_meanVar::dualEval(double point,
                              double minCost,
                              unsigned int t,
                              unsigned int s,
                              unsigned int r)
{
  if(s + 1 >= t){return(-std::numeric_limits<double>::infinity());}
  //if(r + 1 == s){return(-std::numeric_limits<double>::infinity());} // => Vb = 0
  double Mt = (cumsum[t] - cumsum[s]) / (t - s);
  double Mt2 = (cumsum2[t] - cumsum2[s]) / (t - s);
  double Ms = (cumsum[s] - cumsum[r]) / (s - r);
  double Ms2 = (cumsum2[s] - cumsum2[r]) / (s - r);

  // Compute variance terms
  double Va = Mt2 - std::pow(Mt, 2);
  double Vb = Ms2 - std::pow(Ms, 2);

  /// pruning if same mean and same variance
  //if(Mt == Ms && Va == Vb){return(std::numeric_limits<double>::infinity());}

  double u = (Va + Vb) * (1 + std::pow((Mt - Ms) / std::sqrt(Va + Vb), 2));

  if(Vb > 0){point = point * ((u - std::sqrt(std::pow(u, 2) - 4.0 * Va * Vb)) / (2.0 * Vb));}
  else{point = point * (Va / (Va + pow(Mt - Ms, 2)));}

  //std::cout << point << " ";
  double A = (Mt2 - point *  Ms2)/(1 - point);
  double B = (Mt - point *  Ms)/(1 - point);

  return (costRecord[s] - minCost) / (t - s)
    + point * (costRecord[s] - costRecord[r]) / (s - r)
    + 0.5 * (1 - point) * (1 + std::log(A - B*B));
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

double DUST_meanVar::dualMax(double minCost_t, unsigned int t, unsigned int s, unsigned int r)
{
  return 0;
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////



double DUST_meanVar::muMax(double a, double b, double a2, double b2)
{
  double Va = a2 - std::pow(a, 2);
  double Vb = b2 - std::pow(b, 2);
  double u = (Va + Vb) * (1 + std::pow((a - b) / std::sqrt(Va + Vb), 2));

  if(Vb > 0){return((u - std::sqrt(std::pow(u, 2) - 4.0 * Va * Vb)) / (2.0 * Vb));}
  else{return(Va / (Va + pow(a - b, 2)));}
}


double DUST_meanVar::xMax(double a, double b)
{
  return std::numeric_limits<double>::infinity();
}
