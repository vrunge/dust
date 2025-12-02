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
  if(constr_index == "rand"){indices = new Indices_1D_Rand();}
  else{indices = new Indices_1D_Det;}

  /// /// ///
  /// /// /// dualmax_algo METHOD
  /// /// ///
  if(dualmax_algo == "DUSTr"){current_test = &DUST_meanVar::dualMaxAlgoRand;}
  if(dualmax_algo == "DUST"){current_test = &DUST_meanVar::dualMaxAlgoExact;}
  if(dualmax_algo == "PELT"){current_test = &DUST_meanVar::dualMaxPelt;}
  if(dualmax_algo == "OP"){current_test = &DUST_meanVar::dualOP;}

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
  for (unsigned t = indices->get_first() + 1; t <= n; t++)
  {
    ///////////// OP step /////////////
    ///////////// OP step /////////////
    indices->reset();
    double minCost_t = std::numeric_limits<double>::infinity();
    double minCost_tm1 = std::numeric_limits<double>::infinity();

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
      if ((this->*current_test)(minCost_tm1, t-1, indices->get_current(), indices->get_constraint()))
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

/// SPECIAL CASE BEFORE DUAL EVALUATION / OPTIMIZATION

bool DUST_meanVar::isOnePointOrLinear(double a, double b)
{
  if(isLeftBoundary(a) == true){return true;}
  if(a == b){return true;}
  return false;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


/// cases Bern/binom with right bound, not really considered here

bool DUST_meanVar::specialCasePruning(double a,
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

bool DUST_meanVar::dualMaxAlgoRand(double minCost_t,
                           unsigned int t,
                           unsigned int s,
                           unsigned int r)
{
  return (dualEval(dist(engine), minCost_t, t, s, r) > 0);
}


bool DUST_meanVar::dualOP(double minCost_t,
                          unsigned int t,
                          unsigned int s,
                          unsigned int r)
{
  return false;
}


bool DUST_meanVar::dualMaxPelt(double minCost_t,
                              unsigned int t,
                              unsigned int s,
                              unsigned int r)
{
  if(s + 1 >= t){return false;}
  //if(costRecord[s] == std::numeric_limits<double>::infinity()){return false;}
  //if(costMin(t,s) == std::numeric_limits<double>::infinity()){return false;}
  //Rcout << t << " " << s << " Q[s]: " << costRecord[s] << " c(t,s): " << costMin(t,s) << " Q[t]: " << minCost_t << " PELT: "<< (costRecord[s] + costMin(t,s) > minCost_t) << " val: "<< (costRecord[s] + costMin(t,s)) << std::endl;
  return ((costRecord[s] + costMin(t,s)) > minCost_t);
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// DUST exact with decision function inequality based test

bool DUST_meanVar::dualMaxAlgoExact(double minCost_t,
                           unsigned int t,
                           unsigned int s,
                           unsigned int r)
{
  double a = (cumsum[t] - cumsum[s]) / (t - s);
  double b = (cumsum[s] - cumsum[r]) / (s - r);
  double a2 = (cumsum2[t] - cumsum2[s]) / (t - s);
  double b2 = (cumsum2[s] - cumsum2[r]) / (s - r);
  double c = (minCost_t - costRecord[s]) / (t - s);
  double d = (costRecord[s] - costRecord[r]) / (s - r);
  double mu_max = muMax(a, b, a2, b2);

  //Rcout << " t: "  << t << " s: "  << s  << " t: "<< r << std::endl;
  ///
  /// case dual domain = one point OR dual = linear function
  ///
  if(isOnePointOrLinear(a,b) == true){return specialCasePruning(a,b,c,d,mu_max);}

  double R = -(c-d)/(a-b);
  double x_max = xMax(a, b);

  double xstar = 1.0/(a - b) * (DstarPrimeInv(R) - a);

  //Rcout << " xstar: "  << xstar << std::endl;

  if(xstar >= 0 && xstar < x_max)
  {
    // Rcout << " 00 stdcase: " << std::endl;
    return  (- Dstar(a + xstar*(a-b)) - (c + xstar*(c-d))) > 0;
  }

  if(xstar >= x_max) // case finite only
  {
    //Rcout << " a-b: " <<  a-b << " a: " <<  a <<" b: " <<  b << " xstar: " <<  xstar << " x_max " << x_max << " test " << ( (- Dstar_leftboundary() - (c + x_max*(c-d))) > 0) << " x_max*(c-d) " <<  x_max*(c-d) << std::endl;
    return  (- Dstar_leftboundary() - (c + x_max*(c-d))) > 0;
  }

  if(xstar <= 0 && xstar > -1)
  {
    // Rcout << "00  small neg: " << std::endl;
    //Rcout << " Cas1: " <<  - Dstar(a) -c << " - Dstar(a): " <<  - Dstar(a) << " -c: " << -c << std::endl;
    return (- Dstar(a) -c > 0);
  }
  if(xstar <= -1)
  {
    //Rcout << "00  BIGBIGBIGBIGBIGBIGBIGBIG neg: " << x_max << std::endl;
    //if(x_max == std::numeric_limits<double>::infinity())
    //Rcout << " Cas2: " << Dstar_leftboundary() - (c + x_max*(c-d)) <<" Dbound: " << Dstar_leftboundary()<<" constC: " <<  c << " constD: " <<  d << " RES " << ((-Dstar_leftboundary() - (c + x_max*(c-d))) > 0) << std::endl;

    if (x_max < std::numeric_limits<double>::infinity())
    {
      //Rcout << " TEST_TEST: " << (-Dstar_leftboundary() - (c + x_max*(c-d)));
      //Rcout << "x_maxx_maxx_max: " << (-Dstar_leftboundary() - (c + x_max*(c-d))) << std::endl;
      return (-Dstar_leftboundary() - (c + x_max*(c-d))) > 0;
    }
    else
    {
      return (-Dstar_superLinearLimit() - (c-d) > 0);
    }
  }

  //double xstar = 1/(a - b) * (DstarPrimeInv(R) - a);
  //if(xstar < 0){R = DstarPrime(a);}
  //if(xstar >= x_max){return false;}

  //return ((costEval(R, t, s) + costRecord[s]) > minCost_t);
}













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

double DUST_meanVar::decisionEval1(double point, double minCost_t,
                               unsigned int t, unsigned int s, unsigned int r)
{
  return (costRecord[s] - minCost_t) / (t - s) + point * (costRecord[s] - costRecord[r]) / (s - r)
  - 0.5 * std::pow((cumsum[t] - cumsum[s]) / (t - s) - point * ((cumsum[s] - cumsum[r]) / (s - r)), 2) / (1 - point);
}

double DUST_meanVar::decisionEval2(double point1, double point2, double minCost_t,
                               unsigned int t, unsigned int s, unsigned int r)
{
  return 0;
}


double DUST_meanVar::dualEval(double point, double minCost, unsigned int t, unsigned int s, unsigned int r)
{
  if(s + 1 == t){return(-std::numeric_limits<double>::infinity());}
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

bool DUST_meanVar::isLeftBoundary(double a) {return false;}
double DUST_meanVar::Dstar_leftboundary() {return std::numeric_limits<double>::infinity();}

double DUST_meanVar::Dstar_superLinearLimit() {return std::numeric_limits<double>::infinity();}

double DUST_meanVar::Dstar(double x)
{
  return 0.5 * x * x;
}


double DUST_meanVar::DstarPrime(double x)
{
  return x;
}

double DUST_meanVar::DstarPrimeInv(double x)
{
  return x;
}

double DUST_meanVar::DstarSecond(double x)
{
  return 1.0;
}





