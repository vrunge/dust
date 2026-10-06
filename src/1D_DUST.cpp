#include <Rcpp.h>
#include <cmath>
#include "1D_DUST.h"

using namespace Rcpp;

DUST_1D::DUST_1D(std::string dualmax_algo)
  : dualmax_algo(dualmax_algo),
    n(0),
    penalty(0.0)
{}

DUST_1D::~DUST_1D() {}

List DUST_1D::dust(const Rcpp::NumericVector& inData, Nullable<double> inPenalty)
{
  append_data(inData, inPenalty);
  update_partition();
  return get_partition();
}


std::forward_list<unsigned int> DUST_1D::backtrack_changepoints()
{
  std::forward_list<unsigned int> changepoints {n};
  for (int tau = chptRecord[n]; tau != 0; tau = chptRecord[tau])
  {
    changepoints.push_front(tau);
  }
  return changepoints;
}

List DUST_1D::get_info()
{
  return List::create(
    _["backend"] = "scalar",
    _["data_statistic"] = cumsum,
    _["data_length"] = n,
    _["current_penalty"] = penalty,
    _["model"] = get_model(),
    _["pruning_algo"] = dualmax_algo
  );
}
