// Highway backend for the one-parameter segmentation methods.
// It uses contiguous active indices for SIMD evaluation and falls back to
// the scalar implementation when Highway is unavailable. The DUST pruning
// test follows the dual decision rule of Runge, Truong and Querné (2025),
// doi:10.48550/arXiv.2507.02467.

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdint>
#include <string>
#include <type_traits>

using namespace Rcpp;

// --------------------------------------------------------------------- //
// Delegation to the scalar engine when built without Highway.
// --------------------------------------------------------------------- //
#include "1D_DUST.h"
#include "1D_DUST_Impl.h"
#include "1D_OP_Impl.h"
#include "1D_Indices.h"
#include "1D_DualMaxPolicies.h"
#include "1D_Validation.h"
#include "1D_1_GaussModel.h"
#include "1D_2_PoissonModel.h"
#include "1D_3_ExpModel.h"
#include "1D_4_GeomModel.h"
#include "1D_5_BernModel.h"
#include "1D_6_BinomModel.h"
#include "1D_7_NegbinModel.h"
#include "1D_8_VarianceModel.h"

namespace {

template <class Model, class DualMaxPolicy>
List run_dust_engine(NumericVector data, double penalty)
{
  DUST_1D_T<Model, DualMaxPolicy> obj;
  obj.dust(data, Rcpp::Nullable<double>(Rcpp::wrap(penalty)));
  return obj.get_partition();
}

template <class Model>
List run_dust_by_algo(NumericVector data, double penalty, const std::string& method)
{
  if (method == "DUST")   return run_dust_engine<Model, DualMax_DUST<Model>>(data, penalty);
  if (method == "DUSTib") return run_dust_engine<Model, DualMax_DUSTib<Model>>(data, penalty);
  if (method == "PELT")   return run_dust_engine<Model, DualMax_PELT<Model>>(data, penalty);
  if (method == "OP")
  {
    // Genuinely unpruned OP (1D_OP_Impl.h), not the general engine with a
    // no-op DualMaxPolicy -- see _ModuleAssembly_1D.cpp's newModuleByAlgo
    // for why that distinction matters.
    DUST_1D_OP_T<Model> obj;
    obj.dust(data, Rcpp::Nullable<double>(Rcpp::wrap(penalty)));
    return obj.get_partition();
  }
  stop("Unrecognized method \"" + method + "\". Valid values: \"DUST\", \"DUSTib\", \"PELT\", \"OP\".");
}

[[maybe_unused]] List run_dust_by_model(NumericVector data, double penalty, const std::string& model, const std::string& method)
{
  if (model == "gauss")    return run_dust_by_algo<GaussPolicy>(data, penalty, method);
  if (model == "poisson")  return run_dust_by_algo<PoissonPolicy>(data, penalty, method);
  if (model == "exp")      return run_dust_by_algo<ExpPolicy>(data, penalty, method);
  if (model == "geom")     return run_dust_by_algo<GeomPolicy>(data, penalty, method);
  if (model == "bern")     return run_dust_by_algo<BernPolicy>(data, penalty, method);
  if (model == "binom")    return run_dust_by_algo<BinomPolicy>(data, penalty, method);
  if (model == "negbin")   return run_dust_by_algo<NegbinPolicy>(data, penalty, method);
  if (model == "variance") return run_dust_by_algo<VariancePolicy>(data, penalty, method);
  stop("Unrecognized model \"" + model + "\".");
}

} // namespace

// --------------------------------------------------------------------- //
// Highway-vectorized path.
// --------------------------------------------------------------------- //
#ifdef HAVE_HIGHWAY
#include "DUST_1D_HW_Models.h"
#include "DUST_1D_HW_DUSTib.h"
#ifndef DUSTIB_VERSION
#define DUSTIB_VERSION 2
#endif

namespace {

constexpr double EPS = 1e-14;

// Map Highway policies to their scalar counterparts for the DUST test.
template <class VecModel> struct ScalarModelOf;
template <> struct ScalarModelOf<VecGaussPolicy>    { using type = GaussPolicy; };
template <> struct ScalarModelOf<VecPoissonPolicy>  { using type = PoissonPolicy; };
template <> struct ScalarModelOf<VecExpPolicy>      { using type = ExpPolicy; };
template <> struct ScalarModelOf<VecGeomPolicy>     { using type = GeomPolicy; };
template <> struct ScalarModelOf<VecBernPolicy>     { using type = BernPolicy; };
template <> struct ScalarModelOf<VecNegbinPolicy>   { using type = NegbinPolicy; };
template <> struct ScalarModelOf<VecVariancePolicy> { using type = VariancePolicy; };

// Gaussian DUST has no domain branches and is evaluated in SIMD lanes.
// Other models use the scalar test on the contiguous candidate arrays.
template <class VecModel> struct DustUsesVectorTest : std::false_type {};
template <> struct DustUsesVectorTest<VecGaussPolicy> : std::true_type {};

// Evaluate the Gaussian DUST test in SIMD lanes.
template <class VecModel, class D>
void dust_lane(D d, double cumsum_t, double minCost_t, double t_pos,
               const double* A_cumsum, const double* A_cost, const double* A_pos,
               size_t off, uint8_t* prune)
{
  auto vCumsumT  = hn::Set(d, cumsum_t);
  auto vMinCostT = hn::Set(d, minCost_t);
  auto vTpos     = hn::Set(d, t_pos);
  auto vEps      = hn::Set(d, EPS);
  auto vZero     = hn::Zero(d);
  const double kInf = std::numeric_limits<double>::infinity();

  auto vCs = hn::LoadU(d, A_cumsum + off);
  auto vQs = hn::LoadU(d, A_cost + off);
  auto vPs = hn::LoadU(d, A_pos + off);
  auto vCr = hn::LoadU(d, A_cumsum + off - 1);
  auto vQr = hn::LoadU(d, A_cost + off - 1);
  auto vPr = hn::LoadU(d, A_pos + off - 1);

  auto vA = hn::Div(hn::Sub(vCumsumT, vCs), hn::Sub(vTpos, vPs));
  auto vB = hn::Div(hn::Sub(vCs, vCr), hn::Sub(vPs, vPr));
  auto vC = hn::Div(hn::Sub(vMinCostT, vQs), hn::Sub(vTpos, vPs));
  auto vD = hn::Div(hn::Sub(vQs, vQr), hn::Sub(vPs, vPr));

  auto vMuMax = VecModel::vMuMax(d, vA, vB);
  auto maskLeft  = VecModel::vIsLeftBoundary(d, vA);
  auto maskRight = VecModel::vIsRightBoundary(d, vA);
  auto vDiffAB = hn::Sub(vA, vB);
  auto maskDiffSmall = hn::Lt(hn::Abs(vDiffAB), vEps);
  auto maskLinear = hn::Or(hn::Or(maskLeft, maskRight), maskDiffSmall);

  auto vDstarA_plain = VecModel::vDstar(d, vA);
  auto vDstarLeft  = hn::Set(d, VecModel::DstarLeftBoundary());
  auto vDstarRight = hn::Set(d, VecModel::DstarRightBoundary());
  auto vDstarA_eff = hn::IfThenElse(maskLeft, vDstarLeft, hn::IfThenElse(maskRight, vDstarRight, vDstarA_plain));

  auto test1_left  = hn::Gt(hn::Sub(hn::Neg(vC), vDstarLeft), vZero);
  auto test1_right = hn::Gt(hn::Sub(hn::Neg(vC), vDstarRight), vZero);
  auto test2 = hn::Gt(hn::Sub(hn::Neg(vC), vDstarA_eff), vZero);
  auto test3 = hn::Gt(hn::Sub(hn::Neg(hn::Sub(vC, hn::Mul(vMuMax, vD))),
                               hn::Mul(hn::Sub(hn::Set(d,1.0), vMuMax), vDstarA_eff)), vZero);
  auto maskSpecialPrune = hn::Or(hn::Or(hn::And(maskLeft, test1_left), hn::And(maskRight, test1_right)),
                                  hn::And(maskDiffSmall, hn::Or(test2, test3)));

  auto vDstarPrimeA = VecModel::vDstarPrime(d, vA);
  auto vCD = hn::Sub(vC, vD);
  auto vLhs1 = hn::Sub(hn::Neg(hn::Mul(vDiffAB, vDstarPrimeA)), vCD);
  auto maskMainNeg = hn::Lt(vLhs1, vZero);
  auto candMain1 = hn::Gt(hn::Sub(hn::Neg(vDstarA_plain), vC), vZero);

  auto vXmax = VecModel::vXMax(d, vA, vB);
  auto maskXmaxInf = hn::Eq(vXmax, hn::Set(d, kInf));
  auto vThetaStar = hn::Div(hn::Neg(vCD), vDiffAB);

  auto candInfBranch = hn::Gt(hn::Sub(hn::Set(d, -VecModel::DstarPrimeAtInfinity()), vCD), vZero);

  auto vThetaAtXmax = hn::Add(vA, hn::Mul(vXmax, vDiffAB));
  auto maskAGtB = hn::Gt(vA, vB);
  auto condGe = hn::Ge(vThetaStar, vThetaAtXmax);
  auto condLe = hn::Le(vThetaStar, vThetaAtXmax);
  auto maskOvershoot = hn::Or(hn::And(maskAGtB, condGe), hn::AndNot(maskAGtB, condLe));
  auto vDstarAtXmax = hn::IfThenElse(maskAGtB, vDstarRight, vDstarLeft);
  auto candOvershootVal = hn::Gt(hn::Sub(hn::Neg(vDstarAtXmax), hn::Add(vC, hn::Mul(vXmax, vCD))), vZero);
  auto candFiniteBranch = hn::And(maskOvershoot, candOvershootVal);

  auto candXmaxBranch = hn::Or(hn::And(maskXmaxInf, candInfBranch), hn::AndNot(maskXmaxInf, candFiniteBranch));

  auto vCostEvalVal = VecModel::vCostEval(d, vThetaStar, vA); // costEval's "mean" arg is exactly `a`
  auto candBottom = hn::Gt(vCostEvalVal, vC);

  auto candElseBranch = hn::Or(candXmaxBranch, candBottom);
  auto candGeneral = hn::Or(hn::And(maskMainNeg, candMain1), hn::AndNot(maskMainNeg, candElseBranch));

  auto maskPrune = hn::Or(hn::And(maskLinear, maskSpecialPrune), hn::AndNot(maskLinear, candGeneral));

  double buf[hn::MaxLanes(D())];
  hn::StoreU(hn::IfThenElse(maskPrune, hn::Set(d, 1.0), vZero), d, buf);
  size_t N = hn::Lanes(d);
  for (size_t L = 0; L < N; L++) prune[off + L] = (buf[L] != 0.0);
}

// Other models evaluate the same DUST test on contiguous candidate arrays.
template <class VecModel>
void dust_lane_scalar(double cumsum_t, double minCost_t, double t_pos,
                       const double* A_cumsum, const double* A_cost, const double* A_pos,
                       size_t off, size_t end, uint8_t* prune)
{
  using ScalarModel = typename ScalarModelOf<VecModel>::type;
  for (size_t i = off; i < end; i++)
  {
    double a = (cumsum_t - A_cumsum[i]) / (t_pos - A_pos[i]);
    double b = (A_cumsum[i] - A_cumsum[i - 1]) / (A_pos[i] - A_pos[i - 1]);
    double c = (minCost_t - A_cost[i]) / (t_pos - A_pos[i]);
    double d = (A_cost[i] - A_cost[i - 1]) / (A_pos[i] - A_pos[i - 1]);
    prune[i] = DualMax_DUST<ScalarModel>::test_abcd(a, b, c, d);
  }
}

// All three formulations use the same Highway active-set schedule.
template<class VecModel,int Version,class D>
void dustib_lane(D d,double sum,double qt,double t,const double* sums,
                 const double* costs,const double* pos,size_t off,uint8_t* prune) {
  if constexpr (Version==3) {
    for(size_t i=0;i<hn::Lanes(d);++i) {
      const size_t j=off+i;
      const double a=(sum-sums[j])/(t-pos[j]);
      const double b=(sums[j]-sums[j-1])/(pos[j]-pos[j-1]);
      const double c=(qt-costs[j])/(t-pos[j]);
      const double q=(costs[j]-costs[j-1])/(pos[j]-pos[j-1]);
      prune[j]=dustib::test<dustib::vec_id<VecModel>,2>(a,b,c,q);
    }
    return;
  }
  auto cs=hn::LoadU(d,sums+off),cr=hn::LoadU(d,sums+off-1);
  auto qs=hn::LoadU(d,costs+off),qr=hn::LoadU(d,costs+off-1);
  auto ps=hn::LoadU(d,pos+off),pr=hn::LoadU(d,pos+off-1);
  auto n=hn::Sub(hn::Set(d,t),ps),l=hn::Sub(ps,pr);
  auto a=hn::Div(hn::Sub(hn::Set(d,sum),cs),n);
  auto b=hn::Div(hn::Sub(cs,cr),l);
  auto c=hn::Div(hn::Sub(hn::Set(d,qt),qs),n);
  auto q=hn::Div(hn::Sub(qs,qr),l);
  auto mask=dustib::vector_test<dustib::vec_id<VecModel>,Version>(d,a,b,c,q);
  double out[hn::MaxLanes(D())];
  hn::StoreU(hn::IfThenElse(mask,hn::Set(d,1.),hn::Zero(d)),d,out);
  for(size_t i=0;i<hn::Lanes(d);++i) prune[off+i]=(out[i]!=0);
}

// The PELT test needs only the candidate index s.
template <class VecModel, class D>
void pelt_lane(D d, double cumsum_t, double minCost_t, double t_pos,
               const double* A_cumsum, const double* A_cost, const double* A_pos,
               size_t off, uint8_t* prune)
{
  auto vDt = hn::Sub(hn::Set(d, t_pos), hn::LoadU(d, A_pos + off));
  auto vA  = hn::Div(hn::Sub(hn::Set(d, cumsum_t), hn::LoadU(d, A_cumsum + off)), vDt);
  auto vC  = hn::Div(hn::Sub(hn::Set(d, minCost_t), hn::LoadU(d, A_cost + off)), vDt);
  auto maskPrune = hn::Gt(hn::Sub(hn::Neg(vC), VecModel::vDstar(d, vA)), hn::Zero(d));
  double buf[hn::MaxLanes(D())];
  hn::StoreU(hn::IfThenElse(maskPrune, hn::Set(d, 1.0), hn::Zero(d)), d, buf);
  size_t N = hn::Lanes(d);
  for (size_t L = 0; L < N; L++) prune[off + L] = (buf[L] != 0.0);
}

// Smallest-surviving-index check (dust's unconditional "prune last"):
// applies for every method, computed here via the same vCostMin used by
// the OP-step scan, so there is exactly one costMin formula per model,
// not two.
template <class VecModel>
bool smallest_index_prune(double cumsum_t, double minCost_t, double t_pos,
                           double A_cumsum0, double A_cost0, double A_pos0)
{
  const hn::CappedTag<double,1> d1;
  auto vDiff = hn::Set(d1, cumsum_t - A_cumsum0);
  auto vDt   = hn::Set(d1, t_pos - A_pos0);
  auto vLastCost = hn::Add(hn::Set(d1, A_cost0), VecModel::vCostMin(d1, vDiff, vDt));
  double lastCost; hn::StoreU(vLastCost, d1, &lastCost);
  return lastCost > minCost_t;
}

struct OpResult { double minCost; size_t argPos; };

template <class VecModel>
OpResult op_step_scan(double cumsum_t, double t_pos,
                       const double* A_cumsum, const double* A_cost, const double* A_pos,
                       size_t k)
{
  const hn::ScalableTag<double> d;
  const hn::ScalableTag<int64_t> di;
  const size_t N = hn::Lanes(d);

  auto vBestVal = hn::Set(d, std::numeric_limits<double>::infinity());
  auto vBestIdx = hn::Zero(di);
  const auto vCumsumT = hn::Set(d, cumsum_t);
  const auto vTpos    = hn::Set(d, t_pos);

  size_t i = 0;
  for (; i + N <= k; i += N)
  {
    auto vAc = hn::LoadU(d, A_cumsum + i);
    auto vAq = hn::LoadU(d, A_cost + i);
    auto vAp = hn::LoadU(d, A_pos + i);
    auto vDiff = hn::Sub(vCumsumT, vAc);
    auto vDt   = hn::Sub(vTpos, vAp);
    auto vVal  = hn::Add(vAq, VecModel::vCostMin(d, vDiff, vDt));
    auto vIdx  = hn::Iota(di, (int64_t)i);
    auto mask  = hn::Lt(vVal, vBestVal);
    vBestVal   = hn::IfThenElse(mask, vVal, vBestVal);
    vBestIdx   = hn::IfThenElse(hn::RebindMask(di, mask), vIdx, vBestIdx);
  }

  double valBuf[hn::MaxLanes(d)];
  int64_t idxBuf[hn::MaxLanes(d)];
  hn::StoreU(vBestVal, d, valBuf);
  hn::StoreU(vBestIdx, di, idxBuf);

  double bestVal = valBuf[0];
  int64_t bestIdx = idxBuf[0];
  for (size_t L = 1; L < N; L++) if (valBuf[L] < bestVal) { bestVal = valBuf[L]; bestIdx = idxBuf[L]; }

  const hn::CappedTag<double,1> d1;
  for (; i < k; i++)
  {
    auto vDiff = hn::Set(d1, cumsum_t - A_cumsum[i]);
    auto vDt   = hn::Set(d1, t_pos - A_pos[i]);
    auto vVal  = hn::Add(hn::Set(d1, A_cost[i]), VecModel::vCostMin(d1, vDiff, vDt));
    double val; hn::StoreU(vVal, d1, &val);
    if (val < bestVal) { bestVal = val; bestIdx = (int64_t)i; }
  }
  return {bestVal, (size_t)bestIdx};
}

void compact(std::vector<uint32_t>& active, std::vector<double>& A_cumsum,
             std::vector<double>& A_cost, std::vector<double>& A_pos,
             const std::vector<uint8_t>& drop)
{
  size_t w = 0;
  for (size_t i = 0; i < active.size(); i++)
  {
    if (drop[i]) continue;
    active[w] = active[i]; A_cumsum[w] = A_cumsum[i]; A_cost[w] = A_cost[i]; A_pos[w] = A_pos[i];
    w++;
  }
  active.resize(w); A_cumsum.resize(w); A_cost.resize(w); A_pos.resize(w);
}

// methodCode: 0 = DUST, 1 = PELT, 2 = OP.
// Keep both one-shot and incremental calls to hw_step inlined.

#if defined(__GNUC__) || defined(__clang__)
  #define DUST_HW_FORCE_INLINE inline __attribute__((always_inline))
#else
  #define DUST_HW_FORCE_INLINE inline
#endif

// One segmentation step, shared by one-shot and incremental interfaces.
template <class VecModel, int Version = DUSTIB_VERSION>
DUST_HW_FORCE_INLINE void hw_step(unsigned t, double penalty, int methodCode,
                     const std::vector<double>& cumsum,
                     std::vector<double>& costRecord, std::vector<int>& chptRecord, std::vector<int>& nb,
                     std::vector<uint32_t>& active, std::vector<double>& A_cumsum,
                     std::vector<double>& A_cost, std::vector<double>& A_pos,
                     std::vector<uint8_t>& drop,
                     const hn::ScalableTag<double>& d, const hn::CappedTag<double,1>& d1, size_t N)
{
  double t_pos = (double)t;
  OpResult op = op_step_scan<VecModel>(cumsum[t], t_pos, A_cumsum.data(), A_cost.data(), A_pos.data(), active.size());
  double minCost_t = op.minCost + penalty;
  costRecord.push_back(minCost_t);
  chptRecord.push_back((int)active[op.argPos]);

  for (;;)
  {
    size_t k = active.size();
    drop.assign(k, 0);

    // The oldest candidate has no smaller constraint. Apply its separate
    // domain-guarded pruning test for every method.
    if (k >= 1)
      drop[0] = smallest_index_prune<VecModel>(cumsum[t], minCost_t, t_pos, A_cumsum[0], A_cost[0], A_pos[0]);

    if (methodCode == 1) // PELT: general test at [1,k), no constraint needed
    {
      size_t i = 1;
      for (; i + N <= k; i += N) pelt_lane<VecModel>(d, cumsum[t], minCost_t, t_pos, A_cumsum.data(), A_cost.data(), A_pos.data(), i, drop.data());
      for (; i < k; i++) pelt_lane<VecModel>(d1, cumsum[t], minCost_t, t_pos, A_cumsum.data(), A_cost.data(), A_pos.data(), i, drop.data());
    }
    else if (methodCode == 0) // DUST: see DustUsesVectorTest above for the model split
    {
      if constexpr (DustUsesVectorTest<VecModel>::value)
      {
        size_t i = 1;
        for (; i + N <= k; i += N) dust_lane<VecModel>(d, cumsum[t], minCost_t, t_pos, A_cumsum.data(), A_cost.data(), A_pos.data(), i, drop.data());
        for (; i < k; i++) dust_lane<VecModel>(d1, cumsum[t], minCost_t, t_pos, A_cumsum.data(), A_cost.data(), A_pos.data(), i, drop.data());
      }
      else
      {
        dust_lane_scalar<VecModel>(cumsum[t], minCost_t, t_pos, A_cumsum.data(), A_cost.data(), A_pos.data(), 1, k, drop.data());
      }
    }
    else if (methodCode == 3) {
      size_t i=1;
      for (;i+N<=k;i+=N) dustib_lane<VecModel,Version>(d,cumsum[t],minCost_t,t_pos,A_cumsum.data(),A_cost.data(),A_pos.data(),i,drop.data());
      for (;i<k;++i) dustib_lane<VecModel,Version>(d1,cumsum[t],minCost_t,t_pos,A_cumsum.data(),A_cost.data(),A_pos.data(),i,drop.data());
    }
    // methodCode == 2 (OP): smallest-index check only, already done above.

    bool any = false;
    for (size_t i = 0; i < k; i++) if (drop[i]) { any = true; break; }
    if (!any) break;
    compact(active, A_cumsum, A_cost, A_pos, drop);
    if (active.empty() || methodCode != 0) break; // only DUST iterates to a fixed point
  }

  nb.push_back((int)active.size());
  active.push_back(t);
  A_cumsum.push_back(cumsum[t]);
  A_cost.push_back(minCost_t);
  A_pos.push_back(t_pos);
}

// Optimal Partitioning scans all candidates without removing any.
template <class VecModel>
DUST_HW_FORCE_INLINE void hw_op_step(unsigned t, double penalty,
                     const std::vector<double>& cumsum,
                     std::vector<double>& costRecord, std::vector<int>& chptRecord, std::vector<int>& nb,
                     std::vector<uint32_t>& active, std::vector<double>& A_cumsum,
                     std::vector<double>& A_cost, std::vector<double>& A_pos)
{
  double t_pos = (double)t;
  OpResult op = op_step_scan<VecModel>(cumsum[t], t_pos, A_cumsum.data(), A_cost.data(), A_pos.data(), active.size());
  double minCost_t = op.minCost + penalty;
  costRecord.push_back(minCost_t);
  chptRecord.push_back((int)active[op.argPos]);

  nb.push_back((int)active.size());
  active.push_back(t);
  A_cumsum.push_back(cumsum[t]);
  A_cost.push_back(minCost_t);
  A_pos.push_back(t_pos);
}

// Shared backtracking tail, also previously duplicated between run_HW()
// and DUST_1D_HW_T::get_partition(). Doesn't depend on VecModel.
List hw_get_partition(unsigned n, const std::vector<int>& chptRecord,
                       const std::vector<uint32_t>& active, const std::vector<int>& nb,
                       const std::vector<double>& costRecord)
{
  std::vector<int> changepoints;
  for (int tau = chptRecord[n]; tau != 0; tau = chptRecord[tau]) changepoints.push_back(tau);
  std::reverse(changepoints.begin(), changepoints.end());
  changepoints.push_back((int)n);

  std::vector<unsigned int> lastIndexSet(active.rbegin(), active.rend());

  return List::create(
    _["changepoints"] = changepoints,
    _["lastIndexSet"] = lastIndexSet,
    _["backend"] = "highway",
    _["nb"] = nb,
    _["costQ"] = std::vector<double>(costRecord.begin() + 1, costRecord.end())
  );
}

template <class VecModel, int Version = DUSTIB_VERSION>
List run_HW(NumericVector data, double penalty, int methodCode)
{
  for (double y : data) dustib::Math<dustib::vec_id<VecModel>>::validate(y);
  const size_t n = data.size();
  std::vector<double> cumsum(n + 1, 0.0);
  for (size_t i = 0; i < n; i++) cumsum[i+1] = cumsum[i] + VecModel::statistic(data[i]);

  std::vector<double> costRecord; costRecord.reserve(n + 1);
  std::vector<int> chptRecord; chptRecord.reserve(n + 1);
  std::vector<int> nb; nb.reserve(n);
  costRecord.push_back(-penalty);
  chptRecord.push_back(0);

  std::vector<uint32_t> active; active.reserve(n + 1); active.push_back(0);
  std::vector<double> A_cumsum; A_cumsum.reserve(n + 1); A_cumsum.push_back(0.0);
  std::vector<double> A_cost; A_cost.reserve(n + 1); A_cost.push_back(-penalty);
  std::vector<double> A_pos; A_pos.reserve(n + 1); A_pos.push_back(0.0);
  std::vector<uint8_t> drop; drop.reserve(n + 1);

  const hn::ScalableTag<double> d;
  const hn::CappedTag<double,1> d1;
  const size_t N = hn::Lanes(d);

  for (size_t t = 1; t <= n; t++)
    hw_step<VecModel,Version>((unsigned)t, penalty, methodCode, cumsum, costRecord, chptRecord, nb,
                       active, A_cumsum, A_cost, A_pos, drop, d, d1, N);

  return hw_get_partition((unsigned)n, chptRecord, active, nb, costRecord);
}

// Keep the existing DUST/PELT kernels unchanged. DUSTib uses the
// benchmark-selected hybrid: Highway cost scan and scalar inequality tests.
template<class VecModel>
List run_HW_chosen(NumericVector data,double penalty,int methodCode) {
  if(methodCode==3) return run_HW<VecModel,3>(data,penalty,methodCode);
  return run_HW<VecModel,2>(data,penalty,methodCode);
}

// Standalone, genuinely unpruned one-shot OP kernel (hw_op_step above) --
// used by DUST.1D.HW(method="OP") instead of run_HW(..., methodCode=2).
template <class VecModel>
List run_OP_HW(NumericVector data, double penalty)
{
  for (double y : data) dustib::Math<dustib::vec_id<VecModel>>::validate(y);
  const size_t n = data.size();
  std::vector<double> cumsum(n + 1, 0.0);
  for (size_t i = 0; i < n; i++) cumsum[i+1] = cumsum[i] + VecModel::statistic(data[i]);

  std::vector<double> costRecord; costRecord.reserve(n + 1);
  std::vector<int> chptRecord; chptRecord.reserve(n + 1);
  std::vector<int> nb; nb.reserve(n);
  costRecord.push_back(-penalty);
  chptRecord.push_back(0);

  std::vector<uint32_t> active; active.reserve(n + 1); active.push_back(0);
  std::vector<double> A_cumsum; A_cumsum.reserve(n + 1); A_cumsum.push_back(0.0);
  std::vector<double> A_cost; A_cost.reserve(n + 1); A_cost.push_back(-penalty);
  std::vector<double> A_pos; A_pos.reserve(n + 1); A_pos.push_back(0.0);

  for (size_t t = 1; t <= n; t++)
    hw_op_step<VecModel>((unsigned)t, penalty, cumsum, costRecord, chptRecord, nb, active, A_cumsum, A_cost, A_pos);

  return hw_get_partition((unsigned)n, chptRecord, active, nb, costRecord);
}

// --------------------------------------------------------------------- //
// Incremental Highway object. append_data stores sufficient statistics;
// update_partition resumes the dynamic program at the next observation.
class DUST_1D_HW_Base
{
  public:
    virtual ~DUST_1D_HW_Base() {}
    virtual void append_data(const NumericVector& data, Nullable<double> inPenalty) = 0;
    virtual void update_partition() = 0;
    virtual List get_partition() = 0;
    virtual List get_info() = 0;
    List dust(const NumericVector& data, Nullable<double> inPenalty)
    {
      append_data(data, inPenalty);
      update_partition();
      return get_partition();
    }
};

template <class VecModel>
class DUST_1D_HW_T : public DUST_1D_HW_Base
{
  public:
    DUST_1D_HW_T(int methodCode_, std::string model_name_, std::string method_name_)
      : methodCode(methodCode_), model_name(std::move(model_name_)), method_name(std::move(method_name_)) {}

    void append_data(const NumericVector& data, Nullable<double> inPenalty) override
    {
      const size_t m = (size_t)data.size();
      if (m == 0) return;
      for (double y : data) dustib::Math<dustib::vec_id<VecModel>>::validate(y);
      bool first_execution = (n == 0);
      const double checked_penalty = dust_validation::penalty(
        inPenalty, penalty, m, !first_execution);
      n += (unsigned int)m;

      if (first_execution)
      {
        penalty = checked_penalty;
        cumsum.push_back(0.0);
        costRecord.push_back(-penalty);
        chptRecord.push_back(0);
        active.push_back(0);
        A_cumsum.push_back(0.0);
        A_cost.push_back(-penalty);
        A_pos.push_back(0.0);
      }

      for (size_t i = 0; i < m; i++) cumsum.push_back(cumsum.back() + VecModel::statistic(data[i]));
    }

    void update_partition() override
    {
      if (n == 0 || active.empty()) return; // nothing appended yet

      const hn::ScalableTag<double> d;
      const hn::CappedTag<double,1> d1;
      const size_t N = hn::Lanes(d);
      std::vector<uint8_t> drop;

      for (unsigned t = active.back() + 1; t <= n; t++) {
        if(methodCode==3)
          hw_step<VecModel,3>(t,penalty,methodCode,cumsum,costRecord,chptRecord,nb,
                              active,A_cumsum,A_cost,A_pos,drop,d,d1,N);
        else
          hw_step<VecModel,2>(t,penalty,methodCode,cumsum,costRecord,chptRecord,nb,
                              active,A_cumsum,A_cost,A_pos,drop,d,d1,N);
      }
    }

    List get_partition() override
    {
      if (n == 0) stop("append data before requesting a partition");
      if (chptRecord.size() != static_cast<std::size_t>(n) + 1)
        stop("update_partition before requesting a partition");
      return hw_get_partition(n, chptRecord, active, nb, costRecord);
    }

    List get_info() override
    {
      return List::create(
        _["backend"] = "highway",
        _["data_statistic"] = cumsum,
        _["data_length"] = n,
        _["current_penalty"] = penalty,
        _["model"] = model_name,
        _["pruning_algo"] = method_name
      );
    }

  private:
    std::vector<double> cumsum, costRecord;
    std::vector<int> chptRecord, nb;
    std::vector<uint32_t> active;
    std::vector<double> A_cumsum, A_cost, A_pos;
    double penalty = 0.0;
    unsigned int n = 0;
    int methodCode; // 0=DUST, 1=PELT (OP has its own class, DUST_1D_HW_OP_T, below)
    std::string model_name, method_name;
};

// Standalone, genuinely unpruned incremental OP (hw_op_step above) --
// parallels DUST_1D_HW_T but never drops anything from the active arrays,
// so it needs no methodCode and no drop/compact bookkeeping at all.
template <class VecModel>
class DUST_1D_HW_OP_T : public DUST_1D_HW_Base
{
  public:
    DUST_1D_HW_OP_T(std::string model_name_, std::string method_name_)
      : model_name(std::move(model_name_)), method_name(std::move(method_name_)) {}

    void append_data(const NumericVector& data, Nullable<double> inPenalty) override
    {
      const size_t m = (size_t)data.size();
      if (m == 0) return;
      for (double y : data) dustib::Math<dustib::vec_id<VecModel>>::validate(y);
      bool first_execution = (n == 0);
      const double checked_penalty = dust_validation::penalty(
        inPenalty, penalty, m, !first_execution);
      n += (unsigned int)m;

      if (first_execution)
      {
        penalty = checked_penalty;
        cumsum.push_back(0.0);
        costRecord.push_back(-penalty);
        chptRecord.push_back(0);
        active.push_back(0);
        A_cumsum.push_back(0.0);
        A_cost.push_back(-penalty);
        A_pos.push_back(0.0);
      }

      for (size_t i = 0; i < m; i++) cumsum.push_back(cumsum.back() + VecModel::statistic(data[i]));
    }

    void update_partition() override
    {
      if (n == 0 || active.empty()) return; // nothing appended yet

      for (unsigned t = active.back() + 1; t <= n; t++)
        hw_op_step<VecModel>(t, penalty, cumsum, costRecord, chptRecord, nb, active, A_cumsum, A_cost, A_pos);
    }

    List get_partition() override
    {
      if (n == 0) stop("append data before requesting a partition");
      if (chptRecord.size() != static_cast<std::size_t>(n) + 1)
        stop("update_partition before requesting a partition");
      return hw_get_partition(n, chptRecord, active, nb, costRecord);
    }

    List get_info() override
    {
      return List::create(
        _["backend"] = "highway",
        _["data_statistic"] = cumsum,
        _["data_length"] = n,
        _["current_penalty"] = penalty,
        _["model"] = model_name,
        _["pruning_algo"] = method_name
      );
    }

  private:
    std::vector<double> cumsum, costRecord;
    std::vector<int> chptRecord, nb;
    std::vector<uint32_t> active;
    std::vector<double> A_cumsum, A_cost, A_pos;
    double penalty = 0.0;
    unsigned int n = 0;
    std::string model_name, method_name;
};

template <class VecModel>
DUST_1D_HW_Base* newHWModuleByAlgo(const std::string& method, const std::string& model)
{
  if (method == "DUST")   return new DUST_1D_HW_T<VecModel>(0, model, method);
  if (method == "DUSTib") return new DUST_1D_HW_T<VecModel>(3, model, method);
  if (method == "PELT")   return new DUST_1D_HW_T<VecModel>(1, model, method);
  if (method == "OP")     return new DUST_1D_HW_OP_T<VecModel>(model, method);
  stop("dust.object.1D.HW: method \"" + method + "\" has no Highway kernel "
       "(use \"DUST\", \"DUSTib\", \"PELT\" or \"OP\" here; scalar "
       "fallback are handled via dust.object.1D() at the R level).");
}

DUST_1D_HW_Base* newModule1D_HW(const std::string& model, const std::string& method)
{
  if (model == "gauss")    return newHWModuleByAlgo<VecGaussPolicy>(method, model);
  if (model == "poisson")  return newHWModuleByAlgo<VecPoissonPolicy>(method, model);
  if (model == "exp")      return newHWModuleByAlgo<VecExpPolicy>(method, model);
  if (model == "geom")     return newHWModuleByAlgo<VecGeomPolicy>(method, model);
  if (model == "bern")     return newHWModuleByAlgo<VecBernPolicy>(method, model);
  if (model == "binom")    return newHWModuleByAlgo<VecBinomPolicy>(method, model);
  if (model == "negbin")   return newHWModuleByAlgo<VecNegbinPolicy>(method, model);
  if (model == "variance") return newHWModuleByAlgo<VecVariancePolicy>(method, model);
  stop("Unrecognized model \"" + model + "\".");
}

} // namespace

#else // !HAVE_HIGHWAY
// Keep the Rcpp module registered in scalar-only builds.
namespace {

class DUST_1D_HW_Base
{
  public:
    virtual ~DUST_1D_HW_Base() {}
    virtual void append_data(const NumericVector& data, Nullable<double> inPenalty) = 0;
    virtual void update_partition() = 0;
    virtual List get_partition() = 0;
    virtual List get_info() = 0;
    List dust(const NumericVector& data, Nullable<double> inPenalty)
    {
      append_data(data, inPenalty);
      update_partition();
      return get_partition();
    }
};

class DUST_1D_HW_Stub : public DUST_1D_HW_Base
{
  public:
    DUST_1D_HW_Stub()
    {
      stop("This dust build has no Highway backend. dust.object.1D.HW() "
           "should have redirected to dust.object.1D() automatically -- "
           "seeing this means that redirect was bypassed.");
    }
    void append_data(const NumericVector&, Nullable<double>) override {}
    void update_partition() override {}
    List get_partition() override { return List(); }
    List get_info() override { return List(); }
};

DUST_1D_HW_Base* newModule1D_HW(const std::string&, const std::string&)
{
  return new DUST_1D_HW_Stub();
}

} // namespace
#endif // HAVE_HIGHWAY

//' @title Highway segmentation module
//'
//' @name DUST_1D_HW_Obj
//'
//' @description
//' Rcpp module for the Highway segmentation object used by
//' \code{\link{dust.object.1D}}.
//'
//' @export
RCPP_MODULE(DUSTHWMODULE1D)
{
  class_<DUST_1D_HW_Base>("DUST_1D_HW_Obj")
    .factory<const std::string&, const std::string&>(newModule1D_HW)
    .method("append_data", &DUST_1D_HW_Base::append_data)
    .method("update_partition", &DUST_1D_HW_Base::update_partition)
    .method("get_partition", &DUST_1D_HW_Base::get_partition)
    .method("get_info", &DUST_1D_HW_Base::get_info)
    .method("dust", &DUST_1D_HW_Base::dust)
  ;
}

//' DUST.1D.HW.backend
//'
//' @description Reports whether DUST.1D.HW was built with Google Highway
//' SIMD ("highway") or the scalar engine ("scalar"). DUSTib supports both.
//' @return A length-1 character vector.
//' @keywords internal
// [[Rcpp::export(name = "DUST.1D.HW.backend")]]
CharacterVector DUST_1D_HW_backend()
{
#ifdef HAVE_HIGHWAY
  return CharacterVector::create("highway");
#else
  return CharacterVector::create("scalar");
#endif
}

//' DUST.1D.HW
//'
//' @description
//' Same parameters and options as \code{\link{dust.1D}} -- all 8 cost
//' models, all 4 pruning methods -- implemented using
//' Google Highway SIMD for DUST, DUSTib, PELT and OP, with a scalar
//' fallback when Highway is unavailable. DUSTib implements a tested
//' one-constraint inequality certificate. The pruning schedule can change
//' candidate counts and the choice among tied optimal segmentations.
//'
//' @param data Numeric vector, univariate time series.
//' @param penalty Finite nonnegative penalty per change point. Default \code{2*log(length(data))}.
//' @param model One of "gauss", "poisson", "exp", "geom", "bern", "binom", "negbin", "variance".
//' @param method One of "DUST", "DUSTib", "PELT", "OP".
//' @return list(changepoints, lastIndexSet, backend, nb, costQ)
//' @examples
//' ### Gaussian series, 4 segments of 500 points, 3 mean shifts.
//' ### Compare the two backends through the public interface.
//' set.seed(20)
//' true_chpts <- c(500, 1000, 1500, 2000)
//' y <- dataGenerator_1D(chpts = true_chpts, parameters = c(0, 1, -1, 2),
//'                        sdNoise = 1, type = "gauss")
//' y <- data_normalization_1D(y, type = "gauss")
//' penalty <- 2 * log(length(y))
//'
//' res_hw <- dust.1D(y, penalty, model = "gauss", method = "DUST", backend = "highway")
//' res_hw$changepoints  # close to true_chpts
//'
//' res_ref <- dust.1D(y, penalty, model = "gauss", method = "DUST", backend = "scalar")
//' identical(as.integer(res_hw$changepoints), as.integer(res_ref$changepoints))
//' @keywords internal
// [[Rcpp::export(name = "DUST.1D.HW")]]
List DUST_1D_HW(NumericVector data, Nullable<double> penalty = R_NilValue,
                std::string model = "gauss", std::string method = "DUST")
{
  if (data.size() == 0) stop("data must be a nonempty numeric vector");
  double pen = dust_validation::penalty(penalty, 0.0, data.size(), false);

#ifdef HAVE_HIGHWAY

  if (method == "OP")
  {
    // Genuinely unpruned OP (run_OP_HW), not run_HW's methodCode==2 path
    // -- see 1D_DUST_Impl.h / _ModuleAssembly_1D.cpp for why that
    // distinction matters.
    if (model == "gauss")    return run_OP_HW<VecGaussPolicy>(data, pen);
    if (model == "poisson")  return run_OP_HW<VecPoissonPolicy>(data, pen);
    if (model == "exp")      return run_OP_HW<VecExpPolicy>(data, pen);
    if (model == "geom")     return run_OP_HW<VecGeomPolicy>(data, pen);
    if (model == "bern")     return run_OP_HW<VecBernPolicy>(data, pen);
    if (model == "binom")    return run_OP_HW<VecBinomPolicy>(data, pen);
    if (model == "negbin")   return run_OP_HW<VecNegbinPolicy>(data, pen);
    if (model == "variance") return run_OP_HW<VecVariancePolicy>(data, pen);
    stop("Unrecognized model \"" + model + "\".");
  }

  int methodCode;
  if (method == "DUST") methodCode = 0;
  else if (method == "PELT") methodCode = 1;
  else if (method == "DUSTib") methodCode = 3;
  else stop("Unrecognized method \"" + method + "\". Valid values: \"DUST\", \"DUSTib\", \"PELT\", \"OP\".");

  if (model == "gauss")    return run_HW_chosen<VecGaussPolicy>(data, pen, methodCode);
  if (model == "poisson")  return run_HW_chosen<VecPoissonPolicy>(data, pen, methodCode);
  if (model == "exp")      return run_HW_chosen<VecExpPolicy>(data, pen, methodCode);
  if (model == "geom")     return run_HW_chosen<VecGeomPolicy>(data, pen, methodCode);
  if (model == "bern")     return run_HW_chosen<VecBernPolicy>(data, pen, methodCode);
  if (model == "binom")    return run_HW_chosen<VecBinomPolicy>(data, pen, methodCode);
  if (model == "negbin")   return run_HW_chosen<VecNegbinPolicy>(data, pen, methodCode);
  if (model == "variance") return run_HW_chosen<VecVariancePolicy>(data, pen, methodCode);
  stop("Unrecognized model \"" + model + "\".");
#else
  return run_dust_by_model(data, pen, model, method);
#endif
}
