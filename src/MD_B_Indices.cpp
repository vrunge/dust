#include "MD_B_Indices.h"

using namespace Rcpp;

// --------- // Indices_MD // --------- //
//
// This class allows drawing elements of a forwardlist at random, uniformly, among
// the elements following the current position of the forwardlist iterator created
// by Indices_MD::reset_prune() and incremented by Indices_MD::next_prune()
// This class was tailored to be used in the DUST algorithm
//
// Parameters:


// Expliciter le constructeur, avec initialisation des 2 premiers indices
// Initiate handler and nb_max_ vector
Indices_MD::Indices_MD() : nb_max(1) {}
Indices_MD::Indices_MD(const unsigned int& nb_max_) : nb_max(nb_max_) { }

// Vérifier les fuites de mémoire
Indices_MD::~Indices_MD() {}

// simple reset for OP step
void Indices_MD::reset() { current = list.begin(); }

// simple next for OP step
void Indices_MD::next() { ++current; }

// simple break check for OP step
bool Indices_MD::check() { return current != list.end(); }

// add index
void Indices_MD::add(const unsigned int& value)
{
  list.push_back(value);
}

void Indices_MD::set_size(const unsigned int& size) { list.reserve(size); }

// obtain list for main function output
std::vector<unsigned int> Indices_MD::get_list() { return list; }

void Indices_MD::remove_last() { list.pop_back(); }

VariableIndices_MD::VariableIndices_MD() : Indices_MD() {}
VariableIndices_MD::VariableIndices_MD(const unsigned int& nb_max_) : Indices_MD(nb_max_) {}

// full reset for pruning step, including reinitialization of nb_max_ vector
void VariableIndices_MD::reset_prune() {
  // reset iterators
  if (list.size() > 1)
  {
    begin = list.begin();
    current = begin + 1;
  }
  else current = list.begin();
}

// full next for pruning step
void VariableIndices_MD::next_prune() {
  ++current;
  if (current - begin > nb_max) ++begin;
}

// remove current index and its pointer
void VariableIndices_MD::prune_current() { current = list.erase(current); }

// break check for pruning step
bool VariableIndices_MD::check_prune() { return current != list.end(); }

std::vector<unsigned int> VariableIndices_MD::get_constraints()
{
  return std::vector(begin, current);
}

RandomIndices_MD::RandomIndices_MD() : Indices_MD() {}
RandomIndices_MD::RandomIndices_MD(const unsigned int& nb_max_)
  : Indices_MD(nb_max_)
    , rng(std::random_device{}())
    , dist(std::uniform_real_distribution(0., 1.))
    // , dist2(std::uniform_int_distribution(0, 0))
  {
    // rng.seed(1);
  }

// full reset for pruning step, including reinitialization of nb_max_ vector
void RandomIndices_MD::reset_prune() {
  // reset iterators
  if (list.size() > 1)
  {
    current = list.begin() + 1;
  }
  else current = list.begin();
}

// full next for pruning step
void RandomIndices_MD::next_prune() { ++current; }

// remove current index and its pointer
void RandomIndices_MD::prune_current() { current = list.erase(current); }

// break check for pruning step
bool RandomIndices_MD::check_prune() { return current != list.end(); }

std::vector<unsigned int> RandomIndices_MD::get_constraints()
{
  // std::vector<unsigned int> constraints;
  // std::sample(
  //   list.begin(),
  //   current,
  //   std::back_inserter(constraints),
  //   nb_max,
  //   rng
  // );
  // Rcout << "r: ";
  // for (auto i = constraints.begin(); i != constraints.end(); ++i)
  //   Rcout << " " << *i;
  // Rcout << std::endl;
  // return constraints;

  // std::set<unsigned int> constraints;
  // for (unsigned int i = 0; i < nb_max; i++)
  //   constraints.insert(list[std::floor(dist(rng) * (current - list.begin()))]);
  // // Rcout << "r: ";
  // // for (auto i = constraints.begin(); i != constraints.end(); ++i)
  // //   Rcout << " " << *i;
  // // Rcout << std::endl;
  // return std::vector<unsigned int>(constraints.begin(), constraints.end());

  std::vector<unsigned int> constraints;
  constraints.reserve(nb_max);
  unsigned int nbC = current - list.begin();
  for (unsigned int i = 0; i < nb_max; i++)
    constraints.push_back(list[std::floor(dist(rng) * nbC)]);
  // Rcout << "r: ";
  // for (auto i = constraints.begin(); i != constraints.end(); ++i)
  //   Rcout << " " << *i;
  // Rcout << std::endl;
  return constraints;

  // std::vector<unsigned int> constraints;
  // constraints.reserve(nb_max);
  // dist2.param(std::uniform_int_distribution<int>::param_type(0, current - list.begin()));
  // // auto x = std::uniform_int_distribution<int>(0, current - list.begin());
  // for (unsigned int i = 0; i < nb_max; i++)
  //   constraints.push_back(list[dist2(rng)]);
  // // Rcout << "r: ";
  // // for (auto i = constraints.begin(); i != constraints.end(); ++i)
  // //   Rcout << " " << *i;
  // // Rcout << std::endl;
  // return constraints;
}
