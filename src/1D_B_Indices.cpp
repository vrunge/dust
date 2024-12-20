#include <Rcpp.h>
#include <cmath>

#include "1D_B_Indices.h"

using namespace Rcpp;

// --------------------------- //
// --- /////////////////// --- //
// --- // Parent Module // --- //
// --- /////////////////// --- //
// --------------------------- //

Indices_1D::~Indices_1D() {}


void Indices_1D::reset(){current = list.begin();}
void Indices_1D::next(){++current;}
bool Indices_1D::check(){return current != list.end();}
unsigned int Indices_1D::get_current(){return *current;}


std::forward_list<unsigned int> Indices_1D::get_list(){return list;}
void Indices_1D::remove_first(){list.pop_front();}


// --------------------------- //
// --- /////////////////// --- //
// --- // Random Module // --- //
// --- /////////////////// --- //
// --------------------------- //

RandomIndices_1D::RandomIndices_1D(unsigned int size, double alpha)
{
  // length of the random vector
  double k = std::max(2., ceil(pow(size, .2)));
  int len = std::log(alpha) / std::log(1 - 1/k);

  randomU = std::vector<double>(len, 0.);

  // Random number engine
  std::minstd_rand0  engine(std::random_device{}());
  std::uniform_real_distribution<double> dist(0.0, 1.0);

  for(int i = 0; i < len; ++i) {randomU[i] = dist(engine);}
  u = randomU.begin();
}

/////////////////////

///
/// Simon : possible update, save only the pointer
///
void RandomIndices_1D::add(unsigned int value)
{
  list.push_front(value);
  pointers.push_back(&list.front());
  nb++;
}

void RandomIndices_1D::reset_prune()
{
  current = list.begin();
  before = list.before_begin();
  pointersCurrent = pointers.rbegin(); //rbegin for reverse iterator

  nbC = nb - 1;
}

void RandomIndices_1D::next_prune()
{
  before = current;
  ++current;
  ++pointersCurrent;
  new_constraint();
}

void RandomIndices_1D::prune_current()
{
  current = list.erase_after(before);
  pointersCurrent = std::vector<unsigned int*>::reverse_iterator(pointers.erase(std::next(pointersCurrent).base()));
  nb--;
  new_constraint();
}

// --- // If no constraint can be selected, exit loop // --- //
bool RandomIndices_1D::check_prune()
{
  return nbC > 0;
}

void RandomIndices_1D::prune_last()
{
  list.erase_after(before);
  pointers.erase(std::next(pointersCurrent).base());
  nb--;
}

// --- // Select new constraint // --- //
void RandomIndices_1D::new_constraint()
{
  nbC--;
}

unsigned int RandomIndices_1D::get_constraint()
{
  constraint = pointers[floor(nbC * (*u))];
  ++u;
  if (u == randomU.end()){u = randomU.begin();}
  return *constraint;
}


// ---------------------------------- //
// --- ////////////////////////// --- //
// --- // Deterministic Module // --- //
// --- ////////////////////////// --- //
// ---------------------------------- //

void DeterministicIndices_1D::add(unsigned int value)
{
  list.push_front(value);
}

void DeterministicIndices_1D::reset_prune()
{
  before = list.before_begin();
  current = std::next(before);
  constraint = std::next(current);
}

void DeterministicIndices_1D::next_prune()
{
  before = current;
  current = constraint;
  new_constraint();
}

void DeterministicIndices_1D::prune_current()
{
  current = list.erase_after(before);
  new_constraint();
}

bool DeterministicIndices_1D::check_prune()
{
  return constraint != list.end();
}

void DeterministicIndices_1D::prune_last()
{
  list.erase_after(before);
}

void DeterministicIndices_1D::new_constraint()
{
  ++constraint;
}

unsigned int DeterministicIndices_1D::get_constraint()
{
  return *constraint;
}
