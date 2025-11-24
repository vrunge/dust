#include <Rcpp.h>
#include <cmath>

#include "1D_Indices.h"

using namespace Rcpp;
// --------------------------- //
// --- /////////////////// --- //
// --- // Parent Module // --- //
// --- /////////////////// --- //
// --------------------------- //

Indices_1D::~Indices_1D() {}

void Indices_1D::reset(){current = list.begin();}
void Indices_1D::next(){++current;}
void Indices_1D::remove_first(){list.pop_front();}

bool Indices_1D::is_not_the_last(){return current != list.end();}

unsigned int Indices_1D::get_first(){return list.front();}
unsigned int Indices_1D::get_current(){return *current;}
std::forward_list<unsigned int> Indices_1D::get_list(){return list;}


////////////////////////////////////////////////////////////////////////////////
// ---------------------------------- //
// --- ////////////////////////// --- //
// --- // Deterministic Module // --- //
// --- ////////////////////////// --- //
// ---------------------------------- //

void Indices_1D_Det::add_first(unsigned int value){list.push_front(value);}

void Indices_1D_Det::reset_pruning()
{
  before = list.before_begin(); // -1
  current = std::next(before); // 0
  constraint = std::next(current); // 1
}

void Indices_1D_Det::next_pruning()
{
  before = current; // +1
  current = constraint; // +1
  new_constraint(); // +1
}

void Indices_1D_Det::prune_current()
{
  current = list.erase_after(before); // remove the after before = current
  new_constraint(); // move constraint
}

void Indices_1D_Det::prune_last(){list.erase_after(before);}
bool Indices_1D_Det::is_not_the_last_pruning(){return constraint != list.end();}

void Indices_1D_Det::new_constraint(){++constraint;}
unsigned int Indices_1D_Det::get_constraint(){return *constraint;}



////////////////////////////////////////////////////////////////////////////////
// --------------------------- //
// --- /////////////////// --- //
// --- // Random Module // --- //
// --- /////////////////// --- //
// --------------------------- //

Indices_1D_Rand::Indices_1D_Rand() :
  engine(std::random_device{}()), dist(0.0, 1.0)
{}

/////////////////////

///
/// Simon : possible update, save only the pointer
///
void Indices_1D_Rand::add_first(unsigned int value)
{
  list.push_front(value);
  pointers.push_back(&list.front());
  nb++;
}

void Indices_1D_Rand::reset_pruning()
{
  current = list.begin();
  before = list.before_begin();
  pointersCurrent = pointers.rbegin(); //rbegin for reverse iterator

  nbC = nb - 1;
}

void Indices_1D_Rand::next_pruning()
{
  before = current;
  ++current;
  ++pointersCurrent;
  new_constraint();
}

void Indices_1D_Rand::prune_current()
{
  current = list.erase_after(before);
  pointersCurrent = std::vector<unsigned int*>::reverse_iterator(pointers.erase(std::next(pointersCurrent).base()));
  nb--;
  new_constraint();
}

void Indices_1D_Rand::prune_last()
{
  list.erase_after(before);
  pointers.erase(std::next(pointersCurrent).base());
  nb--;
}

// --- // If no constraint can be selected, exit loop // --- //
bool Indices_1D_Rand::is_not_the_last_pruning()
{
  return nbC > 0;
}


// --- // Select new constraint // --- //
void Indices_1D_Rand::new_constraint()
{
  nbC--;
}

unsigned int Indices_1D_Rand::get_constraint()
{
  return *pointers[floor(nbC * dist(engine))];
}

