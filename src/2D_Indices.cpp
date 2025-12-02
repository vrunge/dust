#include <Rcpp.h>
#include <cmath>

#include "2D_Indices.h"

using namespace Rcpp;
// --------------------------- //
// --- /////////////////// --- //
// --- // Parent Module // --- //
// --- /////////////////// --- //
// --------------------------- //

Indices_2D::~Indices_2D() {}

void Indices_2D::reset(){current = list.begin();}
void Indices_2D::next(){++current;}
void Indices_2D::remove_first(){list.pop_front();}

bool Indices_2D::is_not_the_last(){return current != list.end();}

unsigned int Indices_2D::get_first(){return list.front();}
unsigned int Indices_2D::get_current(){return *current;}
std::forward_list<unsigned int> Indices_2D::get_list(){return list;}


////////////////////////////////////////////////////////////////////////////////
// ------------------------------------ //
// --- //////////////////////////// --- //
// --- // Deterministic Module 1 // --- //
// --- //////////////////////////// --- //
// ------------------------------------ //

void Indices_2D_Det1::add_first(unsigned int value){list.push_front(value);}
void Indices_2D_Det1::reset_pruning()
{
  before = list.before_begin(); // -1
  current = std::next(before); // 0
  constraint = std::next(current); // 1
}

void Indices_2D_Det1::next_pruning()
{
  before = current; // +1
  current = constraint; // +1
  new_constraint(); // +1
}

void Indices_2D_Det1::prune_current()
{
  current = list.erase_after(before); // remove the after before = current
  new_constraint(); // move constraint
}

void Indices_2D_Det1::prune_last(){list.erase_after(before);}
bool Indices_2D_Det1::is_not_the_last_pruning(){return constraint != list.end();}
void Indices_2D_Det1::new_constraint(){++constraint;}
unsigned int Indices_2D_Det1::get_constraint_r1(){return *constraint;}
unsigned int Indices_2D_Det1::get_constraint_r2(){return *constraint;}

////////////////////////////////////////////////////////////////////////////////
// ------------------------------------ //
// --- //////////////////////////// --- //
// --- // Deterministic Module 2 // --- //
// --- //////////////////////////// --- //
// ------------------------------------ //

void Indices_2D_Det2::add_first(unsigned int value){list.push_front(value);}
void Indices_2D_Det2::reset_pruning()
{
  before = list.before_begin(); // -1
  current = std::next(before); // 0
  constraint = std::next(current); // 1
}

void Indices_2D_Det2::next_pruning()
{
  before = current; // +1
  current = constraint; // +1
  new_constraint(); // +1
}

void Indices_2D_Det2::prune_current()
{
  current = list.erase_after(before); // remove the after before = current
  new_constraint(); // move constraint
}

void Indices_2D_Det2::prune_last(){list.erase_after(before);}
bool Indices_2D_Det2::is_not_the_last_pruning(){return constraint != list.end();}
void Indices_2D_Det2::new_constraint(){++constraint;}

unsigned int Indices_2D_Det2::get_constraint_r1(){return *constraint;}
unsigned int Indices_2D_Det2::get_constraint_r2()
{
  // Try to look at the next element
  auto next_it = std::next(constraint);
  if (next_it != list.end())
  {
    return *next_it;      // there is a "next" value
  }

  return *constraint;

}




