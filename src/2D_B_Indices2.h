#ifndef INDICES22C_H
#define INDICES22C_H

#include <Rcpp.h>
#include <forward_list>
#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL
using namespace Rcpp;

class Indices_2D2
{
public:
  virtual ~Indices_2D2();

  // --- // Methods // --- //
  void reset();
  void next();
  bool check();

  unsigned int get_current();
  std::forward_list<unsigned int> get_list();
  void remove_first();

  virtual void add(unsigned int value) = 0;

  virtual void reset_prune() = 0;
  virtual void next_prune() = 0;
  virtual void prune_current() = 0;
  virtual bool check_prune() = 0;

  virtual void prune_last() = 0;

  virtual void new_constraint() = 0;
  virtual unsigned int get_constraint() = 0;


  // --- // Fields // --- //
  std::forward_list<unsigned int> list;
  std::forward_list<unsigned int>::iterator current;
  std::forward_list<unsigned int>::iterator before;
};

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

class RandomIndices_2D2 : public Indices_2D2
{
public:
  RandomIndices_2D2(unsigned int size, double alpha = 1e-9);
  // ~I_Random() override;

  void add(unsigned int value) override;

  void reset_prune() override;
  void next_prune() override;
  void prune_current() override;
  bool check_prune() override;

  void prune_last() override;

  void new_constraint() override;
  unsigned int get_constraint() override;

private:
  unsigned int nb = 0;
  unsigned int nbC = 0;

  unsigned int* constraint;
  std::vector<unsigned int*> pointers;
  std::vector<unsigned int*>::reverse_iterator pointersCurrent;

  std::vector<double> randomU;
  std::vector<double>::iterator u;

};

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

class DeterministicIndices_2D2 : public Indices_2D2
{
public:
  void add(unsigned int value) override;

  void reset_prune() override;
  void next_prune() override;
  void prune_current() override;
  bool check_prune() override;

  void prune_last() override;

  void new_constraint() override;
  unsigned int get_constraint() override;

private:
  std::forward_list<unsigned int>::iterator constraint;
};

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

class Random2Indices_2D2 : public Indices_2D2
{
public:
  Random2Indices_2D2(unsigned int size, double alpha = 1e-9);
  // ~I_Random() override;

  void add(unsigned int value) override;

  void reset_prune() override;
  void next_prune() override;
  void prune_current() override;
  bool check_prune() override;

  void prune_last() override;

  void new_constraint() override;
  unsigned int get_constraint() override;

private:
  unsigned int nb = 0;
  unsigned int nbC = 0;

  unsigned int* constraint;
  std::vector<unsigned int*> pointers;
  std::vector<unsigned int*>::reverse_iterator pointersCurrent;

  std::vector<double> randomU;
  std::vector<double>::iterator u;

};

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

class Deterministic2Indices_2D2 : public Indices_2D2
{
public:
  void add(unsigned int value) override;

  void reset_prune() override;
  void next_prune() override;
  void prune_current() override;
  bool check_prune() override;

  void prune_last() override;

  void new_constraint() override;
  unsigned int get_constraint() override;

private:
  std::forward_list<unsigned int>::iterator constraint;
};

#endif
