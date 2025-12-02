#ifndef INDICES1C_H
#define INDICES1C_H

#include <Rcpp.h>

#include <forward_list>
#include <random> /// FOR RANDOM NUMBER IN DUAL EVAL

using namespace Rcpp;

class Indices_2D
{
  public:
    virtual ~Indices_2D();

    ////////// Generic //////////
    ////////// Generic //////////

    void reset(); //iterator current = start
    void next(); //current++
    void remove_first(); //remove the first element

    bool is_not_the_last();

    unsigned int get_first();
    unsigned int get_current();
    std::forward_list<unsigned int> get_list();

    ////////// Virtual //////////
    ////////// Virtual //////////

    virtual void add_first(unsigned int value) = 0;

    virtual void reset_pruning() = 0;
    virtual void next_pruning() = 0;
    virtual void prune_current() = 0;
    virtual void prune_last() = 0;

    virtual bool is_not_the_last_pruning() = 0;

    virtual void new_constraint() = 0;
    virtual unsigned int get_constraint_r1() = 0;
    virtual unsigned int get_constraint_r2() = 0;


  protected:
    // --- // Fields // --- //
    std::forward_list<unsigned int> list;
    std::forward_list<unsigned int>::iterator current;
    std::forward_list<unsigned int>::iterator before;
};


////////////////////////////////////////////////////////////////////////////////
////////////////
//////////////// derived class:  Indices_2D_Det
////////////////

class Indices_2D_Det1 : public Indices_2D
{
  public:
    void add_first(unsigned int value) override;

    void reset_pruning() override;
    void next_pruning() override;
    void prune_current() override;
    bool is_not_the_last_pruning() override;

    void prune_last() override;

    void new_constraint() override;
    unsigned int get_constraint_r1() override;
    unsigned int get_constraint_r2() override;

  private:
    // iterator for list
    // a kind of "after" or "next" pointer
    std::forward_list<unsigned int>::iterator constraint;
    //////////// RANDOM NUMBER GENERATOR ////////////
    std::minstd_rand0 engine;  // Random number engine
    std::uniform_real_distribution<double> dist;  // Uniform distribution [0, 1)

};

////////////////////////////////////////////////////////////////////////////////
////////////////
//////////////// derived class:  Indices_2D_Rand
////////////////

class Indices_2D_Det2 : public Indices_2D
{
public:
  void add_first(unsigned int value) override;

  void reset_pruning() override;
  void next_pruning() override;
  void prune_current() override;
  bool is_not_the_last_pruning() override;

  void prune_last() override;

  void new_constraint() override;
  unsigned int get_constraint_r1() override;
  unsigned int get_constraint_r2() override;

private:
  // iterator for list
  // a kind of "after" or "next" pointer
  std::forward_list<unsigned int>::iterator constraint;
  //////////// RANDOM NUMBER GENERATOR ////////////
  std::minstd_rand0 engine;  // Random number engine
  std::uniform_real_distribution<double> dist;  // Uniform distribution [0, 1)

};


#endif
