#ifndef INDICES1C_H
#define INDICES1C_H

#include <forward_list>
#include <vector>
#include <limits>

/// active indices in decreasing order, link_[i] = next smaller index

class Indices_1D
{
  public:
    void add_first(unsigned int value)
    {
      link_.push_back(front_);
      front_ = value;
    }

    void reset() { current_ = front_; }
    void next() { current_ = link_[current_]; }
    bool is_not_the_last() { return current_ != END; }
    unsigned int get_first() { return front_; }
    unsigned int get_current() { return current_; }

    std::forward_list<unsigned int> get_list()
    {
      std::forward_list<unsigned int> out;
      auto it = out.before_begin();
      for (unsigned int i = front_; i != END; i = link_[i])
      {
        it = out.insert_after(it, i);
      }
      return out;
    }

    void reset_pruning() { before_ = END; current_ = front_; constraint_ = link_[front_]; }
    void next_pruning() { before_ = current_; current_ = link_[current_]; new_constraint(); }

    // remove current_
    void prune_current()
    {
      unsigned int next_elem = link_[current_];
      if (before_ == END) { front_ = next_elem; } else { link_[before_] = next_elem; }
      current_ = next_elem;
      new_constraint();
    }

    void prune_last()
    {
      unsigned int next_elem = link_[current_];
      if (before_ == END) { front_ = next_elem; } else { link_[before_] = next_elem; }
    }

    bool is_not_the_last_pruning() { return constraint_ != END; }

    void new_constraint() { constraint_ = link_[current_]; }
    unsigned int get_constraint() { return constraint_; }

  private:
    static constexpr unsigned int END = std::numeric_limits<unsigned int>::max();

    std::vector<unsigned int> link_;
    unsigned int front_ = END;    // largest active index
    unsigned int current_ = END;  // OP scan, then index s tested for pruning
    unsigned int before_;         // previous index in the chain (or END)
    unsigned int constraint_;     // index r (next smaller active index)
};

#endif
