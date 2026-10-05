#ifndef INDICES1C_H
#define INDICES1C_H

#include <forward_list>
#include <vector>
#include <limits>

// Active candidates are stored in decreasing order. link_[i] gives the
// next smaller candidate; each insertion and removal takes constant time.
// The constraint for candidate s is its next smaller active index.

class Indices_1D
{
  public:
    void add_first(unsigned int value)
    {
      // `value` is always link_.size() at the point of the call (indices
      // are appended in strictly increasing order starting at 0): plain
      // push_back, no heap allocation beyond the vector's own amortised
      // growth (reserved up-front by append_data).
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
      // Only called once, to build the R-facing `lastIndexSet` output:
      // O(k) in the (small) surviving-set size, materialising the same
      // front-to-back (largest-to-smallest) order as the original
      // forward_list did.
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

    // Relinks around current_, advances current_ to the element that takes
    // its place (mirrors forward_list::erase_after's returned iterator).
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
    unsigned int current_ = END;  // cursor: OP-step scan, then (same variable,
                                   // reused sequentially, never concurrently)
                                   // the index `s` currently tested for pruning
    unsigned int before_;         // predecessor of current_ in the link chain, or END
    unsigned int constraint_;     // index `r` = next-smaller active index
};

#endif
