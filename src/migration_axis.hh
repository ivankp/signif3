#ifndef IVANP_MIGRATION_AXIS_HH
#define IVANP_MIGRATION_AXIS_HH

#include "axis.hh"
#include "utility.hh"
#include <stdexcept>

template <typename EdgeType>
class migration_axis {
public:
  using edge_type = EdgeType;
  using edge_cref = ivanp::const_ref_if_not_scalar_t<edge_type>;
  using size_type = ivanp::axis_size_type;

private:
  const ivanp::abstract_axis<edge_type> *_axis;
  unsigned nchecks;

public:
  migration_axis() = default;
  ~migration_axis() = default;

  migration_axis(const ivanp::abstract_axis<edge_type>* axis, unsigned nchecks)
  : _axis(axis), nchecks(nchecks) { }

  inline size_type nedges() const { return _axis->nedges(); }
  inline size_type nbins() const { return _axis->nbins() + nchecks; }

  inline edge_cref edge(size_type i) const noexcept { return _axis->edge(i); }

  inline edge_cref min() const { return _axis->min(); }
  inline edge_cref max() const { return _axis->max(); }

  inline edge_cref lower() const { return _axis->lower(); }
  inline edge_cref upper() const { return _axis->upper(); }

  template <typename T, typename... B>
  size_type find_bin(const std::tuple<const T&, B...>& x) const {
    using namespace ivanp;
    const unsigned ncheck = call( check<B...>, x,
      index_sequence_tail<1,sizeof...(B)+1>{});
    if (__builtin_expect(sizeof...(B) != nchecks,0)) throw std::length_error(
      "sizeof...(B) != nchecks");
    unsigned bin;
    if (ncheck==sizeof...(B)) {
      bin = _axis->find_bin(std::get<0>(x));
      if (bin) bin += sizeof...(B);
    } else bin = ncheck + 1;
    return bin;
  }
};

#endif
