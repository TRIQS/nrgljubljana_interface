#include <triqs/mesh.hpp>
#include <triqs/mesh/utils.hpp>

#include <vector>
#include <string>
#include <algorithm>
#include <initializer_list>

namespace triqs::mesh {

  /**
   * A generic mesh of points with sortable values 
   */
  template <typename Value> class point_mesh {
    static_assert(std::totally_ordered<Value>);

    public:
    using value_t      = Value;
    using index_t      = long;
    using data_index_t = long;

    // --  data
    private:
    std::vector<value_t> _pts;
    size_t _mesh_hash = 0;

    // -------------------- Constructors -------------------
    public:
    /// Construct from a vector of points
    point_mesh(std::vector<value_t> pts) : _pts(std::move(pts)), _mesh_hash(std::accumulate(_pts.begin(), _pts.end(), value_t{0})) {
      if (not std::is_sorted(_pts.begin(), _pts.end())) TRIQS_RUNTIME_ERROR << "Point mesh must be constructed with a sorted list of points";
    }

    /// Initializer list constructor
    point_mesh(std::initializer_list<value_t> l) : point_mesh(std::vector<value_t>(l)) {}

    /// Default constructor
    point_mesh() = default;

    // -------------------- Comparison -------------------

    /// Mesh comparison
    bool operator==(point_mesh const &) const = default;
    bool operator!=(point_mesh const &) const = default;

    // --------------------  Mesh Point -------------------

    struct mesh_point_t {
      using mesh_t = point_mesh;

      public:
      long _index         = 0;
      long _data_index    = 0;
      uint64_t _mesh_hash = 0;
      value_t _value      = {};

      public:
      mesh_point_t() = default;
      mesh_point_t(long index, long data_index, uint64_t mesh_hash, double value)
         : _index(index), _data_index(data_index), _mesh_hash(mesh_hash), _value(value) {}

      /// The index of the mesh point
      [[nodiscard]] long index() const { return _index; }

      /// The data index of the mesh point
      [[nodiscard]] long data_index() const { return _data_index; }

      /// The value of the mesh point
      [[nodiscard]] value_t value() const { return _value; }

      /// The Hash for the mesh configuration
      [[nodiscard]] uint64_t mesh_hash() const noexcept { return _mesh_hash; }

      operator value_t() const { return _value; }

      // https://godbolt.org/z/xoYP3vTW4
#define IMPL_OP(OP)                                                                                                                                  \
  template <typename T> friend auto operator OP(mesh_point_t const &mp, T &&y) { return mp.value() OP std::forward<T>(y); }                          \
  template <typename T>                                                                                                                              \
  friend auto operator OP(T &&x, mesh_point_t const &mp)                                                                                             \
    requires(not std::is_same_v<std::decay_t<T>, mesh_point_t>)                                                                                      \
  {                                                                                                                                                  \
    return std::forward<T>(x) OP mp.value();                                                                                                         \
  }
      IMPL_OP(+)
      IMPL_OP(-)
      IMPL_OP(*)
      IMPL_OP(/)
#undef IMPL_OP
    };

    // -------------------- checks -------------------

    [[nodiscard]] bool is_index_valid(index_t index) const noexcept { return 0 <= index and index < size(); }

    [[nodiscard]] bool is_value_valid(value_t x) const noexcept { return _pts.front() <= x and x <= _pts.back(); }

    // -------------------- to_data_index ------------------
    [[nodiscard]] data_index_t to_data_index(index_t index) const noexcept {
      EXPECTS(is_index_valid(index));
      return index;
    }

    [[nodiscard]] index_t to_data_index(closest_mesh_point_t<value_t> const &cmp) const noexcept {
      EXPECTS(is_value_valid(cmp.value));
      return to_data_index(to_index(cmp));
    }

    // ------------------ to_index -------------------

    [[nodiscard]] index_t to_index(data_index_t data_index) const noexcept {
      EXPECTS(is_index_valid(data_index));
      return data_index;
    }

    [[nodiscard]] index_t to_index(closest_mesh_point_t<value_t> const &cmp) const noexcept {
      EXPECTS(is_value_valid(cmp.value));

      auto itr_r = std::lower_bound(begin(_pts), end(_pts), cmp.value);
      long i_r   = itr_r - _pts.begin();

      if (i_r == 0) { return 0; }
      if (i_r == size()) { return size() - 1; }

      long i_l = i_r - 1;
      if (fabs(cmp.value - _pts[i_l]) < fabs(cmp.value - _pts[i_r]))
        return i_l;
      else
        return i_r;
    }

    // ------------------ operator [] () -------------------

    [[nodiscard]] mesh_point_t operator[](long data_index) const noexcept { return (*this)(data_index); }

    [[nodiscard]] mesh_point_t operator[](closest_mesh_point_t<value_t> const &cmp) const noexcept { return (*this)[this->to_data_index(cmp)]; }

    [[nodiscard]] mesh_point_t operator()(index_t index) const noexcept {
      EXPECTS(is_index_valid(index));
      return {index, index, _mesh_hash, to_value(index)};
    }

    // -------------------- to_value ------------------

    [[nodiscard]] value_t to_value(index_t index) const noexcept {
      EXPECTS(is_index_valid(index));
      return _pts[index];
    }

    // -------------------- Accessors -------------------

    /// The Hash for the mesh configuration
    [[nodiscard]] uint64_t mesh_hash() const noexcept { return _mesh_hash; }

    /// The total number of points in the mesh
    [[nodiscard]] long size() const noexcept { return _pts.size(); }

    /// First index of the mesh
    static constexpr long first_index() { return 0; }

    /// Last index of the mesh
    [[nodiscard]] long last_index() const { return size() - 1; }

    /// The vector of points
    [[nodiscard]] std::vector<value_t> const &points() const noexcept { return _pts; }

    // -------------------------- Range & Iteration --------------------------

    private:
    [[nodiscard]] auto r_() const {
      return itertools::transform(range(size()), [this](long i) { return (*this)[i]; });
    }

    public:
    [[nodiscard]] auto begin() const { return r_().begin(); }
    [[nodiscard]] auto cbegin() const { return r_().cbegin(); }
    [[nodiscard]] auto end() const { return r_().end(); }
    [[nodiscard]] auto cend() const { return r_().cend(); }

    // -------------------- print  -------------------

    friend std::ostream &operator<<(std::ostream &sout, point_mesh const &m) { return sout << "Point mesh of size " << m.size(); }

    //  -------------------------- HDF5  --------------------------

    [[nodiscard]] static std::string hdf5_format() { return "PointMesh"; }

    friend void h5_write(h5::group fg, std::string const &subgroup_name, point_mesh const &m) {
      h5::group gr = fg.create_group(subgroup_name);
      write_hdf5_format(gr, m);
      h5::write(gr, "points", m.points());
    }

    /// Read from HDF5
    friend void h5_read(h5::group fg, std::string const &subgroup_name, point_mesh &m) {
      h5::group gr = fg.open_group(subgroup_name);
      assert_hdf5_format(gr, m);

      auto pts = h5::read<std::vector<value_t>>(gr, "points");
      m        = point_mesh(pts);
    }
  };

  // ------------------------- evaluation -----------------------------

  template <typename Value> auto evaluate(point_mesh<Value> const &m, auto const &f, Value x) {
    EXPECTS(m.is_value_valid(x));

    auto itr_r = std::lower_bound(begin(m.points()), end(m.points()), x);
    long i_r   = itr_r - m.points().begin();

    if (i_r == 0) { return f(0); }
    if (i_r == m.size()) { return f(m.size() - 1); }

    long i_l = i_r - 1;

    Value x_l = m.points()[i_l];
    Value x_r = m.points()[i_r];
    auto del    = x_r - x_l;

    // The interpolation weights
    double w_r = (x - x_l) / del;
    double w_l = (x_r - x) / del;

    ASSERT(x_l <= x && x <= x_r);
    ASSERT(w_l + w_r - 1 < 1e-15);

    return f(i_l) * w_l + f(i_r) * w_r;
  }

} // namespace triqs::mesh
