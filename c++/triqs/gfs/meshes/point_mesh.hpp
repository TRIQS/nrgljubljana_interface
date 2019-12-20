#include <triqs/gfs/meshes.hpp>
#include <triqs/gfs/meshes/mesh_tools.hpp>

#include <vector>
#include <string>
#include <algorithm>
#include <initializer_list>

namespace triqs::gfs {

  /**
   * A generic mesh of points on a sortable domain
   */
  template <typename Domain> class point_mesh {

    public:
    /// The domain where the values live
    using domain_t = Domain;

    /// Type of the points in the domain
    using domain_pt_t = typename domain_t::point_t;

    /// Discrete index that uniquely identifies points
    using index_t = long;

    /// Linear index is the position in the flattened mesh
    using linear_index_t = long;

    /// The type of the mesh points
    using mesh_point_t = mesh_point<point_mesh>;

    // -------------------- Constructors -------------------

    /// Construct from a vector of points
    point_mesh(std::vector<domain_pt_t> pts) : _pts(pts) {
      if (not std::is_sorted(pts.begin(), pts.end())) TRIQS_RUNTIME_ERROR << "Point mesh must be constructed with a sorted list of points";
    }

    /// Initializer list constructor
    point_mesh(std::initializer_list<domain_pt_t> l) : point_mesh(std::vector(l)) {}

    /// Default constructor
    point_mesh() = default;

    // -------------------- Accessors -------------------

    /// The corresponding domain
    [[nodiscard]] domain_t const &domain() const noexcept { return _dom; }

    /// The vector of points
    [[nodiscard]] std::vector<domain_pt_t> const &points() const noexcept { return _pts; }

    /// The number of points in the mesh
    [[nodiscard]] size_t size() const noexcept { return _pts.size(); }

    [[nodiscard]] utility::mini_vector<size_t, 1> size_of_components() const noexcept { return {size()}; }

    // -------------------- utility -------------------

    /// Comparison Operator
    [[nodiscard]] bool operator==(point_mesh const &m) const noexcept { return _dom == m._dom && _pts == m._pts; }
    [[nodiscard]] bool operator!=(point_mesh const &m) const noexcept { return !(*this == m); }

    /// Check if points are in the mesh
    [[nodiscard]] static constexpr bool is_within_boundary(all_t) noexcept { return true; }
    [[nodiscard]] bool is_within_boundary(domain_pt_t x) const noexcept { return _pts.front() <= x && x <= _pts.back(); }
    [[nodiscard]] bool is_within_boundary(index_t idx) const noexcept { return idx >= 0 && idx < _pts.size(); }

    /// Return the point in the domain for a given index
    [[nodiscard]] domain_pt_t index_to_point(index_t idx) const noexcept { return _pts[idx]; }

    /// Return the linear index for a given index
    [[nodiscard]] linear_index_t index_to_linear(index_t idx) const noexcept { return idx; }

    // -------------- Iterator Interface --------------------------

    using const_iterator = mesh_pt_generator<point_mesh>;
    mesh_point_t operator[](index_t i) const { return {*this, i}; }
    [[nodiscard]] const_iterator begin() const { return const_iterator(this); }
    [[nodiscard]] const_iterator end() const { return const_iterator(this, true); }
    [[nodiscard]] const_iterator cbegin() const { return const_iterator(this); }
    [[nodiscard]] const_iterator cend() const { return const_iterator(this, true); }

    // -------------- Evaluation of a function on the domain --------------------------

    [[nodiscard]] interpol_data_lin_t<index_t, 2> get_interpolation_data(domain_pt_t x) const noexcept {
      // indices to the left and right
      index_t i_r = std::distance(_pts.begin(), std::lower_bound(_pts.begin(), _pts.end(), x));
      index_t i_l = i_r - 1;

      // Points to the left and right
      domain_pt_t x_l = _pts[i_l];
      domain_pt_t x_r = _pts[i_r];
      auto del        = _pts[i_r] - _pts[i_l];

      // The interpolation weights
      auto w_r = (x - x_l) / del;
      auto w_l = (x_r - x) / del;

      ASSERT(x_l <= x && x <= x_r);
      ASSERT(w_l + w_r - 1 < 1e-15);

      return {{i_l, i_r}, {w_l, w_r}};
    }

    template <typename F>[[nodiscard]] auto evaluate(F const &f, domain_pt_t x) const noexcept {
      EXPECTS(is_within_boundary(x));
      if (x == _pts.front()) return f[0]; // special case
      auto id = get_interpolation_data(x);
      return id.w[0] * f[id.idx[0]] + id.w[1] * f[id.idx[1]];
    }

    // -------------- HDF5  --------------------------

    friend void h5_write_impl(h5::group fg, std::string const &subgroup_name, point_mesh const &m, const char *tag) {
      h5::group gr = fg.create_group(subgroup_name);
      gr.write_hdf5_scheme_as_string(tag);
      h5_write(gr, "domain", m._dom);
      h5_write(gr, "points", m._pts);
    }

    friend void h5_read_impl(h5::group fg, std::string const &subgroup_name, point_mesh &m, const char *tag_expected) {
      h5::group gr = fg.open_group(subgroup_name);
      gr.assert_hdf5_scheme_as_string(tag_expected, true);
      h5_read(gr, "domain", m._dom);
      h5_read(gr, "points", m._pts);
    }

    // -------------------- serialization -------------------

    friend class boost::serialization::access;
    template <class Archive> void serialize(Archive &ar, const unsigned int version) {
      ar &TRIQS_MAKE_NVP("domain", _dom);
      ar &TRIQS_MAKE_NVP("points", _pts);
    }

    // ------------------------------------------------

    private:
    std::vector<domain_pt_t> _pts;
    domain_t _dom;
  };

  /// Print the mesh
  template <typename Domain> std::ostream &operator<<(std::ostream &out, point_mesh<Domain> const &m) {
    return out << "Point mesh of size " << m.size();
  }

  // -------------- Mesh point class --------------------------

  template <typename Domain>
  class mesh_point<point_mesh<Domain>> : public utility::arithmetic_ops_by_cast<mesh_point<point_mesh<Domain>>, typename Domain::point_t> {

    public:
    /// The associated mesh
    using mesh_t = point_mesh<Domain>;

    /// The index type of the associated mesh
    using index_t = typename mesh_t::index_t;

    /// The linear index type of the associated mesh
    using linear_index_t = typename mesh_t::index_t;

    /// The value type of the mesh point
    using value_t = typename mesh_t::domain_pt_t;
    using cast_t  = value_t;

    private:
    /// The mesh that we belong to
    mesh_t const *m = nullptr;

    /// The index in the mesh
    index_t _index;

    // -------------------- Constructors -------------------

    public:
    mesh_point() = default;
    mesh_point(mesh_t const &mesh, index_t const &index) : m(&mesh), _index(index) {}
    mesh_point(mesh_t const &mesh) : mesh_point(mesh, 0) {}

    // -------------------- Utility -------------------

    operator index_t() const noexcept { return _index; }
    operator value_t() const noexcept { return m->index_to_point(_index); }

    [[nodiscard]] linear_index_t linear_index() const { return m->index_to_linear(_index); }
    [[nodiscard]] index_t index() const { return _index; }
    [[nodiscard]] mesh_t const &mesh() const { return *m; }

    void advance() noexcept { ++_index; }
  };

} // namespace triqs::gfs
