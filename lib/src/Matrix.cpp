#include <My_traits.h>
#include <Matrix.h>
#include <iostream>
#include <algorithm>
#include <tuple>

template<typename Field>
void vector_print(vv_matrix<Field> matrix) {
  unsigned int M = matrix.size();
  unsigned int N = matrix[0].size();
  for (unsigned int i = 0; i < M; ++i) {
    for (unsigned int j = 0; j < N; ++j) {
      std::cout << matrix[i][j] << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}

template<unsigned int M, unsigned int N, typename Field>
Matrix<M, N, Field>::Matrix(Field diag_num): _matrix(M, std::vector<Field>(N)) {
  for (unsigned int i = 0; i < std::min(M, N); ++i) {
    _matrix[i][i] = diag_num;
  }
}

template<unsigned int M, unsigned int N, typename Field>
Matrix<M, N, Field>::Matrix(const vv_matrix<Field> &other, unsigned int left_up_x, unsigned int left_up_y) {
  _matrix.resize(M, std::vector<Field>(N));
  _matrix_copy(other, left_up_x, left_up_y, M, N, _matrix);
}

template<unsigned int M, unsigned int N, typename Field>
Matrix<M, N, Field>::Matrix(std::initializer_list<std::initializer_list<Field> > brace) {
  _matrix.resize(M, std::vector<Field>(N));
  typename std::initializer_list<std::initializer_list<Field> >::iterator row = brace.begin();
  typename std::initializer_list<Field>::iterator pos = row->begin();
  for (unsigned int i = 0; i < M; ++i) {
    for (unsigned int j = 0; j < N; ++j) {
      _matrix[i][j] = *pos;
      ++pos;
    }
    ++row;
    pos = row->begin();
  }
}

template<unsigned int M, unsigned int N, typename Field>
Matrix<M, N, Field>::Matrix(const Matrix<M, N, Field> &other) {
  _matrix.resize(M, std::vector<Field>(N));
  for (unsigned int i = 0; i < M; ++i) {
    for (unsigned int j = 0; j < N; ++j) {
      _matrix[i][j] = other[i][j];
    }
  }
}

template<unsigned int M, unsigned int N, typename Field>
std::vector<Field> &Matrix<M, N, Field>::operator[](unsigned int i) {
  return _matrix[i];
}

template<unsigned int M, unsigned int N, typename Field>
const std::vector<Field> &Matrix<M, N, Field>::operator[](unsigned int i) const {
  return _matrix[i];
}

template<unsigned int M, unsigned int N, typename Field>
void Matrix<M, N, Field>::print_matrix() const {
  for (unsigned int i = 0; i < M; ++i) {
    for (unsigned int j = 0; j < N; ++j) {
      std::cout << _matrix[i][j] << " ";
    }
    std::cout << "\n";
  }
}

template<unsigned int M, unsigned int N, typename Field>
const vv_matrix<Field> &Matrix<M, N, Field>::get_matrix() const {
  return _matrix;
}

template<unsigned int M, unsigned int N, typename Field>
std::vector<Field> Matrix<M, N, Field>::getRow(unsigned int i) const {
  return _matrix[i];
}

template<unsigned int M, unsigned int N, typename Field>
std::vector<Field> Matrix<M, N, Field>::getColumn(unsigned int j) const {
  std::vector<Field> ans;
  for (unsigned int i = 0; i < M; ++i) {
    ans.push_back(_matrix[i][j]);
  }
  return ans;
}

template<unsigned int M, unsigned int N, typename Field>
Matrix<M, N, Field> &Matrix<M, N, Field>::operator+=(const Matrix<M, N, Field> &other) {
  this->_matrix_addition(_matrix, other.get_matrix(), 1);
  return *this;
}

template<unsigned int M, unsigned int N, typename Field>
Matrix<M, N, Field> &Matrix<M, N, Field>::operator-=(const Matrix<M, N, Field> &other) {
  this->_matrix_addition(_matrix, other.get_matrix(), -1);
  return *this;
}

template<unsigned int M, unsigned int N, typename Field>
Matrix<M, N, Field> &Matrix<M, N, Field>::operator*=(const Field &number) {
  for (unsigned int i = 0; i < M; ++i) {
    for (unsigned int j = 0; j < N; ++j) {
      _matrix[i][j] *= number;
    }
  }
  return *this;
}

template<unsigned int M, unsigned int N, typename Field>
Matrix<M, N, Field> &Matrix<M, N, Field>::operator*=(const Matrix<M, N, Field> &other) {
  typename my_enable_if<M == N, bool>::type check;
  std::ignore = check;
  return (*this) = (*this).multiplicate(other);
}

template<unsigned int M, unsigned int N, typename Field>
Matrix<M, N, Field> &Matrix<M, N, Field>::operator=(const Matrix<M, N, Field> &other) {
  _matrix_copy(other.get_matrix(), 0, 0, M, N, _matrix);
  return (*this);
}

template<unsigned int M, unsigned int N, typename Field>
Matrix<M, N, Field> operator+(const Matrix<M, N, Field> &left, const Matrix<M, N, Field> &right) {
  Matrix<M, N, Field> ans = left;
  ans += right;
  return ans;
}

template<unsigned int M, unsigned int N, typename Field>
Matrix<M, N, Field> operator-(const Matrix<M, N, Field> &left, const Matrix<M, N, Field> &right) {
  Matrix<M, N, Field> ans = left;
  ans -= right;
  return ans;
}

template<unsigned int M, unsigned int N, unsigned int L, typename Field>
Matrix<M, L, Field> operator*(const Matrix<M, N, Field> &left, const Matrix<N, L, Field> &right) {
  return left.multiplicate(right);
}

template<unsigned int M, unsigned int N, typename Field>
Matrix<M, N, Field> operator*(const Field &left, const Matrix<M, N, Field> &right) {
  Matrix<M, N, Field> ans = right;
  ans *= left;
  return ans;
}

template<unsigned int M, unsigned int N, typename Field>
bool operator==(const Matrix<M, N, Field> &left, const Matrix<M, N, Field> &right) {
  for (unsigned int i = 0; i < M; ++i) {
    for (unsigned int j = 0; j < N; ++j) {
      if (left[i][j] != right[i][j])
        return false;
    }
  }
  return true;
}

template<unsigned int M, unsigned int N, typename Field>
bool operator!=(const Matrix<M, N, Field> &left, const Matrix<M, N, Field> &right) {
  return !(left == right);
}

template<unsigned int M, unsigned int N, typename Field>
template<unsigned int K>
Matrix<M, K, Field> Matrix<M, N, Field>::multiplicate(const Matrix<N, K, Field> &other) const {
  return multiplication_choose<K, (M == N) && (N == K), compile_max<compile_max<M, N>::value, K>::value < _simple_multiplication_limit>::mlt(*this, other);
}

template<unsigned int M, unsigned int N, typename Field>
Field Matrix<M, N, Field>::det() const {
  typename my_enable_if<N == M, vv_matrix<Field> >::type copy = _matrix;
  int sign = _forward_gauss(copy, M);
  Field ans = 1;
  for (unsigned int i = 0; i < M; ++i) {
    ans *= copy[i][i];
  }
  ans *= sign;
  return ans;
}

template<unsigned int M, unsigned int N, typename Field>
Matrix<N, M, Field> Matrix<M, N, Field>::transposed() const {
  Matrix<N, M, Field> ans;
  for (unsigned int i = 0; i < N; ++i) {
    for (unsigned int j = 0; j < M; ++j) {
      ans[i][j] = (*this)[j][i];
    }
  }
  return ans;
}

template<unsigned int M, unsigned int N, typename Field>
unsigned int Matrix<M, N, Field>::rank() const {
  unsigned int work_size = std::max(M, N);
  vv_matrix<Field> copy(work_size, std::vector<Field>(work_size, 0));
  _matrix_copy(_matrix, 0, 0, M, N, copy);
  this->_forward_gauss(copy, work_size);
  unsigned int ans = 0;
  bool is_zero = true;
  for (unsigned int i = 0; i < work_size; ++i) {
    for (unsigned int j = 0; j < work_size; ++j) {
      if (copy[i][j] != 0) {
        is_zero = false;
        break;
      }
    }
    if (is_zero == true) {
      break;
    } else {
      ++ans;
      is_zero = true;
    }
  }
  return ans;
}

template<unsigned int M, unsigned int N, typename Field>
Matrix<M, N, Field> &Matrix<M, N, Field>::invert() {
  typename my_enable_if<N == M, vv_matrix<Field> >::type copy(M, std::vector<Field>(2 * M, 0));
  _matrix_copy(_matrix, 0, 0, M, M, copy);
  for (unsigned int i = 0; i < M; ++i) {
    copy[i][M + i] = 1;
  }
  //vector_print(copy);
  _forward_gauss(copy, M);
  //vector_print(copy);
  for (unsigned int i = 0; i < M; ++i) {
    if (copy[i][i] == 0)
      throw std::logic_error("det = 0");
  }
  _backward_gauss(copy, M);
  //vector_print(copy);
  _matrix_copy(copy, 0, M, M, M, _matrix);
  return *this;
}

template<unsigned int M, unsigned int N, typename Field>
Matrix<M, N, Field> Matrix<M, N, Field>::inverted() const {
  Matrix<M, N, typename my_enable_if<M == N, Field>::type> ans = *this;
  return ans.invert();
}

template<unsigned int M, unsigned int N, typename Field>
Field Matrix<M, N, Field>::trace() const {
  typename my_enable_if<M == N, Field>::type ans = 0;
  for (unsigned int i = 0; i < M; ++i) {
    ans += _matrix[i][i];
  }
  return ans;
}

// private

template<unsigned int M, unsigned int N, typename Field>
void Matrix<M, N, Field>::_matrix_addition(vv_matrix<Field> &left,
                                           const vv_matrix<Field> &right,
                                           const Field &coeff) const {
  unsigned int m_size = get_m_size(left);
  unsigned int n_size = get_n_size(right);
  for (unsigned int i = 0; i < m_size; ++i) {
    for (unsigned int j = 0; j < n_size; ++j) {
      left[i][j] += coeff * right[i][j];
    }
  }
}

template<unsigned int M, unsigned int N, typename Field>
vv_matrix<Field> Matrix<M, N, Field>::_matrix_summation(const vv_matrix<Field> &left,
                                                        const vv_matrix<Field> &right) const {
  vv_matrix<Field> ans = left;
  _matrix_addition(ans, right, 1);
  return ans;
}

template<unsigned int M, unsigned int N, typename Field>
vv_matrix<Field> Matrix<M, N, Field>::_matrix_subtraction(const vv_matrix<Field> &left,
                                                          const vv_matrix<Field> &right) const {
  vv_matrix<Field> ans = left;
  _matrix_addition(ans, right, -1);
  return ans;
}

//Strassen realization
template<unsigned int M, unsigned int N, typename Field>
unsigned int Matrix<M, N, Field>::get_m_size(const vv_matrix<Field> &matrix) const {
  if (matrix.size() == 0)
    throw std::logic_error("There is a mistake somewhere(get_m_size)");
  return matrix.size();
}

template<unsigned int M, unsigned int N, typename Field>
unsigned int Matrix<M, N, Field>::get_n_size(const vv_matrix<Field> &matrix) const {
  if (matrix.size() == 0)
    throw std::logic_error("There is a mistake somewhere(get_n_size1)");
  if (matrix[0].size() == 0)
    throw std::logic_error("There is a mistake somewhere(get_n_size1)");
  return matrix[0].size();
}

template<unsigned int M, unsigned int N, typename Field>
unsigned int Matrix<M, N, Field>::_nearest_pow_2(unsigned int x) {
  unsigned int ans = 1;
  while (x > ans) {
    ans <<= 1;
  }
  return ans;
}

//Gauss realization

template<unsigned int M, unsigned int N, typename Field>
int Matrix<M, N, Field>::_forward_gauss(vv_matrix<Field> &matrix, unsigned int work_size) const {
  int sign = 1;
  for (unsigned int j = 0; j < work_size; ++j) {

    //vector_print(matrix);
    unsigned int i = j;
    bool found_not_zero = false;
    for (unsigned int t = j; t < work_size; ++t) {
      if (matrix[t][j] != 0) {
        i = t;
        found_not_zero = true;
        break;
      }
    }
    if (!found_not_zero)
      continue;

    if (i != j) {
      sign *= -1;
    }
    std::swap(matrix[j], matrix[i]);
    for (unsigned int t = j + 1; t < work_size; ++t) {
      _vector_addition(matrix[t], matrix[j], -matrix[t][j] / matrix[j][j]);
    }
  }
  return sign;
  //vector_print(matrix);
}

template<unsigned int M, unsigned int N, typename Field>
void Matrix<M, N, Field>::_backward_gauss(vv_matrix<Field> &matrix, unsigned int work_size) const {
  for (unsigned int j = 0; j < work_size; ++j) {
    if (matrix[j][j] == 0)
      continue;

    Field tmp = matrix[j][j];
    for (unsigned int t = j; t < matrix[j].size(); ++t) {
      matrix[j][t] /= tmp;
    }

    for (unsigned int t = 0; t < j; ++t) {
      _vector_addition(matrix[t], matrix[j], -matrix[t][j]);
    }
  }
  //vector_print(matrix);
}

template<unsigned int M, unsigned int N, typename Field>
void Matrix<M, N, Field>::_vector_addition(std::vector<Field> &left,
                                           const std::vector<Field> &right,
                                           const Field &coeff) const {
  unsigned int size = left.size();
  for (unsigned int i = 0; i < size; ++i) {
    left[i] += right[i] * coeff;
  }
}
template<unsigned int M, unsigned int N, typename Field>
vv_matrix<Field> &Matrix<M, N, Field>::get_matrix() {
  return _matrix;
}
template<unsigned int M, unsigned int N, typename Field>
Matrix<M, M, Field> Matrix<M, N, Field>::_matrix_collect(Matrix<M / 2, M / 2, Field> sliced[2][2]) {
  vv_matrix<Field> ans(M, std::vector<Field>(M));
  _matrix_copy(sliced[0][0].get_matrix(), 0, 0, M / 2, M / 2, ans);
  _matrix_copy(sliced[0][1].get_matrix(), 0, 0, M / 2, M / 2, ans, 0, M / 2);
  _matrix_copy(sliced[1][0].get_matrix(), 0, 0, M / 2, M / 2, ans, M / 2, 0);
  _matrix_copy(sliced[1][1].get_matrix(), 0, 0, M / 2, M / 2, ans, M / 2, M / 2);
  return Matrix<M, M, Field>(ans);
}
template<unsigned int M, unsigned int N, typename Field>
void Matrix<M, N, Field>::_matrix_split(const Matrix<M, M, Field> &real, Matrix<M / 2, M / 2, Field> sliced[2][2]) {
  _matrix_copy(real.get_matrix(), 0, 0, M / 2, M / 2, sliced[0][0].get_matrix());
  _matrix_copy(real.get_matrix(), 0, M / 2, M / 2, M / 2, sliced[0][1].get_matrix());
  _matrix_copy(real.get_matrix(), M / 2, 0, M / 2, M / 2, sliced[1][0].get_matrix());
  _matrix_copy(real.get_matrix(), M / 2, M / 2, M / 2, M / 2, sliced[1][1].get_matrix());
}

template<unsigned int M, unsigned int N, typename Field>
void Matrix<M, N, Field>::_matrix_copy(const vv_matrix<Field> &from,
                                       unsigned int left_up_x,
                                       unsigned int left_up_y,
                                       unsigned int size_x,
                                       unsigned int size_y,
                                       vv_matrix<Field> &to,
                                       unsigned int to_x,
                                       unsigned int to_y) {
  for (unsigned int i = 0; i < size_x; ++i) {
    for (unsigned int j = 0; j < size_y; ++j) {
      to[to_x + i][to_y + j] = from[left_up_x + i][left_up_y + j];
    }
  }
}