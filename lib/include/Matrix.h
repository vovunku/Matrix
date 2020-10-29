#ifndef MATRIX_MATRIX_H
#define MATRIX_MATRIX_H

#include <Rational.h>
#include <vector>

template<typename Field>
using vv_matrix = std::vector<std::vector<Field> >;

template<unsigned int M, unsigned int N, typename Field = Rational>
class Matrix {
 public:
  explicit Matrix(Field diag_num = 1);
  explicit Matrix(const vv_matrix<Field> &other, unsigned int left_up_x = 0, unsigned int left_up_y = 0);
  Matrix(std::initializer_list<std::initializer_list<Field> > brace);
  Matrix(const Matrix<M, N, Field> &other);

  std::vector<Field> &operator[](unsigned int i);
  const std::vector<Field> &operator[](unsigned int i) const;

  Matrix<M, N, Field> &operator+=(const Matrix<M, N, Field> &other);
  Matrix<M, N, Field> &operator-=(const Matrix<M, N, Field> &other);
  Matrix<M, N, Field> &operator*=(const Field &number);
  Matrix<M, N, Field> &operator*=(const Matrix<M, N, Field> &other);
  Matrix<M, N, Field> &operator=(const Matrix<M, N, Field> &other);

  template<unsigned int L>
  Matrix<M, L, Field> multiplicate(const Matrix<N, L, Field> &other) const;

  Field det() const;
  Matrix<N, M, Field> transposed() const;
  unsigned int rank() const;
  Matrix<M, N, Field> inverted() const;
  Matrix<M, N, Field> &invert();
  Field trace() const;

  void print_matrix() const;
  const vv_matrix<Field> &get_matrix() const;
  vv_matrix<Field> &get_matrix();
  std::vector<Field> getRow(unsigned int i) const;
  std::vector<Field> getColumn(unsigned int j) const;
 private:
  vv_matrix<Field> _matrix;

  static const unsigned int _simple_multiplication_limit = 15;

  template<unsigned int K, bool all_square, bool is_simple>
  struct multiplication_choose {
    static Matrix<M, K, Field> mlt(const Matrix<M, N, Field>& left, const Matrix<N, K, Field>& right);
  };

  template<unsigned int K, bool all_square>
  struct multiplication_choose<K, all_square, true> {
    static Matrix<M, K, Field> mlt(const Matrix<M, N, Field>& left, const Matrix<N, K, Field>& right){
      vv_matrix<Field> ans(M, std::vector<Field>(K));
      for (unsigned int i = 0; i < M; ++i) {
        for (unsigned int j = 0; j < K; ++j) {
          for (unsigned int t = 0; t < N; ++t) {
            ans[i][j] += left[i][t] * right[t][j];
          }
        }
      }
      return Matrix<M, K, Field>(ans);
    };
  };

  template<unsigned int K>
  struct multiplication_choose<K, false, false> {
    static Matrix<M, K, Field> mlt(const Matrix<M, N, Field>& left, const Matrix<N, K, Field>& right){
      const unsigned int new_size = upper_pow_2<compile_max<compile_max<N, M>::value, K>::value>::value;
      vv_matrix<Field> left_v(new_size, std::vector<Field>(new_size, 0));
      vv_matrix<Field> right_v(new_size, std::vector<Field>(new_size, 0));

      Matrix<M, N, Field>::_matrix_copy(left.get_matrix(), 0, 0, M, N, left_v);
      Matrix<M, N, Field>::_matrix_copy(right.get_matrix(), 0, 0, N, K, right_v);

      Matrix<new_size, new_size, Field> left_m(left_v);
      Matrix<new_size, new_size, Field> right_m(right_v);

      return Matrix<M, K, Field>(left_m.multiplicate(right_m).get_matrix());
    };
  };

  template <unsigned int K>
  struct multiplication_choose<K, true, false> {
    static Matrix<M, K, Field> mlt(const Matrix<M, N, Field>& left, const Matrix<N, K, Field>& right){
      Matrix<M / 2, M / 2, Field> A[2][2];
      Matrix<M / 2, M / 2, Field> B[2][2];

      Matrix<M, N, Field>::_matrix_split(left, A);
      Matrix<M, N, Field>::_matrix_split(right, B);

      Matrix<M / 2, M / 2, Field> P1 = (A[0][0] + A[1][1]) * (B[0][0] + B[1][1]);
      Matrix<M / 2, M / 2, Field> P2 = (A[1][0] + A[1][1]) * B[0][0];
      Matrix<M / 2, M / 2, Field> P3 = A[0][0] * (B[0][1] - B[1][1]);
      Matrix<M / 2, M / 2, Field> P4 = A[1][1] * (B[1][0] - B[0][0]);
      Matrix<M / 2, M / 2, Field> P5 = (A[0][0] + A[0][1]) * B[1][1];
      Matrix<M / 2, M / 2, Field> P6 = (A[1][0] - A[0][0]) * (B[0][0] + B[0][1]);
      Matrix<M / 2, M / 2, Field> P7 = (A[0][1] - A[1][1]) * (B[1][0] + B[1][1]);

      Matrix<M / 2, M / 2, Field> C[2][2];

      C[0][0] = P1 + P4 + P7 - P5;
      C[0][1] = P3 + P5;
      C[1][0] = P2 + P4;
      C[1][1] = P1 - P2 + P3 + P6;

      return Matrix<M, N, Field>::_matrix_collect(C);
    }
  };

  void _matrix_addition(vv_matrix<Field> &left, const vv_matrix<Field> &right, const Field &coeff) const;
  vv_matrix<Field> _matrix_summation(const vv_matrix<Field> &left, const vv_matrix<Field> &right) const;
  vv_matrix<Field> _matrix_subtraction(const vv_matrix<Field> &left, const vv_matrix<Field> &right) const;

  //Strassen realization

  static void _matrix_split(const Matrix<M, M, Field> &real, Matrix<M / 2, M / 2, Field> sliced[2][2]);
  static void _matrix_copy(const vv_matrix<Field> &from, unsigned int left_up_x, unsigned int left_up_y, unsigned int size_x,
                           unsigned int size_y, vv_matrix<Field> &to, unsigned int to_x = 0, unsigned int to_y = 0);
  static Matrix<M, M, Field> _matrix_collect(Matrix<M /  2, M / 2, Field> sliced[2][2]);
  unsigned int get_m_size(const vv_matrix<Field> &matrix) const;
  unsigned int get_n_size(const vv_matrix<Field> &matrix) const;

  //Gauss realization
  int _forward_gauss(vv_matrix<Field> &matrix, unsigned int work_size) const; // returns sign of matrix after changing
  void _backward_gauss(vv_matrix<Field> &matrix, unsigned int work_size) const;
  void _vector_addition(std::vector<Field> &left, const std::vector<Field> &right, const Field &coeff = 1) const;

  static unsigned int _nearest_pow_2(unsigned int x);

  static typename my_enable_if<is_prime<what_base<Field>::base>::value || !is_finite<Field>::value>::type checker() {};
};

template<unsigned int M, unsigned int N, typename Field>
Matrix<M, N, Field> operator+(const Matrix<M, N, Field> &left, const Matrix<M, N, Field> &right);

template<unsigned int M, unsigned int N, typename Field>
Matrix<M, N, Field> operator-(const Matrix<M, N, Field> &left, const Matrix<M, N, Field> &right);

template<unsigned int M, unsigned int N, unsigned int L, typename Field>
Matrix<M, L, Field> operator*(const Matrix<M, N, Field> &left, const Matrix<N, L, Field> &right);

template<unsigned int M, unsigned int N, typename Field>
Matrix<M, N, Field> operator*(const Field &left, const Matrix<M, N, Field> &right);

template<unsigned int M, unsigned int N, typename Field>
bool operator==(const Matrix<M, N, Field> &left, const Matrix<M, N, Field> &right);

template<unsigned int M, unsigned int N, typename Field>
bool operator!=(const Matrix<M, N, Field> &left, const Matrix<M, N, Field> &right);

template<unsigned int M, typename Field = Rational>
using SquareMatrix = Matrix<M, M, Field>;

#include "../src/Matrix.cpp"

#endif //MATRIX_MATRIX_H
