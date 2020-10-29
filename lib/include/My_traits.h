#ifndef MATRIX_MY_TRAITS_H
#define MATRIX_MY_TRAITS_H

//for compilation errors

template<bool B, typename T = void>
struct my_enable_if {};

template<typename T>
struct my_enable_if<true, T> { typedef T type; };

template<bool B>
struct my_bool_constant {
  static const bool value = false;
};

template<>
struct my_bool_constant<true> {
  static const bool value = true;
};

//prime checker

template<unsigned int N, unsigned int S>
struct approx_sqrt {
  static const bool too_big = (S * S >= N);
  static const unsigned int value = (too_big && approx_sqrt<N, S / 2>::too_big) ? approx_sqrt<N, S / 2>::value : S;
};

template<unsigned int N>
struct approx_sqrt<N, 1> {
  static const bool too_big = false;
  static const unsigned int value = 1;
};

template<unsigned int N, unsigned int D>
struct is_divided {
  static const bool value = (N % D == 0) || is_divided<N, D - 1>::value;
};

template<unsigned int N>
struct is_divided<N, 1> {
  static const bool value = false;
};

template<unsigned int N>
struct is_prime : public my_bool_constant<!is_divided<N, approx_sqrt<N, N - 1>::value>::value> {};

template<>
struct is_prime<0> : public my_bool_constant<false> {};

template<>
struct is_prime<1> : public my_bool_constant<false> {};

template<unsigned int N>
struct lower_pow_2 {
  static const unsigned int value = lower_pow_2<N / 2>::value * 2;
};

template<>
struct lower_pow_2<1> {
  static const unsigned int value = 1;
};

template <unsigned int N>
struct upper_pow_2 { // or equal!
  static const unsigned int value = lower_pow_2<N>::value != N ? lower_pow_2<N>::value * 2 : N;
};

template <unsigned int N, unsigned int M>
struct compile_max {
  static const unsigned int value = (N > M) ? N : M;
};

#endif //MATRIX_MY_TRAITS_H
