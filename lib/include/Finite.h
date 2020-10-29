#ifndef MATRIX_FINITE_H
#define MATRIX_FINITE_H

template<unsigned int Base>
class Finite {
 public:
  Finite() : _value(0) {}
  Finite(int number) : _value(((number % static_cast<int>(Base)) + static_cast<int>(Base)) % static_cast<int>(Base)) {}
  Finite(const Finite<Base> &other) : _value(other.get_value()) {}

  Finite<Base> &operator=(const Finite<Base> &other);
  Finite<Base> &operator+=(const Finite<Base> &other);
  Finite<Base> &operator-=(const Finite<Base> &other);
  Finite<Base> &operator*=(const Finite<Base> &other);
  Finite<Base> &operator/=(const Finite<Base> &other);

  Finite<Base> operator++(int);
  Finite<Base> &operator++();
  Finite<Base> operator-();

  template<unsigned int B>
  friend std::istream &operator>>(std::istream &in, Finite<B> &finite_num);

  template<unsigned int B>
  friend std::ostream &operator<<(std::ostream &out, const Finite<B> &finite_num);

  Finite<Base> inverse() const;

  unsigned int get_value() const;
  void set_value(unsigned int new_value);
  static const unsigned int base = Base;

 private:
  unsigned int _value;

  template<typename B>
  Finite<Base> _inverse() const;
};

template<>
class Finite<0> {
 private:
  Finite() {}
};

template<>
class Finite<1> {
 private:
  Finite() {}
};

template<unsigned int Base>
Finite<Base> operator+(const Finite<Base> &left, const Finite<Base> &right);

template<unsigned int Base>
Finite<Base> operator*(const Finite<Base> &left, const Finite<Base> &right);

template<unsigned int Base>
Finite<Base> operator-(const Finite<Base> &left, const Finite<Base> &right);

template<unsigned int Base>
Finite<Base> operator/(const Finite<Base> &left, const Finite<Base> &right);

template<unsigned int Base>
bool operator==(const Finite<Base> &left, const Finite<Base> &right);

template<unsigned int Base>
bool operator==(const Finite<Base> &left, int right);

template<unsigned int Base>
bool operator==(int left, const Finite<Base> &right);

template<unsigned int Base>
bool operator!=(const Finite<Base> &left, const Finite<Base> &right);

template<unsigned int Base>
bool operator!=(const Finite<Base> &left, int right);

template<unsigned int Base>
bool operator!=(int left, const Finite<Base> &right);

template<typename T>
struct what_base {
  static const unsigned int base = 0;
};

template<unsigned int T>
struct what_base<Finite<T> > {
  static const unsigned int base = T;
};

template<typename T>
struct is_finite {
  static const bool value = false;
};

template<unsigned int T>
struct is_finite<Finite<T> > {
  static const bool value = true;
};

#include <../src/Finite.cpp>

#endif //MATRIX_FINITE_H
