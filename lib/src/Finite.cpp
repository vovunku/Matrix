#include <exception>

template<unsigned int Base>
unsigned int Finite<Base>::get_value() const {
  return _value;
}

template<unsigned int Base>
void Finite<Base>::set_value(unsigned int new_value) {
  _value = new_value;
}

template<unsigned int Base>
Finite<Base> &Finite<Base>::operator=(const Finite<Base> &other) {
  _value = other.get_value();
  return *this;
}

template<unsigned int Base>
Finite<Base> &Finite<Base>::operator+=(const Finite<Base> &other) {
  _value += other.get_value();
  if (_value >= Base) {
    _value -= Base;
  }
  return *this;
}

template<unsigned int Base>
Finite<Base> &Finite<Base>::operator-=(const Finite<Base> &other) {
  _value += Base - other.get_value();
  if (_value >= Base) {
    _value -= Base;
  }
  return *this;
}

template<unsigned int Base>
Finite<Base> &Finite<Base>::operator*=(const Finite<Base> &other) {
  long long int result = _value;
  result *= other.get_value();
  result %= Base;
  _value = result;
  return *this;
}

template<unsigned int Base>
Finite<Base> &Finite<Base>::operator/=(const Finite<Base> &other) {
  (*this) *= other.inverse();
  return *this;
}

template<unsigned int Base>
bool operator==(const Finite<Base> &left, const Finite<Base> &right) {
  return left.get_value() == right.get_value();
}

template<unsigned int Base>
bool operator==(const Finite<Base> &left, int right) {
  return left == Finite<Base>(right);
}

template<unsigned int Base>
bool operator==(int left, const Finite<Base> &right) {
  return Finite<Base>(left) == right;
}

template<unsigned int Base>
bool operator!=(const Finite<Base> &left, const Finite<Base> &right) {
  return !(left == right);
}

template<unsigned int Base>
bool operator!=(const Finite<Base> &left, int right) {
  return !(left == Finite<Base>(right));
}

template<unsigned int Base>
bool operator!=(int left, const Finite<Base> &right) {
  return !(Finite<Base>(left) == right);
}

template<unsigned int Base>
Finite<Base> &Finite<Base>::operator++() {
  (*this) += 1;
  return *this;
}

template<unsigned int Base>
Finite<Base> Finite<Base>::operator++(int) {
  Finite<Base> ans = *this;
  ++(*this);
  return ans;
}

template<unsigned int Base>
Finite<Base> Finite<Base>::operator-() {
  Finite<Base> ans(Base - _value);
  return ans;
}

template<unsigned int B>
std::istream &operator>>(std::istream &in, Finite<B> &finite_num) {
  in >> finite_num._value;
  return in;
}

template<unsigned int B>
std::ostream &operator<<(std::ostream &out, const Finite<B> &finite_num) {
  out << finite_num._value;
  return out;
}

template<unsigned int Base>
Finite<Base> operator+(const Finite<Base> &left, const Finite<Base> &right) {
  Finite<Base> ans = left;
  ans += right;
  return ans;
}

template<unsigned int Base>
Finite<Base> operator-(const Finite<Base> &left, const Finite<Base> &right) {
  Finite<Base> ans = left;
  ans -= right;
  return ans;
}

template<unsigned int Base>
Finite<Base> operator*(const Finite<Base> &left, const Finite<Base> &right) {
  Finite<Base> ans = left;
  ans *= right;
  return ans;
}

template<unsigned int Base>
Finite<Base> operator/(const Finite<Base> &left, const Finite<Base> &right) {
  Finite<Base> ans = left;
  ans /= right;
  return ans;
}

template<unsigned int Base>
Finite<Base> p_adic_pow(Finite<Base> number, unsigned int power) {
  Finite<Base> result = 1;
  while (power > 0) {
    if (power & 1) {
      result *= number;
    }
    number *= number;
    power >>= 1;
  }
  return result;
}

template<unsigned int Base>
Finite<Base> Finite<Base>::inverse() const {
  return _inverse<typename my_enable_if<is_prime<Base>::value>::type>();
}

template<unsigned int Base>
template<typename B>
Finite<Base> Finite<Base>::_inverse() const {
  if (get_value() == 0) {
    throw std::invalid_argument("inverting zero");
  }
  return p_adic_pow(*this, Base - 2);
}