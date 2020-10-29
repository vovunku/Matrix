#ifndef rational
#define rational

#include<iostream>
#include<vector>
#include<string>

//BigInteger
//==================================================================================================
enum {
  POSITIVE,
  NEGATIVE
};

class BigInteger {
 public:
  BigInteger();
  BigInteger(int number);
  BigInteger(const BigInteger &other);

  BigInteger &operator=(const BigInteger &other);
  BigInteger &operator*=(const BigInteger &other);
  BigInteger &operator+=(const BigInteger &other);
  BigInteger &operator-=(const BigInteger &other);
  BigInteger &operator/=(const BigInteger &other);
  BigInteger &operator%=(const BigInteger &other);
  BigInteger operator-() const;

  BigInteger &operator++();
  BigInteger operator++(int);
  BigInteger &operator--();
  BigInteger operator--(int);

  friend std::istream &operator>>(std::istream &in, BigInteger &bigNum);
  friend std::ostream &operator<<(std::ostream &out, const BigInteger &bigNum);

  std::string toString() const;

  unsigned int operator[](size_t i) const;

  explicit operator bool() const;

  unsigned int size() const;
  unsigned int whatSign() const;
  void swapBigInt(BigInteger &other);
  void changeSign();
  int order();
 private:
  std::vector<unsigned int> storage;
  unsigned int sign;
  const static int limit = 9;
  const static unsigned int base = 1000 * 1000 * 1000;

  void addNumberToDigit(unsigned long long number, unsigned int digit); //POSITIVE + POSITIVE
  void substractNumberToDigit(unsigned long long number, unsigned int digit); //POSITIVE - POSITIVE(SMALLER)
  void substractsmallerBigInt(const BigInteger &other);
  void addBigInt(const BigInteger &other);
  void deleteLeadingZeros();
  void divideAbsByTwo();
  void divideAbsByTen();
};

bool operator==(const BigInteger &left, const BigInteger &right);
bool operator!=(const BigInteger &left, const BigInteger &right);
bool operator<(const BigInteger &left, const BigInteger &right);
bool operator>(const BigInteger &left, const BigInteger &right);
bool operator<=(const BigInteger &left, const BigInteger &right);
bool operator>=(const BigInteger &left, const BigInteger &right);

BigInteger operator+(const BigInteger &left, const BigInteger &right);
BigInteger operator-(const BigInteger &left, const BigInteger &right);
BigInteger operator*(const BigInteger &left, const BigInteger &right);
BigInteger operator/(const BigInteger &left, const BigInteger &right);
BigInteger operator%(const BigInteger &left, const BigInteger &right);

//*************************************************************************//
//myfunctions

unsigned int BigInteger::size() const {
  return storage.size();
}

unsigned int BigInteger::whatSign() const {
  return sign;
}

unsigned int BigInteger::operator[](size_t i) const {
  return storage[i];
}

void BigInteger::addNumberToDigit(unsigned long long number, unsigned int digit) { //digit count from 0
  if (number >= base) {
    addNumberToDigit(number % base, digit);
    addNumberToDigit(number / base, digit + 1);
  } else if (number > 0 && number < base) {
    for (size_t i = size(); i < digit + 1; ++i) {
      storage.push_back(0);
    }
    number += (*this)[digit];
    storage[digit] = number % base;
    addNumberToDigit(number / base, digit + 1);
  }
  return;
}

void BigInteger::substractNumberToDigit(unsigned long long number, unsigned int digit) {
  if (number >= base) {
    substractNumberToDigit(number / base, digit + 1);
    substractNumberToDigit(number % base, digit);
  } else if (number < base && number > 0) {
    if (number < (*this)[digit]) {
      storage[digit] -= number;
    } else if (number > (*this)[digit]) {
      substractNumberToDigit(1, digit + 1);
      storage[digit] += base - number;
    } else {
      //std::cout<<"1" << " " << (*this).size() << " " << digit << "\n";
      if (size() == digit + 1 && size() > 1) {
        storage.pop_back();
        while (size() > 1 && (*this)[size() - 1] == 0) {
          storage.pop_back();
          deleteLeadingZeros();
        }
      } else {
        storage[digit] = 0;
      }
    }
  }
  return;
}

void BigInteger::substractsmallerBigInt(const BigInteger &other) {
  for (size_t i = 0; i < other.size(); ++i) {
    (*this).substractNumberToDigit(other[i], i);
  }
}

void BigInteger::addBigInt(const BigInteger &other) {
  for (size_t i = 0; i < other.size(); ++i) {
    (*this).addNumberToDigit(other[i], i);
  }
}

bool absBigIntMoreThan(const BigInteger &left, const BigInteger &right) { //compare >
  unsigned int lsize = left.size();
  unsigned int rsize = right.size();
  if (lsize > rsize) {
    return true;
  } else if (lsize < rsize) {
    return false;
  } else {
    for (int i = lsize - 1; i >= 0; --i) {
      if (left[i] > right[i]) {
        return true;
      } else if (left[i] < right[i]) {
        return false;
      }
    }
    return false;
  }
}

int tenPowLimit(unsigned int power) {
  int ans = 1;
  for (unsigned int i = 0; i < power; i++) {
    ans *= 10;
  }
  return ans;
}

void reverseString(std::string &s) {
  char tmp;
  unsigned int size = s.size();
  for (unsigned int t = 0; t < (size / 2); ++t) {
    tmp = s[t];
    s[t] = s[size - t - 1];
    s[size - t - 1] = tmp;
  }
}

void BigInteger::changeSign() {
  if (*this != BigInteger(0)) {
    if (sign == POSITIVE) {
      sign = NEGATIVE;
    } else {
      sign = POSITIVE;
    }
  }
}

void BigInteger::divideAbsByTwo() {
  unsigned int tsize = size();
  storage[0] /= 2;
  for (unsigned int i = 1; i < tsize; ++i) {
    addNumberToDigit(((*this)[i] % 2) * base / 2, i - 1);
    storage[i] /= 2;
  }
  deleteLeadingZeros();
  if (size() == 1 && (*this)[0] == 0) {
    sign = POSITIVE;
  }
}

void BigInteger::divideAbsByTen() {
  unsigned int tsize = size();
  storage[0] /= 10;
  for (unsigned int i = 1; i < tsize; ++i) {
    addNumberToDigit(((*this)[i] % 10) * (base / 10), i - 1);
    storage[i] /= 10;
  }
  deleteLeadingZeros();
  if (size() == 1 && (*this)[0] == 0) {
    sign = POSITIVE;
  }
}

void BigInteger::deleteLeadingZeros() {
  while (size() > 1 && (*this)[size() - 1] == 0) {
    storage.pop_back();
  }
}

void BigInteger::swapBigInt(BigInteger &other) {
  this->storage.swap(other.storage);
  std::swap(this->sign, other.sign);
}

//*************************************************************************//

BigInteger::BigInteger() {
  sign = POSITIVE;
  storage.push_back(0);
}

BigInteger::BigInteger(int number) {
  sign = number >= 0 ? POSITIVE : NEGATIVE;
  unsigned int absoluteNumber = abs(number);
  if (absoluteNumber < base) {
    storage.push_back(absoluteNumber);
  } else {
    storage.push_back(absoluteNumber % base);
    storage.push_back(absoluteNumber / base);
  }
}

BigInteger::BigInteger(const BigInteger &other) {
  this->sign = other.sign;
  this->storage = other.storage;
}

BigInteger &BigInteger::operator+=(const BigInteger &other) {
  if (this->sign == other.sign) {
    (*this).addBigInt(other);
  } else {
    if (absBigIntMoreThan((*this), other)) {
      (*this).substractsmallerBigInt(other);
    } else if ((*this) == -other) {
      //std::cout << "whynothere?";
      this->storage.clear();
      this->storage.push_back(0);
      this->sign = POSITIVE;
    } else {
      BigInteger answer = other;
      answer.substractsmallerBigInt(*this);
      (*this) = answer;
    }
  }
  return (*this);
}

BigInteger &BigInteger::operator-=(const BigInteger &other) {
  return (*this) += -other;
}

BigInteger &BigInteger::operator*=(const BigInteger &other) {
  BigInteger ans = 0;
  unsigned int lsize = (*this).size();
  unsigned int rsize = other.size();
  unsigned long long tmp;
  for (unsigned int i = 0; i < lsize; ++i) {
    for (unsigned int j = 0; j < rsize; ++j) {
      tmp = (*this)[i];
      tmp *= other[j];
      ans.addNumberToDigit(tmp, i + j);
    }
  }
  if (this->whatSign() != other.whatSign()) {
    ans.changeSign();
  }
  (*this).swapBigInt(ans);
  return *this;
}

/*BigInteger& BigInteger::operator /= (const BigInteger& other){
    if (*this == BigInteger(0)){
        return *this;
    }
    if (absBigIntMoreThan(other, *this)){
        return (*this) = BigInteger(0);
    }
    BigInteger left = 0;
    BigInteger right = *this;
    BigInteger median;
    if (this->sign == NEGATIVE){
        right.changeSign();
    }
    right += 1;
    while (right - left > 1){
        //std::cout << left << " " << right;
        median = right + left;
        median.divideAbsByTwo();
        //std::cout << " " <<  median << "\n";
        if (absBigIntMoreThan(median * other, (*this))){
            right = median;
        }else{
            left = median;
        }
    }
    if (this->sign != other.sign){
        left.changeSign();
    }
    (*this).swapBigInt(left);
    return *this;
}*/

BigInteger &BigInteger::operator/=(const BigInteger &other) {
  if (*this == BigInteger(0)) {
    return *this;
  }
  if (absBigIntMoreThan(other, *this)) {
    return *this = BigInteger(0);
  }
  BigInteger tmp = *this;
  unsigned int difference = this->size() - other.size() + 1;
  BigInteger big_ten(0);
  BigInteger work_copy = other;
  BigInteger ans(0);
  big_ten.addNumberToDigit(1, difference);
  work_copy *= big_ten;
  for (unsigned int i = 0; i < difference; ++i) {
    for (unsigned int j = 0; j < limit; ++j) {
      //std::cout << *this << " " << work_copy << " " << ans << '\n';
      while (!absBigIntMoreThan(work_copy, *this)) {
        substractsmallerBigInt(work_copy);
        ans += big_ten;
      }
      work_copy.divideAbsByTen();
      big_ten.divideAbsByTen();
      //int a;
      //std::cin >> a;
    }
  }
  while (!absBigIntMoreThan(work_copy, *this)) {
    substractsmallerBigInt(work_copy);
    ans += big_ten;
  }
  //std::cout << *this << " " << work_copy << " " << ans << "finish" <<  '\n';
  if (this->sign != other.whatSign()) {
    ans.changeSign();
  }
  //std::cout << '\\' << tmp << " " << other << " " << ans << "!!!!!!" << "\n";
  swapBigInt(ans);
  return *this;
}

BigInteger &BigInteger::operator%=(const BigInteger &other) {
  BigInteger tmp = (*this) / other;
  //std::cout << (*this) << " " << other << " " << tmp << "!!!!!!" << "\n";
  return *this = (*this) - other * tmp;
}

BigInteger operator+(const BigInteger &left, const BigInteger &right) {
  BigInteger ans = left;
  ans += right;
  return ans;
}

BigInteger operator-(const BigInteger &left, const BigInteger &right) {
  BigInteger ans = left;
  ans -= right;
  return ans;
}

BigInteger operator*(const BigInteger &left, const BigInteger &right) {
  BigInteger ans = left;
  return ans *= right;
}

BigInteger operator/(const BigInteger &left, const BigInteger &right) {
  BigInteger ans = left;
  ans /= right;
  return ans;
}

BigInteger operator%(const BigInteger &left, const BigInteger &right) {
  BigInteger ans = left;
  ans %= right;
  return ans;
}

BigInteger BigInteger::operator-() const {
  BigInteger answer = (*this);
  answer.changeSign();
  return answer;
}

BigInteger &BigInteger::operator++() {
  return (*this) += 1;
}

BigInteger BigInteger::operator++(int) {
  BigInteger ans = *this;
  ++(*this);
  return ans;
}

BigInteger &BigInteger::operator--() {
  return (*this) += -1;
}

BigInteger BigInteger::operator--(int) {
  BigInteger ans = *this;
  --(*this);
  return ans;
}

BigInteger &BigInteger::operator=(const BigInteger &other) {
  if (this == &other) {
    return *this;
  }
  this->sign = other.sign;
  this->storage = other.storage;
  return *this;
}

bool operator==(const BigInteger &left, const BigInteger &right) {
  if (left.size() != right.size()) {
    return false;
  }
  if (left.whatSign() != right.whatSign()) {
    return false;
  }
  for (size_t i = 0; i < left.size(); ++i) {
    if (left[i] != right[i]) {
      return false;
    }
  }
  return true;
}

bool operator!=(const BigInteger &left, const BigInteger &right) {
  return !(left == right);
}

bool operator>(const BigInteger &left, const BigInteger &right) {
  if (left.whatSign() == right.whatSign()) {
    if (left.whatSign() == POSITIVE) {
      return absBigIntMoreThan(left, right);
    } else {
      return absBigIntMoreThan(right, left);
    }
  } else {
    if (left.whatSign() == POSITIVE) {
      return true;
    } else {
      return false;
    }
  }
}

bool operator<(const BigInteger &left, const BigInteger &right) {
  return right > left;
}

bool operator>=(const BigInteger &left, const BigInteger &right) {
  return !(left < right);
}

bool operator<=(const BigInteger &left, const BigInteger &right) {
  return !(left > right);
}

std::istream &operator>>(std::istream &in, BigInteger &bigNumber) {
  std::string input;
  in >> input;
  int it = 0;
  int size = input.size();
  reverseString(input);
  if (input[size - 1] == '-') {
    bigNumber.sign = NEGATIVE;
    --size;
  } else {
    bigNumber.sign = POSITIVE;
  }
  bigNumber.storage.clear();
  bigNumber.storage.push_back(0);
  while (it < size) {
    bigNumber.addNumberToDigit(static_cast<int>(input[it] - '0') * tenPowLimit(it % bigNumber.limit),
                               static_cast<int>(it / bigNumber.limit));
    ++it;
  }
  if (bigNumber.sign == NEGATIVE && bigNumber.size() == 1 && bigNumber.storage[0] == 0) {
    bigNumber.sign = POSITIVE;
  }
  return in;
}

std::ostream &operator<<(std::ostream &out, const BigInteger &bigNumber) {
  std::string ans = bigNumber.toString();
  unsigned int size = ans.size();
  for (unsigned int i = 0; i < size; ++i) {
    out << ans[i];
  }
  return out;
}

std::string BigInteger::toString() const {
  int size = this->size();
  std::string ans = "";
  if (size == 0) {
    return ans;
  }
  if (size == 1 && (*this)[size - 1] == 0) {
    ans.push_back('0');
    return ans;
  }
  --size;
  unsigned int tmp;
  for (int i = 0; i <= size - 1; ++i) {
    tmp = (*this)[i];
    for (int j = 0; j < limit; ++j) {
      ans.push_back(tmp % 10 + '0');
      tmp /= 10;
    }
  }
  tmp = (*this)[size];
  while (tmp > 0) {
    ans.push_back(tmp % 10 + '0');
    tmp /= 10;
  }
  if (sign == NEGATIVE) {
    ans.push_back('-');
  }
  reverseString(ans);
  return ans;
}

BigInteger::operator bool() const {
  return (*this) != BigInteger(0);
}
//==================================================================================================

class Rational {
 public:
  Rational();
  Rational(int number);
  Rational(const BigInteger &bigNumber);
  Rational(const Rational &rationalNum);

  Rational &operator+=(const Rational &other);
  Rational &operator-=(const Rational &other);
  Rational &operator*=(const Rational &other);
  Rational &operator/=(const Rational &other);
  Rational &operator=(const Rational &other);
  Rational operator-() const;

  friend std::ostream &operator<<(std::ostream &out, const Rational &number);
  friend std::istream &operator>>(std::istream &in, Rational &number);

  std::string toString() const;

  std::string asDecimal(size_t precision = 0) const;

  explicit operator double();

  void swapRational(Rational &other);
  BigInteger getNumerator() const;
  BigInteger getDenominator() const;
 private:
  BigInteger numerator; //integer number
  BigInteger denominator; //natural number
  void normalize();
};
//*********************************************** operators
Rational operator+(const Rational &left, const Rational &right);
Rational operator-(const Rational &left, const Rational &right);
Rational operator*(const Rational &left, const Rational &right);
Rational operator/(const Rational &left, const Rational &right);

bool operator==(const Rational &left, const Rational &right);
bool operator!=(const Rational &left, const Rational &right);
bool operator>(const Rational &left, const Rational &right);
bool operator>=(const Rational &left, const Rational &right);
bool operator<(const Rational &left, const Rational &right);
bool operator<=(const Rational &left, const Rational &right);
//*********************************************** my functions
BigInteger gcd(BigInteger left, BigInteger right) {
  if (left.whatSign() == NEGATIVE) {
    left.changeSign();
  }
  if (right.whatSign() == NEGATIVE) {
    right.changeSign();
  }
  while (right > 0) {
    left %= right;
    left.swapBigInt(right);
  }
  return left;
}

void Rational::normalize() {
  if (numerator == 0) {
    denominator = 1;
  } else {
    if (denominator.whatSign() == NEGATIVE) {
      numerator.changeSign();
      denominator.changeSign();
    }
    BigInteger divider = gcd(numerator, denominator);
    numerator /= divider;
    denominator /= divider;
  }
}

void Rational::swapRational(Rational &other) {
  this->numerator.swapBigInt(other.numerator);
  this->denominator.swapBigInt(other.denominator);
}

BigInteger Rational::getNumerator() const {
  return numerator;
}

BigInteger Rational::getDenominator() const {
  return denominator;
}

BigInteger bigTen(unsigned int pow) {
  BigInteger ans = 1;
  for (unsigned int i = 0; i < pow; ++i) {
    ans *= 10;
  }
  return ans;
}

int BigInteger::order() {
  int tmp = (*this)[size() - 1];
  int ans = 0;
  while (tmp > 0) {
    tmp /= 10;
    ++ans;
  }
  return ans += ((size() - 1) * 9);
}

int myMax(int a, int b) {
  return a > b ? a : b;
}

double powTen(int pow) {
  double ans = 1;
  if (pow > 0) {
    for (int i = 0; i < pow; ++i) {
      ans *= 10;
    }
  } else {
    for (int i = 0; i < -pow; ++i) {
      ans /= 10;
    }
  }
  return ans;
}
//***********************************************

Rational::Rational() {
  numerator = 0;
  denominator = 1;
}

Rational::Rational(int number) {
  numerator = number;
  denominator = 1;
}

Rational::Rational(const BigInteger &bigNumber) {
  numerator = bigNumber;
  denominator = 1;
}

Rational::Rational(const Rational &rationalNum) {
  this->numerator = rationalNum.numerator;
  this->denominator = rationalNum.denominator;
}

Rational &Rational::operator=(const Rational &other) {
  this->numerator = other.numerator;
  this->denominator = other.denominator;
  return *this;
}

Rational &Rational::operator+=(const Rational &other) {
  this->numerator = (this->numerator * other.denominator) + (other.numerator * this->denominator);
  this->denominator *= other.denominator;
  (*this).normalize();
  return *this;
}

Rational &Rational::operator-=(const Rational &other) {
  return *this += -other;
}

Rational &Rational::operator*=(const Rational &other) {
  this->numerator *= other.numerator;
  this->denominator *= other.denominator;
  (*this).normalize();
  return *this;
}

Rational &Rational::operator/=(const Rational &other) {
  if (other.numerator == 0) {
    this->numerator = 0;
    this->denominator = 1;
  } else {
    this->numerator *= other.denominator;
    this->denominator *= other.numerator;
    (*this).normalize();
  }
  return (*this);
}

std::ostream &operator<<(std::ostream &out, const Rational &number) {
  out << number.toString();
  return out;
}

std::istream &operator>>(std::istream &in, Rational &number) {
  in >> number.numerator;
  number.denominator = 1;
  return in;
}

Rational operator+(const Rational &left, const Rational &right) {
  Rational ans = left;
  ans += right;
  return ans;
}

Rational operator-(const Rational &left, const Rational &right) {
  Rational ans = left;
  ans -= right;
  return ans;
}

Rational operator*(const Rational &left, const Rational &right) {
  Rational ans = left;
  ans *= right;
  return ans;
}

Rational operator/(const Rational &left, const Rational &right) {
  Rational ans = left;
  ans /= right;
  return ans;
}

Rational Rational::operator-() const {
  Rational ans = *this;
  ans.numerator.changeSign();
  return ans;
}

bool operator==(const Rational &left, const Rational &right) {
  return (left.getNumerator() == right.getNumerator()) && (left.getDenominator() == right.getDenominator());
}

bool operator!=(const Rational &left, const Rational &right) {
  return !(left == right);
}

bool operator>(const Rational &left, const Rational &right) {
  return (left.getNumerator() * right.getDenominator()) > (right.getNumerator() * left.getDenominator());
}

bool operator<(const Rational &left, const Rational &right) {
  return right > left;
}

bool operator>=(const Rational &left, const Rational &right) {
  return !(left < right);
}

bool operator<=(const Rational &left, const Rational &right) {
  return !(left > right);
}

std::string Rational::asDecimal(size_t precision) const {
  BigInteger ans = numerator * bigTen(precision);
  ans /= denominator;
  if (ans == 0) {
    return "0";
  }
  std::string decimalAns = ans.toString();
  std::string zeros(myMax(precision - ans.order() + 1, 0), '0');
  if (decimalAns[0] == '-') {
    decimalAns.insert(1, zeros);
  } else {
    decimalAns.insert(0, zeros);
  }
  if (precision > 0) {
    decimalAns.insert(decimalAns.size() - precision, ".");
  }
  return decimalAns;
}

std::string Rational::toString() const {
  if (denominator == 1) {
    return numerator.toString();
  } else {
    return numerator.toString() + "/" + denominator.toString();
  }
}

Rational::operator double() {//15 - double precision
  if (numerator == 0) {
    return 0;
  }
  const int DOUBLEMAXSIZE = 15;
  int numOrder = numerator.order();
  int denomOrder = denominator.order();
  BigInteger tmp = numerator;
  std::string buffer;
  int stUp = DOUBLEMAXSIZE + myMax(denomOrder - numOrder, 1);
  tmp *= bigTen(stUp);
  tmp /= denominator;
  int curOrder = tmp.order();
  curOrder -= stUp;
  double nextNum;
  double ans = 0;
  buffer = tmp.toString();
  if (buffer[0] == '-') {
    buffer.erase(buffer.begin());
  }
  double tenSt = powTen(curOrder - 1);
  for (int i = 0; i < DOUBLEMAXSIZE; i++) {
    nextNum = buffer[i] - '0';
    nextNum *= tenSt;
    ans += nextNum;
    tenSt /= 10;
  }
  if (tmp.whatSign() == NEGATIVE) {
    return -ans;
  } else {
    return ans;
  }
}

#endif // rational
