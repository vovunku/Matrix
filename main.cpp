#include <gtest/gtest.h>
#include <iostream>
#include <My_traits.h>
#include <Finite.h>
#include <Matrix.h>
#include <test_matrix.h>

TEST(Prime, check) {
bool tmp;
tmp = is_prime<5>::value;
ASSERT_TRUE(tmp);
tmp = is_prime<9>::value;
ASSERT_FALSE(tmp);
tmp = is_prime<1>::value;
ASSERT_FALSE(tmp);
tmp = is_prime<0>::value;
ASSERT_FALSE(tmp);
tmp = is_prime<17>::value;
ASSERT_TRUE(tmp);
tmp = is_prime<24>::value;
ASSERT_FALSE(tmp);
}

TEST(Finite, CE) {
//Finite<0> a = 5;
//Finite<12> b = 5;
//std::cout << b.inverse().get_value();
}

TEST(Finite, add) {
Finite<9> a;
Finite<9> b = 5;
a += b;
ASSERT_EQ(a, 5);
a += 10;
ASSERT_EQ(a.get_value(), 6);
ASSERT_EQ((a + b).get_value(), 2);
}

TEST(Finite, multilication) {
Finite<10> a = 2;
Finite<10> b = 5;
ASSERT_EQ((a * b).get_value(), 0);
}

TEST(Finite, p_adic_pow) {
Finite<13> a = 2;
ASSERT_EQ(p_adic_pow(a, 4).get_value(), 3);
a = p_adic_pow(a, 4);
ASSERT_EQ(p_adic_pow(a, 3).get_value(), 1);
}

TEST(Finite, inverse) {
Finite<23> a = 3;
ASSERT_EQ(a.inverse().get_value(), 8);
a = 5;
ASSERT_EQ(a.inverse().get_value(), 14);
}

TEST(Finite, division) {
Finite<23> a = 5;
Finite<23> b = 3;
ASSERT_EQ((a / b).get_value(), 17);
a = 10;
b = 4;
ASSERT_EQ((a / b).get_value(), 14);
}

TEST(Finte, cout) {
Finite<5> a = 5;
//std::cout << a;
}

TEST(Matrix, print) {
std::vector<std::vector<Rational> > a = {{1, 2, 5}, {3, 15, 27}};
Matrix<2, 3> b(a);
b.print_matrix();
}

TEST(Matrix, add) {
std::vector<std::vector<Rational> > a = {{1, 2, 3, 4, 5, 6},
                                         {7, 8, 9, 10, 11, 12},
                                         {13, 14, 15, 16, 17, 18},
                                         {19, 20, 21, 22, 23, 24}};
Matrix<4, 6> b(a);
Matrix<4, 6> c = b;
b += c;
for (int i = 0; i < 4; ++i) {
for (int j = 0; j < 6; ++j) {
a[i][j] *= 2;
}
}
Matrix<4, 6> ans(a);
ASSERT_EQ(b, ans);
Matrix<4, 6> zero;
b.print_matrix();
b -= b;
//ASSERT_EQ(b, zero);
}

TEST(Matrix, multiplication_field) {
std::vector<std::vector<Rational> > a = {{1, 2, 3, 4, 5, 6},
                                         {7, 8, 9, 10, 11, 12},
                                         {13, 14, 15, 16, 17, 18},
                                         {19, 20, 21, 22, 23, 24}};
Matrix<4, 6> b(a);
Matrix<4, 6> c = b;
Matrix<4, 6> d = b;
d *= 3;
c += b;
c += b;
//d.print_matrix();
std::cout << "\n";
//c.print_matrix();
ASSERT_EQ(d, c);
}

/*TEST(Matrix, my_log){
    Matrix<2, 2> a;
    ASSERT_EQ(a._nearest_pov_2(7), 8);
    ASSERT_EQ(a._nearest_pov_2(8), 8);
    ASSERT_EQ(a._nearest_pov_2(1), 1);
    ASSERT_EQ(a._nearest_pow_2(0), 1);
}*/

TEST(Matrix, multiplication) {
std::vector<std::vector<Rational> > a = {{1, 2, 3, 4, 5, 6},
                                         {7, 8, 9, 10, 11, 12},
                                         {13, 14, 15, 16, 17, 18},
                                         {19, 20, 21, 22, 23, 24},
                                         {25, 26, 27, 28, 29, 30}};
std::vector<std::vector<Rational> > d = {{130, 140, 150, 160},
                                         {370, 404, 438, 472},
                                         {610, 668, 726, 784},
                                         {850, 932, 1014, 1096}};
Matrix<4, 4> b(a);
Matrix<4, 4> c(a);
Matrix<4, 4> e(d);
ASSERT_EQ(b.multiplicate(c), e);
b.multiplicate(c).print_matrix();
Matrix<4, 5> f(a);
Matrix<5, 6> g(a);
(f.multiplicate(g)).print_matrix();
}

/*255	270	285	300	315	330
645	690	735	780	825	870
1035	1110	1185	1260	1335	1410
1425	1530	1635	1740	1845	1950*/

/*TEST(Matrix, Gauss){
    std::vector <std::vector<Rational> > a = {{1, 2, 3, 4, 5, 6},
                                              {12, 11, 10, 9, 8, 7},
                                              {13, 15, 14, 16, 17, 18},
                                              {19, 20, 21, 5, 23, 24},
                                              {25, 26, 27, 28, 0, 30}};
    Matrix<5, 6> g(a);
    Matrix<2, 6> b(a);
    g._forward_gauss(g._matrix, 5);
    g._backward_gauss(g._matrix, 5);
    //g._gauss(g._matrix, 5);
    b._forward_gauss(b._matrix, 2);
    BigInteger e = 111111;
    e *= 1000000;
    e += 111111;
    e *= e;
    std::cout << e;
}*/

TEST(Matrix, det) {
std::vector<std::vector<Rational> > a = {{1, 2, 3, 4, 5, 6},
                                         {12, 11, 10, 9, 8, 7},
                                         {13, 15, 14, 16, 17, 18},
                                         {19, 20, 21, 5, 23, 24},
                                         {25, 26, 27, 28, 0, 30}};
Matrix<5, 5> b(a);
ASSERT_EQ(b.det(), 19227);
}

TEST(Matrix, all) {
std::vector<std::vector<Rational> > a = {{1, 2, 3, 4, 5, 6},
                                         {12, 11, 10, 9, 8, 7},
                                         {13, 15, 14, 16, 17, 18},
                                         {19, 20, 21, 5, 23, 24},
                                         {25, 26, 27, 28, 0, 30}};
std::vector<std::vector<Finite<5>>> d = {{1, 2, 3, 4, 5, 6},
                                         {12, 11, 10, 9, 8, 7},
                                         {13, 15, 14, 16, 17, 18},
                                         {19, 20, 21, 5, 23, 24},
                                         {25, 26, 27, 28, 0, 30}};
SquareMatrix<5> b(a);
SquareMatrix<5, Finite<5>> c(d);
}

TEST(Matrix, bigone) {
std::cout << ((-1 % 17) + 17) % 17;
testMatrix();
}