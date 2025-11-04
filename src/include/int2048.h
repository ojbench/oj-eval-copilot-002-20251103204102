#pragma once
#ifndef SJTU_BIGINTEGER
#define SJTU_BIGINTEGER

// Integer 1:
// Implement a signed big integer class that only needs to support simple addition and subtraction

// Integer 2:
// Implement a signed big integer class that supports addition, subtraction, multiplication, and division, and overload related operators

// Do not use any header files other than the following
#include <complex>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>

// Do not use "using namespace std;"
#include <cstdint>

namespace sjtu {
class int2048 {
public:
  static const uint32_t BASE = 1000000000u; // 1e9
private:
  std::vector<uint32_t> a;                  // little-endian limbs
  bool neg;

  // helpers (implemented in .cpp)
  void trim();
  bool is_zero() const;
  static int cmp_abs(const int2048 &x, const int2048 &y);
  static void add_abs_to(int2048 &x, const int2048 &y);   // x = |x| + |y|
  static void sub_abs_to(int2048 &x, const int2048 &y);   // x = |x| - |y|, assumes |x|>=|y|

  static std::vector<uint32_t> mul_simple(const std::vector<uint32_t> &x,
                                          const std::vector<uint32_t> &y);
  static std::vector<uint32_t> karatsuba(const std::vector<uint32_t> &x,
                                         const std::vector<uint32_t> &y);
  static std::vector<uint32_t> multiply_limbs(const std::vector<uint32_t> &x,
                                              const std::vector<uint32_t> &y);
  static void divmod_abs(const int2048 &u, const int2048 &v, int2048 &q, int2048 &r);

public:
  // Constructors
  int2048();
  int2048(long long);
  int2048(const std::string &);
  int2048(const int2048 &);

  // The parameter types of the following functions are for reference only, you can choose to use constant references or not
  // If needed, you can add other required functions yourself
  // ===================================
  // Integer1
  // ===================================

  // Read a big integer
  void read(const std::string &);
  // Output the stored big integer, no need for newline
  void print();

  // Add a big integer
  int2048 &add(const int2048 &);
  // Return the sum of two big integers
  friend int2048 add(int2048, const int2048 &);

  // Subtract a big integer
  int2048 &minus(const int2048 &);
  // Return the difference of two big integers
  friend int2048 minus(int2048, const int2048 &);

  // ===================================
  // Integer2
  // ===================================

  int2048 operator+() const;
  int2048 operator-() const;

  int2048 &operator=(const int2048 &);

  int2048 &operator+=(const int2048 &);
  friend int2048 operator+(int2048, const int2048 &);

  int2048 &operator-=(const int2048 &);
  friend int2048 operator-(int2048, const int2048 &);

  int2048 &operator*=(const int2048 &);
  friend int2048 operator*(int2048, const int2048 &);

  int2048 &operator/=(const int2048 &);
  friend int2048 operator/(int2048, const int2048 &);

  int2048 &operator%=(const int2048 &);
  friend int2048 operator%(int2048, const int2048 &);

  friend std::istream &operator>>(std::istream &, int2048 &);
  friend std::ostream &operator<<(std::ostream &, const int2048 &);

  friend bool operator==(const int2048 &, const int2048 &);
  friend bool operator!=(const int2048 &, const int2048 &);
  friend bool operator<(const int2048 &, const int2048 &);
  friend bool operator>(const int2048 &, const int2048 &);
  friend bool operator<=(const int2048 &, const int2048 &);
  friend bool operator>=(const int2048 &, const int2048 &);
};
} // namespace sjtu

#endif
