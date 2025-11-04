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
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <cstdint>

namespace sjtu {

// ----- helpers -----
void int2048::trim() {
  while (!a.empty() && a.back() == 0) a.pop_back();
  if (a.empty()) neg = false;
}

bool int2048::is_zero() const { return a.empty(); }

int int2048::cmp_abs(const int2048 &x, const int2048 &y) {
  if (x.a.size() != y.a.size()) return x.a.size() < y.a.size() ? -1 : 1;
  for (int i = (int)x.a.size() - 1; i >= 0; --i) {
    if (x.a[i] != y.a[i]) return x.a[i] < y.a[i] ? -1 : 1;
  }
  return 0;
}

void int2048::add_abs_to(int2048 &x, const int2048 &y) {
  uint64_t carry = 0;
  size_t n = std::max(x.a.size(), y.a.size());
  if (x.a.size() < n) x.a.resize(n, 0);
  for (size_t i = 0; i < n; ++i) {
    uint64_t cur = carry + (uint64_t)x.a[i] + (i < y.a.size() ? (uint64_t)y.a[i] : 0u);
    x.a[i] = (uint32_t)(cur % BASE);
    carry = cur / BASE;
  }
  if (carry) x.a.push_back((uint32_t)carry);
}

void int2048::sub_abs_to(int2048 &x, const int2048 &y) {
  // assumes |x| >= |y|
  int64_t carry = 0;
  for (size_t i = 0; i < x.a.size(); ++i) {
    int64_t cur = (int64_t)x.a[i] - (i < y.a.size() ? (int64_t)y.a[i] : 0) + carry;
    if (cur < 0) {
      cur += BASE;
      carry = -1;
    } else {
      carry = 0;
    }
    x.a[i] = (uint32_t)cur;
  }
  x.trim();
}

static inline std::vector<uint32_t> add_limbs(const std::vector<uint32_t> &x,
                                              const std::vector<uint32_t> &y) {
  std::vector<uint32_t> r;
  size_t n = std::max(x.size(), y.size());
  r.resize(n);
  uint64_t carry = 0;
  for (size_t i = 0; i < n; ++i) {
    uint64_t cur = carry + (i < x.size() ? x[i] : 0u) + (i < y.size() ? y[i] : 0u);
    r[i] = (uint32_t)(cur % sjtu::int2048::BASE);
    carry = cur / sjtu::int2048::BASE;
  }
  if (carry) r.push_back((uint32_t)carry);
  return r;
}

std::vector<uint32_t> int2048::mul_simple(const std::vector<uint32_t> &x,
                                          const std::vector<uint32_t> &y) {
  if (x.empty() || y.empty()) return {};
  std::vector<uint64_t> tmp(x.size() + y.size(), 0);
  for (size_t i = 0; i < x.size(); ++i) {
    uint64_t carry = 0;
    uint64_t xi = x[i];
    for (size_t j = 0; j < y.size(); ++j) {
      unsigned __int128 cur = (unsigned __int128)tmp[i + j] + xi * (uint64_t)y[j] + carry;
      tmp[i + j] = (uint64_t)(cur % BASE);
      carry = (uint64_t)(cur / BASE);
    }
    tmp[i + y.size()] += carry;
  }
  std::vector<uint32_t> r(tmp.size());
  for (size_t i = 0; i < tmp.size(); ++i) r[i] = (uint32_t)tmp[i];
  while (!r.empty() && r.back() == 0) r.pop_back();
  return r;
}

static void normalize_from_ll(std::vector<int64_t> &v, uint32_t BASE) {
  long long carry = 0;
  for (size_t i = 0; i < v.size(); ++i) {
    long long cur = v[i] + carry;
    if (cur >= 0) {
      v[i] = (int64_t)(cur % BASE);
      carry = (long long)(cur / BASE);
    } else {
      long long k = (-(cur + 1)) / BASE + 1;
      cur += k * BASE;
      v[i] = (int64_t)cur;
      carry = -k;
    }
  }
  while (carry > 0) {
    v.push_back((int64_t)(carry % BASE));
    carry /= BASE;
  }
  while (!v.empty() && v.back() == 0) v.pop_back();
}

std::vector<uint32_t> int2048::karatsuba(const std::vector<uint32_t> &x,
                                         const std::vector<uint32_t> &y) {
  size_t n = x.size(), m = y.size();
  if (!n || !m) return {};
  if (std::min(n, m) < 32) return mul_simple(x, y);
  size_t k = std::min(n, m) / 2;

  std::vector<uint32_t> x0(x.begin(), x.begin() + k);
  std::vector<uint32_t> x1(x.begin() + k, x.end());
  std::vector<uint32_t> y0(y.begin(), y.begin() + k);
  std::vector<uint32_t> y1(y.begin() + k, y.end());

  auto z0 = karatsuba(x0, y0);
  auto z2 = karatsuba(x1, y1);

  // (x0 + x1) * (y0 + y1)
  auto sx = add_limbs(x0, x1);
  auto sy = add_limbs(y0, y1);
  auto z1 = karatsuba(sx, sy);

  // z1 = z1 - z0 - z2, use int64 buffer to allow negatives
  std::vector<int64_t> mid(std::max({z1.size(), z0.size(), z2.size()}), 0);
  for (size_t i = 0; i < z1.size(); ++i) mid[i] += z1[i];
  for (size_t i = 0; i < z0.size(); ++i) mid[i] -= z0[i];
  for (size_t i = 0; i < z2.size(); ++i) mid[i] -= z2[i];

  // assemble result: z0 + (mid << k) + (z2 << 2k)
  std::vector<int64_t> res_ll(z0.size() + mid.size() + z2.size(), 0);
  for (size_t i = 0; i < z0.size(); ++i) res_ll[i] += z0[i];
  for (size_t i = 0; i < mid.size(); ++i) res_ll[i + k] += mid[i];
  for (size_t i = 0; i < z2.size(); ++i) res_ll[i + 2 * k] += z2[i];

  // normalize back to base
  normalize_from_ll(res_ll, BASE);
  std::vector<uint32_t> res(res_ll.size());
  for (size_t i = 0; i < res.size(); ++i) res[i] = (uint32_t)res_ll[i];
  while (!res.empty() && res.back() == 0) res.pop_back();
  return res;
}

std::vector<uint32_t> int2048::multiply_limbs(const std::vector<uint32_t> &x,
                                              const std::vector<uint32_t> &y) {
  if (std::min(x.size(), y.size()) < 32) return mul_simple(x, y);
  return karatsuba(x, y);
}

void int2048::divmod_abs(const int2048 &u, const int2048 &v, int2048 &q, int2048 &r) {
  // assumes v != 0 and both are non-negative
  if (v.is_zero()) return; // undefined, but avoid crash
  if (cmp_abs(u, v) < 0) {
    q = int2048(0);
    r = u;
    return;
  }
  if (v.a.size() == 1) {
    // single-limb divisor
    uint64_t div = v.a[0];
    q.a.assign(u.a.size(), 0);
    uint64_t rem = 0;
    for (int i = (int)u.a.size() - 1; i >= 0; --i) {
      unsigned __int128 cur = (unsigned __int128)rem * BASE + u.a[i];
      uint64_t qq = (uint64_t)(cur / div);
      rem = (uint64_t)(cur % div);
      q.a[i] = (uint32_t)qq;
    }
    q.neg = false;
    q.trim();
    r.a.clear();
    if (rem) r.a.push_back((uint32_t)rem);
    r.neg = false;
    r.trim();
    return;
  }
  size_t n = v.a.size();
  size_t m = u.a.size() - n;
  // Normalize
  uint32_t d = (uint32_t)(BASE / (uint64_t)(v.a.back() + 1));
  std::vector<uint32_t> vn(n);
  std::vector<uint32_t> un(u.a.size() + 1, 0);
  uint64_t carry = 0;
  for (size_t i = 0; i < u.a.size(); ++i) {
    uint64_t cur = (uint64_t)u.a[i] * d + carry;
    un[i] = (uint32_t)(cur % BASE);
    carry = cur / BASE;
  }
  un[u.a.size()] = (uint32_t)carry;
  carry = 0;
  for (size_t i = 0; i < n; ++i) {
    uint64_t cur = (uint64_t)v.a[i] * d + carry;
    vn[i] = (uint32_t)(cur % BASE);
    carry = cur / BASE;
  }

  q.a.assign(m + 1, 0);
  const uint32_t vn1 = vn[n - 1];
  const uint32_t vn2 = vn[n - 2];

  for (int j = (int)m; j >= 0; --j) {
    unsigned __int128 numerator = (unsigned __int128)un[j + n] * BASE + un[j + n - 1];
    uint64_t qhat = (uint64_t)(numerator / vn1);
    uint64_t rhat = (uint64_t)(numerator % vn1);
    while (qhat == BASE || (qhat * vn2 > (unsigned __int128)rhat * BASE + un[j + n - 2])) {
      --qhat;
      rhat += vn1;
      if (rhat >= BASE) break;
    }
    // multiply and subtract
    unsigned __int128 borrow = 0, carry2 = 0;
    for (size_t i = 0; i < n; ++i) {
      unsigned __int128 p = (unsigned __int128)vn[i] * qhat + carry2;
      carry2 = p / BASE;
      unsigned __int128 sub = (unsigned __int128)un[j + i] - (uint64_t)(p % BASE) - borrow;
      if ((int64_t)sub < 0) {
        un[j + i] = (uint32_t)(sub + BASE);
        borrow = 1;
      } else {
        un[j + i] = (uint32_t)sub;
        borrow = 0;
      }
    }
    unsigned __int128 sub = (unsigned __int128)un[j + n] - carry2 - borrow;
    if ((int64_t)sub < 0) {
      // qhat too big, decrement and add back vn
      --qhat;
      unsigned __int128 c = 0;
      for (size_t i = 0; i < n; ++i) {
        unsigned __int128 s = (unsigned __int128)un[j + i] + vn[i] + c;
        if (s >= BASE) {
          un[j + i] = (uint32_t)(s - BASE);
          c = 1;
        } else {
          un[j + i] = (uint32_t)s;
          c = 0;
        }
      }
      un[j + n] += (uint32_t)c;
    } else {
      un[j + n] = (uint32_t)sub;
    }
    q.a[j] = (uint32_t)qhat;
  }
  q.neg = false;
  q.trim();
  // remainder: un[0..n-1] / d
  r.a.assign(n, 0);
  uint64_t rem = 0;
  for (int i = (int)n - 1; i >= 0; --i) {
    unsigned __int128 cur = (unsigned __int128)rem * BASE + un[i];
    r.a[i] = (uint32_t)(cur / d);
    rem = (uint64_t)(cur % d);
  }
  r.neg = false;
  r.trim();
}

// ----- constructors -----
int2048::int2048() : neg(false) {}

int2048::int2048(long long v) { *this = int2048(); if (v < 0) { neg = true; v = -v; } else neg = false; while (v) { a.push_back((uint32_t)(v % BASE)); v /= BASE; } trim(); }

int2048::int2048(const std::string &s) { read(s); }

int2048::int2048(const int2048 &o) = default;

// ----- basic ops -----
void int2048::read(const std::string &s) {
  a.clear();
  neg = false;
  size_t i = 0, n = s.size();
  while (i < n && std::isspace((unsigned char)s[i])) ++i;
  bool sign_set = false;
  if (i < n && (s[i] == '+' || s[i] == '-')) { neg = (s[i] == '-'); ++i; sign_set = true; }
  while (i < n && std::isspace((unsigned char)s[i])) ++i;
  // skip leading zeros
  while (i < n && s[i] == '0') ++i;
  std::vector<uint32_t> parts;
  for (long long j = (long long)n - 1; j >= (long long)i; j -= 9) {
    int l = (int)std::max<long long>(i, j - 8);
    int len = (int)(j - l + 1);
    uint32_t val = 0;
    for (int k = 0; k < len; ++k) val = val * 10 + (s[l + k] - '0');
    parts.push_back(val);
  }
  a = parts;
  trim();
  if (is_zero()) neg = false;
}

void int2048::print() {
  if (is_zero()) { std::cout << 0; return; }
  if (neg) std::cout << '-';
  std::cout << a.back();
  for (int i = (int)a.size() - 2; i >= 0; --i) {
    std::cout << std::setw(9) << std::setfill('0') << a[i];
  }
}

int2048 &int2048::add(const int2048 &o) {
  if (o.is_zero()) return *this;
  if (this->is_zero()) { *this = o; return *this; }
  if (neg == o.neg) {
    add_abs_to(*this, o);
  } else {
    int c = cmp_abs(*this, o);
    if (c >= 0) {
      sub_abs_to(*this, o);
      // sign stays the same
    } else {
      int2048 tmp = o;
      sub_abs_to(tmp, *this);
      *this = std::move(tmp);
    }
  }
  trim();
  return *this;
}

int2048 add(int2048 x, const int2048 &y) { return x.add(y); }

int2048 &int2048::minus(const int2048 &o) {
  if (o.is_zero()) return *this;
  if (this->is_zero()) { *this = o; this->neg = !this->neg; return *this; }
  if (neg != o.neg) {
    add_abs_to(*this, o);
  } else {
    int c = cmp_abs(*this, o);
    if (c >= 0) {
      sub_abs_to(*this, o);
      // sign unchanged
    } else {
      int2048 tmp = o;
      sub_abs_to(tmp, *this);
      tmp.neg = !tmp.neg; // result sign flips
      *this = std::move(tmp);
    }
  }
  trim();
  return *this;
}

int2048 minus(int2048 x, const int2048 &y) { return x.minus(y); }

// ----- operators -----
int2048 int2048::operator+() const { return *this; }

int2048 int2048::operator-() const {
  int2048 r(*this);
  if (!r.is_zero()) r.neg = !r.neg;
  return r;
}

int2048 &int2048::operator=(const int2048 &o) = default;

int2048 &int2048::operator+=(const int2048 &o) { return add(o); }

int2048 operator+(int2048 x, const int2048 &y) { return x += y; }

int2048 &int2048::operator-=(const int2048 &o) { return minus(o); }

int2048 operator-(int2048 x, const int2048 &y) { return x -= y; }

int2048 &int2048::operator*=(const int2048 &o) {
  if (this->is_zero() || o.is_zero()) { a.clear(); neg = false; return *this; }
  auto prod = multiply_limbs(a, o.a);
  a.swap(prod);
  neg = neg ^ o.neg;
  trim();
  return *this;
}

int2048 operator*(int2048 x, const int2048 &y) { return x *= y; }

int2048 &int2048::operator/=(const int2048 &o) {
  // floor division
  if (o.is_zero()) return *this; // undefined, keep
  if (this->is_zero()) { neg = false; return *this; }
  int2048 ua = *this; ua.neg = false;
  int2048 va = o; va.neg = false;
  int2048 q, r;
  divmod_abs(ua, va, q, r); // q,r non-negative
  bool signs_diff = (neg != o.neg);
  if (!r.is_zero() && signs_diff) {
    // q = q + 1 towards more negative
    // floor(-a/b) = -(a/b) - 1 when a%b != 0
    // here q is trunc toward zero of |a|/|b|
    // For floor: decrement (make more negative)
    int2048 one(1);
    q.add(one); // q = q + 1
  }
  q.neg = signs_diff && !q.is_zero();
  // set
  *this = q;
  trim();
  return *this;
}

int2048 operator/(int2048 x, const int2048 &y) { return x /= y; }

int2048 &int2048::operator%=(const int2048 &o) {
  if (o.is_zero()) return *this; // undefined
  if (this->is_zero()) { neg = false; return *this; }
  int2048 ua = *this; ua.neg = false;
  int2048 va = o; va.neg = false;
  int2048 q, r;
  divmod_abs(ua, va, q, r); // q,r >= 0
  bool signs_diff = (neg != o.neg);
  if (!r.is_zero() && signs_diff) {
    // r = |b| - r
    int2048 rb = va;
    rb.minus(r); // rb = |b| - r
    r = std::move(rb);
  }
  r.neg = o.neg && !r.is_zero();
  *this = r;
  trim();
  return *this;
}

int2048 operator%(int2048 x, const int2048 &y) { return x %= y; }

std::istream &operator>>(std::istream &is, int2048 &x) {
  std::string s; is >> s; x.read(s); return is;
}

std::ostream &operator<<(std::ostream &os, const int2048 &x) {
  if (x.is_zero()) { os << 0; return os; }
  if (x.neg) os << '-';
  os << x.a.back();
  for (int i = (int)x.a.size() - 2; i >= 0; --i) {
    os << std::setw(9) << std::setfill('0') << x.a[i];
  }
  return os;
}

bool operator==(const int2048 &x, const int2048 &y) {
  return x.neg == y.neg && x.a == y.a;
}

bool operator!=(const int2048 &x, const int2048 &y) { return !(x == y); }

bool operator<(const int2048 &x, const int2048 &y) {
  if (x.neg != y.neg) return x.neg && !x.is_zero();
  int c = int2048::cmp_abs(x, y);
  return x.neg ? (c > 0) : (c < 0);
}

bool operator>(const int2048 &x, const int2048 &y) { return y < x; }

bool operator<=(const int2048 &x, const int2048 &y) { return !(y < x); }

bool operator>=(const int2048 &x, const int2048 &y) { return !(x < y); }

} // namespace sjtu
