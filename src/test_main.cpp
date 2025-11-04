#include "include/int2048.h"
#include <iostream>
using sjtu::int2048;
int main(){
  int2048 a("123456789012345678901234567890"), b("-98765432109876543210");
  std::cout << (a+b) << '\n';
  std::cout << (a-b) << '\n';
  std::cout << (b-a) << '\n';
  int2048 c("-10"), d("3");
  std::cout << (c/d) << ' ' << (c%d) << '\n'; // -4 2
  std::cout << (int2048(10)/int2048(-3)) << ' ' << (int2048(10)%int2048(-3)) << '\n'; // -4 -2
  std::cout << (int2048(-10)/int2048(-3)) << ' ' << (int2048(-10)%int2048(-3)) << '\n'; // 3 -1
  std::cout << (int2048(0)/int2048(5)) << ' ' << (int2048(0)%int2048(5)) << '\n';
  int2048 m("12345678901234567890");
  int2048 n("98765432109876543210");
  std::cout << (m*n) << '\n';
  std::cout << (n/m) << ' ' << (n%m) << '\n';
  return 0;
}
