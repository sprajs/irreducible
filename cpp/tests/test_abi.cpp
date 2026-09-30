#include "irred/abi.h"
#include "irred/arithmetic.hpp"
#include <iostream>
#include <limits>
#include <stdexcept>
#define CHECK(condition)                                                       \
  do {                                                                         \
    if (!(condition))                                                          \
      throw std::runtime_error("ABI fixture check failed: " #condition);       \
  } while (false)
static irred_i64_buffer buffer(const int64_t *p, uint64_t n) {
  return {sizeof(irred_i64_buffer), IRRED_ABI_VERSION, 1, 0, p, n, n * 8};
}
int main() {
  const int64_t a[] = {0, 4, -7, 9223372036854775806LL}, b[] = {0, -2, 3, 1};
  auto x = buffer(a, 4), y = buffer(b, 4);
  irred_result *r = nullptr;
  CHECK(irred_add(&x, &y, 0, &r) == 0);
  const int64_t *p;
  uint64_t n;
  CHECK(irred_result_view(r, &p, &n) == 0 && n == 4);
  CHECK(p[0] == 0 && p[1] == 2 && p[2] == -4 &&
        p[3] == std::numeric_limits<int64_t>::max());
  int64_t direct;
  for (int i = 0; i < 4; ++i) {
    CHECK(irred::checked_add(a[i], b[i], direct));
    CHECK(direct == p[i]);
  }
  std::cout << "[0,2,-4,9223372036854775807]\n";
  CHECK(irred_result_destroy(r) == 0);
  x = buffer(nullptr, 0);
  y = x;
  CHECK(irred_add(&x, &y, 0, &r) == 0);
  CHECK(irred_result_destroy(r) == 0);
  x.abi_version = IRRED_ABI_VERSION + 1;
  CHECK(irred_add(&x, &y, 0, &r) == IRRED_ABI_MISMATCH && !r);
  x = y;
  x.length = 1;
  CHECK(irred_add(&x, &y, 0, &r) == IRRED_INVALID_INPUT && !r);
  x = y;
  CHECK(irred_add(&x, &y, 1, &r) == IRRED_EXCEPTION && !r);
  CHECK(irred_add(&x, &y, 2, &r) == IRRED_ALLOCATION_FAILURE && !r);
  const int64_t max = std::numeric_limits<int64_t>::max(), one = 1;
  x = buffer(&max, 1);
  y = buffer(&one, 1);
  CHECK(irred_add(&x, &y, 0, &r) == IRRED_OVERFLOW && !r);
}
