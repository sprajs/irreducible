// Developer native audit of external original release bytes. Not a default
// test requiring survey data. Decimal rounding intervals are a HYPOTHESIS:
// independently nearest-rounded copies of one symmetric number must differ
// by at most one decimal quantum. Violations disprove that narrow explanation,
// not the published scientific uncertainty model. No symmetrization is
// performed.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
namespace {
long double quantum(const std::string &token) {
  auto e = token.find_first_of("eE");
  auto mantissa = token.substr(0, e);
  auto dot = mantissa.find('.');
  int decimals = dot == std::string::npos
                     ? 0
                     : static_cast<int>(mantissa.size() - dot - 1);
  int exponent = e == std::string::npos ? 0 : std::stoi(token.substr(e + 1));
  return std::pow(10.L, exponent - decimals);
}
std::vector<std::string> words(const std::string &s) {
  std::istringstream in(s);
  std::vector<std::string> w;
  for (std::string x; in >> x;)
    w.push_back(x);
  return w;
}
} // namespace
int main(int argc, char **argv) {
  try {
    if (argc != 3)
      throw std::runtime_error(
          "usage: test_w01_covariance_audit COVARIANCE TABLE");
    std::ifstream cov(argv[1]);
    std::size_t n = 0;
    if (!(cov >> n) || n == 0 || n > 2000)
      throw std::runtime_error("invalid or unsupported dimension");
    std::vector<long double> a(n * n), q(n * n);
    std::map<int, std::size_t> quantization;
    std::string token;
    for (std::size_t i = 0; i < a.size(); ++i) {
      if (!(cov >> token))
        throw std::runtime_error("truncated covariance");
      std::size_t used = 0;
      a[i] = std::stold(token, &used);
      if (used != token.size() || !std::isfinite(a[i]))
        throw std::runtime_error("invalid finite decimal");
      q[i] = quantum(token);
      ++quantization[static_cast<int>(std::llround(std::log10(q[i])))];
    }
    if (cov >> token)
      throw std::runtime_error("trailing covariance token");
    std::ifstream table(argv[2]);
    std::string line;
    if (!std::getline(table, line))
      throw std::runtime_error("missing table");
    auto headers = words(line);
    auto z = std::find(headers.begin(), headers.end(), "zHD");
    if (z == headers.end())
      throw std::runtime_error("missing zHD");
    auto zi = static_cast<std::size_t>(z - headers.begin());
    std::vector<std::size_t> kept;
    std::size_t rows = 0;
    while (std::getline(table, line)) {
      auto row = words(line);
      if (row.empty())
        continue;
      if (row.size() != headers.size())
        throw std::runtime_error("table column count");
      auto value = std::stold(row[zi]);
      if (!std::isfinite(value))
        throw std::runtime_error("nonfinite zHD");
      if (value > .01L)
        kept.push_back(rows);
      ++rows;
    }
    if (rows != n)
      throw std::runtime_error("matrix table rows differ");
    auto audit = [&](const std::vector<std::size_t> &idx, const char *label) {
      long double maximum = 0, delta_inf = 0, frobenius = 0,
                  diagmin = std::numeric_limits<long double>::max(),
                  diagmax = 0, maxratio = 0;
      std::size_t differing = 0, incompatible = 0, maxi = 0, maxj = 0;
      for (auto i : idx) {
        diagmin = std::min(diagmin, a[i * n + i]);
        diagmax = std::max(diagmax, a[i * n + i]);
        long double rowsum = 0;
        for (auto j : idx) {
          auto diff = std::abs(a[i * n + j] - a[j * n + i]);
          rowsum += diff / 2;
          frobenius += diff * diff / 4;
          if (i < j && diff > 0) {
            ++differing;
            const auto intervalsum = (q[i * n + j] + q[j * n + i]) / 2;
            auto ratio = diff / intervalsum;
            maxratio = std::max(maxratio, ratio);
            if (diff > intervalsum * (1 + 1e-8L))
              ++incompatible;
          }
          if (diff > maximum) {
            maximum = diff;
            maxi = i;
            maxj = j;
          }
        }
        delta_inf = std::max(delta_inf, rowsum);
      }
      std::cout << std::setprecision(21) << "{\"scope\":\"" << label
                << "\",\"rows\":" << idx.size()
                << ",\"max_asymmetry_mag2\":" << maximum
                << ",\"max_asymmetry_original_indices\":[" << maxi << ","
                << maxj << "],\"diagonal_min_mag2\":" << diagmin
                << ",\"diagonal_max_mag2\":" << diagmax
                << ",\"differing_unordered_pairs\":" << differing
                << ",\"nearest_decimal_common_value_interval_violations\":"
                << incompatible
                << ",\"maximum_interval_discrepancy_ratio\":" << maxratio
                << ",\"symmetric_average_delta_norm_inf\":" << delta_inf
                << ",\"symmetric_average_delta_norm_frobenius\":"
                << std::sqrt(frobenius) << "}\n";
    };
    std::vector<std::size_t> all(n);
    for (std::size_t i = 0; i < n; ++i)
      all[i] = i;
    audit(all, "original_full");
    audit(kept, "strict_zHD_gt_0.01");
    for (auto [exponent, count] : quantization)
      std::cout << "{\"decimal_quantum_exponent\":" << exponent
                << ",\"tokens\":" << count << "}\n";
    return 0;
  } catch (const std::exception &e) {
    std::fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
