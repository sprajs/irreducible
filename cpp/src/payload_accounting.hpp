#pragma once
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <vector>
namespace irred::detail {
// Requested retained payload, not allocator bookkeeping/RSS. Conservatively
// charges capacity+1 for every owned string, including inline SSO storage.
class PayloadAccounting {
  size_t total_;
  bool valid_ = true;

public:
  explicit PayloadAccounting(size_t owner) : total_(owner) {}
  void add(size_t count, size_t width) noexcept {
    if (!valid_ ||
        count > (std::numeric_limits<size_t>::max() - total_) / width) {
      valid_ = false;
      return;
    }
    total_ += count * width;
  }
  template <class T> void vector(const std::vector<T> &v) noexcept {
    add(v.capacity(), sizeof(T));
  }
  void string(const std::string &s) noexcept {
    if (s.capacity() == std::numeric_limits<size_t>::max()) {
      valid_ = false;
      return;
    }
    add(s.capacity() + 1, 1);
  }
  void strings(const std::vector<std::string> &v) noexcept {
    vector(v);
    for (const auto &s : v)
      string(s);
  }
  void embedded(std::optional<size_t> total, size_t owner) noexcept {
    if (!total || *total < owner) {
      valid_ = false;
      return;
    }
    add(*total - owner, 1);
  }
  std::optional<size_t> result() const noexcept {
    if (!valid_)
      return {};
    return total_;
  }
};
} // namespace irred::detail
