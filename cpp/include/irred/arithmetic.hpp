#pragma once
#include <cstdint>
#include <limits>
namespace irred {
inline bool checked_add(std::int64_t x,std::int64_t y,std::int64_t& out) noexcept {
 if((y>0 && x>std::numeric_limits<std::int64_t>::max()-y)||(y<0 && x<std::numeric_limits<std::int64_t>::min()-y)) return false;
 out=x+y; return true;
}
}
