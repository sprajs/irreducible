#pragma once
#include <cstdint>
#include <vector>
struct cosmo_result { std::vector<int64_t> values; std::vector<double> f64_values; std::vector<uint32_t> statuses; uint32_t kind=1; std::vector<double> error_estimates; std::vector<uint64_t> evaluations; std::vector<uint8_t> selection_mask{}; std::vector<uint64_t> selection_indices{}; };
