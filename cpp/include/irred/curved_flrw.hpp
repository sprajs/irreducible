#pragma once
#include "irred/numerics.hpp"
#include <span>
#include <string_view>
#include <vector>
namespace irred::cosmology {
// Conserved dust, massless radiation and Lambda; Ok is derived from closure.
struct CurvedFLRWSpec { double h0_km_s_mpc=70, omega_m=.3, omega_r=0, omega_lambda=.7; };
struct CurvedFLRWPolicy {
  double relative_tolerance=1e-9, absolute_tolerance_mpc=1e-7;
  std::size_t maximum_points=4096, maximum_native_bytes=16*1024*1024;
  std::size_t maximum_callbacks_per_point=100000, maximum_total_callbacks=1000000;
  unsigned maximum_depth=40;
};
struct CurvedFLRWRow {
  numerics::Status status=numerics::Status::invalid_input;
  double redshift=0, e=0, h_km_s_mpc=0, radial_mpc=0;
  double transverse_mpc=0, angular_mpc=0, luminosity_mpc=0;
  double e_error_estimate=0, h_error_estimate_km_s_mpc=0;
  double radial_error_estimate_mpc=0, transverse_error_estimate_mpc=0;
  std::size_t callbacks=0;
};
struct CurvedFLRWBatch {
  numerics::Status status=numerics::Status::invalid_input;
  std::vector<CurvedFLRWRow> rows;
  std::size_t callbacks=0;
};
class CurvedFLRW {
public:
  CurvedFLRW()=default;
  explicit CurvedFLRW(CurvedFLRWSpec);
  CurvedFLRW(const CurvedFLRW&)=default;
  CurvedFLRW& operator=(const CurvedFLRW&)=default;
  CurvedFLRW(CurvedFLRW&&) noexcept;
  CurvedFLRW& operator=(CurvedFLRW&&) noexcept;
  numerics::Status status() const noexcept { return status_; }
  CurvedFLRWSpec specification() const noexcept { return spec_; }
  double omega_k() const noexcept { return omega_k_; }
  CurvedFLRWBatch evaluate(std::span<const double>, CurvedFLRWPolicy={}) const;
private:
  CurvedFLRWSpec spec_{};
  double omega_k_=0;
  numerics::Status status_=numerics::Status::invalid_input;
};
inline constexpr std::string_view curved_flrw_id="GR/nonflat-conserved-dust-radiation-lambda-FLRW/v1";
}
