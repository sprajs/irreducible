#include "irred/fits_reader.hpp"
#include "../../tests/verification/handcrafted_fits.hpp"
#include <stdexcept>
#define CHECK(...) do { if(!(__VA_ARGS__)) throw std::runtime_error("contract failed: " #__VA_ARGS__); } while(false)
#include <span>
using namespace irred::observations;
int main(){auto bytes=verifier_fits::table();auto view=std::as_bytes(std::span(bytes));auto d=decode_fits_length(view,std::string(64,'a'),3);
#ifdef IRRED_WITH_CFITSIO
 CHECK(d.status==DecodeStatus::ok);CHECK(d.input.values[0]==-2);CHECK(d.input.values[1]==6);CHECK(d.input.missing==std::vector<std::uint8_t>({0,0,1}));CHECK(d.input.quality==std::vector<std::uint64_t>({0,2,128}));CHECK(d.input.event_ids[0]==d.input.event_ids[1]);CHECK(d.input.measurement_ids[0]!=d.input.measurement_ids[1]);auto p=prepare(d.input,{3,0});CHECK(p.status()==Status::ok);CHECK(p.select(Selection::all).status==Status::missing_required_value);
 auto scaled=verifier_fits::table(true);CHECK(decode_fits_length(std::as_bytes(std::span(scaled)),std::string(64,'a'),3).status==DecodeStatus::invalid_format);
 auto wrong=bytes;std::string all(reinterpret_cast<const char*>(wrong.data()),wrong.size());const auto unit=all.find("TUNIT3");CHECK(unit!=std::string::npos);const auto metre=all.find("m",unit);wrong[metre]='s';CHECK(decode_fits_length(std::as_bytes(std::span(wrong)),std::string(64,'a'),3).status==DecodeStatus::invalid_format);
 CHECK(decode_fits_length(view,std::string(64,'a'),2).status==DecodeStatus::resource_limit);bytes.resize(5860);CHECK(decode_fits_length(std::as_bytes(std::span(bytes)),std::string(64,'a'),3).status==DecodeStatus::invalid_format);
#else
 CHECK(d.status==DecodeStatus::codec_unavailable);
#endif
}
