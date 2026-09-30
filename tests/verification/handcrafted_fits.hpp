#pragma once
// Independently assembled FITS v4.0 BINTABLE fixture, no CFITSIO writer ancestry.
// Specification: https://fits.gsfc.nasa.gov/fits_standard.html (V4.0, 2018-08-13).
// Expected LENGTH exact binary values (-6)*0.5+1=-2, 10*0.5+1=6;
// INT32_MIN is missing (not zero). Quality bits are raw 0,2,128. AUX
// 1,NaN,4 is unsupported scientific metadata retained in original bytes.
// No tolerance is needed for these exactly representable decoder values.
// All values are synthetic controls, never observations or fitted summaries.
#include <bit>
#include <cstdint>
#include <string>
#include <vector>
namespace verifier_fits {
using Bytes=std::vector<unsigned char>;
inline void card(Bytes& out,const std::string& key,const std::string& value="") {
 std::string line=key;line.resize(8,' ');if(!value.empty())line+="= "+value;line.resize(80,' ');out.insert(out.end(),line.begin(),line.end());
}
inline void header_end(Bytes& out){card(out,"END");while(out.size()%2880)out.push_back(' ');}
inline void bigendian(Bytes& out,uint64_t value,unsigned width){for(unsigned n=width;n>0;--n)out.push_back(static_cast<unsigned char>((value>>(8*(n-1)))&255));}
inline void text16(Bytes& out,std::string s){s.resize(16,' ');out.insert(out.end(),s.begin(),s.end());}
inline Bytes table(bool scaled_quality=false) {
 Bytes out;card(out,"SIMPLE","                    T");card(out,"BITPIX","                    8");card(out,"NAXIS","                    0");card(out,"EXTEND","                    T");header_end(out);
 card(out,"XTENSION","'BINTABLE'");card(out,"BITPIX","                    8");card(out,"NAXIS","                    2");card(out,"NAXIS1","                   48");card(out,"NAXIS2","                    3");card(out,"PCOUNT","                    0");card(out,"GCOUNT","                    1");card(out,"TFIELDS","                    5");card(out,"EXTNAME","'VERIFIER' ");
 card(out,"TTYPE1","'ROW_ID'  ");card(out,"TFORM1","'16A'     ");card(out,"TTYPE2","'EVENT_ID'");card(out,"TFORM2","'16A'     ");
 card(out,"TTYPE3","'LENGTH'  ");card(out,"TFORM3","'1J'      ");card(out,"TUNIT3","'m'       ");card(out,"TSCAL3","                  0.5");card(out,"TZERO3","                  1.0");card(out,"TNULL3","          -2147483648");
 card(out,"TTYPE4","'QUALITY' ");card(out,"TFORM4","'1J'      ");if(scaled_quality)card(out,"TSCAL4","                  2.0");card(out,"TTYPE5","'AUX'     ");card(out,"TFORM5","'1D'      ");header_end(out);
 const int32_t raw[3]={-6,10,INT32_MIN};const uint32_t quality[3]={0,2,128};const uint64_t auxiliary[3]={std::bit_cast<uint64_t>(1.0),UINT64_C(0x7ff8000000000000),std::bit_cast<uint64_t>(4.0)};
 for(unsigned i=0;i<3;++i){text16(out,i==0?"A:exposure1":i==1?"A:exposure2":"B:exposure1");text16(out,i<2?"event-A":"event-B");bigendian(out,static_cast<uint32_t>(raw[i]),4);bigendian(out,quality[i],4);bigendian(out,auxiliary[i],8);}
 while(out.size()%2880) out.push_back(0);
 return out;
}
}
