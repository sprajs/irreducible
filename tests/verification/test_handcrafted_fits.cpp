#include "handcrafted_fits.hpp"
#include <fitsio.h>
#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <chrono>
static void expect(bool v,const char* what){if(!v)throw std::runtime_error(what);}
struct Scratch{std::filesystem::path path;~Scratch(){std::error_code ignored;std::filesystem::remove_all(path,ignored);}};
struct Fits{fitsfile* file=nullptr;~Fits(){if(file){int status=0;fits_close_file(file,&status);}}};
static void write(const std::filesystem::path& p,const verifier_fits::Bytes& bytes){std::ofstream f(p,std::ios::binary);f.write(reinterpret_cast<const char*>(bytes.data()),bytes.size());expect(f.good(),"fixture write");}
int main(){
 Scratch tmp{std::filesystem::temp_directory_path()/("cosmology-fits-verifier-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()))};std::filesystem::create_directory(tmp.path);
 auto bytes=verifier_fits::table();expect(bytes.size()==8640,"independent 3-block fixture shape");const auto path=tmp.path/"fixture.fits";write(path,bytes);
 Fits f;int status=0;fits_open_file(&f.file,path.c_str(),READONLY,&status);expect(status==0,"open handcrafted FITS");char ext[]="VERIFIER";fits_movnam_hdu(f.file,BINARY_TBL,ext,0,&status);expect(status==0,"extension select");
 LONGLONG n=0;fits_get_num_rowsll(f.file,&n,&status);expect(status==0&&n==3,"row count");
 std::array<double,3> length{},aux{};std::array<char,3> masks{},auxmasks{};int any=0;
 fits_read_colnull(f.file,TDOUBLE,3,1,1,3,length.data(),masks.data(),&any,&status);expect(status==0&&any==1,"integer null decoded");expect(length[0]==-2&&length[1]==6,"big endian plus scaling/offset and signed value");expect(!masks[0]&&!masks[1]&&masks[2],"integer missingness preserved");
 fits_read_colnull(f.file,TDOUBLE,5,1,1,3,aux.data(),auxmasks.data(),&any,&status);expect(status==0&&any==1,"floating NaN decoded");expect(aux[0]==1&&aux[2]==4&&!auxmasks[0]&&auxmasks[1]&&!auxmasks[2],"floating null separated from valid values");
 std::array<int,3> quality{};int nullvalue=-1;fits_read_col(f.file,TINT,4,1,1,3,&nullvalue,quality.data(),&any,&status);expect(status==0&&quality==std::array<int,3>{0,2,128},"raw quality bits retained");
 std::array<std::array<char,32>,3> ids{},events{};std::array<char*,3> idp{ids[0].data(),ids[1].data(),ids[2].data()},ep{events[0].data(),events[1].data(),events[2].data()};
 fits_read_col_str(f.file,1,1,1,3,nullptr,idp.data(),&any,&status);fits_read_col_str(f.file,2,1,1,3,nullptr,ep.data(),&any,&status);expect(status==0,"IDs decoded");expect(std::string(idp[0])=="A:exposure1"&&std::string(idp[1])=="A:exposure2"&&std::string(ep[0])==std::string(ep[1]),"repeated event distinct exposure identities retained");
 char unit[FLEN_VALUE]{};fits_read_key(f.file,TSTRING,"TUNIT3",unit,nullptr,&status);expect(status==0&&std::string(unit)=="m","unit metadata retained");
 bytes.resize(5760+100);const auto truncated=tmp.path/"truncated.fits";write(truncated,bytes);Fits bad;status=0;fits_open_file(&bad.file,truncated.c_str(),READONLY,&status);if(!status){fits_movnam_hdu(bad.file,BINARY_TBL,ext,0,&status);if(!status)fits_read_colnull(bad.file,TDOUBLE,5,1,1,3,aux.data(),auxmasks.data(),&any,&status);}expect(status!=0,"truncated row rejected");
 std::cout<<"{\"suite\":\"independent_handcrafted_FITS_codec\",\"passed\":true,\"scope\":\"synthetic BINTABLE primitive decoding only; no production adapter qualification\"}\n";
}
