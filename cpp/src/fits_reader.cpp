#include "irred/fits_reader.hpp"
#include <limits>
#include <cstdio>
#ifdef COSMOLOGY_WITH_CFITSIO
#include <fitsio.h>
#endif
namespace irred::observations {
Decoded decode_fits_length(std::span<const std::byte> bytes,const std::string& hash,std::size_t maximum_rows){
 Decoded out;
#ifndef COSMOLOGY_WITH_CFITSIO
 (void)bytes;(void)hash;(void)maximum_rows;out.status=DecodeStatus::codec_unavailable;return out;
#else
 if(bytes.empty()||bytes.size()%2880!=0)return out;
 fitsfile* file=nullptr;int status=0;void* ptr=const_cast<std::byte*>(bytes.data());std::size_t size=bytes.size();
 fits_open_memfile(&file,"memory",READONLY,&ptr,&size,0,nullptr,&status);
 struct Close{fitsfile* f;~Close(){if(f){int x=0;fits_close_file(f,&x);}}} close{file};if(status)return out;
 char ext[]="VERIFIER";fits_movnam_hdu(file,BINARY_TBL,ext,0,&status);LONGLONG nr=0;fits_get_num_rowsll(file,&nr,&status);
 if(status||nr<=0)return out;
 if(static_cast<unsigned long long>(nr)>maximum_rows){out.status=DecodeStatus::resource_limit;return out;}
 int nc=0;fits_get_num_cols(file,&nc,&status);if(status||nc!=5)return out;
 const char* names[]={"ROW_ID","EVENT_ID","LENGTH","QUALITY","AUX"};const int types[]={TSTRING,TSTRING,TLONG,TLONG,TDOUBLE};
 for(int i=0;i<5;++i){char key[32],name[FLEN_VALUE]{};std::snprintf(key,sizeof key,"TTYPE%d",i+1);fits_read_key(file,TSTRING,key,name,nullptr,&status);int type=0;long repeat=0,width=0;fits_get_coltype(file,i+1,&type,&repeat,&width,&status);if(status||std::string(name)!=names[i]||type!=types[i]||repeat!=(i<2?16:1))return out;}
 char unit[FLEN_VALUE]{};fits_read_key(file,TSTRING,"TUNIT3",unit,nullptr,&status);if(status||std::string(unit)!="m")return out;
 double quality_scale=1,quality_zero=0;
 fits_read_key(file,TDOUBLE,"TSCAL4",&quality_scale,nullptr,&status);if(status==KEY_NO_EXIST)status=0;
 fits_read_key(file,TDOUBLE,"TZERO4",&quality_zero,nullptr,&status);if(status==KEY_NO_EXIST)status=0;
 if(status||quality_scale!=1||quality_zero!=0)return out; // preserve raw bit flags, never numeric rescaling
 const auto n=static_cast<std::size_t>(nr);auto& x=out.input;x.profile=Profile::fits_length_fixture_v1;x.role=Role::synthetic_control;x.unit=Unit::metre;x.calibration=Calibration::not_applicable;x.uncertainty=Uncertainty::none;x.table_sha256=hash;x.values.resize(n);x.missing.resize(n);x.quality.resize(n);
 std::vector<char> nulls(n);int any=0;fits_read_colnull(file,TDOUBLE,3,1,1,n,x.values.data(),nulls.data(),&any,&status);
 for(std::size_t i=0;i<n;++i)x.missing[i]=nulls[i]?1:0;
 std::vector<long long> flags(n);fits_read_colnull(file,TLONGLONG,4,1,1,n,flags.data(),nulls.data(),&any,&status);for(std::size_t i=0;i<n;++i){if(nulls[i]||flags[i]<0)return out;x.quality[i]=static_cast<std::uint64_t>(flags[i]);}
 std::vector<char*> pointers(n);std::vector<std::vector<char>> strings(n,std::vector<char>(17));for(std::size_t i=0;i<n;++i)pointers[i]=strings[i].data();
 for(int col=1;col<=2;++col){fits_read_col_str(file,col,1,1,n,nullptr,pointers.data(),&any,&status);if(status||any)return out;auto& ids=col==1?x.measurement_ids:x.event_ids;for(auto p:pointers)ids.emplace_back(p);}
 // Read optional auxiliary column too so corrupt/truncated final rows cannot
 // produce a successful partial table. Its raw bytes remain in source asset.
 std::vector<double> aux(n);fits_read_colnull(file,TDOUBLE,5,1,1,n,aux.data(),nulls.data(),&any,&status);if(status)return out;
 out.status=DecodeStatus::ok;return out;
#endif
}
}
