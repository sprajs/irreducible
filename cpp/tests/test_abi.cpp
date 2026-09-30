#include "cosmology/abi.h"
#include "cosmology/arithmetic.hpp"
#include <cassert>
#include <iostream>
#include <limits>
static cosmo_i64_buffer buffer(const int64_t* p,uint64_t n) {return {sizeof(cosmo_i64_buffer),1,1,0,p,n,n*8};}
int main() {
 const int64_t a[]={0,4,-7,9223372036854775806LL}, b[]={0,-2,3,1};
 auto x=buffer(a,4), y=buffer(b,4); cosmo_result* r=nullptr;
 assert(cosmo_add(&x,&y,0,&r)==0);
 const int64_t* p; uint64_t n; assert(cosmo_result_view(r,&p,&n)==0 && n==4);
 assert(p[0]==0 && p[1]==2 && p[2]==-4 && p[3]==std::numeric_limits<int64_t>::max());
 int64_t direct; for(int i=0;i<4;++i){assert(cosmology::checked_add(a[i],b[i],direct));assert(direct==p[i]);}
 std::cout<<"[0,2,-4,9223372036854775807]\n"; assert(cosmo_result_destroy(r)==0);
 x=buffer(nullptr,0); y=x; assert(cosmo_add(&x,&y,0,&r)==0); assert(cosmo_result_destroy(r)==0);
 x.abi_version=2; assert(cosmo_add(&x,&y,0,&r)==COSMO_ABI_MISMATCH && !r);
 x=y; x.length=1; assert(cosmo_add(&x,&y,0,&r)==COSMO_INVALID_INPUT && !r);
 x=y; assert(cosmo_add(&x,&y,1,&r)==COSMO_EXCEPTION && !r); assert(cosmo_add(&x,&y,2,&r)==COSMO_ALLOCATION_FAILURE && !r);
 const int64_t max=std::numeric_limits<int64_t>::max(),one=1; x=buffer(&max,1); y=buffer(&one,1); assert(cosmo_add(&x,&y,0,&r)==COSMO_OVERFLOW && !r);
}
