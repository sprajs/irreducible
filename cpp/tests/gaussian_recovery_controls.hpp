// Frozen root Fraction/CDF facts independently matched by worker Fraction route.
// Original exact Fraction covariance and precision routes; refined CDF facts.
#pragma once
#include <array>
namespace recovery_peer_facts {
inline constexpr std::array<long double,4> generic_K{0.49101796407185628742514970059880239521L,0.19161676646706586826347305389221556886L,-0.10778443113772455089820359281437125749L,0.29940119760479041916167664670658682635L};
inline constexpr std::array<long double,4> generic_V{0.53892215568862275449101796407185628743L,-0.032934131736526946107784431137724550898L,-0.032934131736526946107784431137724550898L,0.60479041916167664670658682634730538922L};
inline constexpr std::array<long double,4> generic_fixed_sampling{0.36157624870020438165584997669331994693L,0.093406002366524436157624870020438165585L,0.093406002366524436157624870020438165585L,0.17476424396715550934060023665244361576L};
inline constexpr std::array<long double,2> generic_fixed_bias{-0.50898203592814371257485029940119760479L,0.89221556886227544910179640718562874251L};
inline constexpr std::array<long double,4> generic_W{3.3263473053892215568862275449101796407L,-0.57859281437125748502994011976047904192L,-0.57859281437125748502994011976047904192L,1.9595808383233532934131736526946107784L};
inline constexpr std::array<long double,4> generic_fixed_future_sampling{1.9342572340349241636487504033848470723L,0.23045376313241779913227437340887088099L,0.23045376313241779913227437340887088099L,1.0995284879343110186812004733048872315L};
inline constexpr std::array<long double,2> generic_fixed_future_bias{-1.2754491017964071856287425149700598802L,1.4011976047904191616766467065868263473L};
inline constexpr std::array<long double,2> generic_posterior_coverage{0.93839369156036818171072300498831060163L,0.93471047931290815431783941583919071359L};
inline constexpr std::array<long double,2> generic_future_coverage{0.9506078251109300530201510000049473759L,0.89973424842486257198010818121235397739L};
inline constexpr std::array<long double,88> ladder_K{0.13026748185977419840173802269637435407L,0.003859528857733266744517156586767554373L,0.11758165790328821697808714279804836081L,0.095055928687555992345975686104799099161L,-0.035558129709941769293821337726489586676L,-0.013649984034239638011461438778240164248L,0.084876478102577389708847410227356886228L,-0.0036832029632069338434636160736027903452L,-0.031149781000626570358884168592699319749L,-0.042317162830173847938532283169725396979L,-0.030769792535807126890722887930542461897L,0.011941528587010482225334687813702141138L,0.13871914415258959485769685146118022399L,-0.047330615450503452028631398245476349602L,-0.019608858957643791035363850924247747412L,0.11112696062403654659310000095685727915L,0.089141726955950449633064792132624128005L,-0.0082332875865931462899572913913789519673L,0.1072703139688018382633566985607065074L,-0.014962430472249435315908666690902749055L,0.0073625010937218629800088089097432703477L,0.018597514755405964342437717301974687944L,-0.027410534995460282756659551208995726672L,-0.035403447698972386229021179725313972442L,0.11478609300575932232376651028673501975L,0.11207922811049419919397052158641210576L,0.10388540916460494509431650058450247771L,0.1052739402674323197188828217377639988L,-0.013856017747584973997919566569784678176L,-0.026126977936599159470189746626376103659L,0.0095811915556798407712679482838568097167L,0.0060938427572271307434726008668876943466L,0.0011199574116123452064011121614830396197L,0.00578222560454663897642665573502784183L,-0.0024831799991165383526042350368355370736L,-0.1146734204271223592785730234935148496L,0.084325832686327844795706474890638899521L,-0.05622611166165748499083979422464522098L,0.065861944148311003749535641417105665107L,-0.017875669729891762368699758988332803414L,0.0066176622994606377141287945966814452151L,0.0058586384671994823358826941093338253557L,-0.0065481018980462048304105299041766272971L,-0.0074309556783510976716623259175816295388L,0.01335333281083400243781425245621305517L,0.0039216031939823397844671774692093090388L,-0.030073067832160502457329858786377467245L,0.092095963492284738786906547221351385914L,0.050448323198909613038684778386934371945L,-0.093373562974157443512075293062175802114L,-0.010902852378204603577804366057547594069L,0.020122265344095275535955891268032110352L,0.0013369306280793676671461748135337045245L,-0.02963609743960916214689890722732724168L,-0.026855145453957380261050244072444560469L,-0.019475723019414687729181913334838361333L,-0.023462439284833345089439021645252485026L,0.00044425507270096444270191873218764557835L,-0.031102403091101462940759377984869343259L,-0.018551900751764771552949334914999315207L,0.0035076118794840031492946057735838490908L,0.10645252117576797682590984963278849664L,0.13703859005557337212711787247269782362L,0.10263651658016632719490280017074629245L,0.069886618102339474011233855250455438258L,-0.023156550688829661667360933117969727309L,-0.0031279529132490108069714561037575845432L,-0.0025654622347714912768074519307318115791L,-0.0097764333007210446622693995223271562398L,0.018720443951000654373013730081649367019L,0.0035108070614614027155825686607727635965L,-0.0168845358976651859758969982913778544L,0.062579458495287614159801057308238711848L,0.026593514947146014859670141313853847319L,-0.14021379826842037410104101820239609878L,-0.11133315046409774914139386157167892911L,0.023359626952563203857674719123807039838L,-0.025685721354601022048190812172031978073L,0.016915623275526348083299859092988756223L,-0.067621679372212745913179317712778519427L,0.07257114602538791768377902480413447946L,0.04339869070129648231217414323307428876L,-0.057998541215908173227565989660751281405L,0.068647696896514877972893108268672476964L,-0.064385422255238393242561765157412096086L,-0.017703389868880922180919785862582581535L,0.12193510995761227491339157808793042704L,0.11534614820625263424350378568231595691L};
inline constexpr std::array<long double,64> ladder_V{0.0020996045261701254641350381027809230877L,-0.000070920144483550548625359309755564655907L,-0.00037594581253601088736897981683253191498L,-0.00023984967440821654436271666281026274841L,0.00019968490382755350094367100382267341396L,-0.0002926315474453673922667521940637457733L,-0.000024529683567798450957403251849356163747L,-0.00034955272548635351615910524249091827415L,-0.000070920144483550548625359309755564655907L,0.0021719191364144343289405982609406637908L,-0.00052968225701294285436789469465599757472L,-0.000023809082835682377736771381024913845L,0.000048794743416751965323320210217478365358L,-0.0003483111421101855550058127010780859677L,-0.00003855224219531714801527873317709979569L,0.00028615365902299628356150437598363170813L,-0.00037594581253601088736897981683253191498L,-0.00052968225701294285436789469465599757472L,0.0021560354041037377549736573893200991277L,0.000050924084660461196502202219533152726577L,-0.000049639605321854272757768365771120856457L,0.00014194620067355277480146921165647382442L,0.000014237376217791496558332793016029229625L,-0.0000059972787270577866105438610288607506285L,-0.00023984967440821654436271666281026274841L,-0.000023809082835682377736771381024913845L,0.000050924084660461196502202219533152726577L,0.0034283341760618972859272395512008352026L,-0.00012312094156811550201025647368715793872L,0.000076651017352893658471506771079775656601L,-0.000032674425827772077789588941424887538635L,-0.00013109928712474943514239363388785438332L,0.00019968490382755350094367100382267341396L,0.000048794743416751965323320210217478365358L,-0.000049639605321854272757768365771120856457L,-0.00012312094156811550201025647368715793872L,0.0034903312731698948906511831745595077947L,0.000049569675481083634732045840160091686038L,-0.000085583007379094072068575235513546827768L,-0.00040713134122886197276993062589302916895L,-0.0002926315474453673922667521940637457733L,-0.0003483111421101855550058127010780859677L,0.00014194620067355277480146921165647382442L,0.000076651017352893658471506771079775656601L,0.000049569675481083634732045840160091686038L,0.0023177732953297054937241081745443430226L,0.0006373395514228673619281272859811010389L,-0.00038011057622829892556918659209726110004L,-0.000024529683567798450957403251849356163747L,-0.00003855224219531714801527873317709979569L,0.000014237376217791496558332793016029229625L,-0.000032674425827772077789588941424887538635L,-0.000085583007379094072068575235513546827768L,0.0006373395514228673619281272859811010389L,0.0029205785210676262747145633444013499018L,0.00036346106591081265709132978306890023724L,-0.00034955272548635351615910524249091827415L,0.00028615365902299628356150437598363170813L,-0.0000059972787270577866105438610288607506285L,-0.00013109928712474943514239363388785438332L,-0.00040713134122886197276993062589302916895L,-0.00038011057622829892556918659209726110004L,0.00036346106591081265709132978306890023724L,0.0017804365203798003152948025736305044346L};
inline constexpr std::array<long double,64> ladder_fixed_sampling{0.00083294758551584902173417021505253003319L,-0.00022161809444699656459279297567176299053L,0.00021842504048116573476597111843566501688L,-0.00011965895841886424649114531919482199161L,0.00011672867431481056477954080400400696378L,-0.000028297359103945992036516956895358263411L,0.00022013243142662465503006976261383503448L,-0.00010686027619802350429186552962309970386L,-0.00022161809444699656459279297567176299053L,0.00084972502230173357512672426636884140064L,0.000089530756609627979974330006814521206933L,0.000057891201146117891029734686737287787248L,-0.000017237297594895427471676031372344181123L,0.00012565848990287567784071374449693037255L,0.00004190708369734077309524991033429422666L,-0.0000044653296786867157496953239901889969625L,0.00021842504048116573476597111843566501688L,0.000089530756609627979974330006814521206933L,0.00078827330042200487179040798037759969905L,-0.000032489392563736558814114825527510608448L,0.000019722807114846106159618727716812256174L,-0.00011537143969584089537073812466357394179L,-0.000060403395148710467555954657822267000676L,0.000002304228227432884537476470660626426895L,-0.00011965895841886424649114531919482199161L,0.000057891201146117891029734686737287787248L,-0.000032489392563736558814114825527510608448L,0.00044527444471349391772089781279556518987L,0.00004299952960379502870618744061627285303L,-0.000027095082247695641277781041530310452441L,0.0000033132409353337607288676081102921222327L,0.000070357841834197674969910261157351916428L,0.00011672867431481056477954080400400696378L,-0.000017237297594895427471676031372344181123L,0.000019722807114846106159618727716812256174L,0.00004299952960379502870618744061627285303L,0.00035864807571973587064205579724543417547L,-0.000065433410829608037197973269174996920312L,0.00009498702262012770850156418955150106154L,0.00011442865073852800243967133125227011995L,-0.000028297359103945992036516956895358263411L,0.00012565848990287567784071374449693037255L,-0.00011537143969584089537073812466357394179L,-0.000027095082247695641277781041530310452441L,-0.000065433410829608037197973269174996920312L,0.00075455636147510820570249949729654533366L,-0.0002035618669662727963828809686000018378L,-0.000011806437589558314193223478165226011522L,0.00022013243142662465503006976261383503448L,0.00004190708369734077309524991033429422666L,-0.000060403395148710467555954657822267000676L,0.0000033132409353337607288676081102921222327L,0.00009498702262012770850156418955150106154L,-0.0002035618669662727963828809686000018378L,0.00059338357221086506115082653002520545808L,-0.000041749784263526882004205005957487659346L,-0.00010686027619802350429186552962309970386L,-0.0000044653296786867157496953239901889969625L,0.000002304228227432884537476470660626426895L,0.000070357841834197674969910261157351916428L,0.00011442865073852800243967133125227011995L,-0.000011806437589558314193223478165226011522L,-0.000041749784263526882004205005957487659346L,0.0008315228194762254142582979650960777171L};
inline constexpr std::array<long double,8> ladder_fixed_bias{0.018366735957700104574354112940606179182L,0.0018999823446318930839749644452898696769L,-0.0039187585782267428786851372876861260158L,0.00017236585771168062363195894871070048536L,0.0070747917028824486199381831034483842889L,-0.043130783670626450301534163855300432676L,-0.18711326281687046918939128005648124303L,-0.026057930022182838183117948056336961376L};
inline constexpr std::array<long double,4> ladder_W{0.019481336101254921881894269791914932544L,0.0037552736201849119085437118246502389722L,0.0037552736201849119085437118246502389722L,0.020451373293056567347944968831901407585L};
inline constexpr std::array<long double,4> ladder_fixed_future_sampling{0.017880363910892911689397008986919251304L,0.0038381681377841195521429405417198337638L,0.0038381681377841195521429405417198337638L,0.018869272773584848506855852861576690655L};
inline constexpr std::array<long double,2> ladder_fixed_future_bias{-0.0032308017585304749150025371377913436714L,-0.12443903162960689025051865515892808932L};
inline constexpr std::array<long double,8> ladder_posterior_coverage{0.99325628312632365722329365038395003548L,0.99823342316894392707213382794838496443L,0.99867736501594945512405233967369997867L,0.99999994620160111120407191644712020065L,0.99999999524299660653866629203917450796L,0.96890450303013704156583709959443092305L,0.0004294447488094294097454631417419916822L,0.97517196336222120809171595541413533966L};
inline constexpr std::array<long double,2> ladder_future_coverage{0.95916800495723381414035337281055109759L,0.87011666500850589982006866114540126812L};
}
#include "irred/gaussian_simulation.hpp"
#include "irred/gaussian_predictive.hpp"
#include <chrono>
#include <bit>
#include <limits>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
namespace recovery_controls {
namespace s=irred::statistics;namespace n=irred::numerics;using W=long double;
inline constexpr size_t attempts=32768,chunk=1024;
inline constexpr std::array<std::uint64_t,2> seeds{0x123456789abcdef0ULL,0xfedcba9876543210ULL};
inline constexpr W z95=1.95996398454005423552L;
inline void need(bool b,const char*w){if(!b)throw std::runtime_error(w);}
template<class A>std::vector<double> doubles(const A&a){return {a.begin(),a.end()};}
template<class A>std::vector<W> wide(const A&a){return {a.begin(),a.end()};}
inline void near(W a,W b){need(std::abs(a-b)<=2e-12L*(1+std::abs(b)),"frozen analytic allocation");}
inline s::Metadata metadata(std::vector<std::string> ids,std::string measure="product d(mag)"){
 s::Metadata m;m.ordered_ids=std::move(ids);m.measure=std::move(measure);m.source_semantics="synthetic controls";
 m.table_identity="frozen synthetic generating object";m.uncertainty_identity="explicit conditional SPD noise";
 m.ordering_provenance="exact frozen source axes";m.calibration_provenance="one shared delta";
 m.dependence_provenance="explicit independent original prior/training/future noise";return m;}
inline s::Gaussian gaussian(const std::vector<double>&C,s::Metadata m){return s::prepare_gaussian(C,s::MatrixKind::covariance,std::move(m),1000000,1e-10,n::Arithmetic::longdouble_cpu_v1);}
inline s::GeneratingMean generating(const s::Gaussian&g,std::vector<double>value,std::vector<std::string>units,std::string role){
 return {std::move(value),g.metadata().ordered_ids,std::move(units),role+" mean",g.metadata().measure,role+" original ideal Gaussian / emitted halfbin approximation"};}
inline s::GaussianSimulationBatch draw(const s::Gaussian&g,const s::GeneratingMean&m,std::uint64_t seed,std::uint64_t stream,size_t start,size_t count,bool record=true){
 auto admission=[&](std::string_view stage,n::Status status){std::cout<<"{\"kind\":\"batch_admission_refusal\",\"seed\":\""<<seed<<"\",\"stream_base\":\""<<stream<<"\",\"start\":"<<start<<",\"planned_vectors\":"<<count<<",\"attempted_vectors\":0,\"dimension\":"<<g.metadata().ordered_ids.size()<<",\"stage\":\""<<stage<<"\",\"numerical_status\":"<<int(status)<<"}\n";};
 s::GaussianSimulationBatch b;
 try {std::vector<irred::random::Address> addresses;addresses.reserve(count);
  if(count&&count-1>SIZE_MAX-start){admission("address_overflow",n::Status::invalid_input);throw std::runtime_error("recorded sample address overflow");}
  for(size_t i=0;i<count;++i)addresses.push_back({stream,start+i});
  s::GaussianSimulationPolicy p;p.outputs=7;b=s::GaussianSimulation::simulate(g,m,addresses,seed,p);
 }catch(const std::bad_alloc&){admission("address_or_batch_allocation",n::Status::work_limit);throw;}
  catch(const std::length_error&){admission("address_or_batch_length",n::Status::work_limit);throw;}
 if(b.status!=n::Status::ok||b.rows.size()!=count){admission("generation_admission",b.status);throw std::runtime_error("recorded whole batch admission refusal");}
 if(record)std::cout<<"{\"kind\":\"generation_stage\",\"seed\":\""<<seed<<"\",\"stream_base\":\""<<stream<<"\",\"start\":"<<start<<",\"attempted_vectors\":"<<count<<",\"dimension\":"<<b.dimension<<",\"status\":"<<int(b.status)<<",\"generating_identity\":\""<<m.generating_law_identity<<"\"}\n";
 return b;}
inline void generation_benchmark(const s::Gaussian&g,const s::GeneratingMean&m){
 std::vector<s::GaussianSimulationBatch>single;single.reserve(1024);
 const auto begin=std::chrono::steady_clock::now();auto coarse=draw(g,m,seeds[0],0x310000,0,1024,false);
 const auto mid=std::chrono::steady_clock::now();
 for(size_t i=0;i<1024;++i)single.push_back(draw(g,m,seeds[0],0x310000,i,1,false));
 const auto end=std::chrono::steady_clock::now();
 for(size_t i=0;i<1024;++i){const auto &one=single[i];
  need(one.rows[0].status==coarse.rows[i].status,"benchmark matched row status");
  for(size_t j=0;j<coarse.dimension;++j)need(one.values[j]==coarse.values[i*coarse.dimension+j]&&
   one.absolute_error_estimates[j]==coarse.absolute_error_estimates[i*coarse.dimension+j]&&one.words[j]==coarse.words[i*coarse.dimension+j],"benchmark matched words/arithmetic/output quality");}
 std::cout<<"{\"kind\":\"retained_generation_comparison\",\"vectors\":1024,\"outputs\":7,\"threads\":1,\"coarse_seconds\":"<<std::chrono::duration<W>(mid-begin).count()
 <<",\"one_vector_seconds\":"<<std::chrono::duration<W>(end-mid).count()<<",\"coarse_payload_bound\":"<<*s::GaussianSimulation::payload_bound(g,m,1024,7)
 <<",\"one_vector_payload_bound\":"<<*s::GaussianSimulation::payload_bound(g,m,1,7)<<",\"comparison_retained_outputs_payload_envelope\":"<<1024*(*s::GaussianSimulation::payload_bound(g,m,1,7)-*g.retained_payload_bound())+*s::GaussianSimulation::payload_bound(g,m,1024,7)-*g.retained_payload_bound()<<",\"borrowed_factor_retained_bytes\":"<<*g.retained_payload_bound()<<",\"verification_timed\":false,\"word_value_error_bits_equal\":true}\n";
}
inline std::vector<double> observation(const std::vector<double>&design,const std::vector<double>&offset,
 const std::vector<double>&beta,std::span<const double>noise){std::vector<double>y(offset.size());
 for(size_t i=0;i<y.size();++i){W v=offset[i]+W(noise[i]);for(size_t j=0;j<beta.size();++j)v+=W(design[i*beta.size()+j])*beta[j];y[i]=double(v);
 }return y;}
inline n::Status observation_status(std::span<const double>y){for(double x:y){if(!std::isfinite(x))return n::Status::overflow;if(x!=0&&!std::isnormal(x))return n::Status::outside_domain;}return n::Status::ok;}
inline void vector_json(std::span<const double>v){std::cout<<'[';for(size_t j=0;j<v.size();++j){if(j)std::cout<<',';
 if(std::isfinite(v[j]))std::cout<<v[j];else std::cout<<"{\"ieee754_hex\":\""<<std::hex<<std::bit_cast<std::uint64_t>(v[j])<<std::dec<<"\"}";}std::cout<<']';}
inline void failed(std::string_view name,std::uint64_t seed,size_t sample,std::string_view stage,
 std::span<const double>beta,std::span<const double>y,std::span<const double>fy,const s::GaussianSimulationBatch&prior,
 const s::GaussianSimulationBatch&train,const s::GaussianSimulationBatch&future,size_t index,
 s::DensityStatus status=s::DensityStatus::numerical_failure,n::Status numerical=n::Status::invalid_input){
 std::cout<<std::setprecision(17)<<"{\"kind\":\"failed_attempt\",\"campaign\":\""<<name<<"\",\"seed\":\""<<seed<<"\",\"sample\":"<<sample<<",\"stage\":\""<<stage<<"\",\"beta\":";
 vector_json(beta);std::cout<<",\"training\":";vector_json(y);std::cout<<",\"future\":";vector_json(fy);std::cout<<",\"status\":"<<int(status)<<",\"numerical_status\":"<<int(numerical)<<",\"attempted_vectors\":1,\"random_blocks\":[";bool first=true;
 for(const auto*b:{&prior,&train,&future})if(index<b->rows.size())for(size_t j=0;j<b->dimension;++j){if(!first)std::cout<<',';first=false;
  const auto a=b->rows[index].address;std::cout<<"{\"stream\":\""<<a.stream+j<<"\",\"sample\":"<<a.sample<<",\"numerical_status\":"<<int(b->rows[index].status)<<",\"words\":[";
  for(size_t k=0;k<4;++k){if(k)std::cout<<',';std::cout<<b->words[index*b->dimension+j][k];}std::cout<<"]}";}
 std::cout<<"]}\n";
}
struct Moments {size_t count=0;std::vector<W>sum,cross;
 explicit Moments(size_t d):sum(d),cross(d*d){}
 void add(std::span<const W>v){++count;for(size_t i=0;i<v.size();++i){sum[i]+=v[i];for(size_t j=0;j<v.size();++j)cross[i*v.size()+j]+=v[i]*v[j];}}
 void report(std::string_view campaign,std::string_view ensemble,std::string_view role)const{
  std::cout<<std::setprecision(17)<<"{\"kind\":\"empirical_full_moments\",\"campaign\":\""<<campaign<<"\",\"ensemble\":\""<<ensemble<<"\",\"role\":\""<<role<<"\",\"admitted\":"<<count<<",\"selected_diagnostic\":"<<(count!=attempts?"true":"false")<<",\"mean\":[";
  for(size_t i=0;i<sum.size();++i){if(i)std::cout<<',';std::cout<<sum[i]/count;}std::cout<<"],\"full_covariance\":[";
  for(size_t i=0;i<sum.size();++i){for(size_t j=0;j<sum.size();++j){if(i||j)std::cout<<',';std::cout<<(cross[i*sum.size()+j]-sum[i]*sum[j]/count)/(count-1);}}
  std::cout<<"]}\n";
 }
 void check(std::span<const W>mean,std::span<const W>C)const{
  need(count==attempts,"refused attempts prevent empirical moment pass");const size_t d=sum.size();
  for(size_t i=0;i<d;++i){const W observed=sum[i]/count;
   const W allowance=6*std::sqrt(C[i*d+i]/count)+1e-6L*(1+std::abs(mean[i]))+1e-10L;
   need(std::abs(observed-mean[i])<=allowance,"empirical mean discrepancy");
   for(size_t j=0;j<d;++j){const W covariance=(cross[i*d+j]-sum[i]*sum[j]/count)/(count-1);
    const W allocation=6*std::sqrt((C[i*d+i]*C[j*d+j]+C[i*d+j]*C[i*d+j])/(count-1))+1e-6L*(1+std::abs(C[i*d+j]))+1e-10L;
    need(std::abs(covariance-C[i*d+j])<=allocation,"full empirical covariance discrepancy");}}
 }
};
struct Reference {std::vector<W>K,V,T,bias,Wcov,Tfuture,bfuture,posterior_probability,future_probability;};
struct FutureEvaluation {s::DensityStatus status=s::DensityStatus::invalid_input; n::Status numerical_status=n::Status::invalid_input;
 std::vector<double>mean;double quadratic=0;size_t work=0;};
inline void coverage(size_t hits,W probability){const W actual=W(hits)/attempts;
 const W allowance=6*std::sqrt(probability*(1-probability)/attempts)+6.L/attempts+5e-4L+1e-10L;
 need(std::abs(actual-probability)<=allowance,"empirical ideal-model coverage discrepancy");}
inline W chi_square_even(unsigned dimensions,W q){need(dimensions%2==0,"even named chi-square domain");W sum=1,term=1;
 for(unsigned i=1;i<dimensions/2;++i){term*=q/(2*i);sum+=term;}return 1-std::exp(-q/2)*sum;}
// No fixed-truth joint ellipsoid target: that law is biased and anisotropic.
struct ProjectionEvaluation {
 s::DensityStatus status=s::DensityStatus::finite;
 n::Status numerical_status=n::Status::ok;
 bool check=true;
};
struct ReservedRange {std::uint64_t seed,stream;size_t start,count,dimension;};
inline std::vector<ReservedRange> reserved_ranges;
inline void reserve_range(std::uint64_t seed,std::uint64_t stream,size_t start,size_t count,size_t d){
 need(d&&d-1<=UINT64_MAX-stream&&count&&count-1<=SIZE_MAX-start,"campaign range overflow");
 for(const auto&r:reserved_ranges)if(seed==r.seed){
  const bool samples_overlap=start<=r.start+r.count-1&&r.start<=start+count-1;
  const bool streams_overlap=stream<=r.stream+r.dimension-1&&r.stream<=stream+d-1;
  need(!samples_overlap||!streams_overlap,"campaign global address overlap");}
 reserved_ranges.push_back({seed,stream,start,count,d});
}
inline void strings_json(const std::vector<std::string>&v){std::cout<<'[';for(size_t i=0;i<v.size();++i){if(i)std::cout<<',';std::cout<<'"'<<v[i]<<'"';}std::cout<<']';}
inline void campaign_identity(std::string_view name,bool fixed,std::uint64_t stream,size_t p,size_t nt,size_t k,
 const s::GeneratingMean&gp,const s::GeneratingMean&gt,const s::GeneratingMean&gf,
 const s::Gaussian&prior,const s::Gaussian&training,const s::Gaussian&future){
 std::cout<<"\"campaign\":\""<<name<<"\",\"ensemble\":\""<<(fixed?"fixed truth":"prior predictive")<<"\",\"seeds\":[\""<<seeds[0]<<"\",\""<<seeds[1]<<"\"],\"sample_first\":0,\"sample_last_per_seed\":"<<attempts/2-1<<",\"planned_count\":"<<attempts<<",\"chunk\":"<<chunk<<",\"stream_family\":\""<<stream<<"\",\"roles\":[";
 bool first=true;for(size_t role=0;role<3;++role){if(fixed&&role==0)continue;if(!first)std::cout<<',';first=false;
 const auto&g=role==0?gp:role==1?gt:gf;const auto&source=(role==0?prior:role==1?training:future).metadata();const auto d=role==0?p:role==1?nt:k;const auto base=stream+role*0x100;
 std::cout<<"{\"role\":\""<<(role==0?"parameter prior":role==1?"training noise":"future noise")<<"\",\"stream_first\":\""<<base<<"\",\"stream_last\":\""<<base+d-1<<"\",\"dimension\":"<<d<<",\"generating_identity\":\""<<g.generating_law_identity<<"\",\"mean_identity\":\""<<g.identity<<"\",\"measure\":\""<<g.coordinate_measure<<"\",\"ordered_ids\":";strings_json(g.ordered_ids);std::cout<<",\"source_table_identity\":\""<<source.table_identity<<"\",\"source_uncertainty_identity\":\""<<source.uncertainty_identity<<"\",\"source_arithmetic\":\""<<source.arithmetic_id<<"\",\"source_semantics\":\""<<source.source_semantics<<"\"";std::cout<<'}';}
 std::cout<<']';
}
template<class Condition,class PosteriorDensity,class Predict,class Projection>
void campaign(std::string_view name,const s::Gaussian&source,const s::Gaussian&prior_noise,const s::Gaussian&future_noise,
 const std::vector<double>&X,const std::vector<double>&A,const std::vector<double>&o,const std::vector<double>&of,
 const std::vector<double>&m,const std::vector<double>&truth,const std::vector<std::string>&units,
 const Reference&ref,std::uint64_t family,Condition condition,PosteriorDensity posterior_density,Predict predict,Projection projection,bool h0_requested=false){
 const size_t p=m.size(),nt=o.size(),k=of.size();need(p<256&&nt<256&&k<256,"global disjoint role ranges");
 // Fixed before draws; verify K/V via separate exact-input finite differences.
 auto y0=observation(X,o,truth,std::vector<double>(nt));auto c0=condition(y0);need(c0.status==s::DensityStatus::finite,"analytic truth mean");
 for(size_t j=0;j<p;++j)near(c0.value[j]-truth[j],ref.bias[j]);
 for(size_t i=0;i<nt;++i){auto y=y0;y[i]+=.125;auto c=condition(y);need(c.status==s::DensityStatus::finite,"analytic K direction");
  for(size_t j=0;j<p;++j)near((c.value[j]-c0.value[j])*8,ref.K[j*nt+i]);}
 const auto generating_prior=generating(prior_noise,m,units,"proper parameter prior");
 const auto generating_train=generating(source,std::vector<double>(nt),std::vector<std::string>(nt,"mag"),"conditional training noise");
 const auto generating_future=generating(future_noise,std::vector<double>(k),std::vector<std::string>(k,"mag"),"independent conditional future noise");
 for(bool fixed:{false,true}){
  const auto start_time=std::chrono::steady_clock::now();W predictive_seconds=0;size_t refused=0,joint_hits=0,posterior_joint_hits=0,declared_work=0,actually_attempted=0;
  std::vector<size_t>posterior_hits(p),future_hits(k);Moments posterior_moments(p),future_moments(k);
  const auto stream=family+(fixed?0x10000:0);
  std::cout<<"{\"kind\":\"campaign_address_identity\",";campaign_identity(name,fixed,stream,p,nt,k,generating_prior,generating_train,generating_future,prior_noise,source,future_noise);std::cout<<"}\n";
  for(auto seed:seeds){if(!fixed)reserve_range(seed,stream,0,attempts/2,p);reserve_range(seed,stream+0x100,0,attempts/2,nt);reserve_range(seed,stream+0x200,0,attempts/2,k);}
  try {
  for(auto seed:seeds)for(size_t start=0;start<attempts/2;start+=chunk){
   s::GaussianSimulationBatch bp;if(!fixed)bp=draw(prior_noise,generating_prior,seed,stream,start,chunk);
   auto bt=draw(source,generating_train,seed,stream+0x100,start,chunk),bf=draw(future_noise,generating_future,seed,stream+0x200,start,chunk);
   declared_work+=bp.work_units+bt.work_units+bf.work_units;
   for(size_t i=0;i<chunk;++i){++actually_attempted;std::vector<double>beta,y,fy;std::string_view failure="parameter_copy";
    s::DensityStatus failed_status=s::DensityStatus::numerical_failure;n::Status failed_numerical=n::Status::work_limit;
    try {
    if(fixed)beta=truth;else if(bp.rows[i].status==n::Status::ok)beta.assign(bp.values.begin()+i*p,bp.values.begin()+(i+1)*p);
    failure="generation";
    failed_numerical=!fixed&&bp.rows[i].status!=n::Status::ok?bp.rows[i].status:bt.rows[i].status!=n::Status::ok?bt.rows[i].status:bf.rows[i].status;
    if(failed_numerical==n::Status::ok){
     failure="training_observation";y=observation(X,o,beta,std::span(bt.values.data()+i*nt,nt));
     failure="future_observation";fy=observation(A,of,beta,std::span(bf.values.data()+i*k,k));
     const auto ys=observation_status(y),fs=observation_status(fy);
     if(ys!=n::Status::ok||fs!=n::Status::ok){failure=ys!=n::Status::ok?"training_observation":"future_observation";failed_numerical=ys!=n::Status::ok?ys:fs;}
     else {failure="conditioning";auto c=condition(y);failed_status=c.status;failed_numerical=c.numerical_status;
      if(c.status==s::DensityStatus::finite){failure="posterior_density";auto pd=posterior_density(y,beta);failed_status=pd.density.status;failed_numerical=pd.density.numerical_status;
       if(pd.density.status==s::DensityStatus::finite){failure="future_prediction_density";const auto begin=std::chrono::steady_clock::now();auto f=predict(y,fy);predictive_seconds+=std::chrono::duration<W>(std::chrono::steady_clock::now()-begin).count();declared_work+=f.work;failed_status=f.status;failed_numerical=f.numerical_status;
        if(f.status==s::DensityStatus::finite){failure="H0_projection";auto h=projection(beta,c,y);failed_status=h.status;failed_numerical=h.numerical_status;
         if(h.status==s::DensityStatus::finite&&h.check){
          std::vector<W>e(p),ef(k);std::vector<bool>ph(p),fh(k);
          for(size_t j=0;j<p;++j){e[j]=W(c.value[j])-beta[j];ph[j]=std::abs(e[j])<=z95*std::sqrt(ref.V[j*p+j]);posterior_hits[j]+=ph[j];}
          for(size_t j=0;j<k;++j){ef[j]=W(fy[j])-f.mean[j];fh[j]=std::abs(ef[j])<=z95*std::sqrt(ref.Wcov[j*k+j]);future_hits[j]+=fh[j];}
          posterior_moments.add(e);future_moments.add(ef);if(!fixed){joint_hits+=f.quadratic<=6;posterior_joint_hits+=pd.quadratic<=16;}
          std::cout<<std::setprecision(17)<<"{\"kind\":\"successful_attempt\",\"campaign\":\""<<name<<"\",\"ensemble\":\""<<(fixed?"fixed truth":"prior predictive")<<"\",\"seed\":\""<<seed<<"\",\"sample\":"<<start+i<<",\"stream_family\":\""<<stream<<"\",\"beta\":";vector_json(beta);
          std::cout<<",\"training\":";vector_json(y);std::cout<<",\"future\":";vector_json(fy);std::cout<<",\"posterior_mean\":";vector_json(c.value);std::cout<<",\"future_mean\":";vector_json(f.mean);
          std::cout<<",\"posterior_quadratic\":"<<pd.quadratic<<",\"future_quadratic\":"<<f.quadratic<<",\"status\":"<<int(s::DensityStatus::finite)<<",\"numerical_status\":"<<int(n::Status::ok)<<",\"h0_check\":"<<(h0_requested?"true":"null")<<",\"posterior_coverage_bits\":\"";for(bool bit:ph)std::cout<<bit;std::cout<<"\",\"future_coverage_bits\":\"";for(bool bit:fh)std::cout<<bit;std::cout<<"\",\"posterior_joint_bit\":"<<(fixed?"null":pd.quadratic<=16?"true":"false")<<",\"future_joint_bit\":"<<(fixed?"null":f.quadratic<=6?"true":"false")<<"}\n";
          failure={};
         }else if(h.status==s::DensityStatus::finite){failed_status=s::DensityStatus::numerical_failure;failed_numerical=n::Status::conditioning_budget_exceeded;}
        }
       }
      }
     }
    }
    }catch(const std::bad_alloc&){failed_status=s::DensityStatus::numerical_failure;failed_numerical=n::Status::work_limit;}
     catch(const std::exception&){++refused;failed(name,seed,start+i,failure,beta,y,fy,bp,bt,bf,i,failed_status,failed_numerical);throw;}
    if(!failure.empty()){++refused;failed(name,seed,start+i,failure,beta,y,fy,bp,bt,bf,i,failed_status,failed_numerical);}
   }
  }
  }catch(const std::exception&){std::cout<<"{\"kind\":\"campaign_aborted\",\"campaign\":\""<<name<<"\",\"planned_count\":"<<attempts<<",\"actually_attempted\":"<<actually_attempted<<",\"refused_attempts\":"<<refused<<"}\n";throw;}
  std::cout.flush();const W elapsed=std::chrono::duration<W>(std::chrono::steady_clock::now()-start_time).count();
  std::cout<<std::setprecision(17)<<"{\"kind\":\"empirical_receipt\",\"operation\":\"Gaussian-recovery/ideal-model-comparison/v1\",";campaign_identity(name,fixed,stream,p,nt,k,generating_prior,generating_train,generating_future,prior_noise,source,future_noise);
  std::cout<<",\"actually_attempted\":"<<actually_attempted<<",\"attempts\":"<<attempts<<",\"refused\":"<<refused<<",\"generator\":\""<<irred::random::generator_id<<"\",\"generation_plus_predictive_preparation_work_units\":"<<declared_work<<",\"total_seconds_including_stage_and_attempt_logging\":"<<elapsed<<",\"prediction_preparation_density_seconds\":"<<predictive_seconds<<",\"posterior_coverage\":[";
  for(size_t j=0;j<p;++j){if(j)std::cout<<',';std::cout<<"{\"hits\":"<<posterior_hits[j]<<",\"lower\":"<<W(posterior_hits[j])/attempts<<",\"upper\":"<<W(posterior_hits[j]+refused)/attempts<<",\"ideal_reference\":"<<(fixed?ref.posterior_probability[j]:.95L)<<'}';}
  std::cout<<"],\"future_coverage\":[";for(size_t j=0;j<k;++j){if(j)std::cout<<',';std::cout<<"{\"hits\":"<<future_hits[j]<<",\"lower\":"<<W(future_hits[j])/attempts<<",\"upper\":"<<W(future_hits[j]+refused)/attempts<<",\"ideal_reference\":"<<(fixed?ref.future_probability[j]:.95L)<<'}';}
  std::cout<<"],\"posterior_joint_hits\":"<<(fixed?"null":std::to_string(posterior_joint_hits))<<",\"future_joint_hits\":"<<(fixed?"null":std::to_string(joint_hits))<<"}\n";
  if(posterior_moments.count>1&&future_moments.count>1){posterior_moments.report(name,fixed?"fixed truth":"prior predictive","posterior mean error");future_moments.report(name,fixed?"fixed truth":"prior predictive","future prediction error");}
  need(actually_attempted==attempts&&refused==0,"refused or unattempted main vectors withhold empirical qualification");
  const std::vector<W>zero(p),zf(k);posterior_moments.check(fixed?std::span(ref.bias):std::span(zero),fixed?std::span(ref.T):std::span(ref.V));
  future_moments.check(fixed?std::span(ref.bfuture):std::span(zf),fixed?std::span(ref.Tfuture):std::span(ref.Wcov));
  for(size_t j=0;j<p;++j)coverage(posterior_hits[j],fixed?ref.posterior_probability[j]:.95L);
  for(size_t j=0;j<k;++j)coverage(future_hits[j],fixed?ref.future_probability[j]:.95L);
  if(!fixed){coverage(joint_hits,1-std::exp(-3.L));coverage(posterior_joint_hits,chi_square_even(p,16));}
 }
}
} // namespace recovery_controls
