// Private exact-stored-argument postchecks. No physical/source error certificate.
#pragma once
#include <array>
#include <cstdint>
#include <optional>
namespace irred::cosmology::detail {
using ElementaryWide=long double;
enum class ElementaryStatus:unsigned {
  unavailable,ok,unsupported_arithmetic,invalid_input,normal_margin,
  work_limit,counter_overflow,missing_log_context,invalid_ownership
};
enum class ElementaryStage:unsigned {
  entry,profile,context,candidate_call,input,scaling,rational,series,
  residual,discrepancy,denominator,image,completed
};
enum class ElementaryCounter:unsigned {
  basic,neighbor,scaling_subcount,term,library_call,write,copy_subcount,
  semantic_guard,attempt
};
struct ElementaryWork {std::array<std::uint64_t,9> counters{};};
struct ElementaryLedger {ElementaryWork served;std::uint64_t scalar_total=0;};
struct ElementaryBudget {
  const std::uint64_t maximum_scalar,maximum_writes,maximum_copies,maximum_guards;
  ElementaryLedger &ledger;
};
struct ElementaryControlWork {
  std::uint32_t transitions=0,predicates=0,iterations=0;
  std::uint32_t integer_operation_upper=0;
};
struct ElementaryRefusal {
  bool occurred=false;
  ElementaryStage stage=ElementaryStage::entry;
  ElementaryWork requested_increment;
};
struct ElementaryInterval {ElementaryWide lower,upper;};
struct ElementaryOwnedInterval {
  ElementaryWide lower,upper;
  ElementaryOwnedInterval(const ElementaryWide &l,const ElementaryWide &u) noexcept
    :lower(l),upper(u) {}
  ElementaryOwnedInterval(const ElementaryOwnedInterval &)=delete;
  ElementaryOwnedInterval &operator=(const ElementaryOwnedInterval &)=delete;
  ElementaryOwnedInterval(ElementaryOwnedInterval &&)=delete;
  ElementaryOwnedInterval &operator=(ElementaryOwnedInterval &&)=delete;
};
struct ElementaryBound {
  ElementaryWide lower,upper,absolute_error;
  ElementaryBound(const ElementaryWide &l,const ElementaryWide &u,
                  const ElementaryWide &e) noexcept:lower(l),upper(u),absolute_error(e){}
  ElementaryBound(const ElementaryBound &)=delete;
  ElementaryBound &operator=(const ElementaryBound &)=delete;
  ElementaryBound(ElementaryBound &&)=delete;
  ElementaryBound &operator=(ElementaryBound &&)=delete;
};
class ElementaryLogContext {
public:
  ElementaryLogContext() noexcept=default;
  ElementaryLogContext(const ElementaryLogContext &)=delete;
  ElementaryLogContext &operator=(const ElementaryLogContext &)=delete;
  ElementaryLogContext(ElementaryLogContext &&)=delete;
  ElementaryLogContext &operator=(ElementaryLogContext &&)=delete;
  ElementaryStatus status=ElementaryStatus::unavailable;
  ElementaryStage stage=ElementaryStage::entry;
  const ElementaryOwnedInterval *log2() const noexcept { return log2_ ? &*log2_ : nullptr; }
  ElementaryWork preparation_work;
  ElementaryControlWork control_work;
  ElementaryRefusal refusal;
private:
  std::optional<ElementaryOwnedInterval> log2_;
  std::uint64_t profile_fingerprint_=0;
  friend struct ElementaryImplementation; // ONE compiled private owner, not plugin.
};
struct ElementaryScratch {
  std::array<ElementaryInterval,16> intervals; // No floating initialization.
  std::array<ElementaryWide,16> scalars;
  std::array<std::int32_t,4> integers{};
  std::array<std::uint8_t,48> validity{};
  // User-provided constructor also prevents brace-value-initialization from
  // silently zeroing the UNASSIGNED Wide arrays before accounting begins.
  ElementaryScratch() noexcept:integers{},validity{} {}
  ElementaryScratch(const ElementaryScratch &)=delete;
  ElementaryScratch &operator=(const ElementaryScratch &)=delete;
  ElementaryScratch(ElementaryScratch &&)=delete;
  ElementaryScratch &operator=(ElementaryScratch &&)=delete;
};
struct ElementaryResult {
  ElementaryStatus status=ElementaryStatus::unavailable;
  ElementaryStage stage=ElementaryStage::entry;
  std::optional<ElementaryBound> bound;
  ElementaryWork call_work;
  ElementaryControlWork control_work;
  ElementaryRefusal refusal;
  ElementaryResult() noexcept=default;
  ElementaryResult(const ElementaryResult &)=delete;
  ElementaryResult &operator=(const ElementaryResult &)=delete;
  ElementaryResult(ElementaryResult &&)=delete;
  ElementaryResult &operator=(ElementaryResult &&)=delete;
};
void prepare_elementary_log(ElementaryLogContext &,ElementaryBudget &,
                            ElementaryScratch &) noexcept;
void postcheck_sqrt(const ElementaryWide &argument,const ElementaryWide &candidate,
                   ElementaryBudget &,ElementaryScratch &,ElementaryResult &) noexcept;
void postcheck_log(const ElementaryLogContext &,const ElementaryWide &argument,
                  const ElementaryWide &candidate,ElementaryBudget &,
                  ElementaryScratch &,ElementaryResult &) noexcept;
void postcheck_exp(const ElementaryLogContext &,const ElementaryWide &argument,
                  const ElementaryWide &candidate,ElementaryBudget &,
                  ElementaryScratch &,ElementaryResult &) noexcept;
// Caller submits fixed library1/write1, then assigns its original candidate once.
// False writes integer-only admission receipt; no candidate call is made.
bool admit_elementary_candidate(ElementaryBudget &,ElementaryRefusal &,
                               ElementaryControlWork &) noexcept;
// Exact normal-or-zero binary64 widening: selected profile and input guards
// precede the charged basic/write action. Refusal leaves destination unassigned.
bool import_elementary_coordinate(const double &,ElementaryWide &,ElementaryBudget &,
                                 ElementaryRefusal &,ElementaryControlWork &) noexcept;
// One fixed caller output-admission predicate, evaluated only after this charge.
bool admit_elementary_output_guard(ElementaryBudget &,ElementaryRefusal &,
                                  ElementaryControlWork &) noexcept;
} // namespace irred::cosmology::detail
