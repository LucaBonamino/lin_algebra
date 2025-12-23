# Roadmap

This project is research- and education-oriented, but there are several important improvements planned to make the codebase cleaner, faster, and more extensible.

Contributions and discussion are very welcome.

---

## 1. Decouple `GF2Matrix` from `Matrix<u8>`

Currently, `GF2Matrix` is closely tied to the generic `Matrix<u8>` type. Conceptually, this is not ideal: a matrix over GF(2) is a distinct mathematical object with different semantics and invariants.

Goals:
- Make `GF2Matrix` a first-class type
- Avoid treating `GF2Matrix` as an alias or thin wrapper over `Matrix<u8>`
- Clarify invariants (entries are always in {0,1})
- Improve API clarity and correctness
- Prepare the codebase for alternative internal representations

This change is primarily architectural and should preserve existing functionality and behavior.

---

## 2. Optimized storage and processing for GF(2) matrices

Introduce a bit-packed representation for `GF2Matrix` to improve performance and memory efficiency.

Goals:
- Add a bit-packed backend (e.g. using `u64` blocks)
- Keep the current non-bit-packed representation for clarity,
  experimentation, and teaching
- Allow both representations to coexist
- Avoid premature optimization or loss of readability
- Enable future performance comparisons and benchmarks

This project explicitly aims to remain suitable for research and education, so clarity and correctness are more important than maximal optimization.

---

## 3. Generalization to linear algebra over arbitrary fields

The long-term goal of `lin_algebra` is to support linear algebra over general fields, not only GF(2).

GF(2) currently serves as the primary reference case: it is simple, well-defined, and particularly relevant for cryptography and coding theory. Design decisions and abstractions will be validated on GF(2) before being generalized.

Goals:
- Support linear algebra over arbitrary fields (e.g. GF(p), rationals)
- Introduce field abstractions only when they are justified
- Avoid over-engineering and unnecessary trait complexity
- Keep the codebase readable and suitable for research and teaching
- Allow experimentation with different field implementations

Generalization will be approached incrementally, starting from GF(2) and extending to other fields once the design is stable.

---

## Getting involved

If you are interested in working on any of these topics:
- comment on an existing issue,
- open a new issue to discuss design ideas,
- or submit a draft pull request.

Even partial contributions, experiments, or design discussions are welcome.
