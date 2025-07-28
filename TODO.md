# TODO.md - Primeron Sieve Optimizations and Future Improvements

This code is not optimized for CPU or memory and there are many opportunities
to dramatically improve both.  This code primarily intends to show how large
layered wheels can be built fully sieving all primes and how holes are
canceled only once.

This document captures all in-code commentary related to performance,
efficiency, refactoring, and future enhancements. Grouped by subsystem and
lightly edited for clarity.

---

## Inner Wheel Optimizations

* Consider avoiding memory waste by eliminating storage of both `next` and
  `delta_to_next` in `PseudoprimeNode`. Only the delta should be necessary.
* In the inner wheel cancellation loop, multiplication is used for `p * q`.
  This could potentially be replaced with delta-stepping as used in the outer
  wheel.

---

## Outer Slice Optimizations

* The first slice of the outer wheel contains the same information that is
  already in the inner wheel.  The code should avoid building it to save
  memory.
* If a large prime fails to cancel any holes in a slice, it may be faster to
  compute how many subsequent primes will also fail to cancel holes and skip
  them.

---

## Memory and Storage Efficiency

* This algorithm could be made highly parallel.
* Slices and even the inner wheel could be written to disk and paged into memory
  dynamically to allow very large wheels.
* When saved to disk, hole positions can be saved as **relative distances**
  between holes instead of absolute values — enabling fast compression and
  reconstruction.  Each file would store the start and end of the wheel,
  storing wheel start/end values in multiple bases (e.g. base-7!, base-19!)
  enables constant-time divisibility checks for small primes—without needing
  to walk the inner wheel.  The holes themselves would be stored in relative
  distance to the previous prime, typically a small number of the form 2*n.
  All linked lists prev/next values could work like the pseudoprime_map and
  save prev/next as the distance to the next prime to avoid storing larger
  numbers
* After construction, `prev` links in the pseudoprime and prime lists are not
  used. These could be dropped from memory to save space, or perhaps the
  linked list could be built in a way that avoids needing back links.
* The `pseudoprime_map` duplicates information already in
  `holes_structural_for_outer`. Logic could be unified to reduce storage.
* Gaps between holes are typically small and often of the form `2n`. Saving
  these as deltas instead of numbers in next enables better compression.
* The prime wheel could be built incrementally—one prime at a time—by extending
  previous wheels and copying existing patterns. For example, holes created by
  2 and 3 repeat every 2×3 = 6; when adding 5, the pattern repeats every
  2×3×5 = 30. Instead of recomputing, the sieve could reuse earlier wheel
  structures as lightweight references or templates.

---

## Validation and Reporting

* Validation for the inner wheel is exhaustive and may need a smarter fallback
  at scale.
* The current validation loop tests every number up to the next outer wheel. At
  larger sizes, switch to random sampling or batched testing.

---

