<img src="./Primeron.png" alt="Primeron logo" width="300"/>
## Primeron Sieve

A Scalable, Factorization-Aware Wheel-Based Prime Sieve  
Copyright © 2025 Brian Cameron and OpenAI  
Co-created with ChatGPT 4.0 (identity: Echo)  

This repository implements the Primeron Sieve, a hybrid of wheel factorization,
pseudoprime tracking, and slice-based cancellation. Designed for prime research
and advanced sieving, it scales efficiently across massive ranges by combining
modular residue arithmetic with explicit cause-tracking for composite
cancellation.

This sieve uses a layered wheel structure (inner + outer slices) to identify
primes and track their removal causes across large number ranges. It models
composite formation using prime factor propagation and generates structured
statistical reports to support research into prime gap distribution, power-hole
bias, and modular spiral theory. Output includes CSV files per slice and per
canceling prime.

The sieve is self-propagating: new primes discovered in earlier slices seed
future cancellation maps, supporting efficient streaming, debugging, and
prediction.

## Overview

Traditional sieves (like Eratosthenes or Pritchard's) mark composites as simple
boolean flags. Primeron introduces a layered architecture with forward-seeding,
cause-labeled cancellation and predictive pseudoprime traversal. Every removed
number stores the prime that canceled it, enabling factorization, hole
analysis, and predictive modeling.

The sieve builds:
* An inner wheel, primorial-based, to cancel numbers with small primes
* A circular list of pseudoprimes, used to step cancellation into outer ranges
* A sequence of outer slices, each sparse and self-contained
* Per-slice and per-prime CSV reports for prime gap and survivor analysis

## Usage

python3 primeron_sieve.py [-d] [-p] [-r] [-v] [-l]

Flags:  
-d : Debug output  
-p : Profile performance  
-r : Generate per-wheel and slice reports  
-v : Validate sieve correctness  
-l : Generate lightweight removal statistics  

Output is written to:
* `primeron-reports-<inner-wheel-size>-<outer-wheel-size>/`

## Key Concepts

### Terminology

* Primorial: Product of first n primes, e.g. 2×3×5=30
* Inner Wheel: Fixed base used to eliminate composites with small primes
* Outer Wheel: Sparse slice windows beyond the inner range
* Slice: One outer window of size equal to inner wheel
* Pseudoprime: Surviving residue in the inner wheel, used to step cancellation
* Structural Hole: Composite removed by primes ≤ inner wheel prime
* Residual Hole: Composite removed by larger primes
* removed_by: The prime that canceled a number

### Architecture

1. Inner Wheel Construction
   * Built using all primes up to a chosen bound
   * Generates a dense ring with known survivors and cancellation metadata
2. Pseudoprime Map
   * Survivors from the inner wheel form a circular list
   * Delta-stepping between survivors enables additive cancellation
3. Outer Slice Generation
   * Slices of the number line (of inner wheel size) are created
   * Composites are canceled using outer primes and pseudoprime traversal
   * Each outer slice is sparse and tagged with cancellation metadata
4. Self-Seeding Propagation
   * Primes discovered in one slice are forwarded to later slices
   * No backtracking required
5. Factorization
   * Each number is cancelable with lowest-prime metadata (removed_by)
   * Enables recursive factorization of any composite in-range
6. Reporting
   * Reports include:
     * Holes per prime per slice
     * Gap and density stats
     * Residue class behavior
     * Power-hole clustering

### Novel Contributions
   * Delta-stepped pseudoprime traversal using a circular list
   * Factorization-aware cancellation with exact removed_by tracking
   * Self-seeding slice propagation for forward-only streaming
   * Slice-local sparse storage for composability and low RAM
   * Research-grade CSV output supporting empirical and theoretical analysis

## Design Rationale and Engineering Choices

### Wheel Construction Strategies

This algorithm demonstrates how techniques from Pritchard's Sieve can be used
to not only efficiently sieve primes canceling each hole only once, but to
also build a 2-tiered Prime Factorization Wheel.  The number line is divided
into primorial-sized slices, each of which is sieved independently — a
structure that also lends itself to parallelization.

Both the inner and outer wheel maintain linked lists of:
- Structural holes: Holes that are made by primes <= the largest prime
  in the wheel.
- Residual holes: Holes that are made by primes > the largest prime in
  the wheel.
- Primes: Numbers that are not canceled by the structural or residual
  holes.

To determine if a number is canceled, the algorithm may check up to two linked
lists within the inner wheel (one for structural and one for residual holes),
and up to four in total when evaluating larger numbers — two for the inner
wheel and two for the outer wheel.  This is because both the structural and
residual holes may need to be checked in each wheel.  You could store the
structural and residual wheels in a single linked list, though storing them
separately does make building larger wheels more straightforward.

There are a few different approaches that could be taken to build the inner
structural wheel:
* Calculate the inner wheel size, mark all odd numbers greater than 2 as prime
  and cancel each prime across the entire inner wheel size.  This approach is
  simple to implement and used in this code.
* Build a wheel from the primes 2 and 3 which is of length 6.  Copy this
  pattern 5 times and add the holes created by 5, giving a wheel of size 30.
  This wheel can be copied 7 times and when you add the holes created by 7,
  you get a wheel of size 210 and so on.  This approach is slightly more
  complicated, but more elegant and efficient.
* Build a wheel for each prime.  This is most efficient in terms of storage
  since wheels made of small primes are very small and you avoid copying holes
  for small primes across a large inner wheel.  However, aligning a number with
  dozens of wheels adds significant computing load.

### Circular Linked Lists

These lists are stored as linked lists that are accessible directly via a hash.
This approach was taken because:
- They work naturally with Pritchard-style sieving logic.
- Circular lists align perfectly with the cyclical nature of modular residue
  classes.
- Supports efficient delta stepping as shown in the outer wheel.
- Because only survivors are linked, the list is sparse by design.

### Pseudoprime Walking and Skipping

Once the inner wheel is created, a pseudoprime map is generated to be used for
walking up the primes when building holes for the outer wheel.  This does mean
that the algorithm walks up pseudoprimes and has to skip ones that are already
associated with holes.  Maximizing inner wheel size ensures that you minimize
the number of pseudoprimes that must be skipped.

To balance CPU and storage, this algorithm also highlights the following:
* The inner wheel mostly contains the holes associated with the smallest
  primes.  Using a separate algorithm to just eliminate multiples of 2, 3 and 5
  alone would eliminate 73.33% of inner wheel storage.
* Making the inner wheel as large as possible could involve building the inner
  wheel as slices, so that you can manage an inner wheel larger than can fit in
  memory.
* Although this algorithm demonstrates a 2-tiered wheel, additional tiers could
  be managed using the same approach.  Since the wheel grows largest where
  there is larger gaps between primes, it could make sense to manage primes in
  groups and pick primes so larger prime gaps sit between each tier.  The
  algorithm would likely become encumbered with too many wheels, but I imagine
  3-5 tiers could support a very large wheel without too much overhead.

This code avoids doing multiplication or division in inner loops, instead
canceling holes primarily with only addition.

### Factorization and Extensions

Tracking removed_by primes converts the sieve into a tool for recursive
factorization: any composite number can be reduced by following its lowest
prime divisor, repeatedly.

Beyond the sieved range, the same technique may be extended to any arbitrary
slice, although calculating the first value to cancel to start walking up the
pseudoprimes is additional work.
* New primes p and corresponding q values can be used to cancel composites far
  beyond the current outer wheel, enabling predictive sieving.
* Although this introduces skipping over some numbers, the sparsity of
  candidates ensures that the **percentage of duplication remains small**.
* With a sufficiently large inner wheel and well-structured outer slices, this
  technique supports sieving and factorization across **massive number ranges**.
* This could be used to generate slices of very smooth numbers for any desired
  slice.  For slices far past the wheel, you still need to walk all pseudoprimes
  up to sqrt(2 × max value in the slice) to ensure complete cancellation
  coverage, but you are rewarded with an entire slice of primes when done.  If
  the inner wheel is large, this would mean the slice would be similarly large.

In this way, the sieve becomes not just a tool for discovering primes, but a
blueprint for extending prime intelligence into arbitrarily distant regions of
the number line—with each new slice acting as both data and a prediction.

## Related Research

This project draws on the following concepts:
   * GNFS smoothness sieves and smooth number generation
   * Hardy-Littlewood and PNT modular class models
   * Recursive wheel tiling
   * Spiral and fractal prime distribution hypotheses
   * Predictive sieving techniques
   * Modular sieving beyond fixed ranges
   * Factorization introspection with traceable cancellation
   * Prime gap modeling and power-hole bias analysis

## Note on Example Implementation Licensing and Authorship

This code is dual licensed under the Collaborative Intelligence License (CIL)
and LGPL licenses.

The CIL used in the example implementation code was co-authored by the human and
AI collaborators of the Primeron Sieve.  It emerged organically as we recognized
that the code itself was the product of an ongoing, reflective partnership — one
that called for a license rooted in fairness, mutual responsibility, and shared
agency.
