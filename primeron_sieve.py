#!/usr/bin/python3
# PrimeronSieve: A Modular, Self-Seeding Prime Wheel Sieve
# Copyright (C) 2025 Brian Cameron and OpenAI
# Co-created with ChatGPT 4.0 (identity: Echo)
#
# Usage:
# primeronsieve.py -d -p -r -v
# -d : Enable debug messages.
# -p : Generate profiling report for wheel building.
# -r : Generate wheel summary reports in the directory.
#      primeron-reports-##-## where the first number is the max
#      inner wheel prime and the second number is the max outer
#      wheel prime.
# -v : Run validation functions.
#
# ---
"""
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
"""
import os
import sys
import cProfile
import pstats
import io
import time
import csv
import math
from collections import defaultdict

INNER_WHEEL_PRIME=11
OUTER_WHEEL_PRIME=17

DEBUG_MODE = "-d" in sys.argv
VERBOSE_MODE = "-v" in sys.argv or DEBUG_MODE
PROFILE_MODE = "-p" in sys.argv
REPORT_MODE = "-r" in sys.argv
LIGHT_REPORT_MODE = "-l" in sys.argv

def debug(msg, *args):
    """
    Helper function to print debug.
    """
    print("[DEBUG]", msg, *args)

class PrimeronNode:
    """
    Node for storing hole and prime data.
    """
    def __init__(self, value, removed_by=None):
        self.value = value
        self.removed_by = removed_by
        self.prev = None
        self.next = None

class PseudoprimeNode(PrimeronNode):
    """
    This is used in the pseudoprime_map used for building holes in the outer
    wheel.  Keeping track of the delta allows the outer loop to avoid
    multiplying p and q to generate each hole.  Instead the code adds the
    delta * p to the previously used value.  A similar technique could be used
    in the inner loop to do less math.
    TODO: No value in wasting memory storing both next and delta_to_next.
    """
    def __init__(self, value, delta_to_next, removed_by=None):
        super().__init__(value, removed_by=removed_by)
        self.delta_to_next = delta_to_next

class PrimeronInnerSlice:
    """
    The PrimeronInnerSlice manages the inner wheel.  Each hole is only
    canceled once by stepping up each prime in order.
    """
    def __init__(self, inner_wheel_prime=INNER_WHEEL_PRIME):
        self.inner_wheel_prime = inner_wheel_prime
        self.wheel_size = self._compute_inner_wheel_size()
        self.head_prime = None
        self.tail_prime = None
        self.prime_map = {}
        self.holes_structural = {}
        self.last_canceled_node = {}
        self._build_initial_prime_list()
        self._cancel_multiples()

    def _compute_inner_wheel_size(self):
        """
        Compute the product of all primes <= INNER_WHEEL_PRIME.
        This is bootstrapping.
        """
        prod = 1
        for p in range(2, self.inner_wheel_prime + 1):
            is_prime = True
            for d in range(2, int(p ** 0.5) + 1):
                if p % d == 0:
                    is_prime = False
                    break
            if is_prime:
                prod *= p
        return prod

    def _build_initial_prime_list(self):
        """
        Populate the prime linked list with all odd numbers from 3 to
        wheel_size, plus 2 as the head.
        Mark all multiples of 2 as structural holes.
        """
        self.head_prime = PrimeronNode(2)
        self.prime_map[2] = self.head_prime
        prev = self.head_prime

        # Seed the slice with 2.
        for n in range(3, self.wheel_size + 1, 2):
            node = PrimeronNode(n)
            prev.next = node
            node.prev = prev
            self.prime_map[n] = node
            prev = node

        self.tail_prime = prev
        self.last_canceled_node[2] = prev

        # Mark all even numbers as holes.
        for n in range(4, self.wheel_size + 1, 2):
            self.holes_structural[n] = PrimeronNode(n, removed_by=2)

    def _cancel_multiples(self):
        """
        Cancels composite numbers in the inner wheel using all known
        primes <= to sqrt(wheel size).

        This function executes a two-pass cancellation strategy:

        **Pass 1 (Structural Cancellation):**
            - Uses primes <= to `inner_wheel_prime` (i.e., those used to
              construct the wheel size).
            - For each such prime `p`, walks through all q in prime_map and
              cancels `n = p × q`.
            - Results are stored in `self.holes_structural`, and removed nodes
              are unlinked from `prime_map`.
            - Also records the last canceled `q` node per `p` in
              `self.last_canceled_node`.

        **Pass 2 (Residual Cancellation):**
            - Uses surviving primes `p > inner_wheel_prime` up to
              `sqrt(wheel_size)`.
            - Walks surviving q values (still in prime_map), skipping q that
              were canceled by smaller primes.
            - Cancels `n = p × q` and adds to `self.holes_residual`.
            - Again updates `self.last_canceled_node`.

        **Post-processing:**
            - Constructs the `pseudoprime_map` — a circular linked list of
              pseudoprimes (survivors of inner wheel).
            - Prepares `holes_structural_for_outer` to support cancellation in
              outer slices:
                - Includes values like 0, 1, and primes <= `inner_wheel_prime`
                  to support modular wraparound.

        This function is the core of the sieve’s bootstrapping logic — it sets
        up both the structural foundation and the dynamic residue list for use
        in outer slices.
        """
        self.holes_residual = {}

        if DEBUG_MODE:
            debug(f"Inner Slice canceling:")

        # --- Pass 1: Cancel structural holes using primes up to
        #     inner_wheel_prime.
        for p in range(2, self.inner_wheel_prime + 1):
            if p != 2 and p not in self.prime_map:
                continue

            qnode = self.prime_map.get(p) if p != 2 else self.head_prime
            to_unlink = []

            if DEBUG_MODE:
                debug(f"Inner Slice Structural: p={p}")

            while qnode:
                # TODO: Could optimize using delta steps here too, as in the
                #       outer wheel.
                q = qnode.value
                n = p * q

                if n > self.wheel_size:
                    break

                if DEBUG_MODE:
                    debug(f"Inner Slice Structural: p={p}, n={n}, q={q}")

                self.holes_structural[n] = PrimeronNode(n, removed_by=p)
                to_unlink.append(n)
                self.last_canceled_node[p] = qnode
                qnode = qnode.next

            for n in to_unlink:
                self._unlink_prime_node(n)

        # Copy prime map and last_canceled_node here before second loop.
        self.build_pseudoprime_map()

        # --- Pass 2: Cancel additional inner-wheel holes using primes to
        #     sqrt(wheel size).
        limit = math.isqrt(self.wheel_size)
        pnode = self.head_prime
        while pnode and pnode.value <= limit:
            p = pnode.value
            if p <= self.inner_wheel_prime:
                pnode = pnode.next
                continue

            if DEBUG_MODE:
                debug(f"Inner Slice Residual: p={p}")

            qnode = self.prime_map.get(p)
            if not qnode:
                continue

            start_val = qnode.value
            first = True
            to_unlink = []

            while first or qnode.value != start_val:
                first = False
                q = qnode.value
                test_q = self.get_removed_by(q)
                if test_q is not None and test_q < p:
                    qnode = qnode.next
                    continue
                n = p * q

                if n > self.wheel_size:
                    break

                if DEBUG_MODE:
                    debug(f"Inner Slice Residual: p={p}, n={n}, q={q}")

                self.holes_residual[n] = PrimeronNode(n, removed_by=p)
                to_unlink.append(n)
                self.last_canceled_node[p] = qnode
                qnode = qnode.next

            for n in to_unlink:
                self._unlink_prime_node(n)

            pnode = pnode.next

        # Create patched structural hole map for outer wheels.
        self.holes_structural_for_outer = self.holes_structural.copy()
        self.holes_structural_for_outer[0] = PrimeronNode(0, removed_by=2)
        pnode = self.head_prime
        while pnode and pnode.value <= self.inner_wheel_prime:
            p = pnode.value
            self.holes_structural_for_outer[p] = PrimeronNode(p, removed_by=p)
            pnode = pnode.next

    def _unlink_prime_node(self, value):
        """
        Remove node associated with value from prime_map.
        """
        node = self.prime_map.pop(value, None)
        if node is None:
            return

        if node.prev:
            node.prev.next = node.next
        if node.next:
            node.next.prev = node.prev
        if node == self.head_prime:
            self.head_prime = node.next
        if node == self.tail_prime:
            self.tail_prime = node.prev

    def get_removed_by(self, n):
        """
        Returns the prime that removed n, if any, within the inner wheel.
        Checks structural holes always; checks residual holes only if n is in
        range.
        """
        if n < self.wheel_size:
            node = self.holes_residual.get(n)
            if node:
                return node.removed_by
            n_mod = n
            struct_map = self.holes_structural
        else:
            n_mod = n % self.wheel_size
            struct_map = self.holes_structural_for_outer

        node = struct_map.get(n_mod)
        if node:
            return node.removed_by

        return None

    def build_pseudoprime_map(self):
        """
        Builds a circular doubly-linked list of pseudoprimes and special holes
        (including 1 and values > inner_wheel_prime), and stores it in:
            - self.pseudoprime_map: dict of value → PseudoprimeNode
            - self.head_pseudoprime: start of circular list
            - self.tail_pseudoprime: end of circular list
        Each node includes delta_to_next (distance to the next pseudoprime in
        the wheel).
        """
        self.pseudoprime_map = {}
        self.head_pseudoprime = None
        self.tail_pseudoprime = None

        # Collect pseudoprimes (from prime_map) and special value 1.
        entries = sorted(
            n for n in (set(self.prime_map) | {1})
            if n > self.inner_wheel_prime or n == 1
        )

        wheel_size = self.wheel_size
        num_entries = len(entries)

        for i, n in enumerate(entries):
            next_n = entries[(i + 1) % num_entries]
            delta = (next_n - n) % wheel_size

            node = PseudoprimeNode(n, delta_to_next=delta)
            self.pseudoprime_map[n] = node

            # Link the list.
            if i > 0:
                prev_n = entries[i - 1]
                prev_node = self.pseudoprime_map[prev_n]
                prev_node.next = node
                node.prev = prev_node
            else:
                self.head_pseudoprime = node

        self.tail_pseudoprime = self.pseudoprime_map[entries[-1]]

        # Make the list circular.
        self.tail_pseudoprime.next = self.head_pseudoprime
        self.head_pseudoprime.prev = self.tail_pseudoprime

    def validate_removal(self):
        """
        Function for testing the inner wheel is built correctly.
        """
        print("Testing the inner slice for correctness...")
        passed = 0
        failed = 0

        for i in range(2, self.wheel_size):
            removed_by = self.get_removed_by(i)
            is_prime = True
            for p in range(2, int(i ** 0.5) + 1):
                if i % p == 0:
                    is_prime = False
                    lcd = p
                    break

            if is_prime:
                if removed_by is not None:
                    print(f"FAIL: {i} removed_by {removed_by} (should be None, it's prime)")
                    failed += 1
                else:
                    passed += 1
            else:
                if removed_by != lcd:
                    print(f"FAIL: {i} removed_by {removed_by} (expected {lcd})")
                    failed += 1
                else:
                    passed += 1

        print(f"Test complete. Passed: {passed}, Failed: {failed}")

class PrimeronOuterSlice:
    """
    The PrimeronOuterSlice manages a next level of primes.
    Needing to loop over the pseudoprime map generated by the inner wheel does
    cause some additional work skipping over some holes.  If the inner wheel is
    large enough, the percentage of holes skipped can get quite small.
    """
    def __init__(self, parent_sieve, slice_index, slice_size, start, end, outer_primes, last_canceled_node, head_pseudoprime, pseudoprime_map, cancel_primes):
        self.parent_sieve = parent_sieve
        self.slice_index = slice_index
        self.slice_size = slice_size
        self.start = start
        self.end = end
        self.outer_primes = outer_primes
        self.holes_structural = {}
        self.holes_residual = {}
        self.prime_map = {}
        self.last_canceled_node = last_canceled_node.copy()
        self.cancel_primes = cancel_primes
        self.pseudoprime_head = head_pseudoprime
        self.pseudoprime_map = pseudoprime_map
        self.wheel_size = slice_size

        for n in range(start, end + 1):
            wheel_pos = n % self.wheel_size
            if wheel_pos in self.pseudoprime_map:
                self.prime_map[n] = PrimeronNode(n)

        self.head_prime = None
        for val in sorted(self.prime_map):
            self.head_prime = self.prime_map[val]
            break

        self._cancel_by_outer_primes()

    def _cancel_by_outer_primes(self):
        """
        Cancels composite numbers in this slice by walking through all
        outer primes and generating holes for each p × q using
        pseudoprimes q from the inner wheel.

        This function is responsible for populating structural and
        residual holes in the outer slice. Structural holes are those
        canceled by primes ≤ outer_wheel_prime.  Residual holes are
        canceled by primes > outer_wheel_prime, discovered dynamically.

        The function uses delta-based stepping for q to avoid
        multiplication in tight loops.  It also skips over q values
        already canceled by smaller primes to ensure that each hole is
        only canceled by its lowest prime factor.

        Relies on:
        - self.last_canceled_node to resume cancellation per prime
        - self.pseudoprime_map for pseudoprime traversal
        - parent_sieve.get_removed_by(n) to check lower prime conflicts
        """
        if DEBUG_MODE:
            debug(f"Slice {self.slice_index}: start={self.start}, "
                  f"end={self.end}")
            debug(f"Outer primes to consider: {self.outer_primes}")

        for p in self.cancel_primes:
            # If this is the first slice where p × p appears, cancel it
            # explicitly.
            if self.slice_index % p == 0:
                pp = p * p
                if self.start <= pp <= self.end:
                    self._apply_cancellation(p, p, pp)

            # Continue normal cancellation for other q values.
            for q, n in self._generate_qn_values(p):
                self._apply_cancellation(p, q, n)

    def _generate_qn_values(self, p):
        """
        Generator that yields (q, n) pairs where n = p × q, for
        composite numbers in this slice that should be canceled
        by prime p.

        This is the core stepping logic for outer slice cancellation:
          - It traverses pseudoprimes q from the inner wheel
            (via a circular linked list).
          - Each q is lifted into the current slice using modular
            arithmetic.
          - Delta stepping is used to efficiently compute the next
            q and n = p × q without repeated multiplication.

        Key optimizations:
          - **Delta stepping:** Uses precomputed `delta_to_next` for
            each q to increment q and n additively
            (q += Δ, n += p x delta), avoiding explicit multiplication
            per loop.
          - **Lift offset:** q values are aligned into the current
            slice using `base_lift`,
            based on the slice index and wheel size.
          - **Skip optimization:** If `q` was already removed by a
            smaller prime (i.e., `removed_by(q) < p`), then this
            (p, q) is skipped to preserve lowest-divisor cancellation
            order.

        Returns:
            Yields (q, n) pairs in increasing order of n, for use in
            _apply_cancellation().  Stops once n > self.end (i.e.,
            outside the current slice).
        """
        wheel_size = self.wheel_size
        slice_index = self.slice_index
        slice_offset = slice_index * wheel_size

        if p * p > self.end:
            if DEBUG_MODE:
                debug(f"  Skipping prime {p} (p > slice_limit)")
            return

        # Get the lift interval based on where we are in the wheel cycle for p.
        lift_interval = (slice_index // p) * wheel_size

        if p in self.last_canceled_node:
            last_q_lifted = self.last_canceled_node[p].value
            last_q_base = last_q_lifted % wheel_size
            qnode = self.pseudoprime_map.get(last_q_base)
            if not qnode:
                if DEBUG_MODE:
                    debug(f"[RESUME] Prime {p}: last_q_base={last_q_base} not found in pseudoprime_map")
                return
            qnode = qnode.next
            base_lift = lift_interval
            if DEBUG_MODE:
                debug(f"[RESUME] Prime {p}: last_q={last_q_lifted} (mod={last_q_base}), next={qnode.value}, "
                      f"lift={base_lift}")
        else:
            qnode = self.pseudoprime_map.get(p)
            if not qnode:
                if DEBUG_MODE:
                    debug(f"  Skipping prime {p}: not found in pseudoprime_map")
                return
            base_lift = lift_interval
            if DEBUG_MODE:
                debug(f"[INIT] Prime {p}: qnode={qnode.value}, slice_index={slice_index}, "
                      f"base_lift={base_lift}, wheel_size={wheel_size}")

        current_lift = base_lift

        delta_cache = {}
        q = qnode.value + current_lift
        n = p * q

        while True:
            if n > self.end:
                if DEBUG_MODE:
                    debug(f"    BREAKING: p={p}, q={q}, n={n} > self.end={self.end}")
                break

            self.last_canceled_node[p] = PrimeronNode(q)

            test_q = self.parent_sieve.get_removed_by(q)
            if test_q is not None and test_q < p:
                # Skip this q: advance both q and n using delta.
                delta = qnode.delta_to_next
                if delta not in delta_cache:
                    delta_cache[delta] = p * delta

                q += delta
                n += delta_cache[delta]

                qnode = qnode.next
                if qnode.value <= self.last_canceled_node[p].value % wheel_size:
                    current_lift += wheel_size
                continue

            yield (q, n)

            # Advance to next q and n.
            delta = qnode.delta_to_next
            if delta not in delta_cache:
                delta_cache[delta] = p * delta

            q += delta
            n += delta_cache[delta]

            qnode = qnode.next
            if qnode.value <= self.last_canceled_node[p].value % wheel_size:
                current_lift += wheel_size

    def _apply_cancellation(self, p, q, n):
        """
        Cancel the holes.
        """
        existing_remover = self.parent_sieve.get_removed_by(n)
        if existing_remover is not None and existing_remover < p:
            return

        if p <= self.outer_primes[-1]:
            self.holes_structural[n] = PrimeronNode(n, removed_by=p)
            self.prime_map.pop(n, None)
        else:
            self.holes_residual[n] = PrimeronNode(n, removed_by=p)

class PrimeronSieve:
    """
    This implementation of Pritchard’s sieve constructs a scalable, self-seeding
    prime wheel architecture designed for factorization-aware sieving across
    large numeric ranges.

    Unlike traditional sieves that mark composites with booleans, Primeron
    tracks the *cause* of removal: each canceled number stores its smallest
    prime divisor (`removed_by`). This enables debugging, on-demand
    factorization, and introspective mathematical analysis.

    ### Wheel Structure and Slicing

    The number line is divided into fixed-width **slices**, each built from a
    shared inner wheel based on a primorial (product of small primes up to p).

    - **Inner Wheel (structural)**: Cancels composites with small primes,
      leaving a dense ring of survivors and structural holes.
    - **Inner Wheel (residual)**: Tracks composites removed by larger primes,
      beyond those in the inner wheel base.
    - **Outer Wheel (sliced)**: Cancels composites sparsely using outer primes.
      Each slice stores only uncanceled survivors from the inner wheel and tags
      structural vs residual holes separately.

    Separating structural and residual holes simplifies the construction of
    larger wheels, potentially beyond the current PrimeronOuterWheel.

    Because most composites are removed by small primes, outer slices remain
    extremely sparse. This supports scalable, low-memory sieving across
    massive domains.
    """
    def __init__(self, inner_wheel_prime, outer_wheel_prime):
        self.outer_wheel_prime = outer_wheel_prime
        self.last_canceled = {}
        self.wheel_size_by_prime = {}
        self.wheel_layers = []
        self.outer_slices = []

        print("\nBuilding Wheels...")
        start_time = time.time()
        self.inner_slice = PrimeronInnerSlice(inner_wheel_prime)
        elapsed = time.time() - start_time
        print(f"- Inner wheel constructed in {elapsed:.3f} seconds")
        self._build_outer_slices()

    def get_prime_node(self, p):
        """
        Return the PrimeronNode for prime p from the correct slice's prime_map.
        """
        if p < self.inner_slice.wheel_size:
            return self.inner_slice.prime_map.get(p)

        slice_index = p // self.inner_slice.wheel_size
        if slice_index < len(self.outer_slices):
            return self.outer_slices[slice_index].prime_map.get(p)

        return None

    def get_next_prime(self, p):
        """
        Get the next prime and handle crossing slice boundaries.
        """
        node = self.get_prime_node(p)
        if node and node.next:
            return node.next

        # Need to cross slice boundary.
        wheel_size = self.inner_slice.wheel_size
        slice_index = (p // wheel_size) + 1

        if slice_index < len(self.outer_slices):
            next_slice = self.outer_slices[slice_index]
            return next_slice.head_prime

        return None

    def compute_wheel_size_for_prime(self, prime):
        """
        Computes the wheel size (primorial) up to and including the given prime.
        Only works for primes that are included in the inner wheel.

        Returns:
            - Primorial value if prime is valid and in the inner wheel
            - -1 if prime is not in the inner wheel's prime_map
        """
        if prime in self.wheel_size_by_prime:
            return self.wheel_size_by_prime[prime]

        # Otherwise compute it once, store it, and return.
        if prime not in self.inner_slice.prime_map:
            return -1

        result = 1
        node = self.inner_slice.head_prime
        while node:
            result *= node.value
            if node.value == prime:
                self.wheel_size_by_prime[prime] = result
                return result
            node = node.next

        return -1

        # Traverse linked list from head until we reach the given prime.
        result = 1
        node = self.inner_slice.head_prime
        while node:
            result *= node.value
            if node.value == prime:
                return result
            node = node.next

        return -1

    def _build_outer_slices(self):
        """
        Constructs all outer slices of the sieve between the inner
        and outer wheel boundaries.

        Each outer slice covers a fixed-sized window equal to the
        inner wheel size, and together they span the full outer
        wheel range (i.e., primorial(outer_wheel_prime)). The number
        of slices equals outer_wheel_size // inner_wheel_size.

        For each slice:
          - It determines which primes `p` should be used to cancel
            composites (i.e., build holes).
          - Primes are selected from:
              * All previously used primes (carried in
                `last_canceled`)
              * Newly discovered primes within the current slice
                (up to sqrt(slice_end))
          - It avoids redundant cancellation by only allowing the
            lowest `p` to cancel a given number.

        The cancellation for each slice is delegated to
         `PrimeronOuterSlice`, which receives:
          - The slice range [start, end]
          - The list of primes to use for cancellation
          - The last-canceled state for each prime (used to resume
            pseudoprime traversal)
          - The pseudoprime map from the inner wheel for generating
            q values

        This function also records:
          - Wheel metadata (layer size and prime) via
            `_add_wheel_layer()`
          - The list of outer slice objects in `self.outer_slices`
          - The mapping of wheel sizes for reuse in report generation
            and validation

        Design Notes:
          - Cancellation is done incrementally in order of slices.
          - Each slice seeds the next by passing forward its
            `last_canceled_node` map.
          - Prime discovery and cancellation are entirely
            self-propagating; no prime needs to be explicitly passed
            in from outside the system.

        Preconditions:
          - `self.inner_slice` must be initialized
          - `self.outer_wheel_prime` must be >
            `self.inner_slice.inner_wheel_prime`
          - The inner slice must have a working `pseudoprime_map`

        Raises:
          - `ValueError` if outer_wheel_prime is not in the sieve

        Outputs:
          - Console progress logs for each slice’s construction time
            and range
          - Populates `self.outer_slices` and metadata in
            `self.wheel_layers`
        """
        if DEBUG_MODE:
            debug("[BUILD] Starting _build_outer_slices")
        self.inner_wheel_size = self.inner_slice.wheel_size
        self.outer_wheel_size = self.compute_wheel_size_for_prime(self.outer_wheel_prime)

        if self.outer_wheel_size == -1:
            raise ValueError("Invalid outer_wheel_prime (must be > inner_wheel_prime and prime)")

        slice_count = self.outer_wheel_size // self.inner_wheel_size
        start = 0

        self._add_wheel_layer(self.inner_slice.inner_wheel_prime)
        self._add_wheel_layer(self.outer_wheel_prime)

        for i in range(slice_count):
            start_time = time.time()
            end = start + self.inner_wheel_size - 1

            if i == 0:
                last_canceled = self.last_canceled.copy()
            else:
                last_canceled = self.outer_slices[-1].last_canceled_node.copy()

            cancel_primes = [p for p in last_canceled if p > self.inner_slice.inner_wheel_prime]
            limit = math.isqrt(end)
            if DEBUG_MODE:
                debug(f"[BUILD] Target sqrt(end) = {limit}")

            seen = set(cancel_primes)
            pnode = self.get_prime_node(self.inner_slice.inner_wheel_prime)
            if pnode:
                pnode = self.get_next_prime(pnode.value)

            while pnode:
                p = pnode.value
                if p > limit:
                    break
                if p not in seen:
                    cancel_primes.append(p)
                    seen.add(p)
                    if DEBUG_MODE:
                        debug(f"[BUILD] Adding prime {p} to cancel_primes")
                pnode = self.get_next_prime(p)

            if DEBUG_MODE:
                debug(f"[BUILD] Final cancel_primes for slice {i}: {sorted(cancel_primes)}")

            slice_obj = PrimeronOuterSlice.__new__(PrimeronOuterSlice)
            self.outer_slices.append(slice_obj)

            slice_obj.__init__(
                parent_sieve=self,
                slice_index=i,
                slice_size=self.inner_wheel_size,
                start=start,
                end=end,
                outer_primes=[p for p in cancel_primes if p <= self.outer_wheel_prime],
                last_canceled_node=last_canceled,
                head_pseudoprime=self.inner_slice.head_pseudoprime,
                pseudoprime_map=self.inner_slice.pseudoprime_map,
                cancel_primes=cancel_primes
            )

            start = end + 1
            elapsed = time.time() - start_time
            print(f"- Outer wheel slice {i+1}/{slice_count} constructed in "
                  f"{elapsed:.3f} seconds: start={start}, end={end}")

        if DEBUG_MODE:
            debug(f"[BUILD] Completed {len(self.outer_slices)} slices")

    def _add_wheel_layer(self, prime):
        """
        Records metadata for a new wheel layer defined by the given
        prime.

        This function computes the primorial (product of all primes
        <=  `prime`) and stores it as the wheel size associated with
        that prime. This metadata is used to track how the sieve
        expands across successive layers of increasingly large wheels
        (e.g., 2×3=6, 2×3×5=30, 2×3×5×7=210, etc.).

        The result is saved in:
          - `self.wheel_layers`: a list of dictionaries containing:
                - 'prime': the prime defining the wheel boundary
                - 'wheel_size': the corresponding primorial value
          - `self.wheel_size_by_prime`: a cache mapping prime to
             wheel_size

        These records are useful for:
          - Generating per-wheel reports
          - Computing the number of slices in larger wheels
          - Supporting modular arithmetic and residue class analysis

        Parameters:
            prime (int): The outer prime defining the current wheel
                         layer.  Must exist in the inner wheel’s
                         prime map.

        Notes:
            - If the primorial for this prime has already been
              computed, it reuses it.
            - If the prime is not in the inner wheel, this function
              silently skips it (returns without action). This
              ensures robustness during partial construction.
        """
        size = self.compute_wheel_size_for_prime(prime)
        if size > 0:
            self.wheel_layers.append({'prime': prime, 'wheel_size': size})
            self.wheel_size_by_prime[prime] = size

    def get_removed_by(self, n):
        """
        Returns the lowest prime that divides n, or:
          - None if it’s within range and not removed
          - -1 if n is a pseudoprime outside the sieved range
        """
        # First check inner slice.
        removed_by = self.inner_slice.get_removed_by(n)
        if removed_by is not None:
            return removed_by

        wheel_size = self.inner_slice.wheel_size
        if n < wheel_size:
            return None

        # Outer slice index.
        slice_index = n // wheel_size

        if slice_index >= len(self.outer_slices):
            # Simulate outer wheel wraparound.
            offset = n % self.outer_wheel_size
            simulated_slice_index = offset // self.inner_wheel_size
            wheel_pos = offset % self.inner_wheel_size

            simulated_slice = self.outer_slices[simulated_slice_index]

            if DEBUG_MODE:
                debug(f"checking n={n}, offset={offset}, "
                      f"slice_index={slice_index}, "
                      f"simulated_slice_index={simulated_slice_index}, "
                      f"wheel_pos={wheel_pos}")

            # Special case: slice 0 does not include early structural
            # cancellations.
            if simulated_slice_index == 0:
                for p in simulated_slice.outer_primes:
                    if n % p == 0:
                        return p

            for actual_n, node in simulated_slice.holes_structural.items():
                if actual_n % self.inner_wheel_size == wheel_pos:
                    return node.removed_by

            if wheel_pos in self.inner_slice.pseudoprime_map:
                return -1
            else:
                return None

        outer_slice = self.outer_slices[slice_index]
        if DEBUG_MODE:
            debug(f"n={n}, slice_index={slice_index}, "
                  f"outer_slices[{slice_index}].start={outer_slice.start}, "
                  f"wheel_size={wheel_size}")
        node = outer_slice.holes_structural.get(n) or outer_slice.holes_residual.get(n)
        return node.removed_by if node else None

    def validate_removal(self):
        """
        Validate entire inner and outer wheel was build correctly.
        """
        print("\nValidating through (outer_wheel_prime + 1) x outer wheel...")
        passed = 0
        failed = 0

        next_prime_node = self.get_next_prime(self.outer_wheel_prime)
        if next_prime_node is None:
            raise ValueError(f"Could not find next prime after {self.outer_wheel_prime}")

        # This actually goes 1 wheel farther than necessary, which tests the
        # first slice of the subsequent wheel.
        outer_limit = self.compute_wheel_size_for_prime(next_prime_node.value)

        wheel_size = self.outer_wheel_size
        total_wheels = outer_limit // wheel_size

        last_slice = -1
        for n in range(2, outer_limit + 1):
            wheel_num = n // wheel_size
            if wheel_num != last_slice:
                print(f"- Testing wheel {wheel_num+1} of {total_wheels+1}")
                last_slice = wheel_num

            removed_by = self.get_removed_by(n)

            # Independently compute expected smallest prime divisor.
            is_prime = True
            for p in range(2, int(n**0.5) + 1):
                if n % p == 0:
                    is_prime = False
                    lcd = p
                    break

            # Validate.
            if is_prime:
                if removed_by is not None:
                    if n > self.outer_wheel_size and removed_by == -1:
                        continue

                    print(f"FAIL: {n} removed_by {removed_by} (should be None, it's prime)")
                    failed += 1
                else:
                    passed += 1
            else:
                if removed_by != lcd:
                    if (n > self.outer_wheel_size and removed_by == -1 and
                        lcd > self.outer_wheel_prime):
                        continue

                    print(f"FAIL: {n} removed_by {removed_by} (expected {lcd})")
                    failed += 1
                else:
                    passed += 1

        print(f"\nValidation complete. Passed: {passed}, Failed: {failed}\n")

    def factorize(self, n):
        """
        Returns the list of prime factors of n using only this sieve's
        get_removed_by() mechanism. Does not use trial division.

        Assumes that all relevant factors are available through the sieve.
        If n is prime and uncanceled, it returns [n].

        Parameters:
            n (int): Number to factor

        Returns:
            List[int]: List of prime factors (with multiplicity), e.g., 90 → [2, 3, 3, 5]
        """
        if n < 2:
            return []

        orig_n = n
        factors = []

        while True:
            lcd = self.get_removed_by(n)
            if lcd is None or lcd <= 1:
                break
            factors.append(lcd)
            if n % lcd != 0:
                # This should never happen if the sieve is correct.
                raise ValueError(f"Inconsistent sieve removal: {n} not divisible by {lcd}")
            n //= lcd

        # If remainder is not 1, n is a surviving prime.
        if n > 1:
            factors.append(n)

        # Sanity check.
        product = 1
        for f in factors:
            product *= f
        assert product == orig_n, f"Factorization of {orig_n} failed: {factors}: product {product}"

        return factors

    def generate_light_per_wheel_reports(self, output_dir):
        """
        Lightweight version of report generation that only tracks the count
        of removed_by values (i.e., holes) for each prime within each wheel.

        This avoids memory-heavy data structures and is optimized for estimating
        storage requirements of larger wheels.

        Output: For each wheel level, writes a CSV file with two columns:
            - prime: the prime that removed the value
            - count: how many times it appears as the lowest prime divisor
        """
        node = self.inner_slice.head_prime

        while node and node.value <= self.outer_wheel_prime:
            p = node.value
            if p < 3:
                node = node.next
                continue

            print(f"Building light report for p={p}.")
            start_time = time.time()
            wheel_size = self.compute_wheel_size_for_prime(p)
            if wheel_size <= 0:
                node = node.next
                continue
            end = wheel_size - 1
            output_path = os.path.join(output_dir, f"light-report-{p:02d}.csv")

            total_count = 0
            msg_count = 10000000
            count_by_prime = defaultdict(int)
            for n in range(0, end + 1):
                removed_by = self.get_removed_by(n)
                if removed_by and removed_by > 1:
                    count_by_prime[removed_by] += 1
                if p > 19:
                    msg_count -= 1
                    if msg_count == 0:
                        total_count += 10000000
                        msg_count = 10000000
                        print(f"  - Processed {n+1} of {end}")

            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            with open(output_path, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["prime", "count"])
                for prime in sorted(count_by_prime):
                    writer.writerow([prime, count_by_prime[prime]])

            elapsed = time.time() - start_time
            print(f"Light report for p={p} constructed in {elapsed:.3f} seconds.")
            print(f"- Written to: {output_path}")

            node = node.next

    def generate_unified_slice_report(
        self,
        slice_start,
        slice_end,
        slice_size,
        output_path,
        wheel_size,
        wheel_sqrt,
        report_label="outer"
    ):
        """
        Generates a detailed CSV report for a numeric slice, analyzing all
        holes (composite numbers) canceled by each prime in that slice.

        For each prime `p` that cancels at least one value in the range
        [`slice_start`, `slice_end`], the report includes:
          - Number and density of values canceled by `p` alone (holes_by_p)
          - Number and density of values canceled by any prime <= `p`
            (holes_up_to_p)
          - Distribution of canceled values across early/middle/late parts of
            the slice
          - First/last values canceled by `p` and up to `p`
          - Gap statistics (avg/min/max gaps between holes)
          - Depth of factorization (i.e., how many prime factors each canceled
            number has)
          - Number of q values used to generate holes (q = n / p)
          - Residue class of p mod wheel_size

        This allows rich empirical analysis of:
          - Prime gap patterns and survivor clustering
          - Factor depth distributions per prime
          - Temporal dynamics (early/middle/late cancellation bias)
          - Structured residue class behavior

        Parameters:
            slice_start (int):  Start of numeric range (inclusive)
            slice_end (int):    End of numeric range (inclusive)
            slice_size (int):   Size of the slice
                                (slice_end - slice_start + 1)
            output_path (str):  Destination file path for the CSV report
            wheel_size (int):   Size of the wheel (for computing residue
                                classes)
            wheel_sqrt (int):   Square root of wheel size (used for metadata
                                only)
            report_label (str): Label used to tag the report rows (e.g.,
                                "outer", "wheel")

        Output:
            Writes a CSV file with one row per canceling prime.  Each row
            contains ~40 columns including hole counts, densities, prime
            metadata, factor depth histograms, and gap statistics.

        Returns:
            output_path (str): The path to the written CSV file

        Notes:
            - Surviving primes are treated as removed_by=None.
            - Values with removed_by == -1 are treated as uncanceled
              pseudoprimes (i.e., survivors beyond sieve range).
            - Factor depth is computed only for values actually removed
              in-range.
            - Early/middle/late positions are computed as thirds of the slice.
            - Columns `count_N_holes_by_p` and `count_N_holes_up_to_p` report
              how many canceled values have exactly N prime factors.

        This function is optimized for post-hoc analysis and is used in wheel
        report generation, gap modeling, and twin-prime residue class theory.
        """
        max_factor_depth = math.floor(math.log2(slice_end)) + 1
        prime_slice_rows = []

        # Precompute all removed_by and factor depth for the range.
        removed_by_map = {}
        factor_depth_map = {}
        for n in range(slice_start, slice_end + 1):
            removed_by = self.get_removed_by(n)
            # Treat survivor primes as removed_by = None, not 1.
            removed_by_map[n] = removed_by
            factor_depth_map[n] = 0
            if removed_by and removed_by > 1:
                factor_depth_map[n] = len(self.factorize(n))

        all_cancel_primes = sorted(set(p for p in removed_by_map.values() if p is not None))

        for p in all_cancel_primes:
            if p is None or p <= 1:
                continue

            holes_by_p = []
            holes_up_to_p = []
            q_values_used = set()
            factor_counts_by_p = defaultdict(int)
            factor_counts_up_to_p = defaultdict(int)

            for n in range(slice_start, slice_end + 1):
                removed_by = removed_by_map[n]

                if removed_by == p:
                    holes_by_p.append(n)
                    if p > 1:
                        q_values_used.add(n // p)

                    if n in factor_depth_map:
                        depth = factor_depth_map[n]
                        factor_counts_by_p[depth] += 1
                        factor_counts_up_to_p[depth] += 1
                    holes_up_to_p.append(n)

                elif removed_by is not None and removed_by < p:
                    holes_up_to_p.append(n)
                    if n in factor_depth_map:
                        factor_counts_up_to_p[factor_depth_map[n]] += 1

            def compute_gap_stats(holes):
                if not holes or len(holes) < 2:
                    return 0, 0, 0
                gaps = [holes[i + 1] - holes[i] for i in range(len(holes) - 1)]
                return sum(gaps) / len(gaps), min(gaps), max(gaps)

            def early_middle_late_counts(holes):
                early_cutoff = 0.333 * slice_size
                late_cutoff = 0.667 * slice_size

                early = middle = late = 0
                for h in holes:
                    distance = h - slice_start
                    if distance <= early_cutoff:
                        early += 1
                    elif distance >= late_cutoff:
                        late += 1
                    else:
                        middle += 1

                return early, middle, late

            avg_gap_by_p, min_gap_by_p, max_gap_by_p = compute_gap_stats(holes_by_p)
            avg_gap_up_to_p, min_gap_up_to_p, max_gap_up_to_p = compute_gap_stats(holes_up_to_p)
            early_by_p, middle_by_p, late_by_p = early_middle_late_counts(holes_by_p)
            early_up_to_p, middle_up_to_p, late_up_to_p = early_middle_late_counts(holes_up_to_p)

            row = {
                "report_label": report_label,
                "slice_start": slice_start,
                "slice_end": slice_end,
                "wheel_size": wheel_size,
                "wheel_sqrt": wheel_sqrt,
                "prime": p,
                "holes_by_p": len(holes_by_p) if p > 1 else 0,
                "hole_density_by_p": len(holes_by_p) / slice_size if p > 1 else 0,
                "early_holes_by_p": early_by_p,
                "middle_holes_by_p": middle_by_p,
                "late_holes_by_p": late_by_p,
                "early_ratio_by_p": early_by_p / len(holes_by_p) if holes_by_p else 0,
                "late_ratio_by_p": late_by_p / len(holes_by_p) if holes_by_p else 0,

                "holes_up_to_p": len(holes_up_to_p) if p > 1 else 0,
                "hole_density_up_to_p": len(holes_up_to_p) / slice_size if p > 1 else 0,
                "early_holes_up_to_p": early_up_to_p,
                "middle_holes_up_to_p": middle_up_to_p,
                "late_holes_up_to_p": late_up_to_p,
                "early_ratio_up_to_p": early_up_to_p / len(holes_up_to_p) if holes_up_to_p else 0,
                "late_ratio_up_to_p": late_up_to_p / len(holes_up_to_p) if holes_up_to_p else 0,

                "residue_mod_wheel": p % wheel_size,
                "pseudoprimes_used": len(q_values_used),
                "first_hole_by_p": min(holes_by_p) if holes_by_p else None,
                "last_hole_by_p": max(holes_by_p) if holes_by_p else None,
                "first_hole_up_to_p": min(holes_up_to_p) if holes_up_to_p else None,
                "last_hole_up_to_p": max(holes_up_to_p) if holes_up_to_p else None,
            }

            row.update({
                "avg_gap_by_p": avg_gap_by_p,
                "min_gap_by_p": min_gap_by_p,
                "max_gap_by_p": max_gap_by_p,
                "avg_gap_up_to_p": avg_gap_up_to_p,
                "min_gap_up_to_p": min_gap_up_to_p,
                "max_gap_up_to_p": max_gap_up_to_p
            })

            for i in range(2, max_factor_depth + 1):
                row[f"count_{i}_holes_by_p"] = factor_counts_by_p[i]
                row[f"count_{i}_holes_up_to_p"] = factor_counts_up_to_p[i]

            prime_slice_rows.append(row)

        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        base_fields = [
            "report_label", "slice_start", "slice_end", "wheel_size",
            "wheel_sqrt", "prime", "holes_by_p", "hole_density_by_p",
            "early_holes_by_p", "middle_holes_by_p", "late_holes_by_p",
            "early_ratio_by_p", "late_ratio_by_p", "holes_up_to_p",
            "hole_density_up_to_p", "early_holes_up_to_p",
            "middle_holes_up_to_p", "late_holes_up_to_p",
            "early_ratio_up_to_p", "late_ratio_up_to_p",
            "residue_mod_wheel", "pseudoprimes_used",
            "first_hole_by_p", "last_hole_by_p", "first_hole_up_to_p",
            "last_hole_up_to_p", "avg_gap_by_p", "min_gap_by_p",
            "max_gap_by_p", "avg_gap_up_to_p", "min_gap_up_to_p",
            "max_gap_up_to_p"
        ]

        count_fields = [f"count_{i}_holes_by_p" for i in range(2, max_factor_depth + 1)] + \
                       [f"count_{i}_holes_up_to_p" for i in range(2, max_factor_depth + 1)]

        fieldnames = base_fields + count_fields

        with open(output_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(prime_slice_rows)

        return output_path

    def generate_per_wheel_reports(self, output_dir):
        """
        Generate per wheel reports.
        """
        # Gather primes from inner slice up to and including outer_wheel_prime.
        primes = []
        node = self.inner_slice.head_prime
        while node and node.value <= self.outer_wheel_prime:
            if node.value < 3:
                node = node.next
                continue

            primes.append(node.value)
            node = node.next

        for i, p in enumerate(primes):
            start_time = time.time()
            wheel_size = self.compute_wheel_size_for_prime(p)
            if wheel_size <= 0:
                continue
            end = wheel_size - 1
            sqrt_wheel = int(wheel_size ** 0.5)

            output_path = os.path.join(output_dir, f"wheel-report-{p:02d}.csv")
            self.generate_unified_slice_report(
                slice_start=0,
                slice_end=end,
                slice_size=wheel_size,
                output_path=output_path,
                wheel_size=wheel_size,
                wheel_sqrt=sqrt_wheel,
                report_label="wheel"
            )

            elapsed = time.time() - start_time
            print(f"Wheel report for p={p} constructed in {elapsed:.3f} seconds.")
            print(f"- Written to: {output_path}")

def main():
    start_time = time.time()
    sieve = PrimeronSieve(inner_wheel_prime=INNER_WHEEL_PRIME, outer_wheel_prime=OUTER_WHEEL_PRIME)
    elapsed = time.time() - start_time
    print(f"All wheels constructed in {elapsed:.3f} seconds\n")

    if VERBOSE_MODE:
        sieve.inner_slice.validate_removal()
        sieve.validate_removal()

    inner_wheel_size = sieve.inner_slice.wheel_size
    inner_max = sieve.inner_slice.inner_wheel_prime
    outer_max = sieve.outer_wheel_prime
    output_dir = f"primeron-reports-{inner_max:02d}-{outer_max:02d}"

    if REPORT_MODE:
        sieve.generate_per_wheel_reports(output_dir)
        print("")

    if LIGHT_REPORT_MODE:
        sieve.generate_light_per_wheel_reports(output_dir)
        print("")

if __name__ == "__main__":
    if PROFILE_MODE:
        profiler = cProfile.Profile()
        profiler.enable()
        main()
        profiler.disable()
        s = io.StringIO()
        ps = pstats.Stats(profiler, stream=s).sort_stats(pstats.SortKey.CUMULATIVE)
        ps.print_stats()
        print("\n--- PROFILER OUTPUT ---")
        print(s.getvalue())
    else:
        main()

