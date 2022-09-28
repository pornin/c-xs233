# xsb233 and xsk233 Groups

This is a demonstration C implementation of the xsb233 and xsk233
groups. These are prime order groups isomorphic to subgroups of the
standard NIST elliptic curves B-233 and K-233, as defined in
[FIPS 186-4](https://csrc.nist.gov/publications/detail/fips/186/4/final)
(section D.1.3.2), also known as sect233r1 and sect233k1, respectively,
in [SEC 2](https://www.secg.org/sec2-v2.pdf) (section 3.3). The
curve B-233 has order 2\*r for some prime integer r, while K-233 has
order 4\*r' for another prime integer r'; the groups xsb233 and xsk233
have order r and r', respectively. Groups xsb233 and xsk233 are deemed
to be safe for constructing cryptographic protocols, with a security
level of (at least) 112 bits.

For each group, functions for encoding/decoding, comparison, addition,
subtraction, negation, doubling and multiplication by a scalar (both for
the conventional generator, and for a dynamically obtained element) are
provided. Some highlights:

  - All code is strictly **constant-time** (execution time and memory
    access pattern do not vary in relation to the data values).

  - Encoding and decoding are **canonical**. The encoding format uses a
    fixed size (30 bytes). Each group element (including the neutral)
    has a unique valid encoding, and only that exact encoding can be
    successfully decoded into the group element.

  - All group operations use **complete formulas**, thus with no special
    case when dealing with the neutral or when adding a point to itself.

  - Formulas are **efficient**; in fact, the general point addition
    formuals are not only complete, but also faster than previously
    known (incomplete) formulas (point addition is computed in cost
    8M+2S+2mb, or 7M+2S+2mb for about half of the possible curves,
    including xsk233).

## Compilation

See the [Makefile](./Makefile) for options. Compilation produces for
each group a stand-alone, linkable object file (`xsb233.o` and
`xsk233.o`), which implements the API described in
[`xs233.h`](./xs233.h). Test and benchmark utilities are also compiled.

On x86 platforms that support the `pclmul` instructions (i.e. basically
all x86 CPUs from the last 10 years), special code is used that
dramatically improves performance. If `pclmul` is not enabled, then
generic portable code is used (and is about 10 times slower).

The code assumes availability of a 128-bit unsigned integer type under
the name `unsigned __int128`, which normally means that 64-bit mode is
used (32-bit mode is not supported). This has been tested with GCC 11.2.0
and Clang-14.0.0 on a Linux x86 system (Ubuntu 22.04), as well as
Clang-10.0.0 on a Linux ARMv8 system (Ubuntu 20.04).

The benchmark code (`speed_gf`, `speed_xs233`) will compile only on x86
and ARMv8 platforms, because it uses intrinsics and/or inline assembly
to access the cycle counter. Moreover, on ARMv8 systems, the cycle
counter is *usually* inaccessible unless explicitly enabled from kernel
code; on Linux, one can use a custom module that enables such access,
e.g. [this
one](https://github.com/jerinjacobk/armv8_pmu_cycle_counter_el0).

## Benchmarks

On an Intel i5-8259U (Coffee Lake), running at 2.3 GHz (TurboBoost
disabled), the following performance is achieved (in clock cycles):

| Group operation           |     xsb233 |     xsk233 |
| :------------------------ | ---------: | ---------: |
| mul                       |      60621 |      49378 |
| mul (ladder)              |      51537 |      46365 |
| mul (endomorphism)        |          - |  **29602** |
| mulgen                    |      25754 |      20084 |
| mulgen (endomorphism)     |          - |      16546 |

The "mul" operation is multiplication of a point by a scalar; the point
is obtained dynamically, so no precomputed tables apply. This would
correspond to, for instance, an ECDH key exchange (when receiving the
public element from the peer). The "ladder" version of the "mul"
operation implements the same functionality using internally an X-only
Montgomery ladder, as described by [López and
Dahab](https://link.springer.com/content/pdf/10.1007/3-540-48059-5_27.pdf)
in 1999; the ladder version still outputs the complete result, not just
the X coordinate. The "endomorphism" version leverages the Frobenius
endomorphism on the underlying Koblitz curve, for improved performance;
this applies only to curve K-233, hence group xsk233.

"mulgen" is similar to "mul" except that the multiplication applies to
the conventional generator element, for which some precomputed tables
were generated, hence the much improved performance (compared to "mul").
This operation would be the core of, for instance, a new key pair
generation, or a signature generation.
