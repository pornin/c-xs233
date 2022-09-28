/*
 * This file implements operations on the finite field GF(2^233), as
 * used by NIST curves B-233 and K-233, and groups xsb233 and xsk233.
 * It is meant to be included by the curve implementation; all
 * relevant types and functions are 'static' and inlinable.
 */

#include <stddef.h>
#include <stdint.h>
#include <string.h>

/*
 * If GF233_INV_FLT is non-zero, then inversions are computed with
 * Fermat's little theorem, with some optimization of sequences of
 * squarings (aka Itoh-Tsujii's method).
 *
 * If GF233_INV_FLT is zero, then an alternate implementation using
 * a binary GCD variant (Brunner-Curiger-Hofstetter algorithm) is
 * used. On systems where multiplications are fast, the FLT/Itoh-Tsujii
 * method is faster, but requires includes of about 15 kB worth of
 * static data; the binary GCD variant thus yields a smaller binary.
 * It might also deliver better performance on architectures where
 * multiplications in the field are slow.
 *
 * Default setting of GF233_INV_FLT is 1.
 */
#ifndef GF233_INV_FLT
#define GF233_INV_FLT   1
#endif

/*
 * If GF233_PCLMUL is non-zero, then an implementation leveraging
 * the pclmul opcode and SSE4.1 instrinsics will be used. It is
 * substantially faster than the non-pclmul code.
 *
 * If GF233_PCLMUL is zero, then a portable 64-bit implementation will
 * be used.
 *
 * Default setting is to use the pclmul code whenever the compiler
 * advertises support for both pclmul and SSE4.1 in the target system.
 */
#ifndef GF233_PCLMUL
#if (defined __PCLMUL__ && __PCLMUL__) \
	&& (defined __SSE4_1__ && __SSE4_1__)
#define GF233_PCLMUL   1
#else
#define GF233_PCLMUL   0
#endif
#endif

#if GF233_PCLMUL
/* ======================================================================= */
/*
 * SSE4.1 + PCLMUL implementation.
 */

#include <immintrin.h>

/*
 * Field: GF(2^233)
 * Modulus: z^233 + z^74 + 1
 */

typedef union {
	uint64_t w[4];
	__m128i v[2];
} gf;

/*
 * Macro for an initializer for a field element. Values are provided as
 * four 64-bit limbs that encode the value in big-endian order.
 */
#define GFw64be(w3, w2, w1, w0)  { { (w0), (w1), (w2), (w3) } }

/*
 * Helper macro for defining initializers for constant point that use
 * the external representation (with values being opaque sequences of
 * 64-bit integers).
 */
#define GFext(w3, w2, w1, w0)  (w0), (w1), (w2), (w3)

/*
 * Helper macros for tests: masks to apply to random 64-bit limbs
 * in order to generate a random field element.
 */
#define GF_TEST_RMASK0    0xFFFFFFFFFFFFFFFF
#define GF_TEST_RMASK1    0xFFFFFFFFFFFFFFFF
#define GF_TEST_RMASK2    0xFFFFFFFFFFFFFFFF
#define GF_TEST_RMASK3    0xFFFFFFFFFFFFFFFF

/*
 * Multiplications are relatively fast (i.e. they cost less than two
 * squarings).
 */
#define GF_MUL_FAST   1

static const gf GF_ZERO = GFw64be(0, 0, 0, 0);
static const gf GF_ONE = GFw64be(0, 0, 0, 1);

/* d <- a + b */
static inline void
gf_add(gf *d, const gf *a, const gf *b)
{
	d->v[0] = _mm_xor_si128(a->v[0], b->v[0]);
	d->v[1] = _mm_xor_si128(a->v[1], b->v[1]);
}

/*
 * Partial reduction of a 512-bit value (in d0:d1:d2:d3).
 * Output is 256 bits, and is written back in d0:d1; d2:d3 are
 * unmodified.
 */
#define REDUCE_PARTIAL_512   do { \
		__m128i gsl23, gsr41, gsl33, gsr31; \
		__m128i fsl23, fsr41, fsl33, fsr31; \
		__m128i t1, v2; \
 \
		/* Fold d3. */ \
		gsl23 = _mm_slli_epi64(d3, 23); \
		gsr41 = _mm_srli_epi64(d3, 41); \
		gsl33 = _mm_slli_epi64(d3, 33); \
		gsr31 = _mm_srli_epi64(d3, 31); \
		v2 = _mm_xor_si128(_mm_xor_si128(d2, gsr31), \
			_mm_bsrli_si128(_mm_xor_si128(gsr41, gsl33), 8)); \
 \
		/* Fold v2. */ \
		fsl23 = _mm_slli_epi64(v2, 23); \
		fsr41 = _mm_srli_epi64(v2, 41); \
		fsl33 = _mm_slli_epi64(v2, 33); \
		fsr31 = _mm_srli_epi64(v2, 31); \
		d0 = _mm_xor_si128(_mm_xor_si128(d0, fsl23), \
			_mm_bslli_si128(_mm_xor_si128(fsr41, fsl33), 8)); \
		t1 = _mm_alignr_epi8(_mm_xor_si128(gsr41, gsl33), \
			_mm_xor_si128(fsr41, fsl33), 8); \
		d1 = _mm_xor_si128(_mm_xor_si128(d1, fsr31), \
			_mm_xor_si128(gsl23, t1)); \
	} while (0)

/* d <- a * b */
static inline void
gf_mul(gf *d, const gf *a, const gf *b)
{
	/*
	 * We apply one level of Karatsuba:
	 *   E <- a0*b0
	 *   F <- a1*b1
	 *   G <- (a0 + a1)*(b0 + b1)
	 *   D <- E + (G + E + F)*z^128 + F*z^256
	 * with D being the unreduced result (over 511 bits)
	 *
	 * Down to individual elements (denoting xl and xh the low and
	 * high 64-bit halves of x, respectively):
	 *   (e0, e1, e2) <- (a0l*b0l, a0h*b0h, a0l*b0h + a0h*b0l)
	 *   (f0, f1, f2) <- (a1l*b1l, a1h*b1h, a1l*b1h + a1h*b1l)
	 *   (u, v) <- (a0 + a1, b0 + b1)
	 *   (g0, g1, g2) <- (ul*vl, uh*vh, ul*vh + uh*vl)
	 *   d0l <- e0l
	 *   d0h <- e0h ^ e2l
	 *   d1l <- e1l ^ e2h ^ g0l       ^ e0l       ^ f0l
	 *   d1h <- e1h       ^ g0h ^ g2l ^ e0h ^ e2l ^ f0h ^ f2l
	 *   d2l <- f0l       ^ g1l ^ g2h ^ e1l ^ e2h ^ f1l ^ f2h
	 *   d2h <- f0h ^ f2l ^ g1h       ^ e1h       ^ f1h
	 *   d3l <- f1l ^ f2h
	 *   d3l <- f1h
	 */
	__m128i a0, a1, b0, b1, aa, bb, d0, d1, d2, d3;
	__m128i e0, e1, e2, f0, f1, f2, g0, g1, g2, hh, jj;

	a0 = a->v[0];
	a1 = a->v[1];
	b0 = b->v[0];
	b1 = b->v[1];

	/* (e0, e1, e2) <- (a0l*b0l, a0h*b0h, a0l*b0h + a0h*b0l) */
	e0 = _mm_clmulepi64_si128(a0, b0, 0x00);
	e1 = _mm_clmulepi64_si128(a0, b0, 0x11);
	e2 = _mm_xor_si128(
		_mm_clmulepi64_si128(a0, b0, 0x01),
		_mm_clmulepi64_si128(a0, b0, 0x10));

	/* (f0, f1, f2) <- (a1l*b1l, a1h*b1h, a1l*b1h + a1h*b1l) */
	f0 = _mm_clmulepi64_si128(a1, b1, 0x00);
	f1 = _mm_clmulepi64_si128(a1, b1, 0x11);
	f2 = _mm_xor_si128(
		_mm_clmulepi64_si128(a1, b1, 0x01),
		_mm_clmulepi64_si128(a1, b1, 0x10));

	/* (u, v) <- (a0 + a1, b0 + b1)
	   (g0, g1, g2) <- (ul*vl, uh*vh, ul*vh + uh*vl) */
	aa = _mm_xor_si128(a0, a1);
	bb = _mm_xor_si128(b0, b1);
	g0 = _mm_clmulepi64_si128(aa, bb, 0x00);
	g1 = _mm_clmulepi64_si128(aa, bb, 0x11);
	g2 = _mm_xor_si128(
		_mm_clmulepi64_si128(aa, bb, 0x01),
		_mm_clmulepi64_si128(aa, bb, 0x10));

	/* unreduced: d0:d1:d2:d3 */
	hh = _mm_xor_si128(e1, f0);
	jj = _mm_xor_si128(g2, _mm_xor_si128(e2, f2));
	d0 = _mm_xor_si128(e0, _mm_bslli_si128(e2, 8));
	d1 = _mm_xor_si128(
		_mm_xor_si128(e0, hh),
		_mm_xor_si128(g0, _mm_alignr_epi8(jj, e2, 8)));
	d2 = _mm_xor_si128(
		_mm_xor_si128(hh, f1),
		_mm_xor_si128(g1, _mm_alignr_epi8(f2, jj, 8)));
	d3 = _mm_xor_si128(f1, _mm_bsrli_si128(f2, 8));

	REDUCE_PARTIAL_512;

	d->v[0] = d0;
	d->v[1] = d1;
}

/* d <- a^2 */
static inline void
gf_sqr(gf *d, const gf *a)
{
	__m128i a0, a1, d0, d1, d2, d3;

	a0 = a->v[0];
	a1 = a->v[1];

	/*
	 * Squaring is a field automorphism, we can just handle each
	 * word separately.
	 */
	d0 = _mm_clmulepi64_si128(a0, a0, 0x00);
	d1 = _mm_clmulepi64_si128(a0, a0, 0x11);
	d2 = _mm_clmulepi64_si128(a1, a1, 0x00);
	d3 = _mm_clmulepi64_si128(a1, a1, 0x11);

	/*
	 * Alternate code using pshufb to expand 4-bit chunks. On an
	 * Intel i5-i8259U (Coffee Lake), it seems slightly slower than
	 * using pclmul.

	__m128i m16 = _mm_set1_epi8(0x0F);
	__m128i shk = _mm_setr_epi8(
		0x00, 0x01, 0x04, 0x05, 0x10, 0x11, 0x14, 0x15,
		0x40, 0x41, 0x44, 0x45, 0x50, 0x51, 0x54, 0x55);
	__m128i t0, t1, t2, t3;

	t0 = _mm_shuffle_epi8(shk, _mm_and_si128(a0, m16));
	t1 = _mm_shuffle_epi8(shk, _mm_and_si128(_mm_srli_epi16(a0, 4), m16));
	t2 = _mm_shuffle_epi8(shk, _mm_and_si128(a1, m16));
	t3 = _mm_shuffle_epi8(shk, _mm_and_si128(_mm_srli_epi16(a1, 4), m16));
	d0 = _mm_unpacklo_epi8(t0, t1);
	d1 = _mm_unpackhi_epi8(t0, t1);
	d2 = _mm_unpacklo_epi8(t2, t3);
	d3 = _mm_unpackhi_epi8(t2, t3);
	 */

	REDUCE_PARTIAL_512;

	d->v[0] = d0;
	d->v[1] = d1;
}

/* Normalize a into d (i.e. full reduction in the field). */
static inline void
gf_norm(gf *d, const gf *a)
{
	__m128i a0, a1, u, v;

	a0 = a->v[0];
	a1 = a->v[1];
	u = _mm_srli_epi64(a1, 41);
	v = _mm_slli_epi64(u, 10);
	d->v[0] = _mm_xor_si128(a0, _mm_unpackhi_epi64(u, v));
	d->v[1] = _mm_and_si128(a1, _mm_set_epi32(
		0x000001FF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF));
}

/* d <- sqrt(a) */
static inline void
gf_sqrt(gf *d, const gf *a)
{
	/*
	 * pclmul allows computing products and thus squares in
	 * GF(2)[z], but there is no equivalent for square roots,
	 * where the shrinking must be done "manually".
	 *
	 * We first expand the source polynomial by spliting it into
	 * "odd" and "even" parts:
	 *   a = ae + z*ao
	 * with:
	 *   ae = \sum_{i=0}^{127} a_{2*i)*z^{2i}
	 *   ao = \sum_{i=0}^{127} a_{2*i+1)*z^{2i}
	 * Then, we have:
	 *   sqrt(a) = sqrt(ae) + sqrt(z)*sqrt(ao)
	 *
	 * Since ae and ao have non-zero coefficients only for even
	 * powers, their respective square roots can be computed with
	 * "shrinking" (moving bits to squeeze out the zeros), in a
	 * relatively low number of operations (logarithmic in the word
	 * size, with shuffling opcodes taking over once we get down to
	 * full bytes).
	 *
	 * In the field:
	 *   sqrt(z) = z^228 + z^191 + z^154 + z^117 + z^69 + z^32
	 * Note that sqrt(ae) and sqrt(ao) both fit on 128 bits each;
	 * we can thus optimize the multiplication of sqrt(ao) with
	 * sqrt(z).
	 */
	__m128i m1 = _mm_set_epi32(
		0x55555555, 0x55555555, 0x55555555, 0x55555555);
	__m128i m2 = _mm_set_epi32(
		0x33333333, 0x33333333, 0x33333333, 0x33333333);
	__m128i m3 = _mm_set_epi32(
		0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F);
	__m128i sklo = _mm_set_epi8(
		0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0x0E, 0x0C, 0x0A, 0x08, 0x06, 0x04, 0x02, 0x00);
	__m128i skhi = _mm_set_epi8(
		0x0E, 0x0C, 0x0A, 0x08, 0x06, 0x04, 0x02, 0x00,
		0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF);
	__m128i sz0 = _mm_set_epi32(
		0x00200000, 0x00000020, 0x00000001, 0x00000000);
	__m128i sz1 = _mm_set_epi32(
		0x00000010, 0x00000000, 0x80000000, 0x04000000);
	__m128i a0, a1, ae0, ae1, ao0, ao1;
	__m128i d0, d1, d2, e0, e1, e2;
	__m128i sl23, sr41, sl33, sr31, v0, v1, w0, w1;

	a0 = a->v[0];
	a1 = a->v[1];
	ae0 = _mm_and_si128(a0, m1);
	ae1 = _mm_and_si128(a1, m1);
	ao0 = _mm_and_si128(_mm_srli_epi64(a0, 1), m1);
	ao1 = _mm_and_si128(_mm_srli_epi64(a1, 1), m1);

	ae0 = _mm_and_si128(_mm_xor_si128(ae0, _mm_srli_epi64(ae0, 1)), m2);
	ae1 = _mm_and_si128(_mm_xor_si128(ae1, _mm_srli_epi64(ae1, 1)), m2);
	ao0 = _mm_and_si128(_mm_xor_si128(ao0, _mm_srli_epi64(ao0, 1)), m2);
	ao1 = _mm_and_si128(_mm_xor_si128(ao1, _mm_srli_epi64(ao1, 1)), m2);

	ae0 = _mm_and_si128(_mm_xor_si128(ae0, _mm_srli_epi64(ae0, 2)), m3);
	ae1 = _mm_and_si128(_mm_xor_si128(ae1, _mm_srli_epi64(ae1, 2)), m3);
	ao0 = _mm_and_si128(_mm_xor_si128(ao0, _mm_srli_epi64(ao0, 2)), m3);
	ao1 = _mm_and_si128(_mm_xor_si128(ao1, _mm_srli_epi64(ao1, 2)), m3);

	ae0 = _mm_xor_si128(ae0, _mm_srli_epi64(ae0, 4));
	ae1 = _mm_xor_si128(ae1, _mm_srli_epi64(ae1, 4));
	ao0 = _mm_xor_si128(ao0, _mm_srli_epi64(ao0, 4));
	ao1 = _mm_xor_si128(ao1, _mm_srli_epi64(ao1, 4));

	/*
	 * Now, each of ae0, ae1, ao0 and ao1 consists of full bytes (X*)
	 * with interleaved junk (--):
	 *   -- X7 -- X6 -- X5 -- X4 -- X3 -- X2 -- X1 -- X0
	 * We can finish with some pshufb opcodes.
	 */
	ae0 = _mm_xor_si128(
		_mm_shuffle_epi8(ae0, sklo),
		_mm_shuffle_epi8(ae1, skhi));
	ao0 = _mm_xor_si128(
		_mm_shuffle_epi8(ao0, sklo),
		_mm_shuffle_epi8(ao1, skhi));

	/*
	 * We have sqrt(ae) (in ae0) and sqrt(ao) (in ao0); we multiply
	 * the latter by sqrt(z), and add the values together.
	 * Since sqrt(z) has Hamming weight 6, doing the shifts and XORs
	 * "manually" would probably not yield any performance advantage,
	 * and it would be more complicated.
	 */

	/* d0:d1:d2 <- sqrt(z)*sqrt(ao), unreduced */
	d0 = _mm_clmulepi64_si128(ao0, sz0, 0x00);
	d1 = _mm_clmulepi64_si128(ao0, sz0, 0x11);
	d2 = _mm_xor_si128(
		_mm_clmulepi64_si128(ao0, sz0, 0x01),
		_mm_clmulepi64_si128(ao0, sz0, 0x10));
	d0 = _mm_xor_si128(d0, _mm_bslli_si128(d2, 8));
	d1 = _mm_xor_si128(d1, _mm_bsrli_si128(d2, 8));

	e0 = _mm_clmulepi64_si128(ao0, sz1, 0x00);
	e1 = _mm_clmulepi64_si128(ao0, sz1, 0x11);
	e2 = _mm_xor_si128(
		_mm_clmulepi64_si128(ao0, sz1, 0x01),
		_mm_clmulepi64_si128(ao0, sz1, 0x10));
	e0 = _mm_xor_si128(e0, _mm_bslli_si128(e2, 8));
	e1 = _mm_xor_si128(e1, _mm_bsrli_si128(e2, 8));

	d1 = _mm_xor_si128(d1, e0);
	d2 = e1;

	/*
	 * Reduction.
	 * This is made simpler by the fact that the intermediate
	 * product has size at most 356 bits, so that one folding is
	 * enough.
	 */
	sl23 = _mm_slli_epi64(d2, 23);
	sr41 = _mm_srli_epi64(d2, 41);
	v0 = _mm_xor_si128(sl23, _mm_bslli_si128(sr41, 8));
	v1 = _mm_bsrli_si128(sr41, 8);
	sl33 = _mm_slli_epi64(d2, 33);
	sr31 = _mm_srli_epi64(d2, 31);
	w0 = _mm_bslli_si128(sl33, 8);
	w1 = _mm_xor_si128(sr31, _mm_bsrli_si128(sl33, 8));

	d->v[0] = _mm_xor_si128(_mm_xor_si128(ae0, d0), _mm_xor_si128(v0, w0));
	d->v[1] = _mm_xor_si128(d1, _mm_xor_si128(v1, w1));
}

/*
 * Trace (always 0 or 1).
 */
static inline uint64_t
gf_trace(const gf *a)
{
	/*
	 * In the field, Tr(a) = a_0 + a_159.
	 * We have a non-normalized value, with 23 extra bits; normalization
	 * would XOR these 23 bits into our value at offsets 0 and 74; the
	 * net result is that bit a_233 must also be taken into account.
	 */
	return (a->w[0] ^ (a->w[2] >> 31) ^ (a->w[3] >> 41)) & 1;
}

/*
 * Precomputed half-traces of z^(2*i+1) for i = 0 to 115.
 */
static const gf HALFTRACE[] = {
	GFw64be(0x00000197FD722105, 0x5AABAEF4E0DC1096,
	        0x20B7FC13C12E8D69, 0xFD01D0208E170521),
	GFw64be(0x0000019395A41BE5, 0x109128B2E0851116,
	        0x9795725B45048071, 0x5C8585D4DA648BAD),
	GFw64be(0x0000010C2A9B2FEA, 0x1A2BEE86B8843672,
	        0x38D2053841A24074, 0x1902CF5E9217F50D),
	GFw64be(0x000000338A8713C4, 0x254C7AE360522FC4,
	        0x61A5035A3AC98C5D, 0x86126EAC1C0A1166),
	GFw64be(0x000001725848920E, 0xB0518410CE9A5281,
	        0xA2BCC8424D24A625, 0xB619963440E20137),
	GFw64be(0x0000010C68C348AA, 0x0A692C21421930C9,
	        0x2E92898009218751, 0xF596806D30877774),
	GFw64be(0x000001D6F5179967, 0x760DEE3DF1DD38A2,
	        0x323EF75BE3A7DCED, 0xACD2A1364347E213),
	GFw64be(0x000001825B64FFE5, 0x262772A288061B2A,
	        0x4C50DAFB23C00154, 0x1A04E7C6670CFD0D),
	GFw64be(0x0000012B062F72CE, 0x4F6A4CB8300C3A3E,
	        0x2D9123CAB8A44010, 0x3946A716EA5479F0),
	GFw64be(0x00000034C83C8D45, 0x00E1E8839B1E46D3,
	        0x34C6863B098156C4, 0xFA5C0B69829ACC79),
	GFw64be(0x0000006DC7D31CA5, 0xF71249CE72F775AB,
	        0xB82BBD63F49ACF2A, 0xCEDFDDF67383E042),
	GFw64be(0x000000B0F1628609, 0x42190624AC380C95,
	        0xA024D82185222546, 0xE0043C318D731517),
	GFw64be(0x000000B0B1608609, 0x41190624AC280C15,
	        0xA024582195222146, 0x60007C118BF31517),
	GFw64be(0x000001130E386A29, 0xAE073201BA7A15D8,
	        0x5DF5248123414E5A, 0xE351CDFDAD5CB8EB),
	GFw64be(0x000000C9AB45864C, 0xBDC8C6AC5AAF7402,
	        0x0F491B2A78A2C273, 0x6FDBCC0127156FD8),
	GFw64be(0x000000FBB1A6E508, 0x104F346FBA9D43BF,
	        0x4ACD52B04B6B4660, 0xE88953A73AC92769),
	GFw64be(0x000001E13FBF06E3, 0x76E0A4F906E37771,
	        0x2E597729E82D2A3F, 0x528F8ED9D0D63E0C),
	GFw64be(0x00000016CDB3CB21, 0x1BD9A49184B028EC,
	        0xA866B5915D252427, 0x911265623863754B),
	GFw64be(0x00000165AD0EFACF, 0xCE93CE63CB721128,
	        0xD0FB32CBA5A99E9B, 0x9754C4F3247DC02F),
	GFw64be(0x00000123A6B0E0EC, 0xC1FB0822A02C180A,
	        0xB511248A9D000007, 0x3841B1D17336D8C5),
	GFw64be(0x000000E85F8DE429, 0x783A080144BB6E9C,
	        0x8AE8F3A1C401A476, 0xB7CA6960FD60266A),
	GFw64be(0x000000B3181963A8, 0xBA9E305E3ABF298A,
	        0xE2C54590464E4676, 0xBB866DE5633A023C),
	GFw64be(0x000001E47516C9C3, 0x6201E697FC3F58EE,
	        0x323AF69981A7E402, 0xAEC8803E6F57D347),
	GFw64be(0x000000826248024D, 0x410040135A1512A3,
	        0x1240880B9085C645, 0x8C91472303959758),
	GFw64be(0x0000013BA26A572E, 0x0A0462E535277F3A,
	        0x32F508F2028B7183, 0x0BD8AF82A616836A),
	GFw64be(0x0000016C0E09236D, 0x52265C3D20DC75A2,
	        0x743A211BC2E70C28, 0xF807CDF6035E8C12),
	GFw64be(0x000001C279C85167, 0x79CB29E68E72CBE4,
	        0x9718D8DBD91A2F1E, 0x875D76EA1874CB85),
	GFw64be(0x00000069B4970E05, 0x849430341E5C6BD4,
	        0xF189672326466E4D, 0xB89F36ADDC2A90A0),
	GFw64be(0x000000578EC33383, 0x052F1EA0F8516F95,
	        0x536F29513360CD18, 0xC69B3FB1DD99C69F),
	GFw64be(0x000000F0FE1C6ACB, 0x758BA8B8C6B7439F,
	        0x06ECE689F104A676, 0xD69907A1EAD04A2F),
	GFw64be(0x0000010C68C348AA, 0x0A692C21C21930C9,
	        0x2A92898009218751, 0xD596806D30877775),
	GFw64be(0x000000E4DDDE14CB, 0x8DE55C8CAF4E0912,
	        0x180AFE693BE23BDD, 0xFB5571D19200A151),
	GFw64be(0x000001C39BCAF3A7, 0xA703860E8E0A5A19,
	        0x313958D331222342, 0x6258B744B5C6D5C3),
	GFw64be(0x000000981724FB4E, 0xD44AA4184E85097B,
	        0x26A472DAE824A22D, 0x0C85748AA2924B26),
	GFw64be(0x0000003B63169906, 0x8AB90AF4E65C00C6,
	        0xA4C59652050EAE7C, 0xED41002D0F325D69),
	GFw64be(0x00000121CD9AFF4D, 0xA59E1CD35A2A5029,
	        0xCC51B4FB366DC282, 0x3658951270AD6C58),
	GFw64be(0x000001786534E22D, 0x87FF3875DE0C5F1A,
	        0xDD3CB6833F4FE004, 0x3C08BFD5E62CEC83),
	GFw64be(0x0000011F867DDBAB, 0x6A8C945E3C3F3DDA,
	        0x66972FD1826E6543, 0xABC6BDE8E21F5A74),
	GFw64be(0x000001636BB8712E, 0x651A800A26325289,
	        0x913994D2B4002607, 0x821C827530B4C182),
	GFw64be(0x000001572982BF0B, 0xD0BF9A4BD49C53D6,
	        0xCCFF1071C749E470, 0xFC48C7B99A2D687F),
	GFw64be(0x000001B1D1E3EC21, 0x036B0841E85434EE,
	        0x2CF5D9A119098D09, 0x8C16896A3B066C6B),
	GFw64be(0x000001E6F7B35B62, 0xEA886092103F78CE,
	        0x1A9AF5D880844407, 0xAB9AB47D7B14E324),
	GFw64be(0x0000000000000000, 0x0000000080000000,
	        0x0000000000000000, 0x0000000000000001),
	GFw64be(0x0000016AC624EEAC, 0xC2CFBC9EB6005339,
	        0x48F8A2A28B666200, 0x0049C7D6A5DC702B),
	GFw64be(0x00000179EAC21AED, 0x366846CD329B5895,
	        0x169D884B68AB4734, 0xB28DB461DDC5DA64),
	GFw64be(0x00000166C5BCAEEE, 0xB8D20C564AFA5FB5,
	        0xB5FAB62A4C2E8E3A, 0xF719FFF299A799BA),
	GFw64be(0x000000D015E2998E, 0x8D04B40F34504FC1,
	        0x6CCC785232636D19, 0x81097AB8109A6D28),
	GFw64be(0x000001F7B967ADA2, 0x5826184E789C73F6,
	        0x605F5B30C24AC531, 0xFD1BD3EA8E4E501C),
	GFw64be(0x000001270C10A70C, 0x18B2A8D04CE41C02,
	        0x80932432440CA82E, 0x1D05EC1003651120),
	GFw64be(0x000001AB5DBF1A02, 0x75DCB2007CEF3542,
	        0xF311F740FE40E82E, 0x6E86DC98173F8680),
	GFw64be(0x000000D095F2998E, 0x0904B40F34704FC5,
	        0x6CCC7C5212636D1A, 0x81097AB9189A6D28),
	GFw64be(0x0000018B2020C7C7, 0x1D982C5C62C80D65,
	        0xA4B100BB742E8B2D, 0x65446D8B08B30966),
	GFw64be(0x0000011A8E19E3EB, 0x66AF9096BC7E34DE,
	        0x64942599A3466C5B, 0xBA06DD7CEA1F1965),
	GFw64be(0x00000069A4978E05, 0x849434341E5C6BD4,
	        0xF589272326666E4D, 0xB81E36ADDC2A98A0),
	GFw64be(0x0000008EDA97306F, 0xCF930089CA3521D7,
	        0x8EC2C74BB5018602, 0x9D9711E89AF13E6D),
	GFw64be(0x000001652FCB0FAF, 0x5241485EB1DF35A7,
	        0x5EBB3933C98E5DB9, 0xBAC38DA31BCDEB63),
	GFw64be(0x00000008A09F9222, 0x8E91449A90202334,
	        0x9880074025A44007, 0x010716C39960A135),
	GFw64be(0x000001692F5D0A2D, 0x3B4C9AEE209774A6,
	        0x42193F035A4A0530, 0xDBCBDC624B1D4714),
	GFw64be(0x000000887B0FFAAE, 0x30273AAADCE32D4C,
	        0x4600D3C24340E82A, 0x56C76DC86C085F55),
	GFw64be(0x0000005B031990E1, 0xCA9474CD6C1865A4,
	        0xF4AD154986EBA404, 0xE54A08E7083A9D72),
	GFw64be(0x00000094AFFBC96C, 0x5AD9E482DEA52A4F,
	        0xBA263D9ACDA0E222, 0x1D8322493BF6E213),
	GFw64be(0x000000BC333377CB, 0xB4C252A347DD6B9B,
	        0x72C655F968C1BEED, 0xAC9B66E0E69B9238),
	GFw64be(0x0000004DC393DC20, 0xF7B840EE7AF6658B,
	        0xB86B95E0F48ACE2A, 0xCE5F58E473F3E04A),
	GFw64be(0x0000019A7EB4B107, 0xAF3451800A6B065D,
	        0x5874E65336D00A1B, 0x33C40E1CFDCDB05E),
	GFw64be(0x0000011BECDF9C0B, 0x34C5B424D0C73955,
	        0x5335AF616B62C968, 0x0E86A0899CCC9397),
	GFw64be(0x000001428B65418B, 0xA050EAF03867796F,
	        0xABF81B910C8C495E, 0x0A9FA18A2FA666BE),
	GFw64be(0x000001889B2888E6, 0xF27490CD18BF1C85,
	        0xDBD0500ACE4B4472, 0xEA91B8740CFCA6AC),
	GFw64be(0x0000000082C5F560, 0x9CD2863D7AAD240F,
	        0x8E400BF86C27C322, 0x2D9749503FB02B4C),
	GFw64be(0x0000018170D0EE04, 0x437692EF642C183A,
	        0xC8F1CCA29E4BA146, 0x6C01B502B67D617E),
	GFw64be(0x000001BF760F5922, 0xA895925A9027354A,
	        0xC697E3D0074C4006, 0x4B8288D9676D4E35),
	GFw64be(0x00000187915920E5, 0xBEC6304672FA34A1,
	        0x40335D0B6A4ACF6E, 0xA743C962408D1152),
	GFw64be(0x00000058077CE060, 0xD6FF786684194DB9,
	        0xCF6C3E88EFCA2505, 0xB49978F6E5B97BDB),
	GFw64be(0x0000011505A67208, 0xA57FF450C03416FA,
	        0xF89732C03FEC8446, 0xCC049E3BF23FA021),
	GFw64be(0x0000010CEA06BDCA, 0x96BBAA1CB8BC14C6,
	        0xA4D2827865064473, 0xF801C93D0F375C39),
	GFw64be(0x0000007CEDC5AA85, 0x0F9614541E0D6F53,
	        0xF7AEBB03366E6344, 0x29DA3F9DD2AB9BF6),
	GFw64be(0x00000063646A7966, 0x5C46729262C54D7F,
	        0x6BA9A8DAEAC48B79, 0x59D878DBEA8B26B2),
	GFw64be(0x000000748F325FE7, 0x1144926964B04BDC,
	        0x79EE34FB5A49A477, 0xC40976E8E91BB0BA),
	GFw64be(0x000001000FD20848, 0x395CF2984CF71FF6,
	        0xFA903C085EC4AD7E, 0xDBC4FAAEDA3FE730),
	GFw64be(0x000000B1DE31E469, 0x3A37FAD984A32B18,
	        0xD665E5A947CD2077, 0x57D273C4E1398B5B),
	GFw64be(0x0000005E0F64ADF1, 0xC497AC86FE556E12,
	        0xE20E3A3DA722EF18, 0x9C8E6B51826B4344),
	GFw64be(0x00000066ECBFC146, 0xD8E452190E84654B,
	        0x702AA79ACAC52274, 0x194B098837DB8106),
	GFw64be(0x0000011F067DDBAB, 0xEE8C14D63C1F3DDE,
	        0x66B72FD1A26E6540, 0xABC6BDE9EA1E5B62),
	GFw64be(0x000000D7145F42CF, 0xBBB90413A2AE6B5B,
	        0x990F6F8B55250333, 0x7B0F7698A2A1F5D1),
	GFw64be(0x000001CF139F4225, 0x6388A0D1084F721F,
	        0x0A7B5783900D0859, 0x2ACFC751FA84364A),
	GFw64be(0x0000011FE6EA498B, 0x201C105F3C271D6A,
	        0xEAB7A891064F6147, 0x4AD1B88A732E7676),
	GFw64be(0x000001F0184FE5C3, 0xA30BACE3EE4A7759,
	        0x199C4BB91129AB1D, 0x660A9A88A1C4E4F1),
	GFw64be(0x0000004DAA0C4160, 0x9F54D2D26CD04D8C,
	        0xE98B02987ECCAC29, 0xC51C68E56C2B65A4),
	GFw64be(0x0000011F866DDBAB, 0x6A8C145E3C3F3DDA,
	        0x66B72BD1826E6543, 0xABC6BDE8E21E5B62),
	GFw64be(0x000001E6E7F3DB62, 0xEA8A6492103F78CE,
	        0x1C9ABDD880A44407, 0xAB9BB47D7B10E924),
	GFw64be(0x0000000568F3AA60, 0x46B382C9805929B4,
	        0x8E238D08A5090D1C, 0xF0D765F699316F13),
	GFw64be(0x00000069A4978E05, 0x849430349E5C6BD4,
	        0xF589272326466E4D, 0xB81E36ADDC2A98A1),
	GFw64be(0x0000017EC2D8D5AC, 0x6C98A46A247F568E,
	        0x975E8CF2A4282D1F, 0xFB9D8B603B35DA9C),
	GFw64be(0x000001D637D26C07, 0xE8DF68008B601C2D,
	        0xBC7E7CA38F801BCF, 0x0145E84678F7E95F),
	GFw64be(0x000000450A93D342, 0x11C5B648FCF06EB8,
	        0x710B05D85B68EC2E, 0xC41B7E26F54BC491),
	GFw64be(0x00000003E979E481, 0xF98366709697279F,
	        0x0AE19DA1D1AC6730, 0xDB821BB4EEC1373F),
	GFw64be(0x00000197FD622105, 0x5AABAEF4E0DC1096,
	        0x20B7F813C12E8D69, 0xFD01D0208E160437),
	GFw64be(0x000001576BDAD84B, 0xC0FC58EC2E01556D,
	        0xDA3F9CC98ECA2355, 0x10DC988A38FDE206),
	GFw64be(0x000001724808120E, 0xB0518010CE9A5281,
	        0xA6BC80424D04A625, 0xB618963440E60B37),
	GFw64be(0x000001C5D9F16D26, 0x82AE146E7546330E,
	        0x645BDDB2026AF9D9, 0x1E16D2813B5E485C),
	GFw64be(0x00000048AB607640, 0xB103C227FAAF4C3F,
	        0x364818E85183C336, 0x3E883D12EAC28F59),
	GFw64be(0x0000007E0BA278C4, 0xCD25B717F6545EB3,
	        0x6A4E10CAB377EE19, 0x9D0CAE22929A6319),
	GFw64be(0x000001724848120E, 0xB07184104E9B5281,
	        0xA6BC88424D24A725, 0xB698962440E60B36),
	GFw64be(0x000000244ED9EF65, 0x18C2E6A6729D26DA,
	        0x1E02ADBB48A2C764, 0xFD865E78B601FF04),
	GFw64be(0x000000F659F05D4A, 0x9E6C34A73CC04A0B,
	        0x4C8EDCF86A63696D, 0x110C331437C97C30),
	GFw64be(0x000000B073A77369, 0xDFCB8019D695289A,
	        0x2E64D3D9F905E664, 0xCD933561B2433E5B),
	GFw64be(0x0000005747028203, 0xA9A9FAEFEF4C0417,
	        0x748FB00111CBBA89, 0x6D445D44CE8A8821),
	GFw64be(0x0000016F63EDC0AC, 0x21912655A4777A8B,
	        0xAA5B9B82152F2D1F, 0xCA8BF26177F7231D),
	GFw64be(0x000001C29BCAF3A7, 0xA703860E8E0A5A19,
	        0x313858D331222340, 0x6358B744B5C6D5C3),
	GFw64be(0x000000DA6534E22D, 0x87FF3875DE0C5F1A,
	        0xDD2CB6833F4FE340, 0x3808BFD5E62CEC83),
	GFw64be(0x00000116ECDF9689, 0x7CFC922056A33165,
	        0xF276AF61EE40E376, 0x1786858F1CFE971E),
	GFw64be(0x000000FB53611068, 0x8E9DB252C0206730,
	        0xC48DD948274C8142, 0x051A1AD781790C25),
	GFw64be(0x000001EB10A78B80, 0xD53DBAA0FED17FC5,
	        0xD7194310F740EE2D, 0xD4DBEFBD5CBD9FC1),
	GFw64be(0x000001F0F1CE6282, 0x4CD75A2082405C69,
	        0xFC3CDA80AFC00B08, 0x055CFD0F31EFBD57),
	GFw64be(0x00000123E4E887AE, 0xD1B8CA854AB11EB1,
	        0xA3D1A832D4838722, 0xD4D5EE62D1E652BC),
	GFw64be(0x0000002703C2AF40, 0x21EE5A48001303F4,
	        0x7A0318381AC80550, 0xC6C116BED95AF610),
        GFw64be(0x000001C39639F66F, 0x8A667AC94A397999,
                0x70D965EB0AC98652, 0xE1CAE5F0A1DF806C)
};

/*
 * d <- half-trace(a)
 */
static void
gf_halftrace(gf *d, const gf *a)
{
	/*
	 * H(a) = \sum_{i=0}^{116} a^(2^(2*i)).
	 * This is a linear transform. Moreover, for any field element x:
	 *   H(x^2) = H(x)^2
	 *   H(x)^2 + H(x) = x + Tr(x)
	 *   H(x^2 + x) = x
	 *   H(x) = H(sqrt(x)) + sqrt(x)
	 *
	 * The generic method is a basic matrix multiplication. To
	 * reduce the size of that matrix, we apply the method described
	 * below to remove half of the input bits; this is basically the
	 * same method as the one described by Knudsen in "Elliptic Scalar
	 * Multiplication Using Point Halving" (ASIACRYPT'99).
	 *
	 * Split a into its even and odd halves (as in gf_sqrt()):
	 *   a = ae + z*ao
	 * Then:
	 *   H(a) = H(ae) + H(z*ao)
	 *        = H(sqrt(ae)) + H(z*ao) + sqrt(ae)
	 * Note that sqrt(ae) has only half-size (since ae had only
	 * even-powered non-zero coefficients). We can continue
	 * recursively, by spliting sqrt(ae) into an even and odd halves:
	 *   sqrt(ae) = aee + z*aeo
	 *   H(a) = H(sqrt(aee)) + H(z*(ao + aeo)) + sqrt(ae) + sqrt(aee)
	 * with aee only having quarter size. Doing it a few more times,
	 * we end up with H(a) begin the sum of a single bit (0 or 1) and
	 * the half-trace of a value with only odd-powered non-zero
	 * coefficients. For that value, we use a precomputed matrix.
	 */
	__m128i m1 = _mm_set_epi32(
		0x55555555, 0x55555555, 0x55555555, 0x55555555);
	__m128i m2 = _mm_set_epi32(
		0x33333333, 0x33333333, 0x33333333, 0x33333333);
	__m128i m3 = _mm_set_epi32(
		0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F);
	__m128i sklo = _mm_set_epi8(
		0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
		0x0E, 0x0C, 0x0A, 0x08, 0x06, 0x04, 0x02, 0x00);
	__m128i skhi = _mm_set_epi8(
		0x0E, 0x0C, 0x0A, 0x08, 0x06, 0x04, 0x02, 0x00,
		0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF);
	__m128i a0, a1, ae0, ae1, ao0, ao1, d0, d1;
	gf t;

	/*
	 * We first ensure that the input is normalized.
	 */
	gf_norm(&t, a);
	a0 = t.v[0];
	a1 = t.v[1];

	/*
	 * Split the input into the even (ae) and odd (ao) parts.
	 */
	ae0 = _mm_and_si128(a0, m1);
	ae1 = _mm_and_si128(a1, m1);
	ao0 = _mm_and_si128(_mm_srli_epi32(a0, 1), m1);
	ao1 = _mm_and_si128(_mm_srli_epi32(a1, 1), m1);

	/*
	 * Compute sqrt(ae) with some squeezing, as in gf_sqrt().
	 */
	ae0 = _mm_and_si128(_mm_xor_si128(ae0, _mm_srli_epi32(ae0, 1)), m2);
	ae1 = _mm_and_si128(_mm_xor_si128(ae1, _mm_srli_epi32(ae1, 1)), m2);
	ae0 = _mm_and_si128(_mm_xor_si128(ae0, _mm_srli_epi32(ae0, 2)), m3);
	ae1 = _mm_and_si128(_mm_xor_si128(ae1, _mm_srli_epi32(ae1, 2)), m3);
	ae0 = _mm_xor_si128(ae0, _mm_srli_epi32(ae0, 4));
	ae1 = _mm_xor_si128(ae1, _mm_srli_epi32(ae1, 4));
	ae0 = _mm_xor_si128(
		_mm_shuffle_epi8(ae0, sklo),
		_mm_shuffle_epi8(ae1, skhi));

	/*
	 * H(a) = H(sqrt(ae)) + H(z*ao) + sqrt(ae)
	 * We accumulate sqrt(ae) into the result (r0).
	 */
	d0 = ae0;

	/*
	 * Repeat the process on sqrt(ae). At each iteration:
	 *  - The input is in ae0.
	 *  - We split the input into an even part and an odd part.
	 *    The odd part is accumulated in ao, and we compute the
	 *    square root of the even part as the new ae0; the
	 *    square root is also accumulated in d0.
	 *  - Each iteration halves the size of the value. Thus, seven
	 *    iterations are enough to bring down the current value to 1.
	 *
	 * In the last few iterations, the byte shuffling is not needed
	 * (value fits in a single byte) but pshufb cost is only a single
	 * cycle so we do not bother to optimize it away.
	 */
	for (int i = 0; i < 7; i ++) {
		ao0 = _mm_xor_si128(ao0,
			_mm_and_si128(_mm_srli_epi32(ae0, 1), m1));
		ae0 = _mm_and_si128(ae0, m1);
		ae0 = _mm_and_si128(_mm_xor_si128(
			ae0, _mm_srli_epi32(ae0, 1)), m2);
		ae0 = _mm_and_si128(_mm_xor_si128(
			ae0, _mm_srli_epi32(ae0, 2)), m3);
		ae0 = _mm_xor_si128(ae0, _mm_srli_epi32(ae0, 4));
		ae0 = _mm_shuffle_epi8(ae0, sklo);
		d0 = _mm_xor_si128(d0, ae0);
	}

	/*
	 * H(0) = 0, H(1) = 1, so we can just add the current ae0 to d0.
	 */
	d0 = _mm_xor_si128(d0, ae0);
	d1 = _mm_setzero_si128();

	/*
	 * We now have to compute the half-trace of z*ao. The polynomial
	 * z*ao is odd-powered (i.e. ao is even-powered). Since we
	 * normalized the input, ao fits on 233 bits; thus, there are
	 * 116 bits to consider.
	 *
	 * Each inner loop processes a 32-bit word, i.e. 16 bits to
	 * consider.
	 */
	ao0 = _mm_slli_epi32(ao0, 1);
	ao1 = _mm_slli_epi32(ao1, 1);
	for (int i = 0; i < 7; i ++) {
		__m128i x;

		x = _mm_shuffle_epi32(ao0, 0x00);
		for (int j = 15; j >= 0; j --) {
			__m128i xm = _mm_srai_epi32(x, 31);
			x = _mm_slli_epi32(x, 2);
			d0 = _mm_xor_si128(d0, _mm_and_si128(xm,
				HALFTRACE[16 * i + j].v[0]));
			d1 = _mm_xor_si128(d1, _mm_and_si128(xm,
				HALFTRACE[16 * i + j].v[1]));
		}
		if (i == 3) {
			ao0 = ao1;
		} else {
			ao0 = _mm_bsrli_si128(ao0, 4);
		}
	}

	/*
	 * Specialized handling for the last four bits.
	 */
	ao0 = _mm_shuffle_epi32(_mm_slli_epi32(ao0, 24), 0x00);
	for (int j = 3; j >= 0; j --) {
		__m128i xm = _mm_srai_epi32(ao0, 31);
		ao0 = _mm_slli_epi32(ao0, 2);
		d0 = _mm_xor_si128(d0, _mm_and_si128(xm,
			HALFTRACE[112 + j].v[0]));
		d1 = _mm_xor_si128(d1, _mm_and_si128(xm,
			HALFTRACE[112 + j].v[1]));
	}

	d->v[0] = d0;
	d->v[1] = d1;
}

#if GF233_INV_FLT
/*
 * If GF233_INV_FLT is non-zero, then inversions are computed through
 * Fermat's little theorem, with some chains of squarings optimized
 * through the linearity of the Frobenius endomorphism. This method
 * is also known as Itoh-Tsujii inversion. It is fast, but the static
 * tables of constants used for that matrix-based optimization increase
 * the compiled binary size by about 15 kB.
 *
 * Performance measured on an Intel i5-8259U CPU (Coffee Lake):
 * about 2890 cycles for an inversion (+33 for a division).
 */

/* d <- a^(2^n) */
static void
gf_xfrob(gf *d, const gf *a, unsigned n)
{
	if (n == 0) {
		*d = *a;
		return;
	}
	gf_sqr(d, a);
	n --;
	while (n -- > 0) {
		gf_sqr(d, d);
	}
}

static inline void
gf_xfrob_tab(gf *d, const gf *a, const gf *xft)
{
	gf t;
	__m128i a0, a1, d0, d1;

	gf_norm(&t, a);
	a0 = t.v[0];
	a1 = t.v[1];
	d0 = _mm_setzero_si128();
	d1 = _mm_setzero_si128();
	for (int i = 0; i < 7; i ++) {
		__m128i x = _mm_shuffle_epi32(a0, 0x00);
		for (int j = 31; j >= 0; j --) {
			__m128i xm = _mm_srai_epi32(x, 31);
			x = _mm_slli_epi32(x, 1);
			d0 = _mm_xor_si128(d0, _mm_and_si128(xm,
				xft[32 * i + j].v[0]));
			d1 = _mm_xor_si128(d1, _mm_and_si128(xm,
				xft[32 * i + j].v[1]));
		}
		if (i == 3) {
			a0 = a1;
		} else {
			a0 = _mm_bsrli_si128(a0, 4);
		}
	}
	a0 = _mm_shuffle_epi32(_mm_slli_epi32(a0, 23), 0x00);
	for (int j = 8; j >= 0; j --) {
		__m128i xm = _mm_srai_epi32(a0, 31);
		a0 = _mm_slli_epi32(a0, 1);
		d0 = _mm_xor_si128(d0, _mm_and_si128(xm, xft[224 + j].v[0]));
		d1 = _mm_xor_si128(d1, _mm_and_si128(xm, xft[224 + j].v[1]));
	}

	d->v[0] = d0;
	d->v[1] = d1;
}

#include "gf233_frob58.h"
#include "gf233_frob116.h"

#define gf_xfrob58(d, a)    gf_xfrob_tab(d, a, GF233_XFROB_58)
#define gf_xfrob116(d, a)   gf_xfrob_tab(d, a, GF233_XFROB_116)

/* d <- 1/a  (if a == 0 then d is set to zero) */
static void
gf_inv(gf *d, const gf *a)
{
	gf x, y;

	/* a^(2^2-1) */
	gf_sqr(&y, a);
	gf_mul(&x, &y, a);

	/* a^(2^3-1) */
	gf_sqr(&y, &x);
	gf_mul(&x, &y, a);

	/* a^(2^6-1) */
	gf_xfrob(&y, &x, 3);
	gf_mul(&x, &y, &x);

	/* a^(2^7-1) */
	gf_sqr(&y, &x);
	gf_mul(&x, &y, a);

	/* a^(2^14-1) */
	gf_xfrob(&y, &x, 7);
	gf_mul(&x, &y, &x);

	/* a^(2^28-1) */
	gf_xfrob(&y, &x, 14);
	gf_mul(&x, &y, &x);

	/* a^(2^29-1) */
	gf_sqr(&y, &x);
	gf_mul(&x, &y, a);

	/* a^(2^58-1) */
	gf_xfrob(&y, &x, 29);
	gf_mul(&x, &y, &x);

	/* a^(2^116-1) */
	gf_xfrob58(&y, &x);
	gf_mul(&x, &y, &x);

	/* a^(2^232-1) */
	gf_xfrob116(&y, &x);
	gf_mul(&x, &y, &x);

	/* a^(2^233-2) = 1/a */
	gf_sqr(d, &x);
}

/* d <- x/y (returns 0 if y == 0) */
static inline void
gf_div(gf *d, const gf *x, const gf *y)
{
	gf t;

	gf_inv(&t, y);
	gf_mul(d, x, &t);
}

#else
/*
 * If GF233_INV_FLT is zero, then inversions and divisions are computed
 * with a binary GCD variant (constant-time).
 *
 * Performance measured on an Intel i5-8259U CPU (Coffee Lake): about
 * 5190 cycles for a division. This is slower than the FLT-based method,
 * but avoids the need for the tables of constants and thus makes the
 * compiled binary smaller.
 */

/*
 * Given polynomials a and b (held in gf structures) and update
 * coefficients f0, g0, f1 and g1, compute new values for a and b:
 *   (a, b) <- ((a*f0 + b*g0)/z^64, (a*f1 + b*g1)/z^64)
 * Division is computed over plain polynomials and is assumed to be exact
 * (i.e. it is a right shift by 64 bits, and the dropped bits are ignored).
 *
 * f0 is in the high half of fg0, g0 is in the low half of fg0
 * f1/z is in the high half of fg1, g1/z is in the low half of fg1
 * (f1 and g1 are always even).
 */
static inline void
poly_lin_div64(gf *a, gf *b, __m128i fg0, __m128i fg1)
{
	/* Adjust f1 and g1, keeping track of their high bits. */
	__m128i hh = _mm_srai_epi32(fg1, 31);
	__m128i hf1 = _mm_shuffle_epi32(hh, 0xFF);
	__m128i hg1 = _mm_shuffle_epi32(hh, 0x55);
	fg1 = _mm_slli_epi64(fg1, 1);

	/* a*f0 + b*g0 */
	__m128i c01 = _mm_clmulepi64_si128(a->v[0], fg0, 0x10);
	__m128i c12 = _mm_clmulepi64_si128(a->v[0], fg0, 0x11);
	__m128i d01 = _mm_clmulepi64_si128(b->v[0], fg0, 0x00);
	__m128i d12 = _mm_clmulepi64_si128(b->v[0], fg0, 0x01);
	__m128i c23 = _mm_clmulepi64_si128(a->v[1], fg0, 0x10);
	__m128i c34 = _mm_clmulepi64_si128(a->v[1], fg0, 0x11);
	__m128i d23 = _mm_clmulepi64_si128(b->v[1], fg0, 0x00);
	__m128i d34 = _mm_clmulepi64_si128(b->v[1], fg0, 0x01);
	c01 = _mm_xor_si128(c01, d01);
	c12 = _mm_xor_si128(c12, d12);
	c23 = _mm_xor_si128(c23, d23);
	c34 = _mm_xor_si128(c34, d34);

	/* a*f1 + b*g1 */
	__m128i s01 = _mm_clmulepi64_si128(a->v[0], fg1, 0x10);
	__m128i s12 = _mm_clmulepi64_si128(a->v[0], fg1, 0x11);
	__m128i t01 = _mm_clmulepi64_si128(b->v[0], fg1, 0x00);
	__m128i t12 = _mm_clmulepi64_si128(b->v[0], fg1, 0x01);
	__m128i s23 = _mm_clmulepi64_si128(a->v[1], fg1, 0x10);
	__m128i s34 = _mm_clmulepi64_si128(a->v[1], fg1, 0x11);
	__m128i t23 = _mm_clmulepi64_si128(b->v[1], fg1, 0x00);
	__m128i t34 = _mm_clmulepi64_si128(b->v[1], fg1, 0x01);
	s01 = _mm_xor_si128(s01, t01);
	s12 = _mm_xor_si128(s12, t12);
	s23 = _mm_xor_si128(s23, t23);
	s34 = _mm_xor_si128(s34, t34);

	/* Adjust a*f1 + b*g1 to account for the high bits of f1 and g1. */
	s12 = _mm_xor_si128(s12, _mm_xor_si128(
		_mm_and_si128(hf1, a->v[0]),
		_mm_and_si128(hg1, b->v[0])));
	s34 = _mm_xor_si128(s34, _mm_xor_si128(
		_mm_and_si128(hf1, a->v[1]),
		_mm_and_si128(hg1, b->v[1])));

	/* Assemble parts of a, with the right-shift by 64 bits. */
	a->v[0] = _mm_xor_si128(c12,
		_mm_xor_si128(
			_mm_bsrli_si128(c01, 8),
			_mm_bslli_si128(c23, 8)));
	a->v[1] = _mm_xor_si128(c34, _mm_bsrli_si128(c23, 8));
	b->v[0] = _mm_xor_si128(s12,
		_mm_xor_si128(
			_mm_bsrli_si128(s01, 8),
			_mm_bslli_si128(s23, 8)));
	b->v[1] = _mm_xor_si128(s34, _mm_bsrli_si128(s23, 8));
}

/*
 * (u, v) <- (u*f0 + v*g0, u*f1 + v*g1)   (in the field)
 * f0 is in the high half of fg0; g0 is in the low half
 * f1/z is in the high half of fg1; g1/z is in the low half
 * (f1 and g1 are always even)
 */
static void
gf_lin(gf *u, gf *v, __m128i fg0, __m128i fg1)
{
	/* Adjust f1 and g1, keeping track of their high bits. */
	__m128i hh = _mm_srai_epi32(fg1, 31);
	__m128i hf1 = _mm_shuffle_epi32(hh, 0xFF);
	__m128i hg1 = _mm_shuffle_epi32(hh, 0x55);
	fg1 = _mm_slli_epi64(fg1, 1);

	/* u*f0 + v*g0 (unreduced) */
	__m128i c01 = _mm_clmulepi64_si128(u->v[0], fg0, 0x10);
	__m128i c12 = _mm_clmulepi64_si128(u->v[0], fg0, 0x11);
	__m128i d01 = _mm_clmulepi64_si128(v->v[0], fg0, 0x00);
	__m128i d12 = _mm_clmulepi64_si128(v->v[0], fg0, 0x01);
	__m128i c23 = _mm_clmulepi64_si128(u->v[1], fg0, 0x10);
	__m128i c34 = _mm_clmulepi64_si128(u->v[1], fg0, 0x11);
	__m128i d23 = _mm_clmulepi64_si128(v->v[1], fg0, 0x00);
	__m128i d34 = _mm_clmulepi64_si128(v->v[1], fg0, 0x01);
	c01 = _mm_xor_si128(c01, d01);
	c12 = _mm_xor_si128(c12, d12);
	c23 = _mm_xor_si128(c23, d23);
	c34 = _mm_xor_si128(c34, d34);

	/* u*f1 + v*g1 (unreduced) */
	__m128i s01 = _mm_clmulepi64_si128(u->v[0], fg1, 0x10);
	__m128i s12 = _mm_clmulepi64_si128(u->v[0], fg1, 0x11);
	__m128i t01 = _mm_clmulepi64_si128(v->v[0], fg1, 0x00);
	__m128i t12 = _mm_clmulepi64_si128(v->v[0], fg1, 0x01);
	__m128i s23 = _mm_clmulepi64_si128(u->v[1], fg1, 0x10);
	__m128i s34 = _mm_clmulepi64_si128(u->v[1], fg1, 0x11);
	__m128i t23 = _mm_clmulepi64_si128(v->v[1], fg1, 0x00);
	__m128i t34 = _mm_clmulepi64_si128(v->v[1], fg1, 0x01);
	s01 = _mm_xor_si128(s01, t01);
	s12 = _mm_xor_si128(s12, t12);
	s23 = _mm_xor_si128(s23, t23);
	s34 = _mm_xor_si128(s34, t34);

	/* Adjust a*f1 + b*g1 to account for the high bits of f1 and g1. */
	s12 = _mm_xor_si128(s12, _mm_xor_si128(
		_mm_and_si128(hf1, u->v[0]),
		_mm_and_si128(hg1, v->v[0])));
	s34 = _mm_xor_si128(s34, _mm_xor_si128(
		_mm_and_si128(hf1, u->v[1]),
		_mm_and_si128(hg1, v->v[1])));

	/* Assemble and reduce. */
	__m128i u0 = _mm_xor_si128(c01, _mm_bslli_si128(c12, 8));
	__m128i u1 = _mm_xor_si128(c23, _mm_xor_si128(
		_mm_bsrli_si128(c12, 8), _mm_bslli_si128(c34, 8)));
	__m128i h0 = _mm_unpackhi_epi64(
		_mm_slli_epi64(c34, 23),
		_mm_srli_epi64(c34, 41));
	__m128i k0 = _mm_and_si128(
		_mm_slli_epi64(c34, 33),
		_mm_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0, 0));
	__m128i k1 = _mm_bsrli_si128(_mm_srli_epi64(c34, 31), 8);
	u->v[0] = _mm_xor_si128(u0, _mm_xor_si128(h0, k0));
	u->v[1] = _mm_xor_si128(u1, k1);

	__m128i v0 = _mm_xor_si128(s01, _mm_bslli_si128(s12, 8));
	__m128i v1 = _mm_xor_si128(s23, _mm_xor_si128(
		_mm_bsrli_si128(s12, 8), _mm_bslli_si128(s34, 8)));
	__m128i m0 = _mm_unpackhi_epi64(
		_mm_slli_epi64(s34, 23),
		_mm_srli_epi64(s34, 41));
	__m128i n0 = _mm_and_si128(
		_mm_slli_epi64(s34, 33),
		_mm_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0, 0));
	__m128i n1 = _mm_bsrli_si128(_mm_srli_epi64(s34, 31), 8);
	v->v[0] = _mm_xor_si128(v0, _mm_xor_si128(m0, n0));
	v->v[1] = _mm_xor_si128(v1, n1);
}

/* Modulus itself, stored in a structure for an element. */
static const gf GF_MOD = GFw64be((uint64_t)1 << 41, 0, 1 << 10, 1);

/* d <- x/y (returns 0 if y == 0) */
static void
gf_div(gf *d, const gf *x, const gf *y)
{
	/*
	 * Algorithm is the binary GCD variant described first by
	 * Brunner, Curiger and Hofstetter in 1993 ("On computing
	 * multiplicative inverses in GF(2^m)", IEEE Trans. on Computers,
	 * vol. 48, issue 8, pp. 1010-1015). We can perform iterations in
	 * chunks of 64, working only on the low words, and propagating
	 * the changes only once every 64 inner iterations. Compared to
	 * the binary GCD on integers (as in eprint.iacr.org/2020/972),
	 * the carryless nature of multiplications in GF(2)[z] makes it
	 * so that we do not need approximations of the high words; we
	 * can use a simple integer counter that keeps track of the balance
	 * between the two maximum bounds on the value degrees. We thus
	 * avoid the operations that gather the top words, and we can also
	 * do 64 iterations in the inner loop instead of 31. There is also
	 * no corrective step for the sign.
	 *
	 *   a <- y
	 *   b <- m (field modulus = z^233 + z^74 + 1)
	 *   u <- x
	 *   v <- 0
	 *   n <- -1
	 *   invariants:
	 *      a*x = u*y mod m
	 *      b*x = v*y mod m
	 *      n = maxsize(a) - maxsize(b)
	 *   (maxsize(t) is the proven bound on the size of t;
	 *   we use size(t) = 1 + degree(t))
	 *
	 * At each iteration:
	 *   if a_0 != 0:
	 *       if n < 0:
	 *           (a, u, b, v) <- (b, v, a, u)
	 *           n = -n
	 *       a <- a + b
	 *       u <- u + v
	 *   a <- a/z
	 *   n <- n - 1
	 *
	 * b_0 is always equal to 1.
	 *
	 * maxsize(a) + maxsize(b) starts at 467 (since the initial
	 * value of a has size at most 233, while b starts at size 234),
	 * and is decreased by 1 at each iteration; thus, 465 iterations
	 * are sufficient to ensure that size(a) + size(b) <= 2,
	 * assuming that y != 0. It is then not possible, at that point,
	 * that size(b) == 2, since that would mean that a == 0 and
	 * that the GCD of y and b is either z or z+1; we are in a
	 * field, so the GCD can only be 1 or 0. Thus, if y != 0, then
	 * after 465 iterations we have b = 1, and v = y/x in the field.
	 * If y == 0, then a, b, u and v remain unchanged throughout; in
	 * particular, v == 0, which is the result we want to return in
	 * that case. Thus, v always contains the result after 465
	 * iterations.
	 */
	gf a, b, u, v;
	__m128i xn;

	/* 1/z^512 */
	static const gf GF_INVZ512 = GFw64be(
		0x0000000000000000, 0x0800002000000000,
		0x0000004008000000, 0x0000000010020040);

	gf_norm(&a, y);
	b = GF_MOD;
	u = *x;
	v = GF_ZERO;
	xn = _mm_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF);

	for (int i = 0; i < 8; i ++) {
		/*
		 * xa = low 64 bits of a (high half is ignored)
		 * xb = low 64 bits of b (high half is ignored)
		 * fg0 = update factors for a (high = f0, low = g0)
		 * fg1 = update factors for b (high = f1, low = g1)
		 */
		__m128i xa = a.v[0];
		__m128i xb = b.v[0];
		__m128i fg0 = _mm_set_epi32(0, 1, 0, 0);
		__m128i fg1 = _mm_set_epi32(0, 0, 0, 1);
		int num_inner;

		if (i == 7) {
			/* We do only 17 iterations in the last inner loop. */
			num_inner = 17;
			fg0 = _mm_set_epi32(0x00008000, 0, 0, 0);
			fg1 = _mm_set_epi32(0, 0, 0x00008000, 0);
		} else {
			num_inner = 64;
			fg0 = _mm_set_epi32(0, 1, 0, 0);
			fg1 = _mm_set_epi32(0, 0, 0, 1);
		}
		for (int j = 1;; j ++) {
			/* copy xa_0 to high bit of each byte of a_odd */
			__m128i a_odd = _mm_shuffle_epi8(
				_mm_slli_epi32(xa, 7), _mm_setzero_si128());
			/* n_neg = -1 if xn < 0 */
			__m128i n_neg =
				_mm_cmpgt_epi32(_mm_setzero_si128(), xn);
			/* do swap only if xa is odd and xn < 0 */
			__m128i swap = _mm_and_si128(a_odd, n_neg);
			__m128i nxa = _mm_blendv_epi8(xa, xb, swap);
			xb = _mm_blendv_epi8(xb, xa, swap);
			__m128i nfg0 = _mm_blendv_epi8(fg0, fg1, swap);
			fg1 = _mm_blendv_epi8(fg1, fg0, swap);
			__m128i nxn = _mm_blendv_epi8(xn,
				_mm_sub_epi32(_mm_setzero_si128(), xn), swap);
			/* if xa is odd then XOR xb into xa */
			nxa = _mm_blendv_epi8(
				nxa, _mm_xor_si128(nxa, xb), a_odd);
			fg0 = _mm_blendv_epi8(
				nfg0, _mm_xor_si128(nfg0, fg1), a_odd);
			/* xa is even, divide it by 2 */
			xa = _mm_srli_epi64(nxa, 1);
			xn = _mm_sub_epi32(nxn, _mm_set_epi32(1, 1, 1, 1));
			/* We exit on the last iteration _before_ the shift
			   on fg1 because we need to handle the case of
			   f1 and/or g1 overflowing the 64-bit field. */
			if (j == num_inner) {
				break;
			}
			fg1 = _mm_slli_epi64(fg1, 1);

		#if 0
			/*
			 * This is an alternate version of the inner loop,
			 * which uses only SSE2 opcodes, no pshufb or pblendvb.
			 * It is a bit slower.
			 */

			/* copy xa_0 to high bit of each byte of a_odd */
			__m128i a_odd = _mm_shuffle_epi32(
				_mm_srai_epi32(_mm_slli_epi32(xa, 31), 31), 0);
			/* n_neg = -1 if xn < 0 */
			__m128i n_neg =
				_mm_cmpgt_epi32(_mm_setzero_si128(), xn);
			/* do swap only if xa is odd and xn < 0 */
			__m128i swap = _mm_and_si128(a_odd, n_neg);
			__m128i t1 = _mm_and_si128(swap,
				_mm_xor_si128(xa, xb));
			xa = _mm_xor_si128(xa, t1);
			xb = _mm_xor_si128(xb, t1);
			__m128i t2 = _mm_and_si128(swap,
				_mm_xor_si128(fg0, fg1));
			fg0 = _mm_xor_si128(fg0, t2);
			fg1 = _mm_xor_si128(fg1, t2);
			xn = _mm_sub_epi32(xn,
				_mm_and_si128(swap, _mm_add_epi32(xn, xn)));
			/* if xa is odd then XOR xb into xa */
			xa = _mm_xor_si128(xa, _mm_and_si128(a_odd, xb));
			fg0 = _mm_xor_si128(fg0, _mm_and_si128(a_odd, fg1));
			/* xa is even, divide it by 2 */
			xa = _mm_srli_epi64(xa, 1);
			xn = _mm_sub_epi32(xn, _mm_set_epi32(1, 1, 1, 1));
			/* We exit on the last iteration _before_ the shift
			   on fg1 because we need to handle the case of
			   f1 and/or g1 overflowing the 64-bit field. */
			if (j == num_inner) {
				break;
			}
			fg1 = _mm_slli_epi64(fg1, 1);
		#endif
		}

		/*
		 * Propagate changes to a and b:
		 *   (a, b) <- ((a*f0 + b*g0)/z^64, (a*f1 + b*g1)/z^64)
		 *   (u, v) <- (u*f0 + v*g0, u*f1 + v*g1)
		 */
		poly_lin_div64(&a, &b, fg0, fg1);
		gf_lin(&u, &v, fg0, fg1);
	}

	/*
	 * Result is in v, except that it was multiplied by z^512, so we
	 * must divide it by that much; this is a simple multiplication
	 * by a precomputed constant.
	 */
	gf_mul(d, &v, &GF_INVZ512);
}

/* d <- 1/a  (if a == 0 then d is set to zero) */
static inline void
gf_inv(gf *d, const gf *a)
{
	gf_div(d, &GF_ONE, a);
}

#endif

/* Returned value is -1 if a is zero, 0 otherwise. */
static inline uint64_t
gf_iszero(const gf *a)
{
	gf t;
	uint64_t r;

	gf_norm(&t, a);
	r = t.w[0] | t.w[1] | t.w[2] | t.w[3];
	return ((r | -r) >> 63) - 1;
}

/* Returned value is -1 if a == b, 0 otherwise. */
static inline uint64_t
gf_equals(const gf *a, const gf *b)
{
	gf t;
	uint64_t r;

	gf_add(&t, a, b);
	gf_norm(&t, &t);
	r = t.w[0] | t.w[1] | t.w[2] | t.w[3];
	return ((r | -r) >> 63) - 1;
}

/* Decode d from exactly 30 bytes. A failure is reported if the source
   is not canonical (bits 1-7 in the last byte were not all equal to zero).
   On success, -1 is returned, or 0 on failure. On failure, *d is set
   to zero. */
static inline uint64_t
gf_decode(gf *d, const void *src)
{
	/*
	 * We can copy the bytes directly since x86 CPUs use little-endian.
	 */
	memcpy((uint8_t *)(void *)d, src, 30);
	memset((uint8_t *)(void *)d + 30, 0, 2);

	/*
	 * Value is valid as long as the extra bits are zero.
	 */
	uint64_t r = ((const uint8_t *)src)[29];
	r = ((r + 0xFE) >> 8) - 1;
	d->w[0] &= r;
	d->w[1] &= r;
	d->w[2] &= r;
	d->w[3] &= r;
	return r;
}

/* Encode a into exactly 30 bytes. */
static inline void
gf_encode(void *dst, const gf *a)
{
	gf t;

	gf_norm(&t, a);
	memcpy(dst, &t, 30);
}

/* d <- a0 if ctl == 0
   d <- a1 if ctl == -1
   ctl MUST be either 0 or -1 */
static inline void
gf_select(gf *d, const gf *a0, const gf *a1, uint64_t ctl)
{
	__m128i m = _mm_set1_epi64x(ctl);
	d->v[0] = _mm_blendv_epi8(a0->v[0], a1->v[0], m);
	d->v[1] = _mm_blendv_epi8(a0->v[1], a1->v[1], m);
}

/* d <- a if ctl == -1, but unchanged if ctl == 0 */
static inline void
gf_condset(gf *d, const gf *a, uint64_t ctl)
{
	gf_select(d, d, a, ctl);
}

/* d <- d + a if ctl == -1, but unchanged if ctl == 0
   ctl MUST be either 0 or -1 */
static inline void
gf_condadd(gf *d, const gf *a, uint64_t ctl)
{
	__m128i m = _mm_set1_epi64x(ctl);
	d->v[0] = _mm_xor_si128(d->v[0], _mm_and_si128(a->v[0], m));
	d->v[1] = _mm_xor_si128(d->v[1], _mm_and_si128(a->v[1], m));
}

#else
/* ======================================================================= */
/*
 * Generic implementation, uses 64x64 multiplications.
 */

/*
 * We represent values as four limbs:
 *   w[0]   bits 0..57
 *   w[1]   bits 58..115
 *   w[2]   bits 116..173
 *   w[3]   bits 174..232
 * Limbs 0, 1 and 2 contain 58 bits each; limb 3 contains 59 bits. In
 * each limb, the contents is stored in the low bits; extra unused bits
 * are always zero.
 */
typedef struct {
	uint64_t w[4];
} gf;

/*
 * Macro for an initializer for a field element. Values are provided as
 * four 64-bit limbs that encode the value in big-endian order.
 */
#define GFw64be(w3, w2, w1, w0)  { GFext(w3, w2, w1, w0) }

/*
 * Helper macro for defining initializers for constant point that use
 * the external representation (with values being opaque sequences of
 * 64-bit integers).
 */
#define GFext(w3, w2, w1, w0) \
    (uint64_t)(w0) & 0x03FFFFFFFFFFFFFF, \
    (((uint64_t)(w0) >> 58) | ((uint64_t)(w1) << 6)) & 0x03FFFFFFFFFFFFFF, \
    (((uint64_t)(w1) >> 52) | ((uint64_t)(w2) << 12)) & 0x03FFFFFFFFFFFFFF, \
    (((uint64_t)(w2) >> 46) | ((uint64_t)(w3) << 18)) & 0x07FFFFFFFFFFFFFF

/*
 * Helper macros for tests: masks to apply to random 64-bit limbs
 * in order to generate a random field element.
 */
#define GF_TEST_RMASK0    0x03FFFFFFFFFFFFFF
#define GF_TEST_RMASK1    0x03FFFFFFFFFFFFFF
#define GF_TEST_RMASK2    0x03FFFFFFFFFFFFFF
#define GF_TEST_RMASK3    0x07FFFFFFFFFFFFFF

/*
 * Multiplications are slow (even on a Coffee Lake x86, multiplications
 * cost more than 5 squarings).
 */
#define GF_MUL_FAST   0

static const gf GF_ZERO = GFw64be(0, 0, 0, 0);
static const gf GF_ONE = GFw64be(0, 0, 0, 1);

/* d <- a + b */
static inline void
gf_add(gf *d, const gf *a, const gf *b)
{
	d->w[0] = a->w[0] ^ b->w[0];
	d->w[1] = a->w[1] ^ b->w[1];
	d->w[2] = a->w[2] ^ b->w[2];
	d->w[3] = a->w[3] ^ b->w[3];
}

/*
 * Carry-less multiplication of two words.
 * IMPORTANT: both operands must fit in 60 bits each. Their upper 4 bits
 * MUST be zero.
 */
static inline unsigned __int128
clmul(uint64_t x, uint64_t y)
{
	/*
	 * We use integer multiplications with "holes" between the
	 * relevant bits. With input operands up to 60 bits each, we
	 * can use holes of 3 bits; each individual multiplication will
	 * use values with Hamming weight at most 15, and 4 bits are
	 * then enough to accumulate up to 15 ones without spilling into
	 * the next chunk.
	 *
	 * This was used in BearSSL's GHASH implementation on 64-bit
	 * systems with no pclmul support (file src/hash/ghash_ctmul64.c).
	 */
	uint64_t x0, x1, x2, x3;
	uint64_t y0, y1, y2, y3;
	unsigned __int128 z0, z1, z2, z3;

	x0 = x & (uint64_t)0x1111111111111111;
	x1 = x & (uint64_t)0x2222222222222222;
	x2 = x & (uint64_t)0x4444444444444444;
	x3 = x & (uint64_t)0x8888888888888888;
	y0 = y & (uint64_t)0x1111111111111111;
	y1 = y & (uint64_t)0x2222222222222222;
	y2 = y & (uint64_t)0x4444444444444444;
	y3 = y & (uint64_t)0x8888888888888888;
#define M(u, v)   ((unsigned __int128)(u) * (unsigned __int128)(v))
	z0 = M(x0, y0) ^ M(x1, y3) ^ M(x2, y2) ^ M(x3, y1);
	z1 = M(x0, y1) ^ M(x1, y0) ^ M(x2, y3) ^ M(x3, y2);
	z2 = M(x0, y2) ^ M(x1, y1) ^ M(x2, y0) ^ M(x3, y3);
	z3 = M(x0, y3) ^ M(x1, y2) ^ M(x2, y1) ^ M(x3, y0);
#undef M
#define C128(t)   ((unsigned __int128)(t) | ((unsigned __int128)(t) << 64))
	z0 &= C128(0x1111111111111111);
	z1 &= C128(0x2222222222222222);
	z2 &= C128(0x4444444444444444);
	z3 &= C128(0x8888888888888888);
#undef C128
	return z0 | z1 | z2 | z3;
}

/*
 * Given a 465-bit result in r0..r7 (58-bit limbs, r7 uses 59 bits),
 * reduce the value into r0..r3. Values r4 and r5 are modified.
 */
#define REDUCE_465   do { \
		r2 ^= (r7 & 1) << 57; \
		r3 ^= r7 >> 1; \
		r4 ^= (r7 << 15) & 0x03FFFFFFFFFFFFFF; \
		r5 ^= r7 >> 43; \
		r1 ^= (r6 & 1) << 57; \
		r2 ^= r6 >> 1; \
		r3 ^= (r6 << 15) & 0x03FFFFFFFFFFFFFF; \
		r4 ^= r6 >> 43; \
		r0 ^= (r5 & 1) << 57; \
		r1 ^= r5 >> 1; \
		r2 ^= (r5 << 15) & 0x03FFFFFFFFFFFFFF; \
		r3 ^= r5 >> 43; \
		r3 ^= (r4 & 1) << 58; \
		r4 &= ~(uint64_t)1; \
		r0 ^= r4 >> 1; \
		r1 ^= (r4 << 15) & 0x03FFFFFFFFFFFFFF; \
		r2 ^= r4 >> 43; \
	} while (0)

/* d <- a * b */
static inline void
gf_mul(gf *d, const gf *a, const gf *b)
{
	/*
	 * One level of Karatsuba:
	 *   (a0 + a1*z^116)*(b0 + b1*z^116)
	 *   = a0*b0 + (a0*b1 + a1*b0)*z^116 + a1*b1*z^232
	 * With:
	 *   e = a0*b0
	 *   f = a1*b1
	 *   g = (a0 + a1)*(b0 + b1) + e + f
	 * so that g = a0*b1 + a1*b0
	 */
	uint64_t a0, a1, a2, a3;
	uint64_t b0, b1, b2, b3;
	uint64_t c0, c1, d0, d1;
	uint64_t e0, e1, e2, e3;
	uint64_t f0, f1, f2, f3;
	uint64_t g0, g1, g2, g3;
	uint64_t r0, r1, r2, r3, r4, r5, r6, r7;
	unsigned __int128 tl, tm, th;

	a0 = a->w[0];
	a1 = a->w[1];
	a2 = a->w[2];
	a3 = a->w[3];
	b0 = b->w[0];
	b1 = b->w[1];
	b2 = b->w[2];
	b3 = b->w[3];

	/* c0:c1 <- a0:a1 + a2:a3 */
	c0 = a0 ^ a2;
	c1 = a1 ^ a3;

	/* d0:d1 <- b0:b1 + b2:b3 */
	d0 = b0 ^ b2;
	d1 = b1 ^ b3;

	/* e0:e1:e2:e3 <- a0:a1 * b0:b1 */
	tl = clmul(a0, b0);
	th = clmul(a1, b1);
	tm = clmul(a0 ^ a1, b0 ^ b1) ^ tl ^ th;
	e0 = (uint64_t)tl & 0x03FFFFFFFFFFFFFF;
	e1 = (uint64_t)(tl >> 58) ^ ((uint64_t)tm & 0x03FFFFFFFFFFFFFF);
	e2 = (uint64_t)(tm >> 58) ^ ((uint64_t)th & 0x03FFFFFFFFFFFFFF);
	e3 = (uint64_t)(th >> 58);

	/* f0:f1:f2:f3 <- a2:a3 * b2:b3 */
	tl = clmul(a2, b2);
	th = clmul(a3, b3);
	tm = clmul(a2 ^ a3, b2 ^ b3) ^ tl ^ th;
	f0 = (uint64_t)tl & 0x03FFFFFFFFFFFFFF;
	f1 = (uint64_t)(tl >> 58) ^ ((uint64_t)tm & 0x03FFFFFFFFFFFFFF);
	f2 = (uint64_t)(tm >> 58) ^ ((uint64_t)th & 0x03FFFFFFFFFFFFFF);
	f3 = (uint64_t)(th >> 58);

	/* g0:g1:g2:g3 <- c0:c1 * d0:d1 */
	tl = clmul(c0, d0);
	th = clmul(c1, d1);
	tm = clmul(c0 ^ c1, d0 ^ d1) ^ tl ^ th;
	g0 = (uint64_t)tl & 0x03FFFFFFFFFFFFFF;
	g1 = (uint64_t)(tl >> 58) ^ ((uint64_t)tm & 0x03FFFFFFFFFFFFFF);
	g2 = (uint64_t)(tm >> 58) ^ ((uint64_t)th & 0x03FFFFFFFFFFFFFF);
	g3 = (uint64_t)(th >> 58);

	/*
	 * e3 uses 57 bits (e is the product of two 116-bit values)
	 * f3 uses 59 bits (f is the product of two 117-bit values)
	 * g3 uses 59 bits too; however, we have:
	 *     e + f + g = (a0:a1)*(b2:b3) + (a2:a3)*(b0:b1)
	 * and since a0:a1 and b0:b1 both are 116-bit values, the
	 * sum e + f + g must fit on 232 bits in total, hence yielding
	 * a high limb of 58 bits only.
	 *
	 * In total, we can simply add up the limbs, and only the top
	 * limb will exceed 58 bits.
	 */

	/* Assemble unreduced result in r0..r7. */
	r0 = e0;
	r1 = e1;
	r2 = e2 ^ e0 ^ f0 ^ g0;
	r3 = e3 ^ e1 ^ f1 ^ g1;
	r4 = f0 ^ e2 ^ f2 ^ g2;
	r5 = f1 ^ e3 ^ f3 ^ g3;
	r6 = f2;
	r7 = f3;

	/* Reduce the value. */
	REDUCE_465;

	d->w[0] = r0;
	d->w[1] = r1;
	d->w[2] = r2;
	d->w[3] = r3;
}

static inline uint64_t
expand_lo(uint64_t x)
{
	x = (x & 0x000000000000FFFF) | ((x & 0x00000000FFFF0000) << 16);
	x = (x & 0x000000FF000000FF) | ((x & 0x0000FF000000FF00) <<  8);
	x = (x & 0x000F000F000F000F) | ((x & 0x00F000F000F000F0) <<  4);
	x = (x & 0x0303030303030303) | ((x & 0x0C0C0C0C0C0C0C0C) <<  2);
	x = (x & 0x1111111111111111) | ((x & 0x2222222222222222) <<  1);
	return x;
}

/* d <- a^2 */
static inline void
gf_sqr(gf *d, const gf *a)
{
	uint64_t r0, r1, r2, r3, r4, r5, r6, r7;

	/* Polynomial squaring (i.e. bit expansion) */
	r0 = expand_lo(a->w[0] & 0x000000001FFFFFFF);
	r1 = expand_lo(a->w[0] >> 29);
	r2 = expand_lo(a->w[1] & 0x000000001FFFFFFF);
	r3 = expand_lo(a->w[1] >> 29);
	r4 = expand_lo(a->w[2] & 0x000000001FFFFFFF);
	r5 = expand_lo(a->w[2] >> 29);
	r6 = expand_lo(a->w[3] & 0x000000001FFFFFFF);
	r7 = expand_lo(a->w[3] >> 29);

	/* Reduce the value. */
	REDUCE_465;

	d->w[0] = r0;
	d->w[1] = r1;
	d->w[2] = r2;
	d->w[3] = r3;
}

/* Normalize a into d (i.e. full reduction in the field). */
static inline void
gf_norm(gf *d, const gf *a)
{
	/* Values are always kept in normalized format. */
	*d = *a;
}

static inline uint64_t
squeeze(uint64_t x)
{
	x = (x & 0x1111111111111111) | ((x & 0x4444444444444444) >> 1);
	x = (x & 0x0303030303030303) | ((x & 0x3030303030303030) >> 2);
	x = (x & 0x000F000F000F000F) | ((x & 0x0F000F000F000F00) >> 4);
	x = (x & 0x000000FF000000FF) | ((x & 0x00FF000000FF0000) >> 8);
	x = (x & 0x000000000000FFFF) | ((x & 0x0000FFFF00000000) >> 16);
	return x;
}

/* d <- sqrt(a) */
static inline void
gf_sqrt(gf *d, const gf *a)
{
	/*
	 * We split a into "odd" and "even" parts:
	 *   a = ae + z*ao
	 * with:
	 *   ae = \sum_{i=0}^{127} a_{2*i)*z^{2i}
	 *   ao = \sum_{i=0}^{127} a_{2*i+1)*z^{2i}
	 * Then, we have:
	 *   sqrt(a) = sqrt(ae) + sqrt(z)*sqrt(ao)
	 * sqrt(ae) and sqrt(ao) are computed by "squeezing" words
	 * (odd-numbered digits are removed).
	 */
	uint64_t e0, e1, f0, f1;
	uint64_t r0, r1, r2, r3, r4, r5;

	/* e0:e1 <- sqrt(ae)  (e1 has size 59 bits)*/
	e0 = squeeze(a->w[0]) | (squeeze(a->w[1]) << 29);
	e1 = squeeze(a->w[2]) | (squeeze(a->w[3]) << 29);

	/* f0:f1 <- sqrt(ao) */
	f0 = squeeze(a->w[0] >> 1) | (squeeze(a->w[1] >> 1) << 29);
	f1 = squeeze(a->w[2] >> 1) | (squeeze(a->w[3] >> 1) << 29);

	/* sqrt(z) = z^228 + z^191 + z^154 + z^117 + z^69 + z^32 */
	r0 = e0;
	r1 = e1;
	r2 = (e1 >> 58);
	r3 = 0;
	r4 = 0;
	r5 = 0;

	r0 ^= f0 << 32;
	r1 ^= (f0 >> 26) ^ (f1 << 32);
	r2 ^= f1 >> 26;

	r1 ^= f0 << 11;
	r2 ^= (f0 >> 47) ^ (f1 << 11);
	r3 ^= f1 >> 47;

	r2 ^= f0 << 1;
	r3 ^= (f0 >> 57) ^ (f1 << 1);
	r4 ^= f1 >> 57;

	r2 ^= f0 << 38;
	r3 ^= (f0 >> 20) ^ (f1 << 38);
	r4 ^= f1 >> 20;

	r3 ^= f0 << 17;
	r4 ^= (f0 >> 41) ^ (f1 << 17);
	r5 ^= f1 >> 41;

	r3 ^= f0 << 54;
	r4 ^= (f0 >> 4) ^ (f1 << 54);
	r5 ^= f1 >> 4;

	r0 &= 0x03FFFFFFFFFFFFFF;
	r1 &= 0x03FFFFFFFFFFFFFF;
	r2 &= 0x03FFFFFFFFFFFFFF;
	r3 &= 0x03FFFFFFFFFFFFFF;
	r4 &= 0x03FFFFFFFFFFFFFF;
	/* not needed: r5 &= 0x03FFFFFFFFFFFFFF; */

	/*
	 * Maximum size of sqrt(ao)*z^228 is 344 bits, which makes reduction
	 * simpler than the post-multiplication case.
	 */
	r0 ^= (r5 & 1) << 57;
	r1 ^= r5 >> 1;
	r2 ^= (r5 << 15) & 0x03FFFFFFFFFFFFFF;
	r3 ^= r5 >> 43;
	r3 ^= (r4 & 1) << 58; \
	r4 &= ~(uint64_t)1;
	r0 ^= r4 >> 1;
	r1 ^= (r4 << 15) & 0x03FFFFFFFFFFFFFF;
	r2 ^= r4 >> 43;

	d->w[0] = r0;
	d->w[1] = r1;
	d->w[2] = r2;
	d->w[3] = r3;
}

/*
 * Trace (always 0 or 1).
 */
static inline uint64_t
gf_trace(const gf *a)
{
	/* Tr(a) = a_0 + a_159. */
	return (a->w[0] ^ (a->w[2] >> 43)) & 1;
}

/*
 * Precomputed half-traces of z^(2*i+1) for i = 0 to 115.
 */
static const gf HALFTRACE[] = {
	GFw64be(0x00000197FD722105, 0x5AABAEF4E0DC1096,
	        0x20B7FC13C12E8D69, 0xFD01D0208E170521),
	GFw64be(0x0000019395A41BE5, 0x109128B2E0851116,
	        0x9795725B45048071, 0x5C8585D4DA648BAD),
	GFw64be(0x0000010C2A9B2FEA, 0x1A2BEE86B8843672,
	        0x38D2053841A24074, 0x1902CF5E9217F50D),
	GFw64be(0x000000338A8713C4, 0x254C7AE360522FC4,
	        0x61A5035A3AC98C5D, 0x86126EAC1C0A1166),
	GFw64be(0x000001725848920E, 0xB0518410CE9A5281,
	        0xA2BCC8424D24A625, 0xB619963440E20137),
	GFw64be(0x0000010C68C348AA, 0x0A692C21421930C9,
	        0x2E92898009218751, 0xF596806D30877774),
	GFw64be(0x000001D6F5179967, 0x760DEE3DF1DD38A2,
	        0x323EF75BE3A7DCED, 0xACD2A1364347E213),
	GFw64be(0x000001825B64FFE5, 0x262772A288061B2A,
	        0x4C50DAFB23C00154, 0x1A04E7C6670CFD0D),
	GFw64be(0x0000012B062F72CE, 0x4F6A4CB8300C3A3E,
	        0x2D9123CAB8A44010, 0x3946A716EA5479F0),
	GFw64be(0x00000034C83C8D45, 0x00E1E8839B1E46D3,
	        0x34C6863B098156C4, 0xFA5C0B69829ACC79),
	GFw64be(0x0000006DC7D31CA5, 0xF71249CE72F775AB,
	        0xB82BBD63F49ACF2A, 0xCEDFDDF67383E042),
	GFw64be(0x000000B0F1628609, 0x42190624AC380C95,
	        0xA024D82185222546, 0xE0043C318D731517),
	GFw64be(0x000000B0B1608609, 0x41190624AC280C15,
	        0xA024582195222146, 0x60007C118BF31517),
	GFw64be(0x000001130E386A29, 0xAE073201BA7A15D8,
	        0x5DF5248123414E5A, 0xE351CDFDAD5CB8EB),
	GFw64be(0x000000C9AB45864C, 0xBDC8C6AC5AAF7402,
	        0x0F491B2A78A2C273, 0x6FDBCC0127156FD8),
	GFw64be(0x000000FBB1A6E508, 0x104F346FBA9D43BF,
	        0x4ACD52B04B6B4660, 0xE88953A73AC92769),
	GFw64be(0x000001E13FBF06E3, 0x76E0A4F906E37771,
	        0x2E597729E82D2A3F, 0x528F8ED9D0D63E0C),
	GFw64be(0x00000016CDB3CB21, 0x1BD9A49184B028EC,
	        0xA866B5915D252427, 0x911265623863754B),
	GFw64be(0x00000165AD0EFACF, 0xCE93CE63CB721128,
	        0xD0FB32CBA5A99E9B, 0x9754C4F3247DC02F),
	GFw64be(0x00000123A6B0E0EC, 0xC1FB0822A02C180A,
	        0xB511248A9D000007, 0x3841B1D17336D8C5),
	GFw64be(0x000000E85F8DE429, 0x783A080144BB6E9C,
	        0x8AE8F3A1C401A476, 0xB7CA6960FD60266A),
	GFw64be(0x000000B3181963A8, 0xBA9E305E3ABF298A,
	        0xE2C54590464E4676, 0xBB866DE5633A023C),
	GFw64be(0x000001E47516C9C3, 0x6201E697FC3F58EE,
	        0x323AF69981A7E402, 0xAEC8803E6F57D347),
	GFw64be(0x000000826248024D, 0x410040135A1512A3,
	        0x1240880B9085C645, 0x8C91472303959758),
	GFw64be(0x0000013BA26A572E, 0x0A0462E535277F3A,
	        0x32F508F2028B7183, 0x0BD8AF82A616836A),
	GFw64be(0x0000016C0E09236D, 0x52265C3D20DC75A2,
	        0x743A211BC2E70C28, 0xF807CDF6035E8C12),
	GFw64be(0x000001C279C85167, 0x79CB29E68E72CBE4,
	        0x9718D8DBD91A2F1E, 0x875D76EA1874CB85),
	GFw64be(0x00000069B4970E05, 0x849430341E5C6BD4,
	        0xF189672326466E4D, 0xB89F36ADDC2A90A0),
	GFw64be(0x000000578EC33383, 0x052F1EA0F8516F95,
	        0x536F29513360CD18, 0xC69B3FB1DD99C69F),
	GFw64be(0x000000F0FE1C6ACB, 0x758BA8B8C6B7439F,
	        0x06ECE689F104A676, 0xD69907A1EAD04A2F),
	GFw64be(0x0000010C68C348AA, 0x0A692C21C21930C9,
	        0x2A92898009218751, 0xD596806D30877775),
	GFw64be(0x000000E4DDDE14CB, 0x8DE55C8CAF4E0912,
	        0x180AFE693BE23BDD, 0xFB5571D19200A151),
	GFw64be(0x000001C39BCAF3A7, 0xA703860E8E0A5A19,
	        0x313958D331222342, 0x6258B744B5C6D5C3),
	GFw64be(0x000000981724FB4E, 0xD44AA4184E85097B,
	        0x26A472DAE824A22D, 0x0C85748AA2924B26),
	GFw64be(0x0000003B63169906, 0x8AB90AF4E65C00C6,
	        0xA4C59652050EAE7C, 0xED41002D0F325D69),
	GFw64be(0x00000121CD9AFF4D, 0xA59E1CD35A2A5029,
	        0xCC51B4FB366DC282, 0x3658951270AD6C58),
	GFw64be(0x000001786534E22D, 0x87FF3875DE0C5F1A,
	        0xDD3CB6833F4FE004, 0x3C08BFD5E62CEC83),
	GFw64be(0x0000011F867DDBAB, 0x6A8C945E3C3F3DDA,
	        0x66972FD1826E6543, 0xABC6BDE8E21F5A74),
	GFw64be(0x000001636BB8712E, 0x651A800A26325289,
	        0x913994D2B4002607, 0x821C827530B4C182),
	GFw64be(0x000001572982BF0B, 0xD0BF9A4BD49C53D6,
	        0xCCFF1071C749E470, 0xFC48C7B99A2D687F),
	GFw64be(0x000001B1D1E3EC21, 0x036B0841E85434EE,
	        0x2CF5D9A119098D09, 0x8C16896A3B066C6B),
	GFw64be(0x000001E6F7B35B62, 0xEA886092103F78CE,
	        0x1A9AF5D880844407, 0xAB9AB47D7B14E324),
	GFw64be(0x0000000000000000, 0x0000000080000000,
	        0x0000000000000000, 0x0000000000000001),
	GFw64be(0x0000016AC624EEAC, 0xC2CFBC9EB6005339,
	        0x48F8A2A28B666200, 0x0049C7D6A5DC702B),
	GFw64be(0x00000179EAC21AED, 0x366846CD329B5895,
	        0x169D884B68AB4734, 0xB28DB461DDC5DA64),
	GFw64be(0x00000166C5BCAEEE, 0xB8D20C564AFA5FB5,
	        0xB5FAB62A4C2E8E3A, 0xF719FFF299A799BA),
	GFw64be(0x000000D015E2998E, 0x8D04B40F34504FC1,
	        0x6CCC785232636D19, 0x81097AB8109A6D28),
	GFw64be(0x000001F7B967ADA2, 0x5826184E789C73F6,
	        0x605F5B30C24AC531, 0xFD1BD3EA8E4E501C),
	GFw64be(0x000001270C10A70C, 0x18B2A8D04CE41C02,
	        0x80932432440CA82E, 0x1D05EC1003651120),
	GFw64be(0x000001AB5DBF1A02, 0x75DCB2007CEF3542,
	        0xF311F740FE40E82E, 0x6E86DC98173F8680),
	GFw64be(0x000000D095F2998E, 0x0904B40F34704FC5,
	        0x6CCC7C5212636D1A, 0x81097AB9189A6D28),
	GFw64be(0x0000018B2020C7C7, 0x1D982C5C62C80D65,
	        0xA4B100BB742E8B2D, 0x65446D8B08B30966),
	GFw64be(0x0000011A8E19E3EB, 0x66AF9096BC7E34DE,
	        0x64942599A3466C5B, 0xBA06DD7CEA1F1965),
	GFw64be(0x00000069A4978E05, 0x849434341E5C6BD4,
	        0xF589272326666E4D, 0xB81E36ADDC2A98A0),
	GFw64be(0x0000008EDA97306F, 0xCF930089CA3521D7,
	        0x8EC2C74BB5018602, 0x9D9711E89AF13E6D),
	GFw64be(0x000001652FCB0FAF, 0x5241485EB1DF35A7,
	        0x5EBB3933C98E5DB9, 0xBAC38DA31BCDEB63),
	GFw64be(0x00000008A09F9222, 0x8E91449A90202334,
	        0x9880074025A44007, 0x010716C39960A135),
	GFw64be(0x000001692F5D0A2D, 0x3B4C9AEE209774A6,
	        0x42193F035A4A0530, 0xDBCBDC624B1D4714),
	GFw64be(0x000000887B0FFAAE, 0x30273AAADCE32D4C,
	        0x4600D3C24340E82A, 0x56C76DC86C085F55),
	GFw64be(0x0000005B031990E1, 0xCA9474CD6C1865A4,
	        0xF4AD154986EBA404, 0xE54A08E7083A9D72),
	GFw64be(0x00000094AFFBC96C, 0x5AD9E482DEA52A4F,
	        0xBA263D9ACDA0E222, 0x1D8322493BF6E213),
	GFw64be(0x000000BC333377CB, 0xB4C252A347DD6B9B,
	        0x72C655F968C1BEED, 0xAC9B66E0E69B9238),
	GFw64be(0x0000004DC393DC20, 0xF7B840EE7AF6658B,
	        0xB86B95E0F48ACE2A, 0xCE5F58E473F3E04A),
	GFw64be(0x0000019A7EB4B107, 0xAF3451800A6B065D,
	        0x5874E65336D00A1B, 0x33C40E1CFDCDB05E),
	GFw64be(0x0000011BECDF9C0B, 0x34C5B424D0C73955,
	        0x5335AF616B62C968, 0x0E86A0899CCC9397),
	GFw64be(0x000001428B65418B, 0xA050EAF03867796F,
	        0xABF81B910C8C495E, 0x0A9FA18A2FA666BE),
	GFw64be(0x000001889B2888E6, 0xF27490CD18BF1C85,
	        0xDBD0500ACE4B4472, 0xEA91B8740CFCA6AC),
	GFw64be(0x0000000082C5F560, 0x9CD2863D7AAD240F,
	        0x8E400BF86C27C322, 0x2D9749503FB02B4C),
	GFw64be(0x0000018170D0EE04, 0x437692EF642C183A,
	        0xC8F1CCA29E4BA146, 0x6C01B502B67D617E),
	GFw64be(0x000001BF760F5922, 0xA895925A9027354A,
	        0xC697E3D0074C4006, 0x4B8288D9676D4E35),
	GFw64be(0x00000187915920E5, 0xBEC6304672FA34A1,
	        0x40335D0B6A4ACF6E, 0xA743C962408D1152),
	GFw64be(0x00000058077CE060, 0xD6FF786684194DB9,
	        0xCF6C3E88EFCA2505, 0xB49978F6E5B97BDB),
	GFw64be(0x0000011505A67208, 0xA57FF450C03416FA,
	        0xF89732C03FEC8446, 0xCC049E3BF23FA021),
	GFw64be(0x0000010CEA06BDCA, 0x96BBAA1CB8BC14C6,
	        0xA4D2827865064473, 0xF801C93D0F375C39),
	GFw64be(0x0000007CEDC5AA85, 0x0F9614541E0D6F53,
	        0xF7AEBB03366E6344, 0x29DA3F9DD2AB9BF6),
	GFw64be(0x00000063646A7966, 0x5C46729262C54D7F,
	        0x6BA9A8DAEAC48B79, 0x59D878DBEA8B26B2),
	GFw64be(0x000000748F325FE7, 0x1144926964B04BDC,
	        0x79EE34FB5A49A477, 0xC40976E8E91BB0BA),
	GFw64be(0x000001000FD20848, 0x395CF2984CF71FF6,
	        0xFA903C085EC4AD7E, 0xDBC4FAAEDA3FE730),
	GFw64be(0x000000B1DE31E469, 0x3A37FAD984A32B18,
	        0xD665E5A947CD2077, 0x57D273C4E1398B5B),
	GFw64be(0x0000005E0F64ADF1, 0xC497AC86FE556E12,
	        0xE20E3A3DA722EF18, 0x9C8E6B51826B4344),
	GFw64be(0x00000066ECBFC146, 0xD8E452190E84654B,
	        0x702AA79ACAC52274, 0x194B098837DB8106),
	GFw64be(0x0000011F067DDBAB, 0xEE8C14D63C1F3DDE,
	        0x66B72FD1A26E6540, 0xABC6BDE9EA1E5B62),
	GFw64be(0x000000D7145F42CF, 0xBBB90413A2AE6B5B,
	        0x990F6F8B55250333, 0x7B0F7698A2A1F5D1),
	GFw64be(0x000001CF139F4225, 0x6388A0D1084F721F,
	        0x0A7B5783900D0859, 0x2ACFC751FA84364A),
	GFw64be(0x0000011FE6EA498B, 0x201C105F3C271D6A,
	        0xEAB7A891064F6147, 0x4AD1B88A732E7676),
	GFw64be(0x000001F0184FE5C3, 0xA30BACE3EE4A7759,
	        0x199C4BB91129AB1D, 0x660A9A88A1C4E4F1),
	GFw64be(0x0000004DAA0C4160, 0x9F54D2D26CD04D8C,
	        0xE98B02987ECCAC29, 0xC51C68E56C2B65A4),
	GFw64be(0x0000011F866DDBAB, 0x6A8C145E3C3F3DDA,
	        0x66B72BD1826E6543, 0xABC6BDE8E21E5B62),
	GFw64be(0x000001E6E7F3DB62, 0xEA8A6492103F78CE,
	        0x1C9ABDD880A44407, 0xAB9BB47D7B10E924),
	GFw64be(0x0000000568F3AA60, 0x46B382C9805929B4,
	        0x8E238D08A5090D1C, 0xF0D765F699316F13),
	GFw64be(0x00000069A4978E05, 0x849430349E5C6BD4,
	        0xF589272326466E4D, 0xB81E36ADDC2A98A1),
	GFw64be(0x0000017EC2D8D5AC, 0x6C98A46A247F568E,
	        0x975E8CF2A4282D1F, 0xFB9D8B603B35DA9C),
	GFw64be(0x000001D637D26C07, 0xE8DF68008B601C2D,
	        0xBC7E7CA38F801BCF, 0x0145E84678F7E95F),
	GFw64be(0x000000450A93D342, 0x11C5B648FCF06EB8,
	        0x710B05D85B68EC2E, 0xC41B7E26F54BC491),
	GFw64be(0x00000003E979E481, 0xF98366709697279F,
	        0x0AE19DA1D1AC6730, 0xDB821BB4EEC1373F),
	GFw64be(0x00000197FD622105, 0x5AABAEF4E0DC1096,
	        0x20B7F813C12E8D69, 0xFD01D0208E160437),
	GFw64be(0x000001576BDAD84B, 0xC0FC58EC2E01556D,
	        0xDA3F9CC98ECA2355, 0x10DC988A38FDE206),
	GFw64be(0x000001724808120E, 0xB0518010CE9A5281,
	        0xA6BC80424D04A625, 0xB618963440E60B37),
	GFw64be(0x000001C5D9F16D26, 0x82AE146E7546330E,
	        0x645BDDB2026AF9D9, 0x1E16D2813B5E485C),
	GFw64be(0x00000048AB607640, 0xB103C227FAAF4C3F,
	        0x364818E85183C336, 0x3E883D12EAC28F59),
	GFw64be(0x0000007E0BA278C4, 0xCD25B717F6545EB3,
	        0x6A4E10CAB377EE19, 0x9D0CAE22929A6319),
	GFw64be(0x000001724848120E, 0xB07184104E9B5281,
	        0xA6BC88424D24A725, 0xB698962440E60B36),
	GFw64be(0x000000244ED9EF65, 0x18C2E6A6729D26DA,
	        0x1E02ADBB48A2C764, 0xFD865E78B601FF04),
	GFw64be(0x000000F659F05D4A, 0x9E6C34A73CC04A0B,
	        0x4C8EDCF86A63696D, 0x110C331437C97C30),
	GFw64be(0x000000B073A77369, 0xDFCB8019D695289A,
	        0x2E64D3D9F905E664, 0xCD933561B2433E5B),
	GFw64be(0x0000005747028203, 0xA9A9FAEFEF4C0417,
	        0x748FB00111CBBA89, 0x6D445D44CE8A8821),
	GFw64be(0x0000016F63EDC0AC, 0x21912655A4777A8B,
	        0xAA5B9B82152F2D1F, 0xCA8BF26177F7231D),
	GFw64be(0x000001C29BCAF3A7, 0xA703860E8E0A5A19,
	        0x313858D331222340, 0x6358B744B5C6D5C3),
	GFw64be(0x000000DA6534E22D, 0x87FF3875DE0C5F1A,
	        0xDD2CB6833F4FE340, 0x3808BFD5E62CEC83),
	GFw64be(0x00000116ECDF9689, 0x7CFC922056A33165,
	        0xF276AF61EE40E376, 0x1786858F1CFE971E),
	GFw64be(0x000000FB53611068, 0x8E9DB252C0206730,
	        0xC48DD948274C8142, 0x051A1AD781790C25),
	GFw64be(0x000001EB10A78B80, 0xD53DBAA0FED17FC5,
	        0xD7194310F740EE2D, 0xD4DBEFBD5CBD9FC1),
	GFw64be(0x000001F0F1CE6282, 0x4CD75A2082405C69,
	        0xFC3CDA80AFC00B08, 0x055CFD0F31EFBD57),
	GFw64be(0x00000123E4E887AE, 0xD1B8CA854AB11EB1,
	        0xA3D1A832D4838722, 0xD4D5EE62D1E652BC),
	GFw64be(0x0000002703C2AF40, 0x21EE5A48001303F4,
	        0x7A0318381AC80550, 0xC6C116BED95AF610),
        GFw64be(0x000001C39639F66F, 0x8A667AC94A397999,
                0x70D965EB0AC98652, 0xE1CAE5F0A1DF806C)
};

/*
 * d <- half-trace(a)
 */
static void
gf_halftrace(gf *d, const gf *a)
{
	/*
	 * We split a into its even and odd halves (as in gf_sqrt()):
	 *   a = ae + z*ao
	 * Then:
	 *   H(a) = H(ae) + H(z*ao)
	 * Since H(x) = H(sqrt(x)) + sqrt(x) for all x, we can replace H(ae):
	 *   H(a) = H(sqrt(ae)) + H(z*ao) + sqrt(ae)
	 * sqrt(ae) is obtained through simple squeezing, and has half-size,
	 * so it can be split again, recursively. We thus remove all
	 * even-indexed bits from the computation, and thus can use a
	 * half-size table for the matrix (HALFTRACE[] array).
	 */
	uint64_t ao[4], e0, e1;
	uint64_t d0, d1, d2, d3;

	ao[0] = a->w[0];
	ao[1] = a->w[1];
	ao[2] = a->w[2];
	ao[3] = a->w[3];

	/*
	 * We accumulate the odd bits into ao[] (z*ao is already there).
	 * The running value is kept in e0:e1 (only e0 after the first
	 * iteration), and the square roots are accumulated into d0..d3.
	 */

	/* H(a) = H(sqrt(ae)) + H(z*ao) + sqrt(ae) */
	e0 = squeeze(a->w[0]) | (squeeze(a->w[1]) << 29);
	e1 = squeeze(a->w[2]) | (squeeze(a->w[3]) << 29);
	d0 = e0;
	d1 = e1 & 0x03FFFFFFFFFFFFFF;
	d2 = e1 >> 58;
	d3 = 0;

	/* size(e0:e1) = 117 */

	/* H(a) = H(sqrt(aee)) + H(z*(ao + aeo)) + sqrt(ae) + sqrt(aee) */
	ao[0] ^= e0;
	ao[1] ^= e1;
	e0 = squeeze(e0) | (squeeze(e1) << 29);
	d0 ^= e0 & 0x03FFFFFFFFFFFFFF;
	d1 ^= e0 >> 58;

	/* size(e0) = 59 */

	/*
	 * The running value fits on a single word. In 6 iterations, it
	 * is shrunk to a single bit. Since H(1) = 1, we can stop there.
	 */
	for (int i = 0; i < 6; i ++) {
		ao[0] ^= e0;
		e0 = squeeze(e0);
		d0 ^= e0;
	}
	d0 ^= e0;

	/*
	 * We now have to compute the halftrace of all the odd-indexed
	 * bits in ao[] (ignoring the even-indexed bits).
	 */
	for (int i = 0; i < 4; i ++) {
		uint64_t mw = ao[i] >> 1;
		for (int j = 0; j < 29; j ++) {
			uint64_t m;

			m = -(mw & 1);
			mw >>= 2;
			d0 ^= m & HALFTRACE[29 * i + j].w[0];
			d1 ^= m & HALFTRACE[29 * i + j].w[1];
			d2 ^= m & HALFTRACE[29 * i + j].w[2];
			d3 ^= m & HALFTRACE[29 * i + j].w[3];
		}
	}

	d->w[0] = d0;
	d->w[1] = d1;
	d->w[2] = d2;
	d->w[3] = d3;
}

/* d <- a^(2^n) */
static void
gf_xfrob(gf *d, const gf *a, unsigned n)
{
	if (n == 0) {
		*d = *a;
		return;
	}
	gf_sqr(d, a);
	n --;
	while (n -- > 0) {
		gf_sqr(d, d);
	}
}

static inline void
gf_xfrob_tab(gf *d, const gf *a, const gf *xft)
{
	uint64_t d0, d1, d2, d3, m;

	d0 = 0;
	d1 = 0;
	d2 = 0;
	d3 = 0;
	for (int i = 0; i < 4; i ++) {
		uint64_t aw = a->w[i];
		for (int j = 0; j < 58; j ++) {
			m = -(aw & 1);
			aw >>= 1;
			d0 ^= m & xft[58 * i + j].w[0];
			d1 ^= m & xft[58 * i + j].w[1];
			d2 ^= m & xft[58 * i + j].w[2];
			d3 ^= m & xft[58 * i + j].w[3];
		}
	}
	m = -(a->w[3] >> 58);
	d0 ^= m & xft[232].w[0];
	d1 ^= m & xft[232].w[1];
	d2 ^= m & xft[232].w[2];
	d3 ^= m & xft[232].w[3];

	d->w[0] = d0;
	d->w[1] = d1;
	d->w[2] = d2;
	d->w[3] = d3;
}

#include "gf233_frob58.h"
#include "gf233_frob116.h"

#define gf_xfrob58(d, a)    gf_xfrob_tab(d, a, GF233_XFROB_58)
#define gf_xfrob116(d, a)   gf_xfrob_tab(d, a, GF233_XFROB_116)

/* d <- 1/a  (if a == 0 then d is set to zero) */
static void
gf_inv(gf *d, const gf *a)
{
	gf x, y;

	/* a^(2^2-1) */
	gf_sqr(&y, a);
	gf_mul(&x, &y, a);

	/* a^(2^3-1) */
	gf_sqr(&y, &x);
	gf_mul(&x, &y, a);

	/* a^(2^6-1) */
	gf_xfrob(&y, &x, 3);
	gf_mul(&x, &y, &x);

	/* a^(2^7-1) */
	gf_sqr(&y, &x);
	gf_mul(&x, &y, a);

	/* a^(2^14-1) */
	gf_xfrob(&y, &x, 7);
	gf_mul(&x, &y, &x);

	/* a^(2^28-1) */
	gf_xfrob(&y, &x, 14);
	gf_mul(&x, &y, &x);

	/* a^(2^29-1) */
	gf_sqr(&y, &x);
	gf_mul(&x, &y, a);

	/* a^(2^58-1) */
	gf_xfrob(&y, &x, 29);
	gf_mul(&x, &y, &x);

	/* a^(2^116-1) */
	gf_xfrob58(&y, &x);
	gf_mul(&x, &y, &x);

	/* a^(2^232-1) */
	gf_xfrob116(&y, &x);
	gf_mul(&x, &y, &x);

	/* a^(2^233-2) = 1/a */
	gf_sqr(d, &x);
}

/* d <- x/y (returns 0 if y == 0) */
static inline void
gf_div(gf *d, const gf *x, const gf *y)
{
	gf t;

	gf_inv(&t, y);
	gf_mul(d, x, &t);
}

/* Returned value is -1 if a is zero, 0 otherwise. */
static inline uint64_t
gf_iszero(const gf *a)
{
	uint64_t r;

	r = a->w[0] | a->w[1] | a->w[2] | a->w[3];
	return ((r | -r) >> 63) - 1;
}

/* Returned value is -1 if a == b, 0 otherwise. */
static inline uint64_t
gf_equals(const gf *a, const gf *b)
{
	gf t;
	uint64_t r;

	gf_add(&t, a, b);
	r = t.w[0] | t.w[1] | t.w[2] | t.w[3];
	return ((r | -r) >> 63) - 1;
}

/* Decode d from exactly 30 bytes. A failure is reported if the source
   is not canonical (bits 1-7 in the last byte were not all equal to zero).
   On success, -1 is returned, or 0 on failure. On failure, *d is set
   to zero. */
static inline uint64_t
gf_decode(gf *d, const void *src)
{
	const uint8_t *buf;
	uint64_t w[4];

	buf = src;
	for (int i = 0; i < 3; i ++) {
		w[i] = (uint64_t)buf[8 * i + 0]
			| ((uint64_t)buf[8 * i + 1] << 8)
			| ((uint64_t)buf[8 * i + 2] << 16)
			| ((uint64_t)buf[8 * i + 3] << 24)
			| ((uint64_t)buf[8 * i + 4] << 32)
			| ((uint64_t)buf[8 * i + 5] << 40)
			| ((uint64_t)buf[8 * i + 6] << 48)
			| ((uint64_t)buf[8 * i + 7] << 56);
	}
	w[3] = (uint64_t)buf[24 + 0]
		| ((uint64_t)buf[24 + 1] << 8)
		| ((uint64_t)buf[24 + 2] << 16)
		| ((uint64_t)buf[24 + 3] << 24)
		| ((uint64_t)buf[24 + 4] << 32)
		| ((uint64_t)buf[24 + 5] << 40);
	d->w[0] = w[0] & 0x03FFFFFFFFFFFFFF;
	d->w[1] = ((w[0] >> 58) | (w[1] << 6)) & 0x03FFFFFFFFFFFFFF;
	d->w[2] = ((w[1] >> 52) | (w[2] << 12)) & 0x03FFFFFFFFFFFFFF;
	d->w[3] = (w[2] >> 46) | (w[3] << 18);

	/*
	 * Value is valid as long as the extra bits are zero.
	 */
	uint64_t r = ((const uint8_t *)buf)[29];
	r = ((r + 0xFE) >> 8) - 1;
	d->w[0] &= r;
	d->w[1] &= r;
	d->w[2] &= r;
	d->w[3] &= r;
	return r;
}

/* Encode a into exactly 30 bytes. */
static inline void
gf_encode(void *dst, const gf *a)
{
	uint64_t w[4];
	uint8_t tmp[32];

	w[0] = a->w[0] | (a->w[1] << 58);
	w[1] = (a->w[1] >> 6) | (a->w[2] << 52);
	w[2] = (a->w[2] >> 12) | (a->w[3] << 46);
	w[3] = a->w[3] >> 18;
	for (int i = 0; i < 4; i ++) {
		tmp[8 * i + 0] = (uint8_t)w[i];
		tmp[8 * i + 1] = (uint8_t)(w[i] >> 8);
		tmp[8 * i + 2] = (uint8_t)(w[i] >> 16);
		tmp[8 * i + 3] = (uint8_t)(w[i] >> 24);
		tmp[8 * i + 4] = (uint8_t)(w[i] >> 32);
		tmp[8 * i + 5] = (uint8_t)(w[i] >> 40);
		tmp[8 * i + 6] = (uint8_t)(w[i] >> 48);
		tmp[8 * i + 7] = (uint8_t)(w[i] >> 56);
	}
	memcpy(dst, tmp, 30);
}

/* d <- a0 if ctl == 0
   d <- a1 if ctl == -1
   ctl MUST be either 0 or -1 */
static inline void
gf_select(gf *d, const gf *a0, const gf *a1, uint64_t ctl)
{
	d->w[0] = a0->w[0] ^ (ctl & (a0->w[0] ^ a1->w[0]));
	d->w[1] = a0->w[1] ^ (ctl & (a0->w[1] ^ a1->w[1]));
	d->w[2] = a0->w[2] ^ (ctl & (a0->w[2] ^ a1->w[2]));
	d->w[3] = a0->w[3] ^ (ctl & (a0->w[3] ^ a1->w[3]));
}

/* d <- a if ctl == -1, but unchanged if ctl == 0 */
static inline void
gf_condset(gf *d, const gf *a, uint64_t ctl)
{
	gf_select(d, d, a, ctl);
}

/* d <- d + a if ctl == -1, but unchanged if ctl == 0
   ctl MUST be either 0 or -1 */
static inline void
gf_condadd(gf *d, const gf *a, uint64_t ctl)
{
	d->w[0] ^= ctl & a->w[0];
	d->w[1] ^= ctl & a->w[1];
	d->w[2] ^= ctl & a->w[2];
	d->w[3] ^= ctl & a->w[3];
}

#endif
