#ifndef XS233_H__
#define XS233_H__

#include <stddef.h>
#include <stdint.h>

/* ===================================================================== */
/*
 * The xsb233 and xsk233 groups are prime order groups isomorphic to the
 * prime order subgroup of, respectively, the standard NIST curves B-233
 * and K-233 (also known as sect233r1 and sect233k1, respectively). Both
 * groups follow very similar APIs.
 *
 * All group elements (including the neutral) can be encoded into to a
 * fixed-sized 30-byte format. The top 7 seven bits of the last byte
 * are always zero. Encoding is always canonical, and is verified upon
 * decoding (including the zero-bits in the last byte).
 *
 * In the following, group elements are also called "points" since they
 * correspond to points on a curve.
 *
 * The order of xsb233 is slightly above 2^232. The order of xsk233 is
 * slightly above 2^231. Both groups are considered secure for use in
 * cryptographic protocols, with a security level of (at least) "112 bits".
 *
 * Except where explicitly stated, all functions are constant-time, i.e.
 * acceptable for handling secret data.
 */

/*
 * Type for an in-memory representation of a group element.
 * Contents are opaque.
 */
typedef struct {
	uint64_t opaque[16];
} xsb233_point;
typedef struct {
	uint64_t opaque[16];
} xsk233_point;

/*
 * Statically allocated neutral element.
 */
extern const xsb233_point xsb233_neutral;
extern const xsk233_point xsk233_neutral;

/*
 * Statically allocated conventional generator (this generator corresponds
 * to the conventional generator defined by SEC and NIST for B-233 and K-233).
 */
extern const xsb233_point xsb233_generator;
extern const xsk233_point xsk233_generator;

/*
 * Decode a group element from its double-odd representation. This
 * representation always has size exactly 30 bytes. The last byte has
 * value 0x00 or 0x01 (there are thus 7 "free bits" in the last byte,
 * that can be used to encode other data, but this function expects and
 * checks that these bits have value zero, so any reuse of these 7 bits
 * must be handled externally and the value adjusted before calling this
 * function).
 *
 * Returned value is 0xFFFFFFFF on success, 0x00000000 on error. On
 * error, the `*p` structure is set to the group neutral element.
 */
uint32_t xsb233_decode(xsb233_point *p, const void *src);
uint32_t xsk233_decode(xsk233_point *p, const void *src);

/*
 * Encode a group element into bytes. This always writes exactly 30 bytes
 * into the provided `dst` buffer.
 */
void xsb233_encode(void *dst, const xsb233_point *p);
void xsk233_encode(void *dst, const xsk233_point *p);

/*
 * Check whether a group element is the neutral. Returned value is
 * 0xFFFFFFFF for the neutral, 0x00000000 for any other point.
 */
uint32_t xsb233_is_neutral(const xsb233_point *p);
uint32_t xsk233_is_neutral(const xsk233_point *p);

/*
 * Check whether two group elements are equal to each other.
 * Returned value is 0xFFFFFFFF on equality, 0x00000000 otherwise.
 */
uint32_t xsb233_equals(const xsb233_point *p1, const xsb233_point *p2);
uint32_t xsk233_equals(const xsk233_point *p1, const xsk233_point *p2);

/*
 * Add two points together: p3 <- p1 + p2
 */
void xsb233_add(xsb233_point *p3,
	const xsb233_point *p1, const xsb233_point *p2);
void xsk233_add(xsk233_point *p3,
	const xsk233_point *p1, const xsk233_point *p2);

/*
 * Subtract one point from another: p3 <- p1 - p2
 */
void xsb233_sub(xsb233_point *p3,
	const xsb233_point *p1, const xsb233_point *p2);
void xsk233_sub(xsk233_point *p3,
	const xsk233_point *p1, const xsk233_point *p2);

/*
 * Double a point: p3 <- 2*p1
 */
void xsb233_double(xsb233_point *p3, const xsb233_point *p1);
void xsk233_double(xsk233_point *p3, const xsk233_point *p1);

/*
 * Double a point n times: p3 <- (2^n)*p1
 * If n == 0, then this simply copies the source point into the destination.
 * The execution time of this function grows linearly with n. The value
 * n is not considered secret.
 */
void xsb233_xdouble(xsb233_point *p3, const xsb233_point *p1, unsigned n);
void xsk233_xdouble(xsk233_point *p3, const xsk233_point *p1, unsigned n);

/*
 * Negate a point: p3 <- -p1
 */
void xsb233_neg(xsb233_point *p3, const xsb233_point *p1);
void xsk233_neg(xsk233_point *p3, const xsk233_point *p1);

/*
 * Select a point value:
 *   p3 <- p0 if ctl == 0x00000000
 *   p3 <- p1 if ctl == 0xFFFFFFFF
 * ctl MUST be equal to 0x00000000 or 0xFFFFFFFF
 */
void xsb233_select(xsb233_point *p3,
	const xsb233_point *p0, const xsb233_point *p1, uint32_t ctl);
void xsk233_select(xsk233_point *p3,
	const xsk233_point *p0, const xsk233_point *p1, uint32_t ctl);

/*
 * Negate a point conditionally:
 *   p3 <- p1 if ctl == 0x00000000
 *   p3 <- -p1 if ctl == 0xFFFFFFFF
 * ctl MUST be equal to 0x00000000 or 0xFFFFFFFF
 */
void xsb233_condneg(xsb233_point *p3, const xsb233_point *p1, uint32_t ctl);
void xsk233_condneg(xsk233_point *p3, const xsk233_point *p1, uint32_t ctl);

/*
 * Multiply a point by an integer: p3 <- n*p3
 * The integer n is provided in unsigned little-endian encoding, of
 * length n_len bytes (with n_len <= 30); this encoding is not required
 * to be of minimal size (i.e. there may be trailing zeros). The
 * execution time is linear in n_len but does not depend upon the value
 * of the bytes. The size of the scalar MUST NOT exceed 30 bytes.
 */
void xsb233_mul(xsb233_point *p3,
	const xsb233_point *p1, const void *n, size_t n_len);
void xsk233_mul(xsk233_point *p3,
	const xsk233_point *p1, const void *n, size_t n_len);

/*
 * xsb233_mul_ladder() and xsk233_mul_ladder() implement the same
 * functionality as xsb233_mul() and xsk233_mul(), respectively, using
 * internally a Montgomery ladder with López-Dahab formulas. They are
 * provided as externally callable functions for benchmark purposes.
 */
void xsb233_mul_ladder(xsb233_point *p3,
	const xsb233_point *p1, const void *n, size_t n_len);
void xsk233_mul_ladder(xsk233_point *p3,
	const xsk233_point *p1, const void *n, size_t n_len);

/*
 * xsk233_mul_frob() implements the same functionality as xsk233_mul(),
 * but internally uses the Frobenius endomorphism on the underlying
 * Koblitz curve for improved performance.
 */
void xsk233_mul_frob(xsk233_point *p3,
	const xsk233_point *p1, const void *n, size_t n_len);

/*
 * Multiply the conventional generator by an integer: p3 <- n*G
 * The rules are the same as with xs?233_mul(), but this function is
 * typically more efficient than calling the generic multiplication
 * routine on the generator point.
 */
void xsb233_mulgen(xsb233_point *p3, const void *n, size_t n_len);
void xsk233_mulgen(xsk233_point *p3, const void *n, size_t n_len);

/*
 * xsk233_mulgen_frob() implements the same functionality as xsk233_mulgen(),
 * but internally uses the Frobenius endomorphism on the underlying
 * Koblitz curve for improved performance.
 */
void xsk233_mulgen_frob(xsk233_point *p3, const void *n, size_t n_len);

/*
 * TODO: add xs?233_scalar_*() functions for computations modulo the
 * curve order.
 */

/* ===================================================================== */

#endif
