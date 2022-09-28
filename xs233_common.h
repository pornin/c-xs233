/*
 * Common code for xsb233 and xsk233.
 */

/*
 * Field elements.
 */
#include "gf233.h"

/*
 * A non-supersingular curve in a binary field has short Weierstraß
 * equation:
 *   y^2 + x*y = x^3 + A*x^2 + B
 * for some constants A and B. The curve is well-defined as long as B != 0.
 * Curve isomorphism are mappings (x,y) -> (x, y + c*x) for any field
 * element c; this modifies the A constant into A + c^2 + c. This
 * implies that two curves are isomorphic to each other if and only if
 * they use the same B, and their A constants have the same trace. In
 * particular, we can always arrange for A to be either zero (if the
 * original A is such that trace(A) == 0), or a trace-1 element with very
 * low Hamming weight (when working in GF(2^m) for an odd m, as is the
 * case for NIST curves B-233 and K-233, trace(1) == 1, so we can make
 * A == 1 with an isomorphism). Standard curves usually already use this
 * technique. When trace(A) == 0, it is highly recommended to apply an
 * isomorphism to make A == 0, since that saves one multiplication in the
 * general point addition routine. When trace(A) == 1, it is a good idea
 * to set B == 1 through the isomorphism.
 *
 * B-233:
 *   A = 1
 *   B = some pseudo-randomly generated large constant
 *
 * K-233:
 *   A = 0
 *   B = 1
 *
 * We define a = A and b = sqrt(B). The change of variable y -> y + b maps
 * this curve to:
 *   y^2 + x*y = x^3 + a*x^2 + b*x
 * which is the curve that we use thereafter.
 * The point N = (0, 0) has order 2. We define the xs?233 group as the
 * set of points which are _not_ r-torsion points; i.e. the group consists
 * of points P+N for any P such that r*P = 0. The group neutral is N. The
 * addition in the group is: add(P1+N, P2+N) = (P1+P2)+N.
 *
 * We use (x,s) coordinates, with s = y + x^2 + a*x + b. For the point N,
 * the s coordinate is equal to b; for all other group elements, s = y^2/x.
 * In (x,s) coordinates, the curve equation becomes:
 *   s^2 + x*s = x^4 + a*x^2 + b^2
 *
 * It can be shown that:
 *  - x and s cannot be zero simultaneously.
 *  - x == 0 only for a single point in the group, which is the neutral N.
 *  - s can be zero only if trace(b/a^2) == 0; for curves that admit a
 *    point in the group such that s == 0, there can be only a single
 *    point in that case, and its x coordinate is non-zero.
 *
 * On group elements in (x,s) coordinates, the following formulas are complete:
 *
 *   x3 = b*(x1*x2 + s1*x2 + s2*x1) / (x1*x2 + b)^2
 *
 *   s3 = b*(b^2*(s1*s2 + a^2*x1*x2) + (x1*x2)^2*((a^2 + 1)*x1*x2
 *           + s1*x2 + s2*x1 + s1*s2)) / (x1*x2 + b)^4
 *
 * In practice, we use extended coordinates (X:S:Z:T) such that:
 *   Z != 0  (for all points, including the neutral)
 *   x = sqrt(b)*X/Z
 *   s = sqrt(b)*S/Z^2
 *   T = X*Z
 * Denoting by M, S and mb the costs of a multiplication, squaring,
 * and multiplication by the constant sqrt(b) respectively, then
 * we can have the following costs:
 *
 *   Point addition:           8M + 2S + 2mb
 *   with a == 0:              7M + 2S + 2mb
 *
 *   Point doubling:           3M + 4S + 4mb
 *
 * When sqrt(b) is not a "simple" constant (e.g. in the case of curve
 * B-233, sqrt(b) is pseudorandomly generated and multiplying a field
 * element by sqrt(b) is as expensive as generic field multiplication),
 * it is beneficial to switch to some alternate formulas for point
 * doubling:
 *
 *   Single doubling:          3M + 5S + 2mb
 *   Sequence of n doublings:  (3*n)M + (4*n)S + (n+3)*mb
 *
 * If w = y/x = sqrt(s/x), then the mapping from the non-neutral
 * elements to their w coordinate is injective. We thus encode a point
 * by computing its w coordinate, and encoding it. The point N does not
 * have a defined w coordinate, so we (formally) set w(N) to any
 * "free" value. Formally:
 *
 *   Let k be a conventional fixed value. If a == 0 or if trace(b/a^2) == 1,
 *   then we set k = 0. Otherwise, we set k to a value such that
 *   k^2 + k + a != 0 and trace(b/(k^2 + k + a)) == 1.
 *
 *   To encode the group element P = (x,s):
 *      if P == N, then the encoding is that of the field element 0
 *      otherwise, we use the field element sqrt(s/x) + k
 *
 *   To decode a group element:
 *      decode the value as a field element u (reject if non-canonical).
 *      if u == 0, then return N
 *      w <- u + k
 *      d <- w^2 + w + a
 *      if d == 0 then reject  (this cannot happen if trace(a) != 0)
 *      e <- b/d^2
 *      if trace(e) == 1 then reject
 *      x <- d*qsolve(e)       (qsolve(e) returns f such that f^2 + f = e)
 *      if x does not match a group element, then:
 *          x <- x + d
 *          if x still does not match a group element, then reject
 *      s <- x*w^2
 *      return (x, s)
 *
 * On fields GF(2^m) with m odd (this is the usual case), qsolve() is
 * the half-trace:
 *     halftrace(e) = \sum_{i=0}^{(m-1)/2} e^(2^(2*i))
 * Matching x with group elements:
 *
 *   If trace(a) == 1, x (as constructed) matches a group element if and
 *   only if trace(x) == 1. Note that, for such a curve, trace(d) == 1,
 *   and therefore if trace(x) == 0 then trace(x + d) == 1, so x + d will
 *   always work if x does not.
 *
 *   If trace(a) == 0, then the whole curve order is equal to r*2^t for
 *   an odd integer r, and an integer t >= 2. In that case, once x is
 *   obtained with x <- d*qsolve(e), the process is:
 *       x' <- x
 *       w' <- w
 *       for i = 0 to t-2:    (t-1 iterations)
 *           if trace(x' + a) != 0 then reject
 *           lambda <- qsolve(x' + a)
 *           x' <- sqrt((lambda + w' + 1)*x' + b)
 *           w' <- sqrt(lambda + a)
 *       if trace(x' + a) == 0:
 *           x <- x + d
 *       s <- x*w^2
 *           
 * The mechanism used here is point halving; a point P is in the group if
 * and only if P = (2^(t-1))*Q for some point Q, and Q cannot be halved.
 * If we get to Q and Q can still be halved, then we replace P with -P+N,
 * which means adding d to x, but leaving w unchanged.
 */

/*
 * Internal structure for points.
 */
typedef struct {
	gf x, s, z, t;
} inner_point;

/*
 * Internal structure for points in affine coordinates. The coordinates
 * are actually (X,S) such that x = sqrt(b)*X and s = sqrt(b)*S, i.e.
 * the equivalent to extended coordinates with Z = 1 and T = X.
 */
typedef struct {
	gf x, s;
} inner_point_affine;

#if GF233_PCLMUL

/*
 * Point structure received from the external API may not be aligned,
 * since the header does not enforce __m128i alignment. We use the
 * load_?() and store_?() functions to the coordinates.
 */

static inline void
point_load_x(gf *x, const void *p)
{
	const uint8_t *pbuf = p;
	x->v[0] = _mm_loadu_si128((const __m128i *)&pbuf[0]);
	x->v[1] = _mm_loadu_si128((const __m128i *)&pbuf[16]);
}

static inline void
point_load_s(gf *s, const void *p)
{
	const uint8_t *pbuf = p;
	s->v[0] = _mm_loadu_si128((const __m128i *)&pbuf[32]);
	s->v[1] = _mm_loadu_si128((const __m128i *)&pbuf[48]);
}

static inline void
point_load_z(gf *z, const void *p)
{
	const uint8_t *pbuf = p;
	z->v[0] = _mm_loadu_si128((const __m128i *)&pbuf[64]);
	z->v[1] = _mm_loadu_si128((const __m128i *)&pbuf[80]);
}

static inline void
point_load_t(gf *t, const void *p)
{
	const uint8_t *pbuf = p;
	t->v[0] = _mm_loadu_si128((const __m128i *)&pbuf[96]);
	t->v[1] = _mm_loadu_si128((const __m128i *)&pbuf[112]);
}

static inline void
point_store_x(void *p, const gf *x)
{
	uint8_t *pbuf = p;
	_mm_storeu_si128((__m128i *)&pbuf[0], x->v[0]);
	_mm_storeu_si128((__m128i *)&pbuf[16], x->v[1]);
}

static inline void
point_store_s(void *p, const gf *s)
{
	uint8_t *pbuf = p;
	_mm_storeu_si128((__m128i *)&pbuf[32], s->v[0]);
	_mm_storeu_si128((__m128i *)&pbuf[48], s->v[1]);
}

static inline void
point_store_z(void *p, const gf *z)
{
	uint8_t *pbuf = p;
	_mm_storeu_si128((__m128i *)&pbuf[64], z->v[0]);
	_mm_storeu_si128((__m128i *)&pbuf[80], z->v[1]);
}

static inline void
point_store_t(void *p, const gf *t)
{
	uint8_t *pbuf = p;
	_mm_storeu_si128((__m128i *)&pbuf[96], t->v[0]);
	_mm_storeu_si128((__m128i *)&pbuf[112], t->v[1]);
}

static inline void
point_load(inner_point *ip3, const void *p1)
{
	point_load_x(&ip3->x, p1);
	point_load_s(&ip3->s, p1);
	point_load_z(&ip3->z, p1);
	point_load_t(&ip3->t, p1);
}

static inline void
point_store(void *p3, const inner_point *ip1)
{
	point_store_x(p3, &ip1->x);
	point_store_s(p3, &ip1->s);
	point_store_z(p3, &ip1->z);
	point_store_t(p3, &ip1->t);
}

static inline void
point_select(void *p3, const void *p0, const void *p1, uint32_t ctl)
{
	const __m128i *px0 = p0;
	const __m128i *px1 = p1;
	__m128i *px3 = p3;
	__m128i m = _mm_set1_epi32(ctl);

	for (int i = 0; i < 8; i ++) {
		__m128i v0 = _mm_loadu_si128(&px0[i]);
		__m128i v1 = _mm_loadu_si128(&px1[i]);
		_mm_storeu_si128(&px3[i], _mm_blendv_epi8(v0, v1, m));
	}
}

static inline void
inner_select(inner_point *p3,
	const inner_point *p0, const inner_point *p1, uint32_t ctl)
{
	__m128i m = _mm_set1_epi32(ctl);
	p3->x.v[0] = _mm_blendv_epi8(p0->x.v[0], p1->x.v[0], m);
	p3->x.v[1] = _mm_blendv_epi8(p0->x.v[1], p1->x.v[1], m);
	p3->s.v[0] = _mm_blendv_epi8(p0->s.v[0], p1->s.v[0], m);
	p3->s.v[1] = _mm_blendv_epi8(p0->s.v[1], p1->s.v[1], m);
	p3->z.v[0] = _mm_blendv_epi8(p0->z.v[0], p1->z.v[0], m);
	p3->z.v[1] = _mm_blendv_epi8(p0->z.v[1], p1->z.v[1], m);
	p3->t.v[0] = _mm_blendv_epi8(p0->t.v[0], p1->t.v[0], m);
	p3->t.v[1] = _mm_blendv_epi8(p0->t.v[1], p1->t.v[1], m);
}

/*
 * Lookup a point from a window. The digit d has a value between -num
 * and +num (inclusive). If d > 0, this returns win[d - 1]; if d < 0,
 * the point -win[-d - 1] is returned. If d == 0, then the neutral is
 * returned (sqrt(b) is provided as parameter for this purpose).
 * The returned point is written into *q.
 */
static inline void
inner_lookup(inner_point *q, const inner_point *win, int num,
	int8_t d, const gf *sqrt_b)
{
	inner_point r;

	/*
	 * k <- abs(d)
	 * sk <- -1 if d < 0, or 0 otherwise
	 */
	uint32_t k = (uint32_t)(uint8_t)d;
	uint64_t sk = -(uint64_t)(k >> 7);
	k = ((k ^ (uint32_t)sk) - (uint32_t)sk) & 0xFF;

	/*
	 * Get the point corresponding to k. If there is none (i.e. d == 0)
	 * then the values in r remain at zero.
	 */
	__m128i xk = _mm_set1_epi32(k);
	__m128i xi = _mm_setzero_si128();
	memset(&r, 0, sizeof r);
	for (int i = 0; i < num; i ++) {
		xi = _mm_add_epi32(xi, _mm_set1_epi32(1));
		__m128i xm = _mm_cmpeq_epi32(xi, xk);
		r.x.v[0] = _mm_or_si128(r.x.v[0],
			_mm_and_si128(win[i].x.v[0], xm));
		r.x.v[1] = _mm_or_si128(r.x.v[1],
			_mm_and_si128(win[i].x.v[1], xm));
		r.s.v[0] = _mm_or_si128(r.s.v[0],
			_mm_and_si128(win[i].s.v[0], xm));
		r.s.v[1] = _mm_or_si128(r.s.v[1],
			_mm_and_si128(win[i].s.v[1], xm));
		r.z.v[0] = _mm_or_si128(r.z.v[0],
			_mm_and_si128(win[i].z.v[0], xm));
		r.z.v[1] = _mm_or_si128(r.z.v[1],
			_mm_and_si128(win[i].z.v[1], xm));
		r.t.v[0] = _mm_or_si128(r.t.v[0],
			_mm_and_si128(win[i].t.v[0], xm));
		r.t.v[1] = _mm_or_si128(r.t.v[1],
			_mm_and_si128(win[i].t.v[1], xm));
	}

	/*
	 * If d == 0, set the point to the neutral; X and T are already
	 * correct (they are zero), but we should set S to sqrt(b) and Z to 1.
	 */
	uint64_t kz = (uint64_t)((k | -k) >> 31) - 1;
	gf_condset(&r.s, sqrt_b, kz);
	gf_condset(&r.z, &GF_ONE, kz);

	/*
	 * If d < 0 (i.e. sk = 1), negate the point.
	 */
	gf_condadd(&r.s, &r.t, sk);

	*q = r;
}

/*
 * Lookup an affine point from a window. The digit d has a value between -num
 * and +num (inclusive). If d > 0, this returns win[d - 1]; if d < 0,
 * the point -win[-d - 1] is returned. If d == 0, then the neutral is
 * returned (sqrt(b) is provided as parameter for this purpose).
 * The returned point is written into *q.
 */
static void
inner_lookup_affine(inner_point_affine *q,
	const inner_point_affine *win, int num, int8_t d, const gf *sqrt_b)
{
	inner_point_affine r;

	/*
	 * k <- abs(d)
	 * sk <- -1 if d < 0, or 0 otherwise
	 */
	uint32_t k = (uint32_t)(uint8_t)d;
	uint64_t sk = -(uint64_t)(k >> 7);
	k = ((k ^ (uint32_t)sk) - (uint32_t)sk) & 0xFF;

	/*
	 * Get the point corresponding to k. If there is none (i.e. d == 0)
	 * then the values in r remain at zero.
	 */
	__m128i xk = _mm_set1_epi32(k);
	__m128i xi = _mm_setzero_si128();
	r.x.v[0] = _mm_setzero_si128();
	r.x.v[1] = _mm_setzero_si128();
	r.s.v[0] = _mm_setzero_si128();
	r.s.v[1] = _mm_setzero_si128();
	for (int i = 0; i < num; i ++) {
		xi = _mm_add_epi32(xi, _mm_set1_epi32(1));
		__m128i xm = _mm_cmpeq_epi32(xi, xk);
		r.x.v[0] = _mm_or_si128(r.x.v[0],
			_mm_and_si128(win[i].x.v[0], xm));
		r.x.v[1] = _mm_or_si128(r.x.v[1],
			_mm_and_si128(win[i].x.v[1], xm));
		r.s.v[0] = _mm_or_si128(r.s.v[0],
			_mm_and_si128(win[i].s.v[0], xm));
		r.s.v[1] = _mm_or_si128(r.s.v[1],
			_mm_and_si128(win[i].s.v[1], xm));
	}

	/*
	 * If d == 0, set the point to the neutral; X is then correct (zero)
	 * but S should be corrected.
	 */
	uint64_t kz = (uint64_t)((k | -k) >> 31) - 1;
	gf_condset(&r.s, sqrt_b, kz);

	/*
	 * If d < 0 (i.e. sk = 1), negate the point.
	 * Negation is adding T to S; here, T is implicit since it is
	 * equal to X.
	 */
	gf_condadd(&r.s, &r.x, sk);

	*q = r;
}

#else /* GF233_PCLMUL */

static inline void
point_load_x(gf *x, const void *p)
{
	const gf *pv = p;
	*x = pv[0];
}

static inline void
point_load_s(gf *s, const void *p)
{
	const gf *pv = p;
	*s = pv[1];
}

static inline void
point_load_z(gf *z, const void *p)
{
	const gf *pv = p;
	*z = pv[2];
}

static inline void
point_load_t(gf *t, const void *p)
{
	const gf *pv = p;
	*t = pv[3];
}

static inline void
point_store_x(void *p, const gf *x)
{
	gf *pv = p;
	pv[0] = *x;
}

static inline void
point_store_s(void *p, const gf *s)
{
	gf *pv = p;
	pv[1] = *s;
}

static inline void
point_store_z(void *p, const gf *z)
{
	gf *pv = p;
	pv[2] = *z;
}

static inline void
point_store_t(void *p, const gf *t)
{
	gf *pv = p;
	pv[3] = *t;
}

static inline void
point_load(inner_point *ip3, const void *p1)
{
	memcpy(ip3, p1, sizeof(inner_point));
}

static inline void
point_store(void *p3, const inner_point *ip1)
{
	memcpy(p3, ip1, sizeof(inner_point));
}

static inline void
point_select(void *p3, const void *p0, const void *p1, uint32_t ctl)
{
	uint64_t m;

	m = (uint64_t)ctl | ((uint64_t)ctl << 32);
	for (int i = 0; i < 16; i ++) {
		uint64_t w0 = ((const uint64_t *)p0)[i];
		uint64_t w1 = ((const uint64_t *)p1)[i];
		((uint64_t *)p3)[i] = w0 ^ (m & (w0 ^ w1));
	}
}

static inline void
inner_select(inner_point *p3,
	const inner_point *p0, const inner_point *p1, uint32_t ctl)
{
	uint64_t m;

	m = (uint64_t)ctl | ((uint64_t)ctl << 32);
	p3->x.w[0] = p0->x.w[0] ^ (m & (p0->x.w[0] ^ p1->x.w[0]));
	p3->x.w[1] = p0->x.w[1] ^ (m & (p0->x.w[1] ^ p1->x.w[1]));
	p3->x.w[2] = p0->x.w[2] ^ (m & (p0->x.w[2] ^ p1->x.w[2]));
	p3->x.w[3] = p0->x.w[3] ^ (m & (p0->x.w[3] ^ p1->x.w[3]));
	p3->s.w[0] = p0->s.w[0] ^ (m & (p0->s.w[0] ^ p1->s.w[0]));
	p3->s.w[1] = p0->s.w[1] ^ (m & (p0->s.w[1] ^ p1->s.w[1]));
	p3->s.w[2] = p0->s.w[2] ^ (m & (p0->s.w[2] ^ p1->s.w[2]));
	p3->s.w[3] = p0->s.w[3] ^ (m & (p0->s.w[3] ^ p1->s.w[3]));
	p3->z.w[0] = p0->z.w[0] ^ (m & (p0->z.w[0] ^ p1->z.w[0]));
	p3->z.w[1] = p0->z.w[1] ^ (m & (p0->z.w[1] ^ p1->z.w[1]));
	p3->z.w[2] = p0->z.w[2] ^ (m & (p0->z.w[2] ^ p1->z.w[2]));
	p3->z.w[3] = p0->z.w[3] ^ (m & (p0->z.w[3] ^ p1->z.w[3]));
	p3->t.w[0] = p0->t.w[0] ^ (m & (p0->t.w[0] ^ p1->t.w[0]));
	p3->t.w[1] = p0->t.w[1] ^ (m & (p0->t.w[1] ^ p1->t.w[1]));
	p3->t.w[2] = p0->t.w[2] ^ (m & (p0->t.w[2] ^ p1->t.w[2]));
	p3->t.w[3] = p0->t.w[3] ^ (m & (p0->t.w[3] ^ p1->t.w[3]));
}

/*
 * Lookup a point from a window. The digit d has a value between -num
 * and +num (inclusive). If d > 0, this returns win[d - 1]; if d < 0,
 * the point -win[-d - 1] is returned. If d == 0, then the neutral is
 * returned (sqrt(b) is provided as parameter for this purpose).
 * The returned point is written into *q.
 */
static inline void
inner_lookup(inner_point *q, const inner_point *win, int num,
	int8_t d, const gf *sqrt_b)
{
	inner_point r;

	/*
	 * k <- abs(d)
	 * sk <- -1 if d < 0, or 0 otherwise
	 */
	uint64_t k = (uint64_t)(uint8_t)d;
	uint64_t sk = -(k >> 7);
	k = ((k ^ sk) - sk) & 0xFF;

	/*
	 * Get the point corresponding to k. If there is none (i.e. d == 0)
	 * then the values in r remain at zero.
	 */
	memset(&r, 0, sizeof r);
	for (int i = 0; i < num; i ++) {
		uint64_t m;

		m = (-(k ^ (uint64_t)(i + 1)) >> 63) - 1;
		r.x.w[0] |= m & win[i].x.w[0];
		r.x.w[1] |= m & win[i].x.w[1];
		r.x.w[2] |= m & win[i].x.w[2];
		r.x.w[3] |= m & win[i].x.w[3];
		r.s.w[0] |= m & win[i].s.w[0];
		r.s.w[1] |= m & win[i].s.w[1];
		r.s.w[2] |= m & win[i].s.w[2];
		r.s.w[3] |= m & win[i].s.w[3];
		r.z.w[0] |= m & win[i].z.w[0];
		r.z.w[1] |= m & win[i].z.w[1];
		r.z.w[2] |= m & win[i].z.w[2];
		r.z.w[3] |= m & win[i].z.w[3];
		r.t.w[0] |= m & win[i].t.w[0];
		r.t.w[1] |= m & win[i].t.w[1];
		r.t.w[2] |= m & win[i].t.w[2];
		r.t.w[3] |= m & win[i].t.w[3];
	}

	/*
	 * If d == 0, set the point to the neutral; X and T are already
	 * correct (they are zero), but we should set S to sqrt(b) and Z to 1.
	 */
	uint64_t kz = ((k | -k) >> 63) - 1;
	gf_condset(&r.s, sqrt_b, kz);
	gf_condset(&r.z, &GF_ONE, kz);

	/*
	 * If d < 0 (i.e. sk = 1), negate the point.
	 */
	gf_condadd(&r.s, &r.t, sk);

	*q = r;
}

/*
 * Lookup an affine point from a window. The digit d has a value between -num
 * and +num (inclusive). If d > 0, this returns win[d - 1]; if d < 0,
 * the point -win[-d - 1] is returned. If d == 0, then the neutral is
 * returned (sqrt(b) is provided as parameter for this purpose).
 * The returned point is written into *q.
 */
static void
inner_lookup_affine(inner_point_affine *q,
	const inner_point_affine *win, int num, int8_t d, const gf *sqrt_b)
{
	inner_point_affine r;

	/*
	 * k <- abs(d)
	 * sk <- -1 if d < 0, or 0 otherwise
	 */
	uint64_t k = (uint64_t)(uint8_t)d;
	uint64_t sk = -(k >> 7);
	k = ((k ^ sk) - sk) & 0xFF;

	/*
	 * Get the point corresponding to k. If there is none (i.e. d == 0)
	 * then the values in r remain at zero.
	 */
	memset(&r, 0, sizeof r);
	for (int i = 0; i < num; i ++) {
		uint64_t m;

		m = (-(k ^ (uint64_t)(i + 1)) >> 63) - 1;
		r.x.w[0] |= m & win[i].x.w[0];
		r.x.w[1] |= m & win[i].x.w[1];
		r.x.w[2] |= m & win[i].x.w[2];
		r.x.w[3] |= m & win[i].x.w[3];
		r.s.w[0] |= m & win[i].s.w[0];
		r.s.w[1] |= m & win[i].s.w[1];
		r.s.w[2] |= m & win[i].s.w[2];
		r.s.w[3] |= m & win[i].s.w[3];
	}

	/*
	 * If d == 0, set the point to the neutral; X is then correct (zero)
	 * but S should be corrected.
	 */
	uint64_t kz = ((k | -k) >> 63) - 1;
	gf_condset(&r.s, sqrt_b, kz);

	/*
	 * If d < 0 (i.e. sk = 1), negate the point.
	 * Negation is adding T to S; here, T is implicit since it is
	 * equal to X.
	 */
	gf_condadd(&r.s, &r.x, sk);

	*q = r;
}

#endif

/*
 * Recode 40 bits from the provided integer into 8 signed 5-bit digits.
 * For the input x, this stores sd[0], sd[1],... sd[7] and returns y
 * such that:
 *   x = y*2^40 + \sum_{i=0}^{7} sd[i]*2^(5*i)
 *   -15 <= sd[i] <= +16
 * This implies that:
 *   y <= floor((x + 532021755375)/2^40)
 */
static inline uint64_t
u40_recode5(int8_t *sd, uint64_t x)
{
	for (int i = 0; i < 8; i ++) {
		uint64_t v = (x & 0x1F);
		uint64_t cc = (16 - v) >> 63;
		x = (x >> 5) + cc;
		sd[i] = (int8_t)v - ((int8_t)cc << 5);
	}
	return x;
}

/*
 * Recode a 30-byte scalar into 48 signed 5-bit digits. The final carry
 * is returned; it may be equal to 0 or 1 (this is the extra top digit).
 */
static uint64_t
u240_recode5(int8_t *sd, const void *n)
{
	const uint8_t *buf = n;
	uint64_t x = 0;
	for (int i = 0; i < 6; i ++) {
		for (int j = 0; j < 5; j ++) {
			x += (uint64_t)(*buf ++) << (8 * j);
		}
		x = u40_recode5(sd + 8 * i, x);
	}
	return x;
}
