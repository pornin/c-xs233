#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include "xs233.h"

#include "xs233_common.h"

/* see xs233.h */
const xsk233_point xsk233_neutral = { {
	GFext(0, 0, 0, 0),
	GFext(0, 0, 0, 1),
	GFext(0, 0, 0, 1),
	GFext(0, 0, 0, 0)
} };

/* see xs233.h */
const xsk233_point xsk233_generator = { {
	/* G is obtained from mapping the K-233 generator with the
	   y -> y + b change of variable. */
	GFext(0x000001ECB92776D0, 0xFB3DEC476585B906,
	      0x5724EF7E1966BF54, 0xA850E5CBDDAA1BE6),
	GFext(0x000000EDFF3B4D4E, 0x5BAA47FCFDF3669D,
	      0xF7193250076F96C1, 0x66F9E0BF367D9A99),
	GFext(0, 0, 0, 1),
	GFext(0x000001ECB92776D0, 0xFB3DEC476585B906,
	      0x5724EF7E1966BF54, 0xA850E5CBDDAA1BE6)
} };

/*
 * (2^180)*G
 */
static const inner_point XSK233_T180_G = {
	GFw64be(0x0000012A4D27EBCE, 0x40EDDD44BB1C4E85,
	        0x60421AA6BFC092E6, 0x50F67F36F364F4CF),
	GFw64be(0x000001BC2DB6E778, 0xC34A5847B86898E3,
	        0x5EE97661043AA776, 0xA419B83462187C5B),
	GFw64be(0, 0, 0, 1),
	GFw64be(0x0000012A4D27EBCE, 0x40EDDD44BB1C4E85,
	        0x60421AA6BFC092E6, 0x50F67F36F364F4CF)
};

static const inner_point xsk233_inner_neutral = {
	GFw64be(0, 0, 0, 0),
	GFw64be(0, 0, 0, 1),
	GFw64be(0, 0, 0, 1),
	GFw64be(0, 0, 0, 0)
};

static const inner_point_affine XSK233_G0[];
static const inner_point_affine XSK233_G60[];
static const inner_point_affine XSK233_G120[];
static const inner_point_affine XSK233_G180[];

static const inner_point_affine XSK233_TAU_G0[];
static const inner_point_affine XSK233_TAU_G60[];
static const inner_point_affine XSK233_TAU_G120[];
static const inner_point_affine XSK233_TAU_G180[];

/* see xs233.h */
uint32_t
xsk233_decode(xsk233_point *p, const void *src)
{
	/*
	 * Procedure:
	 *   decode w from bytes; reject if invalid (non-canonical).
	 *   if w == 0 then return N
	 *   d <- w^2 + w
	 *   if d == 0 then reject
	 *   e <- 1/(d^2)
	 *   if trace(e) == 1 then reject
	 *   f <- halftrace(e)       (we have f^2 + f = e)
	 *   x <- d*f
	 *   if trace(x) == 1 then reject (point cannot be halved)
	 *   lambda = halftrace(x)
	 *   if trace(x*(lambda + w) + b) == 0:
	 *       x <- x + d
	 *   s <- x*w^2
	 *
	 * Note: trace(x*(lambda + w) + b) == 1 + trace(x*(lambda + w))
	 * (since b = 1 for xsk233).
	 */
	gf w, d, e, f, g, x, s;
	uint64_t ve, rp, wz, r;

	/* decode w; remember if it was zero */
	ve = gf_decode(&w, src);
	wz = gf_iszero(&w);

	/* d <- w^2 + w; must not be zero */
	gf_sqr(&d, &w);
	gf_add(&d, &d, &w);
	rp = ~gf_iszero(&d);

	/* e <- b/(d^2); must have trace 0 */
	gf_sqr(&e, &d);
	gf_inv(&e, &e);
	rp &= gf_trace(&e) - 1;

	/* solve f^2 + f = e */
	gf_halftrace(&f, &e);

	/* x <- d*f; must have trace 0 */
	gf_mul(&x, &d, &f);
	rp &= gf_trace(&x) - 1;

	/* solve lambda^2 + lambda = x */
	gf_halftrace(&g, &x);

	/* If trace(x*(lambda + 1 + w) + b) == 0, then x <- x + d */
	gf_add(&g, &g, &w);
	gf_mul(&g, &g, &x);
	gf_condadd(&x, &d, -gf_trace(&g));

	/* s <- x*w^2 */
	gf_sqr(&s, &w);
	gf_mul(&s, &s, &x);

	/*
	 * ve = -1 if the field element was canonical
	 * wz = -1 if w == 0
	 * rp = -1 if trace(e) != 0 or trace(x) != 0
	 * If ve == -1 then wz == -1.
	 * Decoding succeeds if ve == -1 and either wz == -1 or rp == -1.
	 * We force the point to N on failure, or if wz == -1.
	 */
	r = ve & (rp | wz);
	gf_condset(&x, &GF_ZERO, wz | ~rp);
	gf_condset(&s, &GF_ONE, wz | ~rp);

	/*
	 * In extended coordinates:
	 *   x = sqrt(b)*X/Z
	 *   s = sqrt(b)*S/Z^2
	 *   T = X*Z
	 * We set Z = 1, hence X = x/sqrt(b), S = s/sqrt(b), T = X
	 * On xsk233, sqrt(b) == 1.
	 */
	point_store_x(p, &x);
	point_store_s(p, &s);
	point_store_z(p, &GF_ONE);
	point_store_t(p, &x);
	return (uint32_t)r;
}

/* see xs233.h */
void
xsk233_encode(void *dst, const xsk233_point *p)
{
	/*
	 * w = y/x = sqrt(s/x)
	 * Since:
	 *   x = sqrt(b)*X/Z
	 *   s = sqrt(b)*S/Z^2
	 *   T = X*Z
	 * we have:
	 *   w = sqrt((sqrt(b)*S*Z)/(sqrt(b)*X*Z^2)
	 *     = sqrt(S/T)
	 * For the point N, x == 0, hence T == 0, and the gf_div() call
	 * will yield 0, which is what we want in that case.
	 */
	gf s, t, w;

	point_load_s(&s, p);
	point_load_t(&t, p);
	gf_div(&w, &s, &t);
	gf_sqrt(&w, &w);
	gf_encode(dst, &w);
}

/* see xs233.h */
uint32_t
xsk233_is_neutral(const xsk233_point *p)
{
	/*
	 * The neutral point N is the only point for which x == 0.
	 */
	gf x;

	point_load_x(&x, p);
	return (uint32_t)gf_iszero(&x);
}

/* see xs233.h */
uint32_t
xsk233_equals(const xsk233_point *p1, const xsk233_point *p2)
{
	/*
	 * Mapping (x,s) -> w = sqrt(s/x) is injective for the group.
	 * We can compare the w^2 values (squaring is also injective),
	 * with w^2 = (sqrt(b)*S/Z^2)/(sqrt(b)*X/Z)
	 *          = S/X*Z
	 *          = S/T
	 * Thus, we only need to check whether S1*T2 == S2*T1. There is
	 * no special case to handle because:
	 *    S and T can never be zero simultaneously for a point.
	 *    There is a single group element with T == 0 (the neutral N).
	 *    There is (at most) a single group element with S == 0
	 */
	gf s1, t1, s2, t2;

	point_load_s(&s1, p1);
	point_load_t(&t1, p1);
	point_load_s(&s2, p2);
	point_load_t(&t2, p2);
	gf_mul(&s1, &s1, &t2);
	gf_mul(&s2, &s2, &t1);
	return (uint32_t)gf_equals(&s1, &s2);
}

/*
 * Inner addition function on xsk233.
 */
static inline void
xsk233_inner_add(inner_point *p3, const inner_point *p1, const inner_point *p2)
{
	/*
	 * x1x2 <- X1*X2
	 * s1s2 <- S1*S2
	 * z1z2 <- Z1*Z2
	 * (suppressed: t1t2 <- T1*T2)
	 * d <- (S1 + T1)*(S2 + T2)
	 * (suppressed: e <- (a^2)*t1t2)
	 * f <- x1x2^2
	 * g <- z1z2^2
	 * X3 <- d + s1s2
	 * S3 <- sqrt(b)*(g*s1s2 + f*d)     note: sqrt(b) = 1 for xsk233
	 * Z3 <- sqrt(b)*(f + g)
	 * T3 <- X3*Z3
	 */
	gf x1x2, s1s2, z1z2;
	gf d, f, g, tmp1, tmp2;

	gf_mul(&x1x2, &p1->x, &p2->x);
	gf_mul(&s1s2, &p1->s, &p2->s);
	gf_mul(&z1z2, &p1->z, &p2->z);
	gf_add(&tmp1, &p1->s, &p1->t);
	gf_add(&tmp2, &p2->s, &p2->t);
	gf_mul(&d, &tmp1, &tmp2);
	gf_sqr(&f, &x1x2);
	gf_sqr(&g, &z1z2);
	gf_add(&p3->x, &d, &s1s2);
	gf_mul(&tmp1, &s1s2, &g);
	gf_mul(&tmp2, &d, &f);
	gf_add(&p3->s, &tmp1, &tmp2);
	gf_add(&p3->z, &f, &g);
	gf_mul(&p3->t, &p3->x, &p3->z);
}

/*
 * Inner subtraction function on xsk233.
 */
static inline void
xsk233_inner_sub(inner_point *p3, const inner_point *p1, const inner_point *p2)
{
	inner_point p4;

	p4.x = p2->x;
	gf_add(&p4.s, &p2->s, &p2->t);
	p4.z = p2->z;
	p4.t = p2->t;
	xsk233_inner_add(p3, p1, &p4);
}

/* see xs233.h */
void
xsk233_add(xsk233_point *p3, const xsk233_point *p1, const xsk233_point *p2)
{
	inner_point ip1, ip2, ip3;

	point_load(&ip1, p1);
	point_load(&ip2, p2);
	xsk233_inner_add(&ip3, &ip1, &ip2);
	point_store(p3, &ip3);
}

/*
 * Mixed addition: second operand is in affine coordinates (Z2 == 1, T2 == X2).
 */
static inline void
xsk233_add_mixed(inner_point *p3, const inner_point *p1,
	const inner_point_affine *p2)
{
	/*
	 * When P2 is in affine coordinates, Z2 == 1 and T2 == X2, so
	 * the values Z2 and T2 are implicit, and the multiplication
	 * Z1*Z2 simplifies to using just Z1.
	 */
	gf x1x2, s1s2;
	gf d, f, g, tmp1, tmp2;

	gf_mul(&x1x2, &p1->x, &p2->x);
	gf_mul(&s1s2, &p1->s, &p2->s);
	gf_add(&tmp1, &p1->s, &p1->t);
	gf_add(&tmp2, &p2->s, &p2->x);
	gf_mul(&d, &tmp1, &tmp2);
	gf_sqr(&f, &x1x2);
	gf_sqr(&g, &p1->z);
	gf_add(&p3->x, &d, &s1s2);
	gf_mul(&tmp1, &s1s2, &g);
	gf_mul(&tmp2, &d, &f);
	gf_add(&p3->s, &tmp1, &tmp2);
	gf_add(&p3->z, &f, &g);
	gf_mul(&p3->t, &p3->x, &p3->z);
}

/*
 * Inner doubling function on xsk233.
 */
static inline void
xsk233_inner_double(inner_point *p3, const inner_point *p1)
{
	/*
	 *   zz <- Z^2
	 *   tt <- T^2
	 *   D  <- (S + a*T)*(S + (a + 1)*T)        note: a == 0 on xsk233
	 *   E  <- (S + (a + 1)*T + sqrt(b)*zz)^2   note: b == 1 on xsk233
	 *   X' <- sqrt(b)*tt
	 *   Z' <- D + a*tt
	 *   S' <- sqrt(b)*(E*(E + Z' + tt) + (D + sqrt(b)*X')^2)
	 *   T' <- X'*Z'
	 *
	 * Since a == 0 and b == 1 on xsk233, we get some simplifications:
	 *   tt == X'
	 *   D == S*(S + T)
	 *   D == Z'
	 *   Z' + tt == D + b*tt == Z' + X'
	 */
	gf zz, e, tmp1, tmp2;

	gf_sqr(&zz, &p1->z);
	gf_sqr(&p3->x, &p1->t);
	gf_add(&tmp1, &p1->s, &p1->t);
	gf_mul(&p3->z, &p1->s, &tmp1);
	gf_add(&tmp1, &tmp1, &zz);
	gf_sqr(&e, &tmp1);
	gf_add(&tmp2, &p3->z, &p3->x);
	gf_add(&tmp1, &tmp2, &e);
	gf_mul(&tmp1, &e, &tmp1);
	gf_sqr(&tmp2, &tmp2);
	gf_add(&p3->s, &tmp1, &tmp2);
	gf_mul(&p3->t, &p3->x, &p3->z);
}

/* see xs233.h */
void
xsk233_double(xsk233_point *p3, const xsk233_point *p1)
{
	inner_point ip1, ip3;

	point_load(&ip1, p1);
	xsk233_inner_double(&ip3, &ip1);
	point_store(p3, &ip3);
}

/* see xs233.h */
void
xsk233_xdouble(xsk233_point *p3, const xsk233_point *p1, unsigned n)
{
	inner_point ip;

	if (n == 0) {
		*p3 = *p1;
		return;
	}

	point_load(&ip, p1);
	while (n -- > 0) {
		xsk233_inner_double(&ip, &ip);
	}
	point_store(p3, &ip);
}

/* see xs233.h */
void
xsk233_neg(xsk233_point *p3, const xsk233_point *p1)
{
	/*
	 * Negation of (x,s) is (x,s+x).
	 * We have:
	 *    sqrt(b)*S/Z^2 + sqrt(b)*X/Z = sqrt(b)*(S + X*Z)/Z^2
	 * Since T = X*Z, we only have to add T to S.
	 */
	inner_point ip;

	point_load(&ip, p1);
	gf_add(&ip.s, &ip.s, &ip.t);
	point_store(p3, &ip);
}

/* see xs233.h */
void
xsk233_sub(xsk233_point *p3, const xsk233_point *p1, const xsk233_point *p2)
{
	inner_point ip1, ip2;

	point_load(&ip1, p1);
	point_load(&ip2, p2);
	gf_add(&ip2.s, &ip2.s, &ip2.t);
	xsk233_inner_add(&ip1, &ip1, &ip2);
	point_store(p3, &ip1);
}

/* see xs233.h */
void
xsk233_select(xsk233_point *p3,
	const xsk233_point *p0, const xsk233_point *p1, uint32_t ctl)
{
	point_select(p3, p0, p1, ctl);
}

/* see xs233.h */
void
xsk233_condneg(xsk233_point *p3, const xsk233_point *p1, uint32_t ctl)
{
	inner_point ip;

	point_load(&ip, p1);
	gf_condadd(&ip.s, &ip.t, -(uint64_t)(ctl >> 31));
	point_store(p3, &ip);
}

/* see xs233.h */
void
xsk233_mul(xsk233_point *p3,
	const xsk233_point *p1, const void *n, size_t n_len)
{
	/*
	 * Compute the window: win[i] = (i+1)*P1.
	 */
	inner_point win[16];

	point_load(&win[0], p1);
	xsk233_inner_double(&win[1], &win[0]);
	for (int j = 4; j <= 16; j += 2) {
		xsk233_inner_add(&win[j - 2], &win[j - 3], &win[0]);
		xsk233_inner_double(&win[j - 1], &win[(j >> 1) - 1]);
	}

	/*
	 * Recode the scalar into signed digits. This sets the
	 * following:
	 *
	 *   sd[] contains the digits (at most 48).
	 *
	 *   sd_len is the number of relevant digits. It is computed
	 *   from n_len and does not depend on the scalar byte values
	 *   (i.e. sd_len is considered public). Digits sd[i] for
	 *   i >= sd_len are guaranteed to be zero.
	 *
	 *   *p3 is set to either the input P1, or to the neutral, to
	 *   account for the top digit value (sd[sd_len], which was not
	 *   set).
	 */
	int8_t sd[48];
	size_t sd_len;
	inner_point ip3;

	if (n_len == 30) {
		/*
		 * Normal case when working with full-size scalars.
		 */
		uint32_t td = (uint32_t)u240_recode5(sd, n);
		inner_select(&ip3, &xsk233_inner_neutral, &win[0], -td);
		sd_len = 48;
	} else {
		/*
		 * n_len <= 29.
		 * We get one digit per 5 bits of input.
		 */
		uint8_t tmp[30];

		memcpy(tmp, n, n_len);
		memset(tmp + n_len, 0, (sizeof tmp) - n_len);
		u240_recode5(sd, tmp);
		sd_len = ((n_len * 8) + 4) / 5;
		if (n_len * 8 == 5 * sd_len) {
			uint32_t td = (uint32_t)sd[sd_len];
			inner_select(&ip3, &xsk233_inner_neutral, &win[0], -td);
		} else {
			ip3 = xsk233_inner_neutral;
		}
	}

	/*
	 * Process all signed digits in high-to-low order.
	 */
	while (sd_len -- > 0) {
		inner_point Q;

		for (int i = 0; i < 5; i ++) {
			xsk233_inner_double(&ip3, &ip3);
		}
		inner_lookup(&Q, win, 16, sd[sd_len], &GF_ONE);
		xsk233_inner_add(&ip3, &ip3, &Q);
	}

	point_store(p3, &ip3);
}

/* see xs233.h */
void
xsk233_mul_ladder(xsk233_point *p3,
	const xsk233_point *p1, const void *n, size_t n_len)
{
	/*
	 * We use a Montgomery ladder. This is similar in principle to
	 * the López-Dahab formulas, except that we work with addition
	 * in the group (i.e. (P1+N) + (P2+N) + N = P3+N).
	 *
	 * Given points P0+N = (x0, s0), P1+N and P2+N such that:
	 *   P0+N = (P2+N) - (P1+N) + N
	 * we compute P3+N = (P1+N) + (P2+N) + N with:
	 *
	 *   x0 + x3 = b*x1*x2 / (x1*x2 + b)^2
	 *
	 * We can also double P1+N into P4+N = 2*(P1+N) + N with:
	 *
	 *   x4 = b*x1^2 / (x1^2 + b)^2
	 *
	 * (and similarly for doubling P2+N)
	 *
	 * If x1 = sqrt(b)*X1/Z1, x2 = sqrt(b)*X2/Z2, x3 = sqrt(b)*X3/Z3
	 * and x4 = sqrt(b)*X4/Z4, then we obtain these formulas:
	 *
	 *   X3 = X1*X2*Z1*Z2 + x0*(X1*X2 + Z1*Z2)^2
	 *   Z3 = sqrt(b)*(X1*X2 + Z1*Z2)^2
	 *
	 *   X4 = (X1*Z1)^2
	 *   Z4 = sqrt(b)*(X1 + Z1)^4
	 *
	 * which can be computed in a total of 5M + 4S + 2mb. The
	 * original López-Dahab formulas have cost 5M + 4S + 1mb, which
	 * is better by one multiplication by the constant sqrt(b);
	 * however, on xsk233, sqrt(b) = 1, so that difference
	 * disappears. An advantage of the formulas above is that they
	 * are complete for all source group elements, so there is no
	 * need for any conversion or special handling of the neutral N.
	 *
	 * From P1+N = k*P0 + N, and P2+N = (k+1)*P0 + N, we can compute
	 * P3 and P4 as above, and set (P1, P2) <- (P4, P3), i.e.
	 * ((2*k)*P0 + N, (2*k+1)*P0 + N). But we can also compute the
	 * double of P2 instead of P1, and set (P1, P2) <- (P3, P4),
	 * i.e. ((2*k+1)*P0 + N, (2*k+2)*P0 + N). That way, we can
	 * multiply by an arbitrary scalar in a bit-by-bit way (the scalar
	 * is used in high-to-low order). We start with P1+N = N and
	 * P2+N = P0+N.
	 *
	 * The final s coordinate of the result can be obtained at the
	 * end, by noticing that since P2+N = P1+N + P0+N + N, we have:
	 *
	 *   x2 = b*(x0*x1 + s0*x1 + s1*x0) / (x0*x1 + b)^2
	 *
	 * yielding the following expressions for the coordinates of the
	 * final result:
	 *
	 *   X' = x0*X1*Z2
	 *   S' = x0*Z2*(X1*Z1*Z2*(x0 + s0) + X2*(x0*X1 + sqrt(b)*Z1)^2)
	 *   Z' = x0*Z1*Z2
	 *   T' = X'*Z'
	 *
	 * These last formulas are _not_ complete, though: they have a
	 * special case when x0 == 0, i.e. when the source point is the
	 * neutral N. In that case, we get zero for all coordinates, and
	 * this is invalid for S' and Z', so we must add a correction for
	 * that case at the end. In all other cases, the formulas yield
	 * the correct result (even is P1+N = N or P2+N = N, e.g. when
	 * the scalar k is equal to 0 or -1 mod r).
	 */

	inner_point ip;
	gf x0, s0, X1, Z1, X2, Z2, tt;
	const uint8_t *nbuf;
	uint64_t iz;

	/*
	 * If the scalar is empty then return the neutral.
	 */
	if (n_len == 0) {
		*p3 = xsk233_neutral;
		return;
	}

	point_load(&ip, p1);

	/*
	 * x0 = sqrt(b)*X/Z
	 * s0 = sqrt(b)*S/Z^2
	 * (b == 1 on xsk233)
	 */
	gf_inv(&tt, &ip.z);
	gf_mul(&x0, &ip.x, &tt);
	gf_sqr(&tt, &tt);
	gf_mul(&s0, &ip.s, &tt);

	/*
	 * Initial: P1+N = N, P2+N = P0.
	 */
	X1 = GF_ZERO;
	Z1 = GF_ONE;
	X2 = x0;
	Z2 = GF_ONE;
	nbuf = n;
	for (int i = (int)(n_len << 3) - 1; i >= 0; i --) {
		uint64_t bit;
		gf X3, Z3, X4, Z4, x1x2, z1z2, xt, zt;

		/* Extract next scalar bit. */
		bit = -(uint64_t)((nbuf[i >> 3] >> (i & 7)) & 1);

		/* P3+N = P1+N + P2+N + N */
		gf_mul(&x1x2, &X1, &X2);
		gf_mul(&z1z2, &Z1, &Z2);
		gf_add(&Z3, &x1x2, &z1z2);
		gf_sqr(&Z3, &Z3);
		gf_mul(&X3, &x1x2, &z1z2);
		gf_mul(&tt, &x0, &Z3);
		gf_add(&X3, &X3, &tt);

		/* P4+N = 2*P1 + N (if bit == 0) or 2*P2 + N (if bit == 1) */
		gf_select(&xt, &X1, &X2, bit);
		gf_select(&zt, &Z1, &Z2, bit);
		gf_add(&Z4, &xt, &zt);
		gf_sqr(&Z4, &Z4);
		gf_mul(&X4, &xt, &zt);
		gf_sqr(&Z4, &Z4);
		gf_sqr(&X4, &X4);

		/* (P1, P2) <- (P4, P3) (if bit == 0) or (P3, P4) */
		gf_select(&X1, &X4, &X3, bit);
		gf_select(&Z1, &Z4, &Z3, bit);
		gf_select(&X2, &X3, &X4, bit);
		gf_select(&Z2, &Z3, &Z4, bit);
	}

	/*
	 * Rebuild all coordinates of the output.
	 */
	gf_mul(&tt, &Z1, &Z2);       /* tt <- Z1*Z2 */
	gf_mul(&ip.z, &x0, &tt);

	gf_mul(&ip.x, &X1, &Z2);
	gf_mul(&ip.x, &x0, &ip.x);

	gf_mul(&tt, &X1, &tt);       /* tt <- X1*Z1*Z2 */
	gf_add(&ip.s, &x0, &s0);
	gf_mul(&ip.s, &tt, &ip.s);
	gf_mul(&tt, &x0, &X1);
	gf_add(&tt, &tt, &Z1);
	gf_sqr(&tt, &tt);
	gf_mul(&tt, &X2, &tt);
	gf_add(&ip.s, &ip.s, &tt);
	gf_mul(&tt, &x0, &Z2);
	gf_mul(&ip.s, &ip.s, &tt);

	gf_mul(&ip.t, &ip.x, &ip.z);

	/*
	 *   X' = x0*X1*Z2
	 *   S' = x0*Z2*(X1*Z1*Z2*(x0 + s0) + X2*(x0*X1 + sqrt(b)*Z1)^2)
	 *   Z' = x0*Z1*Z2
	 *   T' = X'*Z'
	 */

	/*
	 * Corrective action in case the input was the neutral: in that
	 * case, the final conversion put 0 in all coordinates, and we
	 * must adjust S' and Z'.
	 */
	iz = gf_iszero(&x0);
	gf_condset(&ip.s, &GF_ONE, iz);
	gf_condset(&ip.z, &GF_ONE, iz);

	point_store(p3, &ip);
}

/* see xs233.h */
void
xsk233_mulgen(xsk233_point *p3, const void *n, size_t n_len)
{
	/*
	 * Recode the scalar into signed digits. We enforce a
	 * full-width scalar by padding with zeros if necessary.
	 */
	int8_t sd[48];
	inner_point ip3;

	if (n_len == 30) {
		/*
		 * Normal case when working with full-size scalars.
		 */
		uint32_t td = (uint32_t)u240_recode5(sd, n);
		inner_select(&ip3, &xsk233_inner_neutral, &XSK233_T180_G, -td);
	} else {
		uint8_t tmp[30];

		memcpy(tmp, n, n_len);
		memset(tmp + n_len, 0, (sizeof tmp) - n_len);
		u240_recode5(sd, tmp);
		ip3 = xsk233_inner_neutral;
	}

	/*
	 * Process all signed digits in high-to-low order.
	 */
	for (int i = 11; i >= 0; i --) {
		inner_point_affine Qa;

		for (int j = 0; j < 5; j ++) {
			xsk233_inner_double(&ip3, &ip3);
		}
		inner_lookup_affine(&Qa, XSK233_G0, 16, sd[i], &GF_ONE);
		xsk233_add_mixed(&ip3, &ip3, &Qa);
		inner_lookup_affine(&Qa, XSK233_G60, 16, sd[i + 12], &GF_ONE);
		xsk233_add_mixed(&ip3, &ip3, &Qa);
		inner_lookup_affine(&Qa, XSK233_G120, 16, sd[i + 24], &GF_ONE);
		xsk233_add_mixed(&ip3, &ip3, &Qa);
		inner_lookup_affine(&Qa, XSK233_G180, 16, sd[i + 36], &GF_ONE);
		xsk233_add_mixed(&ip3, &ip3, &Qa);
	}
	point_store(p3, &ip3);
}

/*
 * Point Multiplication with the Frobenius Operator
 * ------------------------------------------------
 *
 * We use the algorithm in section D.4 from FIPS 186-4 to obtain a
 * (partially reduced) tau-NAF representation of the scalar. The method
 * corresponds to algorithms 3 and 7 in the classic paper from J.
 * Solinas: "Efficient Arithmetic on Koblitz Curves"
 * https://doi.org/10.1007/978-1-4757-6856-5_6
 * (the paper is an expanded version from a prior publication at Crypto'97)
 *
 * Since we want a constant-time implementation, we "un-NAF" the output
 * by grouping successive NAF digits. Given [u0 u1 u2 u3 u4] (with each
 * digit u_i being equal to 0, +1 or -1), we compute the integer
 * u0 + 2*u1 + 4*u2 + 8*u3 + 16*u4. Since the sequence is a non-adjacent
 * form (of any two consecutive digits, at most one is non-zero), this
 * recoding is injective (no information is lost in the process) and
 * can yield only integers in the -21 to +21 range (all values are possible
 * in that range).
 *
 * The scalar multiplication algorithm will use the digits in high-to-low
 * order. For each digit:
 *
 *   1. Apply the Frobenius operator 5 times on the current point
 *      (i.e. raise the coordinates to the power 32).
 *
 *   2. Add the precompute point:
 *        u0*P + u1*frob_1(P) + u2*frob_2(P) + u3*frob_3(P) + u4*frob_4(P)
 *      with frob_i() being i applications of the Frobenius operator.
 *
 * We only need to precompute a window for digits +1 to +21 (since we can
 * negate the window point to get the equivalent of negateing all u_i).
 */

static inline uint64_t
dec64le(const void *src)
{
	const uint8_t *buf;

	buf = src;
	return (uint64_t)buf[0]
		| ((uint64_t)buf[1] << 8)
		| ((uint64_t)buf[2] << 16)
		| ((uint64_t)buf[3] << 24)
		| ((uint64_t)buf[4] << 32)
		| ((uint64_t)buf[5] << 40)
		| ((uint64_t)buf[6] << 48)
		| ((uint64_t)buf[7] << 56);
}

#if (defined _MSC_VER && defined _M_X64) \
	|| (defined __x86_64__ && (defined __GNUC__ || defined __clang__))

#include <immintrin.h>

#define addc(c, a, b, d)   _addcarry_u64(c, a, b, (unsigned long long *)(d))
#define subb(c, a, b, d)   _subborrow_u64(c, a, b, (unsigned long long *)(d))

#else

static inline unsigned char
addc(unsigned char c, uint64_t a, uint64_t b, uint64_t *d)
{
	unsigned __int128 z;

	z = (unsigned __int128)a + (unsigned __int128)b + c;
	*d = (uint64_t)z;
	return (unsigned char)(z >> 64);
}

static inline unsigned char
subb(unsigned char c, uint64_t a, uint64_t b, uint64_t *d)
{
	unsigned __int128 z;

	z = (unsigned __int128)a - (unsigned __int128)b - c;
	*d = (uint64_t)z;
	return (unsigned char)((z >> 64) & 1);
}

#endif

/*
 * Unsigned 64x64->128 multiplication.
 */
static inline void
mull(uint64_t *lo, uint64_t *hi, uint64_t a, uint64_t b)
{
	unsigned __int128 z;

	z = (unsigned __int128)a * (unsigned __int128)b;
	*lo = (uint64_t)z;
	*hi = (uint64_t)(z >> 64);
}

static inline void
sub128(uint64_t *d, const uint64_t *a, const uint64_t *b)
{
	unsigned char cc;

	cc = subb(0, a[0], b[0], &d[0]);
	(void)subb(cc, a[1], b[1], &d[1]);
}

/*
 * 128x128 multiplication, truncated to 128 bits.
 */
static inline void
mul128_trunc(uint64_t *d, const uint64_t *a, const uint64_t *b)
{
	uint64_t lo, hi;

	mull(&lo, &hi, a[0], b[0]);
	hi += a[0] * b[1] + a[1] * b[0];
	d[0] = lo;
	d[1] = hi;
}

/*
 * d <- a*u
 * a[] is signed, ux is unsigned
 */
static inline void
mul64x128(uint64_t *d, const uint64_t *a, uint64_t ux)
{
	uint64_t d0, d1, d2, lo;
	unsigned char cc;

	mull(&d0, &d1, a[0], ux);
	mull(&lo, &d2, a[1], ux);
	cc = addc(0, lo, d1, &d1);
	d2 += cc;
	d2 -= ux & (*(int64_t *)&a[1] >> 63);
	d[0] = d0;
	d[1] = d1;
	d[2] = d2;
}

static inline void
sub192(uint64_t *d, const uint64_t *a, const uint64_t *b)
{
	unsigned char cc;

	cc = subb(0, a[0], b[0], &d[0]);
	cc = subb(cc, a[1], b[1], &d[1]);
	(void)subb(cc, a[2], b[2], &d[2]);
}

static inline void
neg192(uint64_t *d, const uint64_t *a)
{
	unsigned char cc;

	cc = subb(0, 0, a[0], &d[0]);
	cc = subb(cc, 0, a[1], &d[1]);
	(void)subb(cc, 0, a[2], &d[2]);
}

/*
 * d <- a + b
 * a[] has size 192 bits. bx is signed.
 */
static inline void
add192p64(uint64_t *d, const uint64_t *a, const uint64_t bx)
{
	uint64_t sb;
	unsigned char cc;

	sb = (uint64_t)(*(int64_t *)&bx >> 63);
	cc = addc(0, a[0], bx, &d[0]);
	cc = addc(cc, a[1], sb, &d[1]);
	(void)addc(cc, a[2], sb, &d[2]);
}

/*
 * In the notations of FIPS 186-4 (D.4), we use:
 *   C = 5
 *   K = (m + 5)/2 + C = 124
 */

static int64_t
approx_divr(uint64_t *q, const uint64_t *si, uint64_t np0, uint64_t np1)
{
	/*
	 * g' <- si*n'
	 * Note that si may be negative (high bit set). We compute an
	 * unsigned multiplication, then adjust the result by subtracting
	 * n'*2^128 if si < 0.
	 */
	uint64_t g0, g1, g2, g3, lo, hi, t;
	unsigned char cc;

	mull(&g0, &g1, np0, si[0]);
	mull(&lo, &g2, np0, si[1]);
	cc = addc(0, lo, g1, &g1);
	mull(&lo, &g3, np1, si[1]);
	cc = addc(cc, lo, g2, &g2);
	g3 += cc;
	mull(&lo, &hi, np1, si[0]);
	cc = addc(0, g1, lo, &g1);
	cc = addc(cc, g2, hi, &g2);
	g3 += cc;
	t = (uint64_t)(*(int64_t *)&si[1] >> 63);
	cc = subb(0, g2, np0 & t, &g2);
	(void)subb(cc, g3, np1 & t, &g3);

	/*
	 * h' <- floor(g' / 2**m)
	 * h' is signed but cannot be large.
	 * On xsk233, s0 and s1 are both negative, and the possible range
	 * for h' is -69..0.
	 */
	int64_t h = *(int64_t *)&g3 >> (233 - 192);

	/*
	 * j' <- h'*Vm
	 * On xsk233, Vm = -137381546011108235394987299651366779.
	 * Since -69 <= h' <= 0 on xsk233, we know that j' is
	 * nonnegative, and fits on 124 bits.
	 */
	uint64_t j0, j1;
	mull(&j0, &j1, -h, 0xBBEC6B57C5CEAF7B);
	j1 += (uint64_t)-h * 0x001A756EE456F351;

	/*
	 * l' <- round((g' + j') / 2^(K - C))
	 * We always have K - C = 119. To perform the rounding, we add
	 * 2^118 prior to the shift (this is merged with the addition
	 * of j').
	 * Take care that the shift must be arithmetic (g' + j' may be
	 * negative).
	 *
	 * On xsk233, |l'| fits on 121 bits. We keep the high half of
	 * l' in a signed integer.
	 */
	uint64_t jh, l0;
	int64_t l1;

	jh = (uint64_t)(*(int64_t *)&j1 >> 63);
	cc = addc(0, g0, j0, &g0);
	cc = addc(cc, g1, j1 + 0x0040000000000000, &g1);
	cc = addc(cc, g2, jh, &g2);
	(void)addc(cc, g3, jh, &g3);
	l0 = (g1 >> 55) | (g2 << 9);
	t = (g2 >> 55) | (g3 << 9);
	l1 = *(int64_t *)&t;

	/*
	 * q' <- round(l' / 2^C)
	 * We use C = 5. We add 16 prior to the shift to perform the
	 * rounding; we need to keep the original l0 too.
	 */
	cc = addc(0, l0, 16, &t);
	l1 += cc;
	q[0] = (t >> 5) | ((uint64_t)l1 << 59);
	q[1] = (uint64_t)(l1 >> 5);

	/*
	 * We return l' - q*2^C; this is a small integer since q is
	 * the rounded version of l'/2^C.
	 */
	t = l0 - (q[0] << 5);
	return *(int64_t *)&t;
}

/*
 * Recode a scalar in PRTNAF, then signed digits for Frobenius-based
 * multiplications.
 *
 * Input scalar n MUST have size 30 bytes (240 bits); it is interpreted
 * with unsigned little-endian convention. Exactly 48 digits are produced
 * and written in sd[]. Each digit is in the -21 to +21 range.
 *
 * This function is non-static, for test and benchmark purposes.
 */
void
xsk233_tau_recode(int8_t *sd, const void *n)
{
	/*
	 * The s0 and s1 constants for xsk233, as specified in FIPS 186-4.
	 * They are here in two's complement representation and
	 * little-endian order (the two values are negative integers).
	 */
	static const uint64_t TAU_S0[] = {
		0xC388ACB7EF3EFC55, 0xFFFAA26900502B63
	};
	static const uint64_t TAU_S1[] = {
		0xE955EBC334C9411A, 0xFFF77D28D2851C91
	};
	static const uint64_t TAU_S0_MINUS_S1[] = {
		0xDA32C0F4BA75BB3B, 0x000325402DCB0ED1
	};
	static const uint64_t TAU_S1_TIMES_2[] = {
		0xD2ABD78669928234, 0xFFEEFA51A50A3923
	};

	/*
	 * Decode n as a big integer, and reduce it modulo r.
	 *
	 * Since r = 2^231 + r0, it suffices to extract the top 9 bits
	 * (t) and subtract t*r0 from n % 2^231. In curve K-233, we have:
	 *   r0 = 0x00069D5BB915BCD46EFB1AD5F173ABDF
	 * If the subtraction yields a borrow, then we have to add back r
	 * (once will be enough).
	 */
	const uint8_t *nbuf;
	uint64_t n0, n1, n2, n3, t, lo, hi;
	unsigned char cc;

#define R0_HI   0x00069D5BB915BCD4
#define R0_LO   0x6EFB1AD5F173ABDF

	nbuf = n;
	n0 = dec64le(nbuf);
	n1 = dec64le(nbuf + 8);
	n2 = dec64le(nbuf + 16);
	n3 = dec64le(nbuf + 22) >> 16;
	t = n3 >> 39;
	n3 &= 0x0000007FFFFFFFFF;
	mull(&lo, &hi, t, R0_LO);
	cc = subb(0, n0, lo, &n0);
	cc = subb(cc, n1, hi, &n1);
	cc = subb(cc, n2, 0, &n2);
	n3 -= cc;
	mull(&lo, &hi, t, R0_HI);
	cc = subb(0, n1, lo, &n1);
	cc = subb(cc, n2, hi, &n2);
	n3 -= cc;
	t = (uint64_t)(*(int64_t *)&n3 >> 63);
	cc = addc(0, n0, R0_LO & t, &n0);
	cc = addc(cc, n1, R0_HI & t, &n1);
	cc = addc(cc, n2, 0, &n2);
	(void)addc(cc, n3, 0x0000008000000000 & t, &n3);

#undef R0_HI
#undef R0_LO

	/*
	 * approx_divr() works on an approximation of n:
	 *   np = floor(n / 2^(a-C+(m-9)/2)
	 * With a = 0 (on xsk233), C = 5, and m = 233, this boils
	 * down to a right shift by 107 bits. Output fits on 125 bits.
	 */
	uint64_t np0 = (n1 >> 43) | (n2 << 21);
	uint64_t np1 = (n2 >> 43) | (n3 << 21);

	/*
	 * q0 = s0*n/r (rounded, extra bits in eta0)
	 * q1 = s1*n/r (rounded, extra bits in eta1)
	 * The "extra bits" are nominally a fractional value between -0.5
	 * and +0.5, but we get them scaled up by 5 bits, as integers.
	 */
	uint64_t q0[2], q1[2];
	int64_t eta0, eta1;

	eta0 = approx_divr(q0, TAU_S0, np0, np1);
	eta1 = approx_divr(q1, TAU_S1, np0, np1);

	/*
	 * Rounding in the plane.
	 * Depending on the values in eta0 and eta1, there are some
	 * possible adjustments to apply to q0 and q1 (by +1 or -1 at
	 * most for each integer).
	 *
	 *   h0 <- 0
	 *   h1 <- 0
	 *   eta <- 2*eta0 - eta1  (here without the scaling up)
	 *   if eta >= 1:
	 *       if (eta0 + 3*eta1) < -1:
	 *           h1 <- -1
	 *       else:
	 *           h0 <- 1
	 *   else:
	 *       if (eta0 - 4*eta1) >= 2:
	 *           h1 <- -1
	 *   if eta < -1:
	 *       if (eta0 + 3*eta1) >= 1:
	 *           h1 <- 1
	 *       else:
	 *           h0 <- -1
	 *   else:
	 *       if (eta0 - 4*eta1) < -2:
	 *           h1 <- 1
	 */
	int64_t eta = 2 * eta0 - eta1;
	int64_t eta3 = eta0 + 3*eta1;
	int64_t eta4 = eta0 - 4*eta1;
	int64_t em1 = ~((eta - 32) >> 8);    /* -1 if eta >= 1 */
	int64_t em2 = (eta3 + 32) >> 8;      /* -1 if (eta0 + 3*eta1) < -1 */
	int64_t em3 = ~((eta4 - 64) >> 8);   /* -1 if (eta0 - 4*eta1) >= 2 */
	int64_t em4 = (eta + 32) >> 8;       /* -1 if eta < -1 */
	int64_t em5 = ~((eta3 - 32) >> 8);   /* -1 if (eta0 + 3*eta1) >= 1 */
	int64_t em6 = (eta4 + 64) >> 8;      /* -1 if (eta0 - 4*eta1) < -2 */
	int64_t h0 = 0;
	int64_t h1 = 0;
	h1 |= em1 & em2 & -1;
	h0 |= em1 & ~em2 & 1;
	h1 |= ~em1 & em3 & -1;
	h1 ^= (em4 & em5) & (h1 ^ 1);
	h0 |= em4 & ~em5 & -1;
	h1 ^= (~em4 & em6) & (h1 ^ 1);

	cc = addc(0, q0[0], (uint64_t)h0, &q0[0]);
	(void)addc(cc, q0[1], (uint64_t)(h0 >> 63), &q0[1]);
	cc = addc(0, q1[0], (uint64_t)h1, &q1[0]);
	(void)addc(cc, q1[1], (uint64_t)(h1 >> 63), &q1[1]);

	/*
	 * u = n - (s0-s1)*q0 - 2*s1*q1
	 * v = s1*q0 - s0*q1
	 *
	 * These values are called r0 and r1, respectively, in FIPS 186-4.
	 * They fit on 118 bits (see Solinas's paper, section 8.3:
	 * the Euclidian norm of the (u,v) vector, equal to u^2 + v^2,
	 * is lower than 4.7095185*r; this implies that |u| and |v| must
	 * both be lower than sqrt(4.7095185*r), which is approximately
	 * 2^116.62. Thus, 118 bits are enough for each of u and v,
	 * including the sign bit. We can therefore compute these values
	 * modulo 2^128.
	 */
	uint64_t u[2], v[2], tt[2];
	u[0] = n0;
	u[1] = n1;
	mul128_trunc(tt, TAU_S0_MINUS_S1, q0);
	sub128(u, u, tt);
	mul128_trunc(tt, TAU_S1_TIMES_2, q1);
	sub128(u, u, tt);
	mul128_trunc(v, TAU_S1, q0);
	mul128_trunc(tt, TAU_S0, q1);
	sub128(v, v, tt);

	/*
	 * We can now produce the NAF digits. We do so by groups of 5,
	 * each group yielding a digit. After 48 digits, the values u
	 * and v must be zero (the decomposition above guarantees that
	 * the NAF length is at most m + a + 3 = 236 values).
	 *
	 * Each inner iteration is:
	 *   - Inspect the two low bits of u and the low bit of v to
	 *     get the next NAF digit:
	 *       if u is even then the next digit is 0
	 *       if u is odd then the next digit is:
	 *           +1 if bit 1 of u equals bit 0 of v
	 *           -1 if bit 1 of u is different from bit 0 of v
	 *     Let z be the next NAF digit.
	 *   - u <- (u - z)/2
	 *   - (u, v) <- (v - u, -u)
	 *
	 * The last two steps can be rewritten into:
	 *     (u, v) <- (v + (z - u)/2, (z - u)/2)
	 *
	 * The next NAF digit then depends only on the two low bits of u,
	 * and the low bit of v. In fact, the next k NAF digits depend only
	 * on the k+1 low bits of u and k+1 low bits of v. We can thus
	 * run 60 iterations using only the low words of u and v, then
	 * propagate the modifications to the full-length u and v only
	 * once every 60 iterations. 60 iterations correspond to 12 output
	 * digits.
	 *
	 * At any point, values u and v can be expressed as linear
	 * combinations of the original u and v, and of the produced
	 * digits (z_0, z_1...). We can define scaled up version of u and v:
	 *
	 *   su_0 = u_0 (initial value of u)
	 *   sv_0 = v_0 (initial value of v)
	 *   su_{i+1} = 2*sv_i - su_i + (2^i)*z_i
	 *   sv_{i+1} = (2^i)*z_i - su_i
	 *
	 * Then the values of u and v after i iterations (i.e. at the
	 * start of iteration i, since we number iterations from 0) are
	 * su_i/2^i and sv_i/2^i, respectively. We then always have:
	 *
	 *   su_i = A_i*u_0 + B_i*v_0 + Lu_i(z_0, z_1...)
	 *   sv_i = C_i*u_0 + D_i*v_0 + Lv_i(z_0, z_1...)
	 *
	 * with:
	 *
	 *   A_0 = 1
	 *   B_0 = 0
	 *   C_0 = 0
	 *   D_0 = 1
	 *   Lu_0(...) = 0
	 *   Lv_0(...) = 0
	 *   A_{i+1} = 2*C_i - A_i
	 *   B_{i+1} = 2*D_i - B_i
	 *   C_{i+1} = -A_i
	 *   D_{i+1} = -B_i
	 *   Lu_{i+1} = 2*Lv_i - Lu_i + (2^i)*z_i
	 *   Lv_{i+1} = -Lu_i + (2^i)*z_i
	 *
	 * It is then straightfoward to compute the coefficients after 60
	 * iterations:
	 *
	 *   A_60 = -1146312001
	 *   B_60 =   493856550
	 *   C_60 =  -246928275
	 *   D_60 =  -899383726
	 *
	 * The coefficients inside Lu_60 and Lv_60 can be large (up to
	 * 2^59), but the sum of their absolute values is less than
	 * 2^61. Thus, even without taking into account that the z
	 * digits are NAF and thus only half of them (at most) can be
	 * non-zero, it is guaranteed that the applied Lu_60 and Lv_60
	 * can fit in 62 bits (including the sign bits).
	 *
	 * We thus organize the loop as follows:
	 *  - We separate the 240 iterations into 4 sequences of 60.
	 *  - Each sequence of 60 iterations is itself a sub-sequence
	 *    of 12 chunks of 5 iterations.
	 *  - Each chunk yields one output digit.
	 *  - A sequence of 60 iteration works on truncated values of
	 *    u and v (low word only). The Lu and Lv values are updated
	 *    with each iteration.
	 *  - The full-size u and v are updated only after each sequence
	 *    of 60 iterations.
	 */
	for (int i = 0; i < 4; i ++) {
		/*
		 * Extract next 12 digits; these depends only on the
		 * low 61 bits of the two running values.
		 *
		 * Instead of dividing z - xu by 2, we left-shift the
		 * mask xk, which contains the current "bit of interest".
		 */
		uint64_t xu = u[0];
		uint64_t xv = v[0];
		uint64_t Lu = 0;
		uint64_t Lv = 0;
		uint64_t xk = 1;
		for (int j = 0; j < 12; j ++) {
			uint64_t d = 0;
			for (int k = 0; k < 5; k ++) {
				uint64_t t1, t2, t3, t4, z;

				t1 = xk & xu;
				t2 = t1 + t1;
				t3 = xv + xv;
				t2 &= xu ^ t3;
				z = t1 - t2;
				d += z;
				xv = z - xu;
				xu = xv + t3;
				xk += xk;
				t4 = Lv + Lv;
				Lv = z - Lu;
				Lu = Lv + t4;
			}
			d >>= 5 * j;
			sd[12 * i + j] = (int8_t)*(int64_t *)&d;
		}

#if 0
		/*
		 * Alternate version using SSE2 to mutualize updates to
		 * (xu, xv) with updates to (Lu, Lv) (the values Lu and Lv
		 * being kept in the upper halves of the registers). It
		 * does not seem to be faster in practice (on Skylake cores).
		 */
		__m128i xu = _mm_set_epi64x(0, u[0]);
		__m128i xv = _mm_set_epi64x(0, v[0]);
		__m128i xk = _mm_set_epi64x(0, 1);
		for (int j = 0; j < 12; j ++) {
			__m128i xd = _mm_setzero_si128();
			for (int k = 0; k < 5; k ++) {
				__m128i t1, t2, t3, z;

				t1 = _mm_and_si128(xk, xu);
				t2 = _mm_add_epi64(t1, t1);
				t3 = _mm_add_epi64(xv, xv);
				t2 = _mm_and_si128(t2, _mm_xor_si128(xu, t3));
				z = _mm_sub_epi64(t1, t2);
				z = _mm_shuffle_epi32(z, 0x44);
				xd = _mm_add_epi64(xd, z);
				xv = _mm_sub_epi64(z, xu);
				xu = _mm_add_epi64(xv, t3);
				xk = _mm_add_epi64(xk, xk);
			}

			uint64_t d = _mm_cvtsi128_si64(xd);
			sd[12 * i + j] = (int8_t)(*(int64_t *)&d >> (5 * j));
		}
		uint64_t Lu = _mm_cvtsi128_si64(_mm_bsrli_si128(xu, 8));
		uint64_t Lv = _mm_cvtsi128_si64(_mm_bsrli_si128(xv, 8));
#endif

		/*
		 * u' <- (A_60*u + B_60*v + Lu)/2^60
		 * v' <- (C_60*u + D_60*v + Lv)/2^60
		 *
		 * coefficients:
		 *   A_60 = -1146312001
		 *   B_60 =   493856550
		 *   C_60 =  -246928275
		 *   D_60 =  -899383726
		 *
		 * Take care that all values are signed. Intermediate values
		 * necessarily fit on 192 bits (including the sign bits).
		 */
		uint64_t nu[3], nv[3], nt[3];

		mul64x128(nu, v, 493856550);
		mul64x128(nt, u, 1146312001);
		sub192(nu, nu, nt);
		add192p64(nu, nu, Lu);

		mul64x128(nv, u, 246928275);
		neg192(nv, nv);
		mul64x128(nt, v, 899383726);
		sub192(nv, nv, nt);
		add192p64(nv, nv, Lv);

		u[0] = (nu[0] >> 60) | (nu[1] << 4);
		u[1] = (nu[1] >> 60) | (nu[2] << 4);
		v[0] = (nv[0] >> 60) | (nv[1] << 4);
		v[1] = (nv[1] >> 60) | (nv[2] << 4);
	}

	/*
	 * Simpler version (but slower) for the production of the digits
	 * from the computed u and v. In this version, we apply the
	 * transforms on the full 128-bit values u and v directly.
	 *

	uint64_t ul = u[0];
	uint64_t uh = u[1];
	uint64_t vl = v[0];
	uint64_t vh = v[1];
	for (int i = 0; i < 48; i ++) {
		uint64_t d = 0;
		for (int j = 0; j < 5; j ++) {
			uint64_t z;

			z = -(ul & 1) & (1 - ((ul + vl + vl) & 2));
			d += z << j;
			cc = subb(0, ul, z, &ul);
			(void)subb(cc, uh, *(int64_t *)&z >> 1, &uh);
			uint64_t xl = (ul >> 1) | (uh << 63);
			uint64_t xh = (uint64_t)(*(int64_t *)&uh >> 1);
			cc = subb(0, vl, xl, &ul);
			(void)subb(cc, vh, xh, &uh);
			cc = subb(0, 0, xl, &vl);
			(void)subb(cc, 0, xh, &vh);
		}
		sd[i] = (int8_t)*(int64_t *)&d;
	}
	*/
}

/*
 * Apply the Frobenius endomorphism on a point (i.e. square all coordinates).
 */
static inline void
xsk233_inner_frob(inner_point *p3, const inner_point *p1)
{
	gf_sqr(&p3->x, &p1->x);
	gf_sqr(&p3->z, &p1->z);
	gf_sqr(&p3->s, &p1->s);
	gf_sqr(&p3->t, &p1->t);
}

/* see xs233.h */
void
xsk233_mul_frob(xsk233_point *p3,
	const xsk233_point *p1, const void *n, size_t n_len)
{
	/*
	 * If P = (x, s), then let frob(P) = (x^2, s^2) (Frobenius
	 * endomorphism on the curve), and let frob_n(P) = (x^(2^n), s^(2^n))
	 * (n successive applications of frob()).
	 *
	 * Let [d0 d1 d2 d3 d4] be a NAF (i.e. each d_i is either 0, -1
	 * or +1, and no two consecutive d_i are non-zero). The index
	 * d is computed from the NAF as:
	 *   d = d0 + 2*d1 + 4*d2 + 8*d3 + 16*d4
	 * For NAF inputs, this encoding is a bijection between the
	 * 5-digit NAFs, and the integers in the -21 to +21 range.
	 *
	 * For NAFs such that d > 0, we set:
	 *   win[d - 1] = d0*P + d1*frob_1(P) + d2*frob_2(P)
	 *                + d3*frob_3(P) + d4*frob_4(P)
	 */
	inner_point win[21];

	point_load(&win[0], p1);
	xsk233_inner_frob(&win[1], &win[0]);
	xsk233_inner_frob(&win[3], &win[1]);
	xsk233_inner_sub(&win[2], &win[3], &win[0]);
	xsk233_inner_add(&win[4], &win[3], &win[0]);
	xsk233_inner_frob(&win[5], &win[2]);
	xsk233_inner_frob(&win[7], &win[3]);
	xsk233_inner_sub(&win[6], &win[7], &win[0]);
	xsk233_inner_add(&win[8], &win[7], &win[0]);
	xsk233_inner_frob(&win[15], &win[7]);
	for (int i = 1; i <= 5; i += 2) {
		xsk233_inner_sub(&win[15 - i], &win[15], &win[i - 1]);
		xsk233_inner_add(&win[15 + i], &win[15], &win[i - 1]);
	}
	for (int i = 10; i <= 14; i += 2) {
		xsk233_inner_frob(&win[i - 1], &win[(i >> 1) - 1]);
	}
	for (int i = 18; i <= 20; i += 2) {
		xsk233_inner_frob(&win[i - 1], &win[(i >> 1) - 1]);
	}

	/*
	 * Recode the scalar in the proper format for use with the
	 * endomorphism.
	 */
	int8_t sd[48];

	if (n_len == 30) {
		xsk233_tau_recode(sd, n);
	} else {
		uint8_t tmp[30];

		memcpy(tmp, n, n_len);
		memset(tmp + n_len, 0, (sizeof tmp) - n_len);
		xsk233_tau_recode(sd, tmp);
	}

	/*
	 * Set the accumulator to the point corresponding to the
	 * top digit, then process all other digits in high-to-low
	 * order.
	 */
	inner_point ip;

	inner_lookup(&ip, win, 21, sd[47], &GF_ONE);
	for (int i = 46; i >= 0; i --) {
		inner_point iq;

#if GF_MUL_FAST
		for (int j = 0; j < 5; j ++) {
			gf_sqr(&ip.x, &ip.x);
			gf_sqr(&ip.s, &ip.s);
			gf_sqr(&ip.z, &ip.z);
		}
		gf_mul(&ip.t, &ip.x, &ip.z);
#else
		for (int j = 0; j < 5; j ++) {
			xsk233_inner_frob(&ip, &ip);
		}
#endif
		inner_lookup(&iq, win, 21, sd[i], &GF_ONE);
		xsk233_inner_add(&ip, &ip, &iq);
	}

	point_store(p3, &ip);
}

/* see xs233.h */
void
xsk233_mulgen_frob(xsk233_point *p3, const void *n, size_t n_len)
{
	/*
	 * Recode the scalar in the proper format for use with the
	 * endomorphism.
	 */
	int8_t sd[48];

	if (n_len == 30) {
		xsk233_tau_recode(sd, n);
	} else {
		uint8_t tmp[30];

		memcpy(tmp, n, n_len);
		memset(tmp + n_len, 0, (sizeof tmp) - n_len);
		xsk233_tau_recode(sd, tmp);
	}

	/*
	 * Set the accumulator by looking up the top digits of the
	 * four chunks of the scalar.
	 */
	inner_point ip;
	inner_point_affine Qa;

	inner_lookup_affine(&Qa, XSK233_TAU_G0, 21, sd[11], &GF_ONE);
	ip.x = Qa.x;
	ip.s = Qa.s;
	ip.z = GF_ONE;
	ip.t = Qa.x;
	inner_lookup_affine(&Qa, XSK233_TAU_G60, 21, sd[23], &GF_ONE);
	xsk233_add_mixed(&ip, &ip, &Qa);
	inner_lookup_affine(&Qa, XSK233_TAU_G120, 21, sd[35], &GF_ONE);
	xsk233_add_mixed(&ip, &ip, &Qa);
	inner_lookup_affine(&Qa, XSK233_TAU_G180, 21, sd[47], &GF_ONE);
	xsk233_add_mixed(&ip, &ip, &Qa);

	/*
	 * Process all other signed digits in high-to-low order.
	 */
	for (int i = 10; i >= 0; i --) {
#if GF_MUL_FAST
		for (int j = 0; j < 5; j ++) {
			gf_sqr(&ip.x, &ip.x);
			gf_sqr(&ip.s, &ip.s);
			gf_sqr(&ip.z, &ip.z);
		}
		gf_mul(&ip.t, &ip.x, &ip.z);
#else
		for (int j = 0; j < 5; j ++) {
			xsk233_inner_frob(&ip, &ip);
		}
#endif
		inner_lookup_affine(&Qa, XSK233_TAU_G0, 21,
			sd[i], &GF_ONE);
		xsk233_add_mixed(&ip, &ip, &Qa);
		inner_lookup_affine(&Qa, XSK233_TAU_G60, 21,
			sd[i + 12], &GF_ONE);
		xsk233_add_mixed(&ip, &ip, &Qa);
		inner_lookup_affine(&Qa, XSK233_TAU_G120, 21,
			sd[i + 24], &GF_ONE);
		xsk233_add_mixed(&ip, &ip, &Qa);
		inner_lookup_affine(&Qa, XSK233_TAU_G180, 21,
			sd[i + 36], &GF_ONE);
		xsk233_add_mixed(&ip, &ip, &Qa);
	}
	point_store(p3, &ip);
}

/* ===================================================================== */

/* i*(2^0)*G  for i = 1 to 16 */
static const inner_point_affine XSK233_G0[] = {
	{ GFw64be(0x000001ECB92776D0, 0xFB3DEC476585B906,
	          0x5724EF7E1966BF54, 0xA850E5CBDDAA1BE6),
	  GFw64be(0x000000EDFF3B4D4E, 0x5BAA47FCFDF3669D,
	          0xF7193250076F96C1, 0x66F9E0BF367D9A99) },
	{ GFw64be(0x000000A6217325BC, 0x2426B0E995AD7E3F,
	          0xA8BA1439CFCDBFA5, 0x6ED496768224E403),
	  GFw64be(0x00000019814E3587, 0xB208935EFF71E59C,
	          0x9BADA851345A0DF3, 0x4882B58E728B7BE8) },
	{ GFw64be(0x0000016C2788C24E, 0x8C5649B4D39E9E04,
	          0x18C1B299C24A9D7A, 0xB803B0D5EB06A145),
	  GFw64be(0x000000F2CBD8D377, 0x07C9343A50677BDB,
	          0x89134BDD9A1FD022, 0xC15CA3CF2CE45743) },
	{ GFw64be(0x00000010A1B22D06, 0x2EE97AA3EB14B7C1,
	          0xFB84C8189C4DB8DC, 0x1977069CE3A4D0C9),
	  GFw64be(0x000000CBA628F7C1, 0x4186BD1D7385A9ED,
	          0xB36FB8BCB7B59025, 0x22E23B717D36C690) },
	{ GFw64be(0x000001DE7F60DAC2, 0xAE3E4049EDBE3DB2,
	          0xCFF3A40A2B274F93, 0xA5ABD0FF71717DE7),
	  GFw64be(0x000000BE430FDB21, 0xF8406B98D0E4F387,
	          0xDB4E2D7A559EB47D, 0xAD071F950C8D9079) },
	{ GFw64be(0x0000015A8FF292B7, 0xA73D8B0F45AAAD90,
	          0x3770E3D054A1B362, 0x32E749E434F7ED94),
	  GFw64be(0x000001F35579FF15, 0x026DFDB1F0D7BF51,
	          0x5ECCFE93B7E2483C, 0xECD28CF9E43438F6) },
	{ GFw64be(0x0000007EC69A3B8B, 0xB1B3F204B85CCEAC,
	          0x398990AF8CDDB931, 0xAB419C601977FF19),
	  GFw64be(0x0000000DB05CE355, 0xFEE06A6B56BE41DA,
	          0x39E578D1F1FE8A4B, 0x72F0705D39657407) },
	{ GFw64be(0x000001E19E3E2447, 0xB21D461CC8F46360,
	          0x05897282B9DABF79, 0xD25CBED568A3C273),
	  GFw64be(0x0000006A949E2A7A, 0xA7FCF7C456E7731F,
	          0x0BABBAC03580EEDC, 0x3A433E5D596F345B) },
	{ GFw64be(0x000001AD967BE56E, 0xE41729E9F137B2A8,
	          0x77128439256D4B22, 0x09D26AB6E60E8D0B),
	  GFw64be(0x0000013811ECF2B4, 0x8240B8C78E351CDC,
	          0xB0796CBB97B968CA, 0xB8AD4521D9D809FE) },
	{ GFw64be(0x00000033491350DC, 0xE181007EE44B320F,
	          0xC2AB8DF8FEB533AE, 0x74C7F388903A788D),
	  GFw64be(0x000000D69DF6294B, 0x0CF9D252333979D6,
	          0x6BB854C15892C61C, 0xB5C6E24A7FA573F9) },
	{ GFw64be(0x0000005DC7D46A74, 0x88FB936F90CFC137,
	          0x1D766BAF70A28D24, 0x733F8FE92B7BC1BD),
	  GFw64be(0x000000A3F30ECA08, 0xDB5184CE2CAF40DB,
	          0xDD3DA01773FFCCD5, 0xB10409AB34229B5A) },
	{ GFw64be(0x000000FB60FBF891, 0x1D0986ACFE8FCB3C,
	          0x04B434C6D625C0AB, 0x67323257CFBA7BD5),
	  GFw64be(0x0000001A20183F63, 0xCB49C50C7382E4B5,
	          0xC4A3235FA8C2416A, 0x5B6052F95ED36ABE) },
	{ GFw64be(0x00000180588C8EF7, 0xCD07B21D500DA5A0,
	          0x54FE83B6CEF15DE9, 0xEC9EB9DB9CC694A8),
	  GFw64be(0x0000018E68F0D5E6, 0xD4ED1FD888D12380,
	          0x57348CCC63C113E3, 0x29094B2ED6CD42C9) },
	{ GFw64be(0x000000939799CBE4, 0x3A3D0AF58314FC90,
	          0x2BE89A1F7A89C148, 0x738B56D392AC710F),
	  GFw64be(0x000000305CFF88FB, 0xB68A47EF61E9EC5E,
	          0x6111DEB07D1D1E8A, 0xF76978C22AD47A60) },
	{ GFw64be(0x000000E1105954BA, 0xF90E2377CA97EFC7,
	          0x6DC9CDD094AAA26C, 0x1C6C013E358343D9),
	  GFw64be(0x000001105A1F275A, 0xA55D7958AE3F5729,
	          0xFF86602046C88ACE, 0x716640478C8BFF24) },
	{ GFw64be(0x000000F7029F1365, 0x51C442435F6F251E,
	          0x4CAEB88CD30D9443, 0x4079195330A98EA4),
	  GFw64be(0x000000A6096BD15E, 0xF44EEB29C8E84434,
	          0xB9415BF2EE92CA71, 0xEB7E2CC0C3404301) }
};

/* i*(2^60)*G  for i = 1 to 16 */
static const inner_point_affine XSK233_G60[] = {
	{ GFw64be(0x000001B9EB6C50CB, 0x61CD295506528988,
	          0x8029062B1063AB62, 0x404AFF3FAA016846),
	  GFw64be(0x0000005BC887A05F, 0x4EEB4359725072FB,
	          0xB32A96D7206B24F5, 0x665998AE2E9D7F48) },
	{ GFw64be(0x00000026A411FDC2, 0x42E5396E10DCD79D,
	          0xD64C84726E307352, 0x364CEF2F5F97FE70),
	  GFw64be(0x000000079837E135, 0x0DEE28A2364FB55D,
	          0x0116DA4CD10ACDB0, 0x4B144B86F65C236F) },
	{ GFw64be(0x000000C07FEDBB26, 0xD66E58CE515015C2,
	          0xA7938A6A1C519620, 0xB60490E275627CE6),
	  GFw64be(0x00000172383AE009, 0xBAF619CA4F88DB35,
	          0x432B7AAC828BB936, 0x3E8BEF9812E151C4) },
	{ GFw64be(0x0000002257DEC6F7, 0x394C8C1B086B3454,
	          0x1749C5B31AD5E64C, 0x9A7E41576F2E791E),
	  GFw64be(0x000000D53B1CF461, 0xE6B9DD6B3A006A9F,
	          0x2476AD5B710DD80A, 0xBCC2BB867ED6DDA6) },
	{ GFw64be(0x0000004144A89F92, 0x512B1841C0345A6D,
	          0xFAB6D752E7D6739A, 0x03F5F25C8164D74B),
	  GFw64be(0x0000011997E08323, 0xE020ED7149EAEB78,
	          0x6560A1112D64DF80, 0x8055100DA5EC8459) },
	{ GFw64be(0x000001E1350326E0, 0x742FA0E39F558243,
	          0xCDF5A7EDB6E3AA74, 0x6D821F94EB95A1BD),
	  GFw64be(0x000000A7BACB5FD2, 0xE7EEA19B8E6010A3,
	          0x5D40455303E336DC, 0xCCB10521075B9E40) },
	{ GFw64be(0x000000096DE4A7F0, 0x474214A474E10622,
	          0x2D38094366146757, 0x404866E174555296),
	  GFw64be(0x000000E7E94A6836, 0x52865EE6E2EA8C60,
	          0xE1F60163D1004C8F, 0x10BB65DCAB6B554B) },
	{ GFw64be(0x000000A5ECF152B2, 0xE81982B6EC28ECD0,
	          0x8511246E76BC6072, 0xF568804A63BBA671),
	  GFw64be(0x00000178277D6AD6, 0x493B0565B7A3969B,
	          0x5994857A17265C2F, 0x01495CBE73D4706C) },
	{ GFw64be(0x0000014B0CA5CB4A, 0xA39117FC72A18199,
	          0x9B25B417F8AB6B4F, 0x21E51FD6744E1070),
	  GFw64be(0x0000013D4EE3DC93, 0xCA79905FF5247987,
	          0x743F17DAA403BCA0, 0x465FBBD4A2CC223A) },
	{ GFw64be(0x000000830F091CE6, 0x184219FF1C4A2099,
	          0x76888CB5565970DD, 0xC5E163D392E7CD38),
	  GFw64be(0x0000005583E07311, 0xDF8286C4B97328D1,
	          0xC8C61EE12C16D1BA, 0x9255D9057C03DFCB) },
	{ GFw64be(0x0000004198BC7036, 0xB289A591D9415972,
	          0x05F0440F9655BAE4, 0xCFD417291373E94F),
	  GFw64be(0x000000A38811D48F, 0x7439D6FCDB8390BC,
	          0x0A91FF39CE63BEEF, 0x90153B99D75CB0E7) },
	{ GFw64be(0x0000007A29B2D631, 0x40F452EFAA70BA41,
	          0xD4AE2DC3866FE1E3, 0x993DD46D6F927FDB),
	  GFw64be(0x0000000B94D94469, 0xD3EC6DAA5F86E00C,
	          0x6E290CBDD24E1EF4, 0xCF1CBE7F2C900A6B) },
	{ GFw64be(0x000001040B5DF6A9, 0xCA011C583AC8F3F6,
	          0x4CAAA00E26069F8E, 0x488EBA4D49B3D446),
	  GFw64be(0x0000012CA96DF7A6, 0xD1CA259B97CB418F,
	          0x7363C1A1BB2B7D0E, 0x15833C1B0E446F88) },
	{ GFw64be(0x000001DA4AABCAA1, 0x94F7797518E68BD5,
	          0x9F19345D2E328513, 0xB19733268A0A0FE6),
	  GFw64be(0x000000EB44E3D508, 0xC07DAE6506D8EC07,
	          0x4BD09C362A833FF5, 0x9D1A40466FB13CB9) },
	{ GFw64be(0x0000019B3A518B4F, 0x4B18E0DE97FE76B2,
	          0x455CCBB952FAB3A6, 0x65653FD63E43DD75),
	  GFw64be(0x000000E22F856FC2, 0x62B376E23D9F5974,
	          0xDB13C718727F8FFF, 0xE43CDF8EAC57989F) },
	{ GFw64be(0x0000007AF834BAB5, 0xFBBC04C4AD54A16A,
	          0x3B042AA3F7F97954, 0x154E587BE15500E3),
	  GFw64be(0x00000122D95690CA, 0x22E5D7632B4252B8,
	          0xCDB60DA01503B112, 0x97E8EF6B13277F22) }
};

/* i*(2^120)*G  for i = 1 to 16 */
static const inner_point_affine XSK233_G120[] = {
	{ GFw64be(0x0000003D000C5FA0, 0x89E2D2B4ADF2D72A,
	          0xF4DA1C0AED95313F, 0x9628CECE4DF50601),
	  GFw64be(0x000000378848CC75, 0xEAA1FDDC8126ECBF,
	          0x16289C9029DDF11E, 0x858D955BD4966CC6) },
	{ GFw64be(0x00000029CCFD790F, 0x82DE4D6818DC2783,
	          0xCE83B6DA4EDF3D7F, 0x78C5B6D464BC879A),
	  GFw64be(0x000000344CEBE69C, 0x29780EB437B27D54,
	          0xFBF75AC7212571C5, 0x69C27CE1DEDAB525) },
	{ GFw64be(0x0000014FFD1D0E82, 0xA9F1B5BF2464C25C,
	          0x9FD53610E26EEB29, 0x312EE867D2683A16),
	  GFw64be(0x00000002D6D4923E, 0x0A20A90C421C94E1,
	          0x1E2CDA17B86C0D0B, 0x61C69040B574AC2D) },
	{ GFw64be(0x0000006E2D09285B, 0x24D43F6167B96457,
	          0x7223BC586CFAB200, 0x44901667CBCB4896),
	  GFw64be(0x000001DE089678C9, 0x26027D0E9213BDAF,
	          0xAF9DF0E305AAD4F0, 0x578055047E0BCF6A) },
	{ GFw64be(0x000000A33CDFBD0F, 0xEDE74F39E6287473,
	          0x16293DD35E37B980, 0x177EB14B7DA66E11),
	  GFw64be(0x0000017383596476, 0x737DDBB963DB8EE6,
	          0xEF9598267D5EC146, 0x73743BC5EBEE8E40) },
	{ GFw64be(0x000001E8E5ED0966, 0x4F9E819DC2052132,
	          0x1BC608C6895BF5F9, 0xA56E460C8ECD3D7B),
	  GFw64be(0x0000009E902CA42D, 0x5CD6C6B813C90D4E,
	          0xFCAA32EC7830F0D3, 0xF0515F4EEC15F200) },
	{ GFw64be(0x0000016656C4D737, 0xD08F79BF23178EFA,
	          0xED511DBF79BEED32, 0x6F2B789A101DCBFE),
	  GFw64be(0x000001504C689D84, 0x09984C05EA9071D4,
	          0xD013771C87D4CC0B, 0x717E255DB0C6571B) },
	{ GFw64be(0x00000195D1A83711, 0xB7E91E3D8DCD9EEC,
	          0x3FB1988492C990A9, 0xB2732BF22617DE45),
	  GFw64be(0x0000008384A9807E, 0x391BCA166AD1D576,
	          0xE6F962204938333F, 0x88779F19E91F5351) },
	{ GFw64be(0x000001534A10E26A, 0xA0D3D5B05E25CBF1,
	          0x6BECBC21F6B449D9, 0xFED03699431D4026),
	  GFw64be(0x0000002B82266C45, 0x21C3AC8035C6FAB4,
	          0xA2DAFF26EFA92405, 0x3C63F3DAC89FDD09) },
	{ GFw64be(0x000000BBE0DC72EE, 0x283ACC6FCA9CDF29,
	          0x52BC9C79B3DDDC5B, 0x84F4D4C014DD1BB1),
	  GFw64be(0x00000191B09225CF, 0x0BAAE5266662ABDB,
	          0xF2E903468437B394, 0x0F5080E6BE7CBCCF) },
	{ GFw64be(0x00000109ADBFADEF, 0xC37427E220329909,
	          0x3CC3BE355A740AA5, 0x543984280861DED6),
	  GFw64be(0x000000DDA3861451, 0xCB60D919D3BF32C3,
	          0x82F8A5179831AFAE, 0x75D2E50616E195F5) },
	{ GFw64be(0x000001BE334008ED, 0xDCD31DE33E6B354B,
	          0x6E5D2566570F1810, 0x6E1F2A24597F0E28),
	  GFw64be(0x000001BE170AA115, 0xAA4BA2F76398D2D0,
	          0xAC53F3068E9D6120, 0xEE6677C59AF7DE8E) },
	{ GFw64be(0x000001429406C156, 0xA552215022812CE1,
	          0x74365F1B94B46D3E, 0x1BDE0B9A17650E58),
	  GFw64be(0x00000050255935A6, 0x56B620F26BF28421,
	          0x0B0E27F48E7827C0, 0x1DFF59D1BC95F2DC) },
	{ GFw64be(0x0000003C61A5F7C8, 0x4D2C0EDB1CCB2523,
	          0x5AC093B8ED31C2A8, 0x0BC39FA18F904E0E),
	  GFw64be(0x000000FD14545CDE, 0x5FBFA8A0B73A4255,
	          0xA6D87187319C8A4B, 0x676275F15C13EBA5) },
	{ GFw64be(0x0000009E9A94E4CC, 0x7D3338AD7F0B5575,
	          0x97600C8BC18FA68C, 0x7CB28B68C9BC6B5A),
	  GFw64be(0x0000004EF96D0E30, 0xF90235288567A664,
	          0xAA6A17B0598A390C, 0x0C7A9234559E9924) },
	{ GFw64be(0x0000014487CC3FBF, 0x5BAB88C4D32E7E1C,
	          0x5F99DEAAC70F562D, 0x180E1AD4171466F1),
	  GFw64be(0x000001407C3174BA, 0x4EDDD9E21B0C82E8,
	          0x177D4931F629EE81, 0xF031CBE029495193) }
};

/* i*(2^180)*G  for i = 1 to 16 */
static const inner_point_affine XSK233_G180[] = {
	{ GFw64be(0x000000993790BDA2, 0xD507195588FD8FA2,
	          0x477F39D11101B1BF, 0x5F5FBBA45CF4651B),
	  GFw64be(0x000001C9CF2AB803, 0x1D29C897B8E229DE,
	          0xE3A33A4562C505AC, 0xB1DFC5F8F5E3869A) },
	{ GFw64be(0x0000011A1595CA59, 0x942E97D688C5A247,
	          0x3BCD07168351EDCB, 0x3E925A47575D8ED5),
	  GFw64be(0x000001D7A9AA2321, 0x1E67BB8977610095,
	          0xF4BD817EDE7A2E84, 0x461DB8F9FF3821F1) },
	{ GFw64be(0x000001E52EEE00D7, 0x87CEEBEE228AC0E1,
	          0xEB488FEFB786FDB4, 0x6A03B6F794E8351E),
	  GFw64be(0x000000BDF5C8B696, 0xF655FEE9C30838BB,
	          0xA0BBDB64A068FF66, 0xBF234D2712BD9ACF) },
	{ GFw64be(0x000000D4D85A5EA1, 0x81E7112717439FA5,
	          0x5B68FAE0C2B71770, 0xCF7CEE6788843368),
	  GFw64be(0x0000010DC10AF268, 0x3B0D665E85C8BA03,
	          0xD30E6EC59CF2403D, 0xF9E16006647579B5) },
	{ GFw64be(0x000000129B629392, 0x3135359665FCC39C,
	          0x6A53CE8EB3DC2733, 0xAF1517C157CC5A38),
	  GFw64be(0x0000007D3B0E42DD, 0x9360BA86F4052479,
	          0xE17547015ABB0961, 0x2704B6C870081AE7) },
	{ GFw64be(0x000001A277C5407F, 0xB3C080C0215F437A,
	          0xEDF75E75EBAFBAC6, 0xF02662152725C872),
	  GFw64be(0x000000CD31A3D65A, 0x7009A093FD2DF479,
	          0x5E1791A61677589C, 0xB7EC51AF9DAAC862) },
	{ GFw64be(0x000001A317C08239, 0x20B162A784D3220C,
	          0xA90E89FB9AB060FF, 0xDD980EE53BF6A6BF),
	  GFw64be(0x000000DB9C68B9AA, 0x6F2A87103F9F516A,
	          0x74A502ECAA0327CE, 0xABF64CD367E07CA9) },
	{ GFw64be(0x0000001E52E7F4D9, 0xA438A5F4594350EA,
	          0xEEE23B1B908EF0D5, 0xB9AB5865114E514A),
	  GFw64be(0x000000C0E99B1CB0, 0xF404CAECB286E55F,
	          0xCCC5F674232AAEF7, 0x4F35A22D8D88AC0D) },
	{ GFw64be(0x000000A9956E934D, 0xF96AEB61E2794300,
	          0xD6FFC5433208B5D4, 0x7CD23EE4635467C3),
	  GFw64be(0x000001876037ACAD, 0xD4864CD28BBE8AD2,
	          0x431467808DF30510, 0x8C094D84586629C8) },
	{ GFw64be(0x000001904B33964E, 0xC5A76F1504825512,
	          0x766FC85590118FE4, 0x2BD5387855B4557C),
	  GFw64be(0x00000026CD22EA68, 0x766C80E739FC881D,
	          0x675FA0D5757D6F96, 0xCB1A13719B6AC4C3) },
	{ GFw64be(0x000000B66B92E3FC, 0x61F6D7CAF8896D80,
	          0x4C02413F8A7706F1, 0xBE280D119CF310E1),
	  GFw64be(0x00000159E46F3D4F, 0x095C91546F3DF224,
	          0xA166571059FD3C29, 0x58183CBB89629AB9) },
	{ GFw64be(0x0000017CB27FC203, 0x2121DCEA39971A88,
	          0xC622E87091906305, 0x9ABF43FCCAF52B84),
	  GFw64be(0x00000062C8E783F9, 0x80FF5BFED188C0C5,
	          0x56E4458A548EB425, 0xEE8153B4762F4335) },
	{ GFw64be(0x000001C350833149, 0xA789C209AF6E0CA3,
	          0xA47FD929F08327A1, 0x77760C777DA8E341),
	  GFw64be(0x000000102E0CF74B, 0xCAE1D47B2C4BA8DF,
	          0x26C9AE7AD13CE16C, 0x84D507CB56891733) },
	{ GFw64be(0x000000CAD69C5E76, 0x86933C6323F581A8,
	          0x9DB065C05FFDD7A5, 0x796DB258006EE93E),
	  GFw64be(0x000000CE647FFC50, 0x2ADAAF01720D71CC,
	          0x77A530222721B291, 0x5E996E0586FD7EB3) },
	{ GFw64be(0x0000003F2D1A3946, 0x4A8FED56A052F441,
	          0xD2A537B70A8928F0, 0x22D8C62DE2C05605),
	  GFw64be(0x000000D673560492, 0x8964D38FF93B70B4,
	          0x318D9BA58E747A8D, 0xFC23629D954A0D69) },
	{ GFw64be(0x00000173BBA27EE4, 0xBD6FE1671B78943C,
	          0x692B102440ED77B5, 0x1AC52061FED33270),
	  GFw64be(0x000000C646EC42D0, 0x5A573666523CE6F8,
	          0x6B4953E558F6503B, 0xFA23152BC8111BAA) }
};

/* ===================================================================== */
/*
 * Precomputed tables for multiples of G, for use with xsk233_mulgen_frob().
 */

/* tau_i(frob_0(G))  for i = 1 to 21 */
static const inner_point_affine XSK233_TAU_G0[] = {
	{ GFw64be(0x000001ECB92776D0, 0xFB3DEC476585B906,
	          0x5724EF7E1966BF54, 0xA850E5CBDDAA1BE6),
	  GFw64be(0x000000EDFF3B4D4E, 0x5BAA47FCFDF3669D,
	          0xF7193250076F96C1, 0x66F9E0BF367D9A99) },
	{ GFw64be(0x000000BAD6FDBF74, 0xA36BBE16ED5862F1,
	          0x736C297428E015B1, 0x09B6A5C01E58CAC7),
	  GFw64be(0x000001C30F24330A, 0x82B7967DC1BF09EC,
	          0x4B60C85D1744ADAB, 0xD6CB4A74ACFE9AD1) },
	{ GFw64be(0x000000590A3FE1E8, 0x8816FE30E5301A41,
	          0x27048C06260C1412, 0x88A142FCA16F3465),
	  GFw64be(0x000000E7B668F2C8, 0x19ABB446D341FA02,
	          0x1380AE6F1767B578, 0xC1BA2A7D82B2694A) },
	{ GFw64be(0x000000522E63B7B0, 0x264AFCA22190701A,
	          0x4EBA1053F16504CD, 0xE50F073F818B878A),
	  GFw64be(0x00000020D2C8B1F9, 0x8B3FB832C6768CD9,
	          0xDA6262AF8C723CE2, 0xEBD020FA7949768B) },
	{ GFw64be(0x000000FF4A3E695A, 0xBF1FB6845B8F12B2,
	          0x98221D3622497219, 0x901969DC5ACDE434),
	  GFw64be(0x000001F5D7247ECD, 0x2DF0F1659DA5396B,
	          0x6FC51A73405EFE5D, 0xFD4591A2DD008634) },
	{ GFw64be(0x000000986050203E, 0x04160AFA89126E00,
	          0xC9F2DB54C8010DC3, 0x8501F25850C61FC5),
	  GFw64be(0x000000AA4C7416D5, 0xAB1FB615C71AC5B5,
	          0xCC7B9C455877578D, 0x55EBF310561496B4) },
	{ GFw64be(0x0000014DF0AC95C9, 0xCFF1BFDB0EAC0DFB,
	          0xB0E3F81C51F2F886, 0xBB3D32FE5D78C512),
	  GFw64be(0x00000132CFF7B531, 0x0AC749BE037F6806,
	          0x0D1E526372D9395A, 0x94863CDA0DEB1260) },
	{ GFw64be(0x0000014C8B2A1B87, 0xFD0316BB88109C49,
	          0x4D6BFAEA36EA99E4, 0xC40F9451E7587F23),
	  GFw64be(0x0000008494A6CCD7, 0xEAF29FA605D84CAF,
	          0x3A56128130A5234A, 0x9F303C39F5F9E31A) },
	{ GFw64be(0x0000012A12B8D5DA, 0x7E48363C4B5DD0E5,
	          0x7B280864498748FC, 0x40291B46253E90C8),
	  GFw64be(0x000000C646DD4AF4, 0x4B646B4C138EF8F0,
	          0x27D34EF11726CDD1, 0x6AE75E76AB384E2B) },
	{ GFw64be(0x000000A6217325BC, 0x2426B0E995AD7E3F,
	          0xA8BA1439CFCDBFA5, 0x6ED496768224E403),
	  GFw64be(0x00000019814E3587, 0xB208935EFF71E59C,
	          0x9BADA851345A0DF3, 0x4882B58E728B7BE8) },
	{ GFw64be(0x00000010D2D631ED, 0x5B61DECA65B9C1C7,
	          0x2398A4C2287F8F0B, 0xDFCB4348F65A9B30),
	  GFw64be(0x0000013B4B0FECEF, 0x28CE4B483F284DFB,
	          0x2769E0EE0D7C0C56, 0x1DD9BF5F0D8FDD78) },
	{ GFw64be(0x0000018CF16F1330, 0xD0CAAA8B2AD25F0D,
	          0xB2B9662B47A68D81, 0xD31E7A05017970BE),
	  GFw64be(0x00000145EBD23AB1, 0x13689D1F399C7149,
	          0x5AD10027795C67C4, 0xF72C59A7D9BD794B) },
	{ GFw64be(0x0000017CFC9671A4, 0xDB9F93407BFFFF9B,
	          0x5866BA7CE61E0556, 0x0BF95FDBD74CC796),
	  GFw64be(0x00000053518F69B9, 0x43FF8790571478EA,
	          0xD5AEB614F39F5799, 0xEB84299F25CFEA8A) },
	{ GFw64be(0x0000000FDD402B70, 0xB38B7F267544A35B,
	          0x9144143229C19666, 0x68043DAEA2F2907E),
	  GFw64be(0x000001D493A49E85, 0xBF24536B85CF0DD7,
	          0x224D34CBBBB646CB, 0xAE8B21101E412D85) },
	{ GFw64be(0x000001CF6323F42A, 0xB54CFA7B0331AF14,
	          0xA4E51C183E19D176, 0x91001E9C80D802C8),
	  GFw64be(0x0000010EEB3688E6, 0x47069F6514D344F5,
	          0xD96A82DA012EC67F, 0xA327BE8227BEEC6A) },
	{ GFw64be(0x000000E7D7C676CC, 0x253C7C6ECB6FBFD1,
	          0xD124BBE5C6071811, 0x8475FC4C71DC69CF),
	  GFw64be(0x0000011EA1AC622B, 0x8D82E4118E8C90AF,
	          0xE7943759449969A9, 0xE5586D276AD7F31D) },
	{ GFw64be(0x000000A2BC7F90A6, 0xA87E2A363AABA589,
	          0xBD1C3C473713B171, 0x0463C012605D96BC),
	  GFw64be(0x00000066D641C87D, 0x8DDFCE1C2A769051,
	          0x30A1C1B03CE67F84, 0x3AFA09388AD3D57F) },
	{ GFw64be(0x0000004A08E8B6B2, 0xBA41C21D92C4CCEA,
	          0x96133360C7827F98, 0xE8B8813B09DBE393),
	  GFw64be(0x0000010D3AD67709, 0x8117069EDA52B373,
	          0xF69F99D9706F63E4, 0x7A7FB13BC1464AB8) },
	{ GFw64be(0x000000A5347F39FE, 0x8935DF228EFD97B9,
	          0x9D9FBC6427F9D8A8, 0xF2A38A80B50FB783),
	  GFw64be(0x000001E146094F77, 0x4C18D97F1C19F0CC,
	          0xFDF85640BBE94304, 0x07BE4F485A4F969C) },
	{ GFw64be(0x000001C603320DE3, 0x58D7F0716F74D4E1,
	          0x24DA8F84F97101F2, 0x7C8AAF07FBB63162),
	  GFw64be(0x0000007B6EE81109, 0xAD18BB4CA059F802,
	          0x9846017C79F462DF, 0xD4EE49F5E960E7A2) },
	{ GFw64be(0x000000FD194D6786, 0x0A89F55B751D37CC,
	          0xE9D0B8ED38F754F7, 0xFA525685F6A37DBA),
	  GFw64be(0x0000010D57B02E34, 0xF99D6FFBD7E55AC4,
	          0x3BF532C90CE8F86E, 0x2C7F55ADEAFF4B86) }
};

/* tau_i(frob_60(G))  for i = 1 to 21 */
static const inner_point_affine XSK233_TAU_G60[] = {
	{ GFw64be(0x000001A4BEE62594, 0xFD3E12B6608DF4F1,
	          0x8F2488B389E1925C, 0xC46D1823E9376902),
	  GFw64be(0x0000014FEC69F2AD, 0x641F14A4D95136A3,
	          0x35B06F74FDEDA510, 0xD948BAD4AD289199) },
	{ GFw64be(0x00000092E2C0E587, 0x484BD689E12B2931,
	          0x45616C211BA50961, 0x6CAED96B80E53BC7),
	  GFw64be(0x0000012A16DF3FB0, 0x77DB5E53CE1D6A57,
	          0xA0A9416D135F6888, 0x9187DE7602D6CDA7) },
	{ GFw64be(0x000001AB2C7C7E4D, 0x230D1B559FD8B14B,
	          0xD9A99A8048F0B967, 0xF3533D6A85EFC6B9),
	  GFw64be(0x000000E74893357D, 0xAA9F7BFA5BD52911,
	          0xBC202DCBE18A6810, 0x5053DDAE1EE92001) },
	{ GFw64be(0x0000000BB6FA26A3, 0xA347CC33A2CA3892,
	          0x6E6B8477018C3C86, 0x66D670C3859949DF),
	  GFw64be(0x000000C3182B3CF3, 0x892F99FD34E6D95A,
	          0xCF508D634EA9A210, 0xEC8F674D8FF61172) },
	{ GFw64be(0x000000EDA8A002FD, 0xE16AC5C3C1AD66FD,
	          0xEA2C4B646455B1FA, 0x27EE2E1DA629EA9F),
	  GFw64be(0x000001671E893CC4, 0xF0D58DDEB799E251,
	          0xD49E6F37E61EDDB9, 0xAEB86AE33B68D344) },
	{ GFw64be(0x000000416BCC68A8, 0xBAC0DD2AEF6CA501,
	          0xE6C7C6477A2377A8, 0xA53784C822A9B1F3),
	  GFw64be(0x000000888651F2CF, 0xFC0962469E69D156,
	          0xF209B5B96DBDBAA8, 0x9CC6348081A7B916) },
	{ GFw64be(0x0000000B7E4FCDF2, 0x08B31D19B2AC6F3B,
	          0xD860DF7957BBAA6E, 0xAE6960540DE0490B),
	  GFw64be(0x000000A1967BFCEF, 0x6580F73DAC9763CB,
	          0x58203510061F4C45, 0x9A7B42ABA6F5800D) },
	{ GFw64be(0x0000004FEAB8BF15, 0x02A340788FFAE301,
	          0x006DAE0EAF5B6A2D, 0x2756F534D64F6F1B),
	  GFw64be(0x00000122C253BCAD, 0xB05446C9CC85C9BE,
	          0xD972D363097CEC2B, 0x4F7CB598E62C6FD9) },
	{ GFw64be(0x000000D8BD181D37, 0xA08B382E191551B2,
	          0xD0D521FE73173703, 0x4C9FDA5C7D961BE2),
	  GFw64be(0x00000114FC8E76D8, 0x75F2B1B0F423B28A,
	          0x4D8533E749A9A666, 0x93E54814390DAF66) },
	{ GFw64be(0x000000D890E51CBA, 0xB432B3314DA88C29,
	          0x1F14E65EAE655DA7, 0x3DDA4E2FFCEE6702),
	  GFw64be(0x000001FEB4D50D17, 0xF4340BD4D35517BF,
	          0x61AD4EEE09071FB1, 0xB0EE420890A8C860) },
	{ GFw64be(0x000000CCE5B8029F, 0xB3D1846786295655,
	          0x6AC18D72B4E3D625, 0x940DFDB090C7CB19),
	  GFw64be(0x0000014AC264C7F1, 0xAEBE376C93D1AF78,
	          0x70280387038709AF, 0xE208F6D990CE8D61) },
	{ GFw64be(0x0000009DDA1692BD, 0x9FEC2CAF9D9D8ECC,
	          0x9D0EECD59B13CADA, 0x2D32099115BF1F69),
	  GFw64be(0x000000CBCF914FCB, 0xBCF3E751ED47C5C3,
	          0xC3883431F129AA4A, 0x20E9D08CCE3E8744) },
	{ GFw64be(0x000000159CF4088F, 0xC312505A77ABFC5D,
	          0xFA7944C2E47EA609, 0xAE1A869F71F86A6A),
	  GFw64be(0x000000EC84A8D458, 0x7597A3F1C66FA8C0,
	          0xC73038AA142A07EA, 0xBDC4B2C296583CD7) },
	{ GFw64be(0x00000008735F3543, 0x391F6F6D4C4CBDC2,
	          0xBB2EF71E5767EF32, 0x6D4E2A92E7B9F11A),
	  GFw64be(0x000000882D930B0A, 0x0894A3DDB27B8A9E,
	          0x6D346A58CF673421, 0xCFDF43AC27DDF007) },
	{ GFw64be(0x0000001907E51A25, 0x512FB1FB4B96BA0A,
	          0x7C61FE972FDE95AB, 0x42F2CFF7057865DA),
	  GFw64be(0x000000D042586B8D, 0x6C7623599186EB43,
	          0xB3526E27222F8494, 0x45EA19707DD6E331) },
	{ GFw64be(0x0000007B6EF6AA7C, 0xC47791451C6EEDAE,
	          0xFF86041DE06503AA, 0xA66A56905490544E),
	  GFw64be(0x00000124D92D342D, 0xA2CB1DF0DC56147E,
	          0x78CFD1C5D27A5064, 0x58EC20E3BB38A18F) },
	{ GFw64be(0x0000013D4577F603, 0x5976A564785C3AF4,
	          0x9D645676C7362D7F, 0xC4D1EB4B09557A16),
	  GFw64be(0x000001E45BCB5A8F, 0x3B9DDFA8AC0D2878,
	          0x6B1093362F49B39F, 0xF5E092683AA928F4) },
	{ GFw64be(0x000001132C237D5E, 0xBF85A13F8D344706,
	          0x05946E5FDC000581, 0xDD9C90B7C36C8529),
	  GFw64be(0x00000011A70D563D, 0xB04944694E3A90F5,
	          0x1C6E556BE4994C13, 0x17B6D231324828E5) },
	{ GFw64be(0x0000018C1EF58946, 0x1CA6B6055C0DA9CA,
	          0x15727A9429546AE4, 0x02EDC10237494096),
	  GFw64be(0x00000100D91C1B1E, 0x0CE16BED2B31E488,
	          0x74CA2D33049FD35F, 0x0FFED167FD762F9E) },
	{ GFw64be(0x000001B07E149B76, 0xC6DEB63193D00196,
	          0x03D28B00973738A6, 0x35317C022594FAC1),
	  GFw64be(0x000000D11056D65C, 0x884B8AB529F0BA70,
	          0xDB6BB812B16A0F9D, 0x9DC5DEB3EBCA45A3) },
	{ GFw64be(0x0000017E50D50F69, 0xD4C29F3A4D792634,
	          0x1D3B3366F432202D, 0x7E940373B21981FC),
	  GFw64be(0x000001193C064946, 0x4DEB6A4267534FF4,
	          0x7ABD5B72FFBF11A1, 0x0928F8AB867B8C45) }
};

/* tau_i(frob_120(G))  for i = 1 to 21 */
static const inner_point_affine XSK233_TAU_G120[] = {
	{ GFw64be(0x00000120298617E0, 0xA278562B92DAB053,
	          0xA5634AC70DEB22CF, 0xF874E871A9644A15),
	  GFw64be(0x0000001BC2BC743E, 0x9427DB7920D4C44C,
	          0x290401EB85B4C557, 0xA75C8A57F0D945A2) },
	{ GFw64be(0x00000087984E5A35, 0x2879D665862AE0D7,
	          0x6E4DBD48FA97EE2D, 0xE667D55886E71CCF),
	  GFw64be(0x000000B820895E4D, 0x4A91E790F2313E22,
	          0xAE8C863ABC04053C, 0x98785D0C68067481) },
	{ GFw64be(0x000001ABCF519F80, 0xBE0BC67C0FD25210,
	          0x6FF656F1A346C078, 0x1E76FF1761346EA0),
	  GFw64be(0x00000096F1E25920, 0xA96F0C30B6228D42,
	          0x644BE038D11D2188, 0x9F109200D719286E) },
	{ GFw64be(0x0000005347719ACA, 0xDD64CB355CFF2C86,
	          0x59BC7D6534092F02, 0x777F142DDFCB3E58),
	  GFw64be(0x000000D0C09E8766, 0xCDF2281882903872,
	          0xDB4316EC5C834217, 0xD4C6AA5656233E31) },
	{ GFw64be(0x000000094927793A, 0x3951021E061F5440,
	          0x450B045129378BB9, 0x8024CEF9860F9276),
	  GFw64be(0x0000005E252F203A, 0x02BEF9BC8D49C2AE,
	          0x06DF22838D15986E, 0x5F1B31ECD69DAF06) },
	{ GFw64be(0x000001B6B9147F03, 0x648D9014508D0C5D,
	          0x40A0591BD0AA82F9, 0xD359D201951F3BEA),
	  GFw64be(0x0000006556A8AD40, 0x53A9ABD38E0A6C93,
	          0xBF72591AAC02BA43, 0x4760289D025B5D58) },
	{ GFw64be(0x000000144E1FA950, 0x88C7B0D432CA45F6,
	          0xA6B3FFFAA2A37F66, 0xD41D2A85927D16C2),
	  GFw64be(0x0000011D313AB65F, 0x51C75A734D65D498,
	          0xC5403539868C8B7B, 0x6D4D3058124B3B72) },
	{ GFw64be(0x000001D0357B36BB, 0x253AA84B04DD4630,
	          0x55C9E3FAAE0B646A, 0xFEC76D700E1DF5B2),
	  GFw64be(0x00000005219EFE78, 0xB3F0480F9A25EFBE,
	          0x99935F4D46BA6365, 0x9146A44ED33CB29D) },
	{ GFw64be(0x0000011C2E8A4344, 0xDB6134603016DED0,
	          0xB6BA0BAE8BB7A12F, 0x3D28F11790FB86B8),
	  GFw64be(0x000000435F517628, 0xD453CDD52F4F868F,
	          0x35F55575E387BDC4, 0x195F061FAB80257F) },
	{ GFw64be(0x000000CD2A1A1381, 0x0EC3A73FE26FE667,
	          0xA3A436448F724B41, 0xAB8CC95954195C91),
	  GFw64be(0x000001D5AE0EC22F, 0xC0518BB363E273A4,
	          0x0B33381BBD09C8C9, 0x25B852773B510A39) },
	{ GFw64be(0x000001A260B54352, 0x1ED725640DA04143,
	          0x8E1E5518DB726F18, 0xAE5736A7A9D1DCAE),
	  GFw64be(0x00000045348FAFA7, 0x4EF4CA91447D5D1E,
	          0xE253DA3AF4F1E1D8, 0x77A1D8AD91ADE688) },
	{ GFw64be(0x000000A8B3432167, 0x5B0AECE4608BC992,
	          0x6C01B7A3AC580964, 0x7DD4685CBDD80166),
	  GFw64be(0x0000010E3BC303E4, 0xECF88A86C5E4C26B,
	          0x5C46A5A31D0E0345, 0x331F2D21D3B7EADA) },
	{ GFw64be(0x000001F00DDD8438, 0x6D0538A5B5AA4A85,
	          0xAC2FDC0059AC6BC8, 0xAAC1BB0218512BBD),
	  GFw64be(0x00000005C5D7A67A, 0xBB2B85C006A0F048,
	          0x095FEB56800146DD, 0x7D8B40B209406BB8) },
	{ GFw64be(0x000001055DDDDF6E, 0x64A6C6AD9F7D90F4,
	          0x3AC02E04156A0E38, 0x731C49FADE721B37),
	  GFw64be(0x0000002A8F9125E9, 0xC0BEE052CA6BB88B,
	          0xA8FD25CA33FA5CCB, 0x58A99960E128099D) },
	{ GFw64be(0x000001620AC23211, 0x69464C8DF9A020CF,
	          0x6704C227580F467B, 0xE45B7FB767F4328B),
	  GFw64be(0x000000A58A302119, 0xF734E46A7D0ED798,
	          0x4B35A43CF99207B0, 0xAFC344F17345345E) },
	{ GFw64be(0x00000063DC25D746, 0x6E56AAEDB6B7FA64,
	          0xC31A2F0900A2A879, 0xA8D91A46901CD9E9),
	  GFw64be(0x0000018593DDBAF1, 0x9ABEEFE436AF9EDD,
	          0x79A6BFEDB9626656, 0x48AB2EB3BB71B5E4) },
	{ GFw64be(0x000000C864C4A573, 0xDF54A3554D50C72C,
	          0x1DA5925C5A8D21EA, 0x2024FFA55CBF9868),
	  GFw64be(0x000001CFAA5C2A0D, 0x3735E4EE1FB16889,
	          0x9BF3C49F6A474D15, 0x9A6D5467EF90E061) },
	{ GFw64be(0x000001E40ACD4676, 0x686D6F9D462DA1B7,
	          0xB62A2F5A9F3F4814, 0x9B2CAB2D8426CA72),
	  GFw64be(0x00000199B1B99D93, 0x742D68356DD3B8A6,
	          0x4695FD53D6E493CD, 0x6BF41F312F8689FB) },
	{ GFw64be(0x00000028C9FB2908, 0x08702AB65D033098,
	          0x4781016A8E3AFB63, 0x8446859D5EC5D3E8),
	  GFw64be(0x0000019379AD5D66, 0x33B27671810CC040,
	          0x0C662E4324736899, 0x36B2A70994A84B69) },
	{ GFw64be(0x0000003A2D3E921A, 0x6A7F9D0EB8CC53C7,
	          0xFFDF67E876F49B1B, 0x6ABE4F4A5FF346E3),
	  GFw64be(0x000001072DE20947, 0xEFF3884952E7BC04,
	          0xBFEBEE93D3140EEE, 0x160BC7A35C1423C7) },
	{ GFw64be(0x00000190E73C466F, 0x449E8EC9F39DEE78,
	          0x0B68FC349DE064A8, 0x4731CBB3F5C18CF5),
	  GFw64be(0x000001D8A892F738, 0xF9BC5936A24F7DF8,
	          0x9254EBB63E833439, 0xEEF24A61DAF3E770) }
};

/* tau_i(frob_180(G))  for i = 1 to 21 */
static const inner_point_affine XSK233_TAU_G180[] = {
	{ GFw64be(0x000000DB107C602A, 0xFBDB481AD7D4C935,
	          0xB2660BFBA686AC34, 0x3CBFC155D59BC343),
	  GFw64be(0x00000064A0BA9D50, 0xFCD8BB3ED7A57946,
	          0x6006B2568216C8BF, 0x9942700ED834B5A7) },
	{ GFw64be(0x000000B4084F77C5, 0xCC9E429EC451408C,
	          0x035A3A372821A3E8, 0x9C3861D288B3920B),
	  GFw64be(0x0000001EE786B3B6, 0x68AEA3347240351E,
	          0xDBE087BF926CCEE2, 0x498EA5DD0B5A1504) },
	{ GFw64be(0x0000010C954013B6, 0x6FE30D2607AD81BB,
	          0x58CA70EC05598EE7, 0x9AB8B64C9BD9CCA8),
	  GFw64be(0x00000164D306CCF3, 0xB25E52560EF60731,
	          0xD2B547BD1D9DCF46, 0x2FA42227D649BA09) },
	{ GFw64be(0x000001442DE4A7B5, 0x0EC20C01EC8C4A22,
	          0xCA58C694EA763949, 0xD0C805353C11179D),
	  GFw64be(0x00000028C2BFC75D, 0xC32C9CD2F0565177,
	          0x6014402F7A846F58, 0x54138C84BF7CA7FF) },
	{ GFw64be(0x000000F3FC260AC0, 0x9ED7ABF0F2F66695,
	          0x71A80FA7A80BEF0E, 0x7343C3EA0BADBD99),
	  GFw64be(0x000001B276EF4783, 0x344F3BAA8DA55451,
	          0x6CC72438989E5B09, 0x0D656B10689E89E6) },
	{ GFw64be(0x000000EC158A74D2, 0x00B39BCBCA70F75C,
	          0xC782A7486981AFF2, 0x39A440F3F2D9E12D),
	  GFw64be(0x000001919235E579, 0x29714B53D25B6875,
	          0xEA63CC4B6D4DF32B, 0xCB451FC3D1BCD236) },
	{ GFw64be(0x000000743DCEF430, 0x24622A664CEEE143,
	          0x37FED6807B0CD4C3, 0xF86A5C5DB6489F64),
	  GFw64be(0x0000013578214010, 0xD2B7ECEB378B2BF0,
	          0x1032D198C39D0065, 0x3DE44EE6A63FBEBC) },
	{ GFw64be(0x000001EA58BEC9B2, 0x5CC61FB42565FA89,
	          0x62CD2009B0661821, 0x3D5C72571379761C),
	  GFw64be(0x000001B090220ED7, 0x3D4CE89016FD20EB,
	          0x5B239652AA2E841C, 0x8F0D919F9BCF4044) },
	{ GFw64be(0x00000045A9EC7BC2, 0x58DE5C0928431503,
	          0x1768A58B1988584C, 0x7A709271E70A1B66),
	  GFw64be(0x000000DE05FAAD71, 0x1FF36593A4EA5A98,
	          0x377316F0BF0420F2, 0xD6582C32D502EE32) },
	{ GFw64be(0x0000006A00574E3D, 0xC6CA886DFE7C5441,
	          0xE2CFA0DDFB3C1E3B, 0x8F1A1E75D9DB85B3),
	  GFw64be(0x000000BD8C3AAFC0, 0x6142C37C9360147B,
	          0x5CA2C5E55E9088E7, 0x18DDDD0895DB697C) },
	{ GFw64be(0x0000017693CF3216, 0xB88C201341554D13,
	          0x211234A5BFC40413, 0x81E2774F1DF71D50),
	  GFw64be(0x000000CA8B2DCA09, 0xB33F0252656D1D8D,
	          0xFDE18B8D73FC318D, 0x942F0C95EC9EE2CA) },
	{ GFw64be(0x000000266C1512EA, 0xB4C168A366DC2E66,
	          0x298EC508C6C3B35B, 0x852FCA89B83C1EA0),
	  GFw64be(0x00000027D872B845, 0xBC713AD1578A0240,
	          0xFC5D0CFE7B31AF4C, 0xF2CA615FC8FA360B) },
	{ GFw64be(0x000001BA8570BEF7, 0x8DB2A139424E757C,
	          0xE40AC87A2827F6BF, 0x38490B40B4A89E74),
	  GFw64be(0x0000018082B6DB86, 0x4030B2B081F67C8E,
	          0xA5A70739D34FD4F1, 0xF78AD5BE63FFC177) },
	{ GFw64be(0x0000015ED136C0A0, 0x37E7A07AF1B2A4E7,
	          0xEA0557154532856E, 0x7F3B419817828EBF),
	  GFw64be(0x000001867B836B68, 0xF82DE9DB8A0602FD,
	          0xD271C63E70D6D3B1, 0x778367FBC444C553) },
	{ GFw64be(0x00000132BF38FE88, 0x36F138E92688E31B,
	          0x5CCF5F3110A6B838, 0x2B217D45EB1DA342),
	  GFw64be(0x00000110668552C2, 0x330728329967AEC3,
	          0x640312081467F187, 0x77C0FD532C1C4505) },
	{ GFw64be(0x000000F304A8004B, 0xED82BE968B6FD51B,
	          0xE0D82E21F40BACCB, 0x4DBEF72461CE423D),
	  GFw64be(0x0000008DCB96B1A6, 0xEEC4AC7CC8BF177B,
	          0x647F98EBE2DA863A, 0xE986507B3688E6C3) },
	{ GFw64be(0x0000017DAB37CDBE, 0xFA685FCD0E007862,
	          0xAF7CB58035ABA8F8, 0xD940027B881A0016),
	  GFw64be(0x000001AB75CAAEC8, 0xEB0280F4143AF709,
	          0xFB4CDD146AEB2996, 0x672D5217FB583A03) },
	{ GFw64be(0x000000406E33C247, 0x0BEBE26A1BC8D037,
	          0xAC6BEC4A36B95D19, 0x47D1C99053D0DBCE),
	  GFw64be(0x0000012FA9BEF7A0, 0xCF55AA900CA11276,
	          0xF153CE719D9C717F, 0x2719B371B012DE9B) },
	{ GFw64be(0x0000013F9BBF9491, 0xF5F6FE4DB4839172,
	          0xF4E6DA8533916527, 0x7E109573FCA91941),
	  GFw64be(0x000000816051AC8B, 0x51F7132AAA5FE811,
	          0xC854A118492246B4, 0xF5423A57C853527B) },
	{ GFw64be(0x000000FDC628F1D1, 0x5F4D07F0035474DD,
	          0x40554518C0A60B0F, 0xBD49D95C55AB1328),
	  GFw64be(0x00000026D0B35C1B, 0x137C4108C24169D8,
	          0x00086C40CA844A8F, 0x1051894AF28DE822) },
	{ GFw64be(0x00000090A9465985, 0xA64E59E4AE616EB7,
	          0x88FDA8C1618CB063, 0x07454251AC8BDCED),
	  GFw64be(0x000001B1A5180D4B, 0xE90D7F691FBDCB9E,
	          0x54FCB51F933D502D, 0x960399B0FE009ED3) }
};
