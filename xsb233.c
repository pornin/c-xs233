#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include "xs233.h"

#include "xs233_common.h"

/* The b constant for xsb233. */
static const gf GF_B = GFw64be(
	0x00000187F85627B9, 0x7874E747EE31E06D,
	0x71CAAEEA52F21253, 0xE5F946D061DA9138);

/* sqrt(b) for xsb233. */
static const gf GF_SB = GFw64be(
	0x00000007D5EF4389, 0xDFF11ECDBA39C309,
	0x70D3CE35CEBBA584, 0x73F64B4DC0F2686C);

/* 1/sqrt(b) for xsb233. */
static const gf GF_INV_SB = GFw64be(
	0x000001A66ECA2D9C, 0x9A99F70ABB81FA42,
	0xA7802D45255C81CC, 0xB6EE66F66C9A7586);

/* b^2 + 1 for xsb233. */
static const gf GF_BBp1 = GFw64be(
	0x00000066647EDE6C, 0x332C7F8C0923BB58,
	0x213B333B20E9CE42, 0x81FE115F7D8F90AC);

/* see xs233.h */
const xsb233_point xsb233_neutral = { {
	/* Neutral is N = (0, b) in (x, s) coordinates; since
	   s = sqrt(b)*S/Z^2 and Z = 1, we must set S to sqrt(b). */
	GFext(0, 0, 0, 0),
	GFext(0x00000007D5EF4389, 0xDFF11ECDBA39C309,
	      0x70D3CE35CEBBA584, 0x73F64B4DC0F2686C),
	GFext(0, 0, 0, 1),
	GFext(0, 0, 0, 0)
} };

/* see xs233.h */
const xsb233_point xsb233_generator = { {
	/* G is obtained from mapping the B-233 generator with the
	   y -> y + b change of variable. */
	GFext(0x0000017512FF306E, 0xA01C3744FBE1B95C,
	      0x002B20344E602BB0, 0x4E1F03DEF0CC6851),
	GFext(0x00000184608AAFB7, 0x3AEF057D75B25F32,
	      0x1C37AE41375A7126, 0xEF3F941AE66B7661),
	GFext(0, 0, 0, 1),
	GFext(0x0000017512FF306E, 0xA01C3744FBE1B95C,
	      0x002B20344E602BB0, 0x4E1F03DEF0CC6851)
} };

/*
 * (2^180)*G
 */
static const inner_point XSB233_T180_G = {
	GFw64be(0x000000673F17BED8, 0x9F71A5C5DBBAF08A,
	        0xC0B6B19F1916E75A, 0x1BF81499AA2B8734),
	GFw64be(0x000001A45AD8B438, 0xBD985C847030444E,
	        0x5EC9C240E30D1D23, 0x7D95F02930922F27),
	GFw64be(0, 0, 0, 1),
	GFw64be(0x000000673F17BED8, 0x9F71A5C5DBBAF08A,
	        0xC0B6B19F1916E75A, 0x1BF81499AA2B8734)
};

static const inner_point xsb233_inner_neutral = {
	/* Neutral is N = (0, b) in (x, s) coordinates; since
	   s = sqrt(b)*S/Z^2 and Z = 1, we must set S to sqrt(b). */
	GFw64be(0, 0, 0, 0),
	GFw64be(0x00000007D5EF4389, 0xDFF11ECDBA39C309,
	        0x70D3CE35CEBBA584, 0x73F64B4DC0F2686C),
	GFw64be(0, 0, 0, 1),
	GFw64be(0, 0, 0, 0)
};

static const inner_point_affine XSB233_G0[];
static const inner_point_affine XSB233_G60[];
static const inner_point_affine XSB233_G120[];
static const inner_point_affine XSB233_G180[];

/* see xs233.h */
uint32_t
xsb233_decode(xsb233_point *p, const void *src)
{
	/*
	 * Procedure:
	 *   decode w from bytes; reject if invalid (non-canonical).
	 *   if w == 0 then return N
	 *   d <- w^2 + w + a        (note: d is always non-zero)
	 *   e <- b/(d^2)
	 *   if trace(e) == 1 then reject
	 *   f <- halftrace(e)       (we have f^2 + f = e)
	 *   x <- d*f
	 *   if trace(x) == 1 then:
	 *       x <- x + d
	 *   s <- x*w^2
	 */
	gf w, d, e, f, x, s;
	uint64_t ve, rp, wz, r;

	/* decode w; remember if it was zero */
	ve = gf_decode(&w, src);
	wz = gf_iszero(&w);

	/* d <- w^2 + w + a */
	gf_sqr(&d, &w);
	gf_add(&d, &d, &w);
	gf_add(&d, &d, &GF_ONE);

	/* e <- b/(d^2) */
	gf_sqr(&e, &d);
	gf_div(&e, &GF_B, &e);

	/* input is not valid if trace(e) == 1 */
	rp = gf_trace(&e) - 1;

	/* solve f^2 + f = e */
	gf_halftrace(&f, &e);

	/* x <- d*f or d*(f + 1)  (whichever makes trace(x) == 0) */
	gf_mul(&x, &d, &f);
	gf_condadd(&x, &d, -gf_trace(&x));

	/* s <- x*w^2 */
	gf_sqr(&s, &w);
	gf_mul(&s, &s, &x);

	/*
	 * ve = -1 if the field element was canonical
	 * wz = -1 if w == 0
	 * rp = -1 if trace(e) == 0
	 * If ve == -1 then wz == -1.
	 * Decoding succeeds if ve == -1 and either wz == -1 or rp == -1.
	 * We force the point to N on failure, or if wz == -1.
	 */
	r = ve & (rp | wz);
	gf_condset(&x, &GF_ZERO, wz | ~rp);
	gf_condset(&s, &GF_B, wz | ~rp);

	/*
	 * In extended coordinates:
	 *   x = sqrt(b)*X/Z
	 *   s = sqrt(b)*S/Z^2
	 *   T = X*Z
	 * We set Z = 1, hence X = x/sqrt(b), S = s/sqrt(b), T = X
	 */
	gf_mul(&x, &x, &GF_INV_SB);
	gf_mul(&s, &s, &GF_INV_SB);

	point_store_x(p, &x);
	point_store_s(p, &s);
	point_store_z(p, &GF_ONE);
	point_store_t(p, &x);
	return (uint32_t)r;
}

/* see xs233.h */
void
xsb233_encode(void *dst, const xsb233_point *p)
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
xsb233_is_neutral(const xsb233_point *p)
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
xsb233_equals(const xsb233_point *p1, const xsb233_point *p2)
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
 * Inner addition function on xsb233.
 */
static inline void
xsb233_inner_add(inner_point *p3, const inner_point *p1, const inner_point *p2)
{
	/*
	 * x1x2 <- X1*X2
	 * s1s2 <- S1*S2
	 * z1z2 <- Z1*Z2
	 * t1t2 <- T1*T2
	 * d <- (S1 + T1)*(S2 + T2)
	 * e <- (a^2)*t1t2              note: a = 1 for xsb233
	 * f <- x1x2^2
	 * g <- z1z2^2
	 * X3 <- d + s1s2
	 * S3 <- sqrt(b)*(g*(s1s2 + e) + f*(d + e))
	 * Z3 <- sqrt(b)*(f + g)
	 * T3 <- X3*Z3
	 */
	gf x1x2, s1s2, z1z2, t1t2;
	gf d, f, g, tmp1, tmp2;

	gf_mul(&x1x2, &p1->x, &p2->x);
	gf_mul(&s1s2, &p1->s, &p2->s);
	gf_mul(&z1z2, &p1->z, &p2->z);
	gf_mul(&t1t2, &p1->t, &p2->t);
	gf_add(&tmp1, &p1->s, &p1->t);
	gf_add(&tmp2, &p2->s, &p2->t);
	gf_mul(&d, &tmp1, &tmp2);
	/* a == 1, therefore e == t1t2 */
	gf_sqr(&f, &x1x2);
	gf_sqr(&g, &z1z2);
	gf_add(&p3->x, &d, &s1s2);
	gf_add(&tmp1, &s1s2, &t1t2);
	gf_mul(&tmp1, &tmp1, &g);
	gf_add(&tmp2, &d, &t1t2);
	gf_mul(&tmp2, &tmp2, &f);
	gf_add(&tmp1, &tmp1, &tmp2);
	gf_mul(&p3->s, &tmp1, &GF_SB);
	gf_add(&tmp2, &f, &g);
	gf_mul(&p3->z, &tmp2, &GF_SB);
	gf_mul(&p3->t, &p3->x, &p3->z);
}

/* see xs233.h */
void
xsb233_add(xsb233_point *p3, const xsb233_point *p1, const xsb233_point *p2)
{
	inner_point ip1, ip2, ip3;

	point_load(&ip1, p1);
	point_load(&ip2, p2);
	xsb233_inner_add(&ip3, &ip1, &ip2);
	point_store(p3, &ip3);
}

/*
 * Mixed addition: second operand is in affine coordinates (Z2 == 1, T2 == X2).
 */
static inline void
xsb233_add_mixed(inner_point *p3, const inner_point *p1,
	const inner_point_affine *p2)
{
	/*
	 * When P2 is in affine coordinates, Z2 == 1 and T2 == X2, so
	 * the values Z2 and T2 are implicit, and the multiplication
	 * Z1*Z2 simplifies to using just Z1.
	 */
	gf x1x2, s1s2, t1t2;
	gf d, f, g, tmp1, tmp2;

	gf_mul(&x1x2, &p1->x, &p2->x);
	gf_mul(&s1s2, &p1->s, &p2->s);
	gf_mul(&t1t2, &p1->t, &p2->x);
	gf_add(&tmp1, &p1->s, &p1->t);
	gf_add(&tmp2, &p2->s, &p2->x);
	gf_mul(&d, &tmp1, &tmp2);
	/* a == 1, therefore e == t1t2 */
	gf_sqr(&f, &x1x2);
	gf_sqr(&g, &p1->z);
	gf_add(&p3->x, &d, &s1s2);
	gf_add(&tmp1, &s1s2, &t1t2);
	gf_mul(&tmp1, &tmp1, &g);
	gf_add(&tmp2, &d, &t1t2);
	gf_mul(&tmp2, &tmp2, &f);
	gf_add(&tmp1, &tmp1, &tmp2);
	gf_mul(&p3->s, &tmp1, &GF_SB);
	gf_add(&tmp2, &f, &g);
	gf_mul(&p3->z, &tmp2, &GF_SB);
	gf_mul(&p3->t, &p3->x, &p3->z);
}

/*
 * Inner doubling function on xsb233.
 */
static inline void
xsb233_inner_double(inner_point *p3, const inner_point *p1)
{
	/*
	 * We use here formulas with cost 3M+5S+2mb, and not the nominally
	 * more efficient 3M+4S+4mb, because on xsb233, the b constant is
	 * not simple and each multiplication by b or sqrt(b) has about the
	 * same cost as a generic multiplication.
	 *
	 *   xx <- X^2
	 *   zz <- Z^2
	 *   X' <- T^2
	 *   S' <- sqrt(b)*(xx*S + zz*(S + T))^2
	 *   Z' <- sqrt(b)*(xx + zz)^2
	 *   T' <- X'*Z'
	 */
	gf xx, zz, tmp1, tmp2;

	gf_sqr(&xx, &p1->x);
	gf_sqr(&zz, &p1->z);
	gf_sqr(&p3->x, &p1->t);
	gf_mul(&tmp1, &xx, &p1->s);
	gf_add(&tmp2, &p1->s, &p1->t);
	gf_mul(&tmp2, &zz, &tmp2);
	gf_add(&tmp1, &tmp1, &tmp2);
	gf_sqr(&tmp1, &tmp1);
	gf_mul(&p3->s, &tmp1, &GF_SB);
	gf_add(&tmp2, &xx, &zz);
	gf_sqr(&tmp2, &tmp2);
	gf_mul(&p3->z, &tmp2, &GF_SB);
	gf_mul(&p3->t, &p3->x, &p3->z);
}

/* see xs233.h */
void
xsb233_double(xsb233_point *p3, const xsb233_point *p1)
{
	inner_point ip1, ip3;

	point_load(&ip1, p1);
	xsb233_inner_double(&ip3, &ip1);
	point_store(p3, &ip3);
}

/*
 * Inner functions for sequences of doublings on xsb233; it should be
 * faster than calling xsb233_inner_double() n times, if n >= 2.
 * This function ASSUMES that n >= 1.
 */
static void
xsb233_inner_xdouble(inner_point *p3, const inner_point *p1, unsigned n)
{
	/*
	 * When n >= 2, we optimize by switching to lambda coordinates:
	 *   P+N (x,s) -> P (x,l):      1S + 1mb
	 *   P (x,l) -> 2*P (x,l):      3M + 4S + 1mb
	 *   P (x,l) -> 2*P+N (x,s):    3M + 3S + 3mb
	 * Total for n doublings:
	 *   using plain doubling:   (3*n)*M + (5*n)*S + (2*n)*mb
	 *   using internal lambda:  (3*n)*M + (4*n)*S + (n+3)*mb
	 * We thus get a gain when n*S > (3-n)*mb; this is always true
	 * for n >= 3, and usually true for n == 2.
	 */
	gf X, Z, L, D, E, F, U, t1, t2;

	/*
	 * Convert P1 -> P1+N in (x,lambda) coordinates.
	 *
	 *   X' <- sqrt(b)*Z^2
	 *   Z' <- T
	 *   L' <- S + (a + 1)*T   (a == 1 on xsb233)
	 *
	 * If P1 is the neutral, then we get (X:L:Z) = (X:0:X) with X != 0.
	 */
	gf_sqr(&X, &p1->z);
	gf_mul(&X, &X, &GF_SB);
	Z = p1->t;
	L = p1->s;

	/*
	 * Perform n-1 doublings in (x,lambda) coordinates.
	 * Formulas from: https://eprint.iacr.org/2013/131
	 *
	 *   D <- Z^2
	 *   E <- (L + X)^2
	 *   F <- L*(L + Z)
	 *   U <- F + a*D     (a == 1 on xsb233)
	 *   X' <- U^2
	 *   Z' <- D*U
	 *   L' <- E*(E + U + D) + (a^2 + b^2)*D^2 + X' + (a + 1)*Z'
	 */
	while (n -- > 1) {
		gf_sqr(&D, &Z);
		gf_add(&F, &L, &Z);
		gf_mul(&F, &L, &F);
		gf_add(&E, &L, &X);
		gf_sqr(&E, &E);
		gf_add(&U, &F, &D);
		gf_sqr(&X, &U);
		gf_mul(&Z, &D, &U);
		gf_add(&t1, &E, &D);
		gf_add(&t1, &t1, &U);
		gf_sqr(&t2, &D);
		gf_mul(&t1, &t1, &E);
		gf_mul(&t2, &t2, &GF_BBp1);
		gf_add(&L, &t1, &X);
		gf_add(&L, &L, &t2);
	}

	/*
	 * Final doubling with conversion back to (x,s).
	 *
	 *   D <- Z^2
	 *   E <- (L + X)^2
	 *   F <- L*(L + Z)
	 *   X' <- sqrt(b)*D
	 *   Z' <- F + a*D         (a == 1 on xsb233)
	 *   T' <- X'*Z'
	 *   S' <- sqrt(b)*(E*(E + Z' + D) + (F + sqrt(b)*X')^2)
	 */
	gf_sqr(&D, &Z);
	gf_add(&F, &L, &Z);
	gf_mul(&F, &L, &F);
	gf_add(&E, &L, &X);
	gf_sqr(&E, &E);
	gf_mul(&p3->x, &D, &GF_SB);
	gf_add(&p3->z, &F, &D);
	gf_add(&t1, &E, &D);
	gf_add(&t1, &t1, &p3->z);
	gf_mul(&t1, &t1, &E);
	gf_mul(&t2, &p3->x, &GF_SB);
	gf_mul(&p3->t, &p3->x, &p3->z);
	gf_add(&t2, &t2, &F);
	gf_sqr(&t2, &t2);
	gf_add(&t1, &t1, &t2);
	gf_mul(&p3->s, &t1, &GF_SB);
}


/* see xs233.h */
void
xsb233_xdouble(xsb233_point *p3, const xsb233_point *p1, unsigned n)
{
	inner_point ip;

	if (n == 0) {
		*p3 = *p1;
		return;
	}
	if (n == 1) {
		xsb233_double(p3, p1);
		return;
	}

	point_load(&ip, p1);
	xsb233_inner_xdouble(&ip, &ip, n);
	point_store(p3, &ip);
}

/* see xs233.h */
void
xsb233_neg(xsb233_point *p3, const xsb233_point *p1)
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
xsb233_sub(xsb233_point *p3, const xsb233_point *p1, const xsb233_point *p2)
{
	inner_point ip1, ip2;

	point_load(&ip1, p1);
	point_load(&ip2, p2);
	gf_add(&ip2.s, &ip2.s, &ip2.t);
	xsb233_inner_add(&ip1, &ip1, &ip2);
	point_store(p3, &ip1);
}

/* see xs233.h */
void
xsb233_select(xsb233_point *p3,
	const xsb233_point *p0, const xsb233_point *p1, uint32_t ctl)
{
	point_select(p3, p0, p1, ctl);
}

/* see xs233.h */
void
xsb233_condneg(xsb233_point *p3, const xsb233_point *p1, uint32_t ctl)
{
	inner_point ip;

	point_load(&ip, p1);
	gf_condadd(&ip.s, &ip.t, -(uint64_t)(ctl >> 31));
	point_store(p3, &ip);
}

/* see xs233.h */
void
xsb233_mul_ladder(xsb233_point *p3,
	const xsb233_point *p1, const void *n, size_t n_len)
{
	/*
	 * We use a Montgomery ladder as exposed in the CHES'99
	 * paper from López and Dahab: "Fast Multiplication on Elliptic
	 * Curves over GF(2^m) without Precomputation"
	 * https://link.springer.com/content/pdf/10.1007/3-540-48059-5_27.pdf
	 */

	/*
	 * If the scalar is empty then return the neutral.
	 */
	if (n_len == 0) {
		*p3 = xsb233_neutral;
		return;
	}

	inner_point ip;

	point_load(&ip, p1);

	/*
	 * We assume that P != N (a corrective step is applied at the end).
	 *
	 * Convert P -> P0 = P+N on original curve (x,y affine coordinates)
	 *
	 *   x0 = b/x
	 *      = sqrt(b)*Z/X
	 *   y0 = b*(s + x**2 + (a + 1)*x + b)/x**2 + b
	 *      = sqrt(b)*(S + (a + 1)*T + sqrt(b)*(X + Z)^2)/X^2 + b
	 *
	 * A single inversion of X is needed (since we assumed that P != N,
	 * we have X != 0). On xsb233, we have a == 1, hence a + 1 == 0.
	 */
	gf x0, y0, iX;

	gf_inv(&iX, &ip.x);
	gf_mul(&x0, &ip.z, &iX);
	gf_mul(&x0, &x0, &GF_SB);
	gf_add(&y0, &ip.x, &ip.z);
	gf_sqr(&y0, &y0);
	gf_sqr(&iX, &iX);
	gf_mul(&y0, &y0, &GF_SB);
	gf_add(&y0, &y0, &ip.s);
	gf_mul(&y0, &y0, &iX);
	gf_mul(&y0, &y0, &GF_SB);
	gf_add(&y0, &y0, &GF_B);

	/*
	 * We represent x coordinates of points in projective format:
	 *   x = X/Z
	 * For the point-at-infinity, we have Z == 0 and X != 0.
	 * For all other points, we have Z != 0.
	 *
	 * We apply the ladder.
	 * At each step, we have points P1 and P2 (x coordinates only)
	 * such that:
	 *   P1 = k*P0
	 *   P2 = (k + 1)*P0
	 * then we compute P3 = P1 + P2, and double either P1 or P2,
	 * thus doing the following:
	 *   (P1, P2) <- (2*P1, P3)     equivalent to: k <- 2*k
	 *   (P1, P2) <- (P3, 2*P2)     equivalent to: k <- 2*k + 1
	 * We thus explore the scalar bits from top to bottom.
	 * Initial values:
	 *   P1 = inf
	 *   P2 = P0
	 *
	 * Formulas (assuming we double P1 into P4 = 2*P1):
	 *   X3 = X1*X2*Z1*Z2 + x0*(X1*Z2 + X2*Z1)^2
	 *   Z3 = (X1*Z2 + X2*Z1)^2
	 *   X4 = (X1 + sqrt(b)*Z1)^4
	 *   Z4 = (X1*Z1)^2
	 */
	gf X1, Z1, X2, Z2;
	const uint8_t *nbuf;

	X1 = GF_ONE;
	Z1 = GF_ZERO;
	X2 = x0;
	Z2 = GF_ONE;
	nbuf = n;
	for (int i = (int)(n_len << 3) - 1; i >= 0; i --) {
		uint64_t bit;
		gf X3, Z3, X4, Z4, x1z2, x2z1, xt, zt, tt;

		/* Extract next scalar bit. */
		bit = -(uint64_t)((nbuf[i >> 3] >> (i & 7)) & 1);

		/* P3 = P1 + P2 */
		gf_mul(&x1z2, &X1, &Z2);
		gf_mul(&x2z1, &X2, &Z1);
		gf_add(&Z3, &x1z2, &x2z1);
		gf_sqr(&Z3, &Z3);
		gf_mul(&X3, &x1z2, &x2z1);
		gf_mul(&tt, &x0, &Z3);
		gf_add(&X3, &X3, &tt);

		/* P4 = 2*P1 (if bit == 0) or 2*P2 (if bit == 1) */
		gf_select(&xt, &X1, &X2, bit);
		gf_select(&zt, &Z1, &Z2, bit);
		gf_mul(&tt, &zt, &GF_SB);
		gf_mul(&Z4, &xt, &zt);
		gf_add(&tt, &xt, &tt);
		gf_sqr(&tt, &tt);
		gf_sqr(&Z4, &Z4);
		gf_sqr(&X4, &tt);

		/* (P1, P2) <- (P4, P3) (if bit == 0) or (P3, P4) */
		gf_select(&X1, &X4, &X3, bit);
		gf_select(&Z1, &Z4, &Z3, bit);
		gf_select(&X2, &X3, &X4, bit);
		gf_select(&Z2, &Z3, &Z4, bit);
	}

	/*
	 * We have the x coordinates of n*P and (n + 1)*P in X1/Z1
	 * and X2/Z2, respectively. We rebuild the y coordinate of P1:
	 *   y1 = (x1 + x0)*((x1 + x0)*(x2 + x0) + x0^2 + y0)/x0 + y0
	 *
	 * In projective coordinates, and remembering that x0 was computed
	 * as x0 = sqrt(b)*Z/X (for the original input point P):
	 *
	 *   Zr = sqrt(b)*Z*(Z1^2)*Z2
	 *
	 *   Yr = X*(X1 + x0*Z1)*((X1 + x0*Z1)*(X2 + x0*Z2) + (x0^2 + y0)*Z1*Z2)
	 *      + y0*Zr
	 *
	 *   y1 = Yr/Zr
	 *
	 * Once we have the result as P1, we need to convert it back to
	 * the group:
	 *   - add b to y1 to switch to the curve y^2 + x*y = x*(x^2 + a*x + b)
	 *   - add N
	 *   - convert to (x,s) coordinates
	 * Yielding:
	 *
	 *   x' = b/x1
	 *
	 *   s' = b*(y1 + x1^2 + (a + 1)*x)/x1^2
	 *
	 * Putting things together:
	 *
	 *   Z' = X1*Z*Z2
	 *
	 *   X' = sqrt(b)*Z*Z1*Z2
	 *
	 *   S' = Z*Z2*(Yr + sb*Z'*(X1 + (a + 1)*Z1))
	 *
	 * Note that since the initial point P was on the group and we
	 * first converted it by adding N (and moving to the original curve),
	 * then we worked on the r-torsion subgroup of the curve, which does
	 * not contain N. N is the only point whose x coordinate is 0.
	 * Therefore, neither P1 nor P2 may be N, which means that X1 != 0.
	 *
	 * However, it is possible that P1 or P2 is the point-at-infinity,
	 * in which case Z1 == 0 or Z2 == 0, respectively. P1 and P2 cannot
	 * be _both_ the point-at-infinity, under the assumption that the
	 * original source point P was not the neutral in the group.
	 *
	 *   If Z1 == 0 then P2 == P0, thus X2 + x0*Z2 == 0, and we get:
	 *       X' == 0
	 *       Z' == X1*Z*Z2 != 0
	 *       S' == sqrt(b)*Z'^2
	 *   which is a valid representation of the point N in the group.
	 *
	 *   If Z2 == 0 then X', Z' and S' are all equal to zero, which
	 *   is not valid and must be fixed. That situation happens when
	 *   the scalar is n = -1 mod r, in which case the result should
	 *   be the point -P.
	 *
	 * We also have to fix the result if the source point P was in
	 * fact the neutral.
	 */
	uint64_t inf2 = gf_iszero(&Z2);
	uint64_t srcn = gf_iszero(&ip.x);
	gf zz2, X3, Z3, S3, Yr, Zr, v1, v2;

	/* Z3 = X1*Z*Z2 */
	gf_mul(&zz2, &ip.z, &Z2);
	gf_mul(&Z3, &X1, &zz2);

	/* X3 = sqrt(b)*Z*Z1*Z2 */
	gf_mul(&X3, &Z1, &zz2);
	gf_mul(&X3, &X3, &GF_SB);

	/* Zr = sqrt(b)*Z*(Z1^2)*Z2 */
	gf_mul(&Zr, &X3, &Z1);

	/* Yr = X*(X1 + x0*Z1)*((X1 + x0*Z1)*(X2 + x0*Z2) + (x0^2 + y0)*Z1*Z2)
	        + y0*Zr */
	gf_mul(&v1, &x0, &Z1);
	gf_add(&v1, &v1, &X1);
	gf_mul(&v2, &x0, &Z2);
	gf_add(&v2, &v2, &X2);
	gf_mul(&Yr, &v2, &v1);
	gf_mul(&v1, &v1, &ip.x);
	gf_sqr(&v2, &x0);
	gf_add(&v2, &v2, &y0);
	gf_mul(&v2, &v2, &Z1);
	gf_mul(&v2, &v2, &Z2);
	gf_add(&Yr, &Yr, &v2);
	gf_mul(&Yr, &Yr, &v1);
	gf_mul(&v2, &y0, &Zr);
	gf_add(&Yr, &Yr, &v2);

	/* S3 = Z*Z2*(Yr + sb*Z3*(X1 + (a + 1)*Z1))  (on xsb233: a + 1 == 0) */
	gf_mul(&v1, &Z3, &X1);
	gf_mul(&v1, &v1, &GF_SB);
	gf_add(&v1, &Yr, &v1);
	gf_mul(&S3, &v1, &zz2);

	/*
	 * If P2 == inf, set P3 to -P1.
	 */
	gf_add(&v2, &ip.s, &ip.t);
	gf_condset(&X3, &ip.x, inf2);
	gf_condset(&S3, &v2, inf2);
	gf_condset(&Z3, &ip.z, inf2);

	/*
	 * If the source was the neutral, set P3 to the neutral.
	 */
	gf_condset(&X3, &GF_ZERO, srcn);
	gf_condset(&S3, &GF_SB, srcn);
	gf_condset(&Z3, &GF_ONE, srcn);

	/*
	 * Return the result; also compute T3 = X3*Z3.
	 */
	ip.x = X3;
	ip.s = S3;
	ip.z = Z3;
	gf_mul(&ip.t, &X3, &Z3);
	point_store(p3, &ip);
}

/* see xs233.h */
void
xsb233_mul(xsb233_point *p3,
	const xsb233_point *p1, const void *n, size_t n_len)
{
	/*
	 * Compute the window: win[i] = (i+1)*P1.
	 */
	inner_point win[16];

	point_load(&win[0], p1);
	xsb233_inner_double(&win[1], &win[0]);
	for (int j = 4; j <= 16; j += 2) {
		xsb233_inner_add(&win[j - 2], &win[j - 3], &win[0]);
		xsb233_inner_double(&win[j - 1], &win[(j >> 1) - 1]);
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
		inner_select(&ip3, &xsb233_inner_neutral, &win[0], -td);
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
			inner_select(&ip3, &xsb233_inner_neutral, &win[0], -td);
		} else {
			ip3 = xsb233_inner_neutral;
		}
	}

	/*
	 * Process all signed digits in high-to-low order.
	 */
	while (sd_len -- > 0) {
		inner_point Q;

		xsb233_inner_xdouble(&ip3, &ip3, 5);
		inner_lookup(&Q, win, 16, sd[sd_len], &GF_SB);
		xsb233_inner_add(&ip3, &ip3, &Q);
	}

	point_store(p3, &ip3);
}

/* see xs233.h */
void
xsb233_mulgen(xsb233_point *p3, const void *n, size_t n_len)
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
		inner_select(&ip3, &xsb233_inner_neutral, &XSB233_T180_G, -td);
	} else {
		uint8_t tmp[30];

		memcpy(tmp, n, n_len);
		memset(tmp + n_len, 0, (sizeof tmp) - n_len);
		u240_recode5(sd, tmp);
		ip3 = xsb233_inner_neutral;
	}

	/*
	 * Process all signed digits in high-to-low order.
	 */
	for (int i = 11; i >= 0; i --) {
		inner_point_affine Qa;

		xsb233_inner_xdouble(&ip3, &ip3, 5);
		inner_lookup_affine(&Qa, XSB233_G0, 16, sd[i], &GF_SB);
		xsb233_add_mixed(&ip3, &ip3, &Qa);
		inner_lookup_affine(&Qa, XSB233_G60, 16, sd[i + 12], &GF_SB);
		xsb233_add_mixed(&ip3, &ip3, &Qa);
		inner_lookup_affine(&Qa, XSB233_G120, 16, sd[i + 24], &GF_SB);
		xsb233_add_mixed(&ip3, &ip3, &Qa);
		inner_lookup_affine(&Qa, XSB233_G180, 16, sd[i + 36], &GF_SB);
		xsb233_add_mixed(&ip3, &ip3, &Qa);
	}
	point_store(p3, &ip3);
}

/* ===================================================================== */

/* i*(2^0)*G  for i = 1 to 16 */
static const inner_point_affine XSB233_G0[] = {
	{ GFw64be(0x0000017512FF306E, 0xA01C3744FBE1B95C,
	          0x002B20344E602BB0, 0x4E1F03DEF0CC6851),
	  GFw64be(0x00000184608AAFB7, 0x3AEF057D75B25F32,
	          0x1C37AE41375A7126, 0xEF3F941AE66B7661) },
	{ GFw64be(0x000001BA5A1A98AB, 0xD43990DB73110A42,
	          0xF4362148DEF8D94D, 0xF8E8FF0C92CEC2DD),
	  GFw64be(0x00000051F548846C, 0xAEEFA79799752E4B,
	          0xD1E1B03D209F7BEC, 0x7EC487F109C3FA0A) },
	{ GFw64be(0x000001E8F4C28F9F, 0x8143129EC1CA11A3,
	          0x966B51CCD1E92105, 0xBAAE0694CB766225),
	  GFw64be(0x000000527E73A1C1, 0x81EF0D4F3E525175,
	          0x5C79AF224D0BC351, 0xAC60069A18F2A57C) },
	{ GFw64be(0x0000019D74FFF308, 0x7D7E303A258E0543,
	          0xD9123EEF1BAB7ED1, 0x1917E49E4C304CB6),
	  GFw64be(0x00000119D7515401, 0x24B4C5A95B6B7AC0,
	          0x833A5638E26A7AC5, 0x997F2041821DF3F7) },
	{ GFw64be(0x0000013F7B4BFE4D, 0xE28109784F5DC21C,
	          0x2E8675F91688EB4D, 0x7F7E155F54EAF3D0),
	  GFw64be(0x000001637FBD54BE, 0x9371F1CA1AC7BC2B,
	          0x0B00407E22992277, 0xEEEEA4563FD690D6) },
	{ GFw64be(0x000000A571BD296F, 0xE5D3E2ADC41408BC,
	          0xB4134A328E64DBD9, 0x4397BE717FEB59D8),
	  GFw64be(0x00000080ABFB40A5, 0x6DF56A43BD8D7FCE,
	          0xD1CAC1C96DD85D0C, 0xB9C5762E8BB230AF) },
	{ GFw64be(0x00000053F18B912E, 0x30A0A341790849BC,
	          0x6B62835E647D3BA1, 0xF28EDBE9FF2A5050),
	  GFw64be(0x0000010EF89F9CDC, 0x9884F580902E9BE4,
	          0x30B21053F384E52E, 0x426A388A2C24D9F4) },
	{ GFw64be(0x000001635B27A771, 0xF2B1F0C1BEB48706,
	          0x1766A98DBB405C07, 0x7563DE8E9F531ADD),
	  GFw64be(0x000000A9632525C6, 0xBECDD7099CA73844,
	          0x601654D904BFE330, 0x461795D09DCE7A61) },
	{ GFw64be(0x0000016BB51C0770, 0x0B295B3D4C8B945E,
	          0x2C41510335129F55, 0x9CFF428AD9A0D4B3),
	  GFw64be(0x0000007C28FDE3B9, 0xB56FEA66D119D737,
	          0x11EB26D9FDFAF657, 0xC74BD80EEF9E48D5) },
	{ GFw64be(0x000000CBB77E7B28, 0xE9F3809F53D7A2E7,
	          0xC44BF5538318BFE4, 0x036EC8476A029FA9),
	  GFw64be(0x0000009BEAD2A791, 0x5D869091669EA555,
	          0xDAEED326020D5336, 0xA74276FE8B133E58) },
	{ GFw64be(0x000000B6B46FCE37, 0x56EE06C9AF720EDE,
	          0x203E8D1FA072018C, 0xAF66C1C26180D415),
	  GFw64be(0x000000F2C612AC90, 0x4C65620C9E6E28D9,
	          0xCC9709DDB3DA3547, 0x8681771403E4BD4A) },
	{ GFw64be(0x0000010A1918C2B6, 0xFE9132CD7F4FEBC3,
	          0x4A489E2AC5F93CCD, 0x5D694EF0F578152F),
	  GFw64be(0x00000135B6D159DF, 0x886BAE1FE4D952C0,
	          0xFC760A0BCE02BEFA, 0x56A3CB0A03E449D8) },
	{ GFw64be(0x000001EB25BA84E1, 0x1FA1BB91450AC31C,
	          0x4473D2FD96BAA47A, 0xB057CB775B95C790),
	  GFw64be(0x0000019A52F8B5EC, 0x3A6161AFEA76618E,
	          0xA00D97BB0240CEE9, 0x289B2FFA6216E2D0) },
	{ GFw64be(0x000000632F1F04A5, 0xEDCDFF0B6FB0D0F7,
	          0x35718F252BFC95CB, 0xC59453A053D479CB),
	  GFw64be(0x000000294AB057AC, 0xABF002D8034DA119,
	          0x85E862FE0C169D85, 0x279962DE7385D12D) },
	{ GFw64be(0x0000001B76BDA641, 0xAE39DC6853CAFE19,
	          0x824D5476B2F4CC22, 0xDF637AA665E887B6),
	  GFw64be(0x000001B93743B6CC, 0xC0CBA6A9008BD16C,
	          0x4EF2F3C9B7639B84, 0x7D8DE72E9F82A688) },
	{ GFw64be(0x0000016D93A0620C, 0xC5B689638C6C31A5,
	          0xF878A56287663576, 0x67503E785FF01ABD),
	  GFw64be(0x0000004FB7645338, 0xE55A2F92CAE60CF1,
	          0xAB5DA3D096C86434, 0xD3EB025D7557B778) }
};

/* i*(2^60)*G  for i = 1 to 16 */
static const inner_point_affine XSB233_G60[] = {
	{ GFw64be(0x00000182173029E3, 0x3207F40988D25E0A,
	          0x1F986A23B121FC97, 0xBFD07CE9979969B1),
	  GFw64be(0x0000014D499E38F2, 0x2151213F9E22B17F,
	          0x1D0EFD4CA11A415F, 0x6BF702329D6F9EE3) },
	{ GFw64be(0x000000372B8ACC0A, 0x9A7E5822D74F941A,
	          0x442384166AE4481C, 0xA739E65C66CA893E),
	  GFw64be(0x000000BC795A2664, 0xF4FB87D23898AE3A,
	          0x7387D740FBD7C293, 0x0FFD3F3BD703FDAF) },
	{ GFw64be(0x00000053F7B46B4F, 0x1F1EBE465BFA5030,
	          0x79E44D372924A874, 0xC9989D304CCA1530),
	  GFw64be(0x0000017EF3A560D7, 0x3493F9960B150503,
	          0x58F0EF1153F097EC, 0x4E39FA2358086EB0) },
	{ GFw64be(0x0000007178409A1D, 0x6D03F1977A88FC17,
	          0x542848174F565307, 0x36E7A69ACFB8A141),
	  GFw64be(0x0000001AA0B0C544, 0xA2132243DFEEB43B,
	          0x33A1265EE88230B0, 0x01DC4CB6CE0BB283) },
	{ GFw64be(0x0000004EBA9B6554, 0x1D5D352EAAE8356C,
	          0x7781943B6988C1C3, 0x2632D26578637F7B),
	  GFw64be(0x000001E5306759B7, 0x1BCFA839ED79B488,
	          0x6F8481FC2B9C5F90, 0xE6413D774B22BEBF) },
	{ GFw64be(0x0000017FBF8E65A6, 0x05BDCB223CBA013E,
	          0x9EAEDC7514DD44FE, 0xF84CD1007B8B5E13),
	  GFw64be(0x0000003C44E8AA87, 0x0CB11B64B2457E5E,
	          0x4382D1B93560BF41, 0xB5292AAC4F80C14B) },
	{ GFw64be(0x0000013AEBC6EA1D, 0xE4256575430E9624,
	          0xE850CD178844ECDE, 0x51464E66E948B1F3),
	  GFw64be(0x0000007AB5A26669, 0x25AA39337B410448,
	          0x76A14E971C42B1EB, 0xF7C588BA7713919D) },
	{ GFw64be(0x000000AFF8FD0A3E, 0xB341AB31ED4AF57D,
	          0xE6E248444C8A380A, 0x30B5EC05E51CC014),
	  GFw64be(0x00000089448B67DB, 0x78031BAB2194A130,
	          0x5C022F71C5394F8D, 0x52550C88A5EA974B) },
	{ GFw64be(0x00000170ECF68CA0, 0xB047A172C549B1B0,
	          0x5A787F1A45620703, 0x6896B6BCCECEEC94),
	  GFw64be(0x00000061F2FC9387, 0x4137AE95920771B1,
	          0x3E86B2A55C4F9194, 0xBFEB8C857F83CBD4) },
	{ GFw64be(0x0000010FF5535048, 0xD36578345ED1CA83,
	          0x4CAC8AB401699EB3, 0xF2D1E88BDB87478E),
	  GFw64be(0x000001E7563AB20E, 0xFD31F9FFD7212A3E,
	          0x0F3E29FBB2967301, 0x14A34E5E77419E2C) },
	{ GFw64be(0x000000F2E7F2E554, 0xDFD31DB23B10DFCB,
	          0xD961A6B3024A7845, 0x929D2AFC660865BE),
	  GFw64be(0x0000014AA8B682E4, 0xCE1C6785D231DFBF,
	          0xC0570B89EB60B9F1, 0xAC3BF694716F6D10) },
	{ GFw64be(0x000001359953EDEB, 0x3BD70BD7234C72FE,
	          0x67F219185B195F4F, 0xF1BA8329AFC9D557),
	  GFw64be(0x00000169EC28609A, 0x711DFA402BAF850C,
	          0x4F39B23B68F8243C, 0x939554AA029890BA) },
	{ GFw64be(0x0000011EF9E2E290, 0xF14DE4037619B4B8,
	          0x053D0C06232031F9, 0xACD428288908247F),
	  GFw64be(0x000001C3C59C1C13, 0xDF5A9D7DB2F845F4,
	          0x9A4DE5535C6757C1, 0x5F86FDD6824FB123) },
	{ GFw64be(0x000001F6162B6768, 0x3C406076D3BE448C,
	          0x3AC18D1C6ED4368F, 0xC6235F4533E5B26D),
	  GFw64be(0x0000018892E137B8, 0x678A208DA577611D,
	          0x9250354B8449A6F5, 0x15FCACE4F8DA6DF5) },
	{ GFw64be(0x000000A18A916040, 0x2D67AE3F4DD3D509,
	          0x50E48297D79C7D37, 0x059B0CC451E76DED),
	  GFw64be(0x000000BAFA41D6C0, 0x134E707EDAE0154F,
	          0x89C3600E50831606, 0xF5E96F4536546938) },
	{ GFw64be(0x000001C311227D49, 0x3DD9CDE8EB05586E,
	          0xF45B18607FE4CC1A, 0xC7948B54A9423D1F),
	  GFw64be(0x000000087F2DBCDF, 0xE3C7026C8C9E726F,
	          0x60FAA23062E658D9, 0x9E9D99D108F5470C) }
};

/* i*(2^120)*G  for i = 1 to 16 */
static const inner_point_affine XSB233_G120[] = {
	{ GFw64be(0x000001D5899AF530, 0x2439961E6599B17E,
	          0xD05694AC90FC42B2, 0x7BBC9B2DC4E0F799),
	  GFw64be(0x0000014AFF4A70D0, 0xBE8E20D96EBEB8E9,
	          0xADA59DCA8CAD7B79, 0x9871BF2CD9FF8E82) },
	{ GFw64be(0x0000003BEF38B980, 0x2CDA88E1B4A493BB,
	          0xF0563E8621B8A461, 0xABCB666EECF822D1),
	  GFw64be(0x000000F17F26D563, 0x912A87F6BDBCE417,
	          0x666AF7CA600B423C, 0xD4E5326F2F76D968) },
	{ GFw64be(0x000001CD047307D6, 0x2829A326F1344C53,
	          0x4899C6EA4C73D601, 0xDD61E03A8CCB79BB),
	  GFw64be(0x0000011C0636F910, 0x420522519FAD9141,
	          0x4CF01B1D1821A017, 0xC928C93595523D68) },
	{ GFw64be(0x00000145EE0EBB89, 0xBCCFE1CC73911EA0,
	          0xB1DE25990922D1C4, 0xACBF4C2A82A0935B),
	  GFw64be(0x0000001EFB356397, 0x6DB1F70CB3238D10,
	          0x7962D2A09EDCEF68, 0x1B1A39438F0B3AF9) },
	{ GFw64be(0x00000163224648E1, 0xC7F8955711D478F2,
	          0x6B6B2A3FF089D9B5, 0x36B2ACBE3355B2EA),
	  GFw64be(0x0000013FCE531663, 0xAA8DD129492BD3B4,
	          0x4970F872C519FA36, 0x117F3D3A5FC668E7) },
	{ GFw64be(0x00000080662D0C7F, 0x77291D1A0FCE7CE6,
	          0x21536A84BDFCC595, 0x7C9D2D21415D7379),
	  GFw64be(0x00000124F6C6BBEF, 0x8F49BEDCFCBE8100,
	          0x7BE0BC97B6A1349D, 0xDF356C0B7AF663FC) },
	{ GFw64be(0x000000AA886C00DD, 0x0A21569E4C720F07,
	          0xD8DA6A2DBADB1E56, 0x4B905E3A90CBEB13),
	  GFw64be(0x000001F6FA0A71DE, 0x94094A4D475D1C0F,
	          0x70A91016D884B82F, 0xF642894F3F5350B1) },
	{ GFw64be(0x000001CF6A05CBCB, 0x7780FD99F274DA15,
	          0x854507EE383D83F1, 0x13C7B96A99548283),
	  GFw64be(0x0000012D2D0F0F68, 0xD2AF8F1423BB76DA,
	          0xB22352603FAEE448, 0x6B40114D0A67B29C) },
	{ GFw64be(0x0000002F6CDA1861, 0x6F64CBEBA71846C5,
	          0x756D252D54E4AC0A, 0x1A3370AD580C09A7),
	  GFw64be(0x0000011DC42108AD, 0xD1FD0D4B8B1F0AD7,
	          0x95AFB90E15970BA9, 0x554B80EE25B3908D) },
	{ GFw64be(0x00000092050E6E7D, 0x1CB1AC1DCE45225D,
	          0x42800904939698F8, 0xA40309F9A0804DD5),
	  GFw64be(0x000001DF5EB5CB26, 0x89E28D6D8DB25F46,
	          0xFAD9964C126B8A0B, 0xC821A41640D1C174) },
	{ GFw64be(0x000000077DFF31B5, 0x377B004AF5BC56D6,
	          0xC6025FD1F065BDC5, 0x39194CEFA7D4E701),
	  GFw64be(0x0000014023B16017, 0x22C5EA22F073E537,
	          0xE217C153F2F590DE, 0x94CDE91B7EE19EE1) },
	{ GFw64be(0x0000009FE1A598DB, 0x2265583AA535AA74,
	          0x70C26307E7FB1569, 0x2E3EF98A6B2B6351),
	  GFw64be(0x0000003F21233852, 0x602B8FD756CEF166,
	          0x78A5B8E66BE35BFE, 0xEFB0D68AC048DA27) },
	{ GFw64be(0x000001B30D1AE884, 0xE2E2DB25A08EC752,
	          0x3A524338652B72E6, 0x4F5BD3363F703144),
	  GFw64be(0x00000124B73932B1, 0x8280AF936830BBD7,
	          0xBEB937B1F1A7E624, 0x1DE10D903A5867B3) },
	{ GFw64be(0x000000151C04C3E7, 0xCCCC3D98A8082F56,
	          0x495DC71A02A0DB8C, 0x9EEB47D7FDD06230),
	  GFw64be(0x000001B9586B5DD2, 0x639B641CBAC3BAE7,
	          0x9DDF5221F79B1E0C, 0xA27D8C510512E89D) },
	{ GFw64be(0x000001D61BE4FB5F, 0xB08B6852C23433E9,
	          0xBA6AC5D563D8C0A3, 0xD158BE35FE6716D8),
	  GFw64be(0x0000016D35094E28, 0x0682ADFFEE94197C,
	          0xAB1C43EA8C457A85, 0x3F28F2C9FDA2F7DC) },
	{ GFw64be(0x00000107917CF4FE, 0xA5359C9091BE963B,
	          0x36A6A66E787A63ED, 0x5456AE016C1A45C9),
	  GFw64be(0x000000891FE54D04, 0x84521BADA92E3974,
	          0x85771193572865D8, 0x29207F0F1A478208) }
};

/* i*(2^180)*G  for i = 1 to 16 */
static const inner_point_affine XSB233_G180[] = {
	{ GFw64be(0x000000673F17BED8, 0x9F71A5C5DBBAF08A,
	          0xC0B6B19F1916E75A, 0x1BF81499AA2B8734),
	  GFw64be(0x000001A45AD8B438, 0xBD985C847030444E,
	          0x5EC9C240E30D1D23, 0x7D95F02930922F27) },
	{ GFw64be(0x000001A46001E6D4, 0x78D95D9786DEC1E0,
	          0x76CBB5611B8926AA, 0xCF21BCFE22EBE92C),
	  GFw64be(0x000001698B879B23, 0xE7D692732FD65CA7,
	          0x3E4D359D8A413147, 0x2C1C71A890CFDD21) },
	{ GFw64be(0x000000DFE6E6E78C, 0xDBDA4C4B753F9326,
	          0x53BD3FEF28315EAC, 0x1DC64EFD011BA71E),
	  GFw64be(0x000001372A556ACB, 0x13424CD244834ADB,
	          0xDB8B31B577A7F106, 0x3B5563BE15E78A7C) },
	{ GFw64be(0x00000102E5D7D288, 0x5A092937F290BE51,
	          0x6D5E381B214CB62F, 0x84E291C95C286A25),
	  GFw64be(0x000001ABA5458846, 0x93D08200FA1DF533,
	          0xAF90CB8EF6F525D2, 0x3795AC7236EE6542) },
	{ GFw64be(0x00000178705C29EC, 0x07D63AE9092FF255,
	          0x0638A33BE5B006FA, 0x0D75CFB801992B01),
	  GFw64be(0x0000007B85D5E5F5, 0xF10DB31BE36BCE5D,
	          0x73DAE4DFEAA40A8C, 0x5CB4E3325B4485AC) },
	{ GFw64be(0x00000036AED4477B, 0x31964FA7C53BDD48,
	          0x9D4C7FE1F325B6AD, 0xE83ACF215EC3E637),
	  GFw64be(0x0000003FFC46CEF1, 0x204C77D84473F105,
	          0x9B14BFE2709ABA58, 0x6033E466379CC433) },
	{ GFw64be(0x00000011B04BA847, 0xB83802F7E6B39E84,
	          0x5E6350888F7BABAD, 0x0513B7CFB45B877C),
	  GFw64be(0x000001218E623F33, 0xCB4AB8DCED591393,
	          0xA4FA6F8FFD9BE7D4, 0x9ABA14677A72F175) },
	{ GFw64be(0x000000DE621B61D4, 0x8F0DC081DB5F1E37,
	          0xE50A67AA8A56C748, 0x331CE53052DFED0B),
	  GFw64be(0x00000095DD724D86, 0x09C07F786539768D,
	          0x009C5A73DC681DE1, 0xDF46F4FAAC0353B3) },
	{ GFw64be(0x000001839FB0B866, 0x055C10348D464F71,
	          0x38372CF37F843942, 0xCA4D62F1F8A73A9D),
	  GFw64be(0x0000016E6DA2304C, 0x15713C6750FEF82E,
	          0x8B0CA3ABD6E60C0C, 0x990BC5C5B6B32ADB) },
	{ GFw64be(0x000001BCFF839DBF, 0x34258669AE53E097,
	          0xEC922BE9AD4508F5, 0x1C526269EAEE0805),
	  GFw64be(0x0000003DBB5D9269, 0xE0A74D62E509741C,
	          0xA192D033A7F14014, 0xB3AD9E25E24CE93D) },
	{ GFw64be(0x0000016A6601B57A, 0x270C4AFB8BA64A56,
	          0xB1D9CC627C131ED4, 0x88CCF2F146F6B6F1),
	  GFw64be(0x000001C5F8ED525E, 0xFF0C880CF20029E4,
	          0x4E8371DB82819DCA, 0x2D41F1096676A374) },
	{ GFw64be(0x00000072A08A61E9, 0x7BFE2FEAF7B23B56,
	          0x5423013535DFEEE1, 0xDD319460FCBB552E),
	  GFw64be(0x00000114F4EBC129, 0xACF49A1512F986C7,
	          0x1F4549905FF5CDBA, 0x5751A7CF0DE4DD8D) },
	{ GFw64be(0x000000561056D6CC, 0xCF81B5BA65E180A9,
	          0x7DBAB8A2B8F7C3D8, 0x33E650DCE1CF074B),
	  GFw64be(0x000000EB2A0ADE8C, 0x5E4BC8AA1B029E0E,
	          0x8F014E9993859227, 0x546AC7390D5A8AEC) },
	{ GFw64be(0x000000B8AF92907D, 0x9D7751FC8B148BBD,
	          0x4A9FC5A557DF0D39, 0xB874264A4D667886),
	  GFw64be(0x000001BDCE4C4422, 0xC9897C0727498ED4,
	          0x10B21DFA69E23FF3, 0xFF5602052ACE422B) },
	{ GFw64be(0x000001FCC2036CBC, 0xADC798617C5AF05A,
	          0x1B2430EFED3407C0, 0x70F7BB42BECDCC5C),
	  GFw64be(0x000000E19F69375E, 0xA2AACAA059187D3E,
	          0x0FBE05F9FCACAAD6, 0x1C12C3BCC5D3A86E) },
	{ GFw64be(0x0000002D314C6976, 0xA8764B7F0DFE4945,
	          0x7C9256E9EE1C9AB6, 0x560E8693657388AB),
	  GFw64be(0x00000035D1885965, 0x4E399D60C594A499,
	          0x70B69FDB22737EFF, 0x6249AC03D0CE3E7C) }
};
