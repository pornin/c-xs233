#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "gf233.h"

/*
 * Perfunctory BLAKE2s implementation, used for making pseudorandom values
 * for tests.
 */

typedef struct {
	uint8_t buf[64];
	uint32_t h[8];
	uint64_t ctr;
	size_t out_len;
} blake2s_context;

static inline uint32_t
dec32le(const void *src)
{
	const uint8_t *buf = src;

	return (uint32_t)buf[0]
		| ((uint32_t)buf[1] << 8)
		| ((uint32_t)buf[2] << 16)
		| ((uint32_t)buf[3] << 24);
}

static inline void
enc32le(void *dst, uint32_t x)
{
	uint8_t *buf = dst;

	buf[0] = (uint8_t)x;
	buf[1] = (uint8_t)(x >> 8);
	buf[2] = (uint8_t)(x >> 16);
	buf[3] = (uint8_t)(x >> 24);
}

static const uint32_t IV[] = {
	0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A,
	0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19
};

static void
process_block(uint32_t *h, const uint8_t *data, uint64_t t, int f)
{
	uint32_t v[16], m[16];
	int i;

	memcpy(v, h, 8 * sizeof(uint32_t));
	memcpy(v + 8, IV, sizeof IV);
	v[12] ^= (uint32_t)t;
	v[13] ^= (uint32_t)(t >> 32);
	if (f) {
		v[14] = ~v[14];
	}
	for (i = 0; i < 16; i ++) {
		m[i] = dec32le(data + (i << 2));
	}

#define ROR(x, n)   (((x) << (32 - (n))) | ((x) >> (n)))

#define G(a, b, c, d, x, y)   do { \
		v[a] += v[b] + (x); \
		v[d] = ROR(v[d] ^ v[a], 16); \
		v[c] += v[d]; \
		v[b] = ROR(v[b] ^ v[c], 12); \
		v[a] += v[b] + (y); \
		v[d] = ROR(v[d] ^ v[a], 8); \
		v[c] += v[d]; \
		v[b] = ROR(v[b] ^ v[c], 7); \
	} while (0)

#define ROUND(s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, sA, sB, sC, sD, sE, sF) \
	do { \
		G(0, 4,  8, 12, m[s0], m[s1]); \
		G(1, 5,  9, 13, m[s2], m[s3]); \
		G(2, 6, 10, 14, m[s4], m[s5]); \
		G(3, 7, 11, 15, m[s6], m[s7]); \
		G(0, 5, 10, 15, m[s8], m[s9]); \
		G(1, 6, 11, 12, m[sA], m[sB]); \
		G(2, 7,  8, 13, m[sC], m[sD]); \
		G(3, 4,  9, 14, m[sE], m[sF]); \
	} while (0)

	ROUND( 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15);
	ROUND(14, 10,  4,  8,  9, 15, 13,  6,  1, 12,  0,  2, 11,  7,  5,  3);
	ROUND(11,  8, 12,  0,  5,  2, 15, 13, 10, 14,  3,  6,  7,  1,  9,  4);
	ROUND( 7,  9,  3,  1, 13, 12, 11, 14,  2,  6,  5, 10,  4,  0, 15,  8);
	ROUND( 9,  0,  5,  7,  2,  4, 10, 15, 14,  1, 11, 12,  6,  8,  3, 13);
	ROUND( 2, 12,  6, 10,  0, 11,  8,  3,  4, 13,  7,  5, 15, 14,  1,  9);
	ROUND(12,  5,  1, 15, 14, 13,  4, 10,  0,  7,  6,  3,  9,  2,  8, 11);
	ROUND(13, 11,  7, 14, 12,  1,  3,  9,  5,  0, 15,  4,  8,  6,  2, 10);
	ROUND( 6, 15, 14,  9, 11,  3,  0,  8, 12,  2, 13,  7,  1,  4, 10,  5);
	ROUND(10,  2,  8,  4,  7,  6,  1,  5, 15, 11,  9, 14,  3, 12, 13,  0);

#undef ROR
#undef G
#undef ROUND

	for (i = 0; i < 8; i ++) {
		h[i] ^= v[i] ^ v[i + 8];
	}
}

/*
 * State rules:
 *
 *   buf    buffered data
 *   h      current state
 *   ctr    number of bytes injected so far
 *
 * Initially, ctr == 0 and h contains the XOR of IV and parameter block;
 * buf[] is empty. For any ctr > 0, buf[] is non-empty; it might contain
 * a full block worth of data (processing of the block is delayed until
 * we know whether this is the final block or not).
 *
 * If a key is injected, then it counts as a first full block.
 */

static void
blake2s_init(blake2s_context *bc, size_t out_len)
{
	memcpy(bc->h, IV, sizeof bc->h);
	bc->h[0] ^= 0x01010000 ^ (uint32_t)out_len;
	bc->ctr = 0;
	bc->out_len = out_len;
}

/* unused
static void
blake2s_init_key(blake2s_context *bc, size_t out_len,
	const void *key, size_t key_len)
{
	blake2s_init(bc, out_len);
	if (key_len > 0) {
		bc->h[0] ^= (uint32_t)key_len << 8;
		memcpy(bc->buf, key, key_len);
		memset(bc->buf + key_len, 0, (sizeof bc->buf) - key_len);
		bc->ctr = sizeof bc->buf;
	}
}
*/

static void
blake2s_update(blake2s_context *bc, const void *data, size_t len)
{
	uint64_t ctr;
	size_t p;

	/* Special case: if no input data, return immediately. */
	if (len == 0) {
		return;
	}

	ctr = bc->ctr;

	/* First complete the current block, if not already full. */
	p = (size_t)ctr & ((sizeof bc->buf) - 1);
	if (ctr == 0 || p != 0) {
		/* buffer is not full */
		size_t clen;

		clen = sizeof bc->buf - p;
		if (clen >= len) {
			memcpy(bc->buf + p, data, len);
			bc->ctr = ctr + len;
			return;
		}
		memcpy(bc->buf + p, data, clen);
		ctr += clen;
		data = (const uint8_t *)data + clen;
		len -= clen;
	}

	/* Process the buffered block. */
	process_block(bc->h, bc->buf, ctr, 0);

	/* Process all subsequent full blocks, except the last. */
	while (len > sizeof bc->buf) {
		ctr += sizeof bc->buf;
		process_block(bc->h, data, ctr, 0);
		data = (const uint8_t *)data + sizeof bc->buf;
		len -= sizeof bc->buf;
	}

	/* Copy the last block (possibly partial) into the buffer. */
	memcpy(bc->buf, data, len);
	bc->ctr = ctr + len;
}

static void
blake2s_final(blake2s_context *bc, void *dst)
{
	int i;
	uint8_t tmp[32];
	size_t p;

	/* Pad the current block with zeros, if not full. If the
	   buffer is empty (no key, no data) then fill it with zeros
	   as well. */
	p = (size_t)bc->ctr & ((sizeof bc->buf) - 1);
	if (bc->ctr == 0 || p != 0) {
		memset(bc->buf + p, 0, (sizeof bc->buf) - p);
	}

	process_block(bc->h, bc->buf, bc->ctr, 1);
	for (i = 0; i < 8; i ++) {
		enc32le(tmp + (i << 2), bc->h[i]);
	}
	memcpy(dst, tmp, bc->out_len);
}

/* unused
static void
blake2s(void *dst, size_t dst_len, const void *key, size_t key_len,
	const void *src, size_t src_len)
{
	blake2s_context bc;

	blake2s_init_key(&bc, dst_len, key, key_len);
	blake2s_update(&bc, src, src_len);
	blake2s_final(&bc, dst);
}
*/

static void
rndgf(gf *d, const void *seed, size_t seed_len, uint32_t k)
{
	blake2s_context bc;
	uint8_t tmp[32];

	tmp[0] = (uint8_t)k;
	tmp[1] = (uint8_t)(k >> 8);
	tmp[2] = (uint8_t)(k >> 16);
	tmp[3] = (uint8_t)(k >> 24);
	blake2s_init(&bc, 32);
	blake2s_update(&bc, seed, seed_len);
	blake2s_update(&bc, tmp, 4);
	blake2s_final(&bc, tmp);
	for (int i = 0; i < 4; i ++) {
		uint64_t s = 0;
		for (int j = 0; j < 8; j ++) {
			s |= (uint64_t)tmp[8 * i + j] << (8 * j);
		}
		d->w[i] = s;
	}
	d->w[0] &= GF_TEST_RMASK0;
	d->w[1] &= GF_TEST_RMASK1;
	d->w[2] &= GF_TEST_RMASK2;
	d->w[3] &= GF_TEST_RMASK3;
}

static void
ref_gf_mul(gf *d, const gf *a, const gf *b)
{
	uint64_t x[4], y[4], z[4];
	uint8_t tmp[32];

	memset(tmp, 0, sizeof tmp);
	gf_encode(tmp, a);
	for (int i = 0; i < 4; i ++) {
		x[i] = (uint64_t)tmp[8 * i + 0]
			| ((uint64_t)tmp[8 * i + 1] << 8)
			| ((uint64_t)tmp[8 * i + 2] << 16)
			| ((uint64_t)tmp[8 * i + 3] << 24)
			| ((uint64_t)tmp[8 * i + 4] << 32)
			| ((uint64_t)tmp[8 * i + 5] << 40)
			| ((uint64_t)tmp[8 * i + 6] << 48)
			| ((uint64_t)tmp[8 * i + 7] << 56);
	}
	gf_encode(tmp, b);
	for (int i = 0; i < 4; i ++) {
		y[i] = (uint64_t)tmp[8 * i + 0]
			| ((uint64_t)tmp[8 * i + 1] << 8)
			| ((uint64_t)tmp[8 * i + 2] << 16)
			| ((uint64_t)tmp[8 * i + 3] << 24)
			| ((uint64_t)tmp[8 * i + 4] << 32)
			| ((uint64_t)tmp[8 * i + 5] << 40)
			| ((uint64_t)tmp[8 * i + 6] << 48)
			| ((uint64_t)tmp[8 * i + 7] << 56);
	}

	z[0] = 0;
	z[1] = 0;
	z[2] = 0;
	z[3] = 0;
	for (int i = 0; i < 233; i ++) {
		uint64_t m = -((y[i >> 6] >> (i & 63)) & 1);
		z[0] ^= m & x[0];
		z[1] ^= m & x[1];
		z[2] ^= m & x[2];
		z[3] ^= m & x[3];

		uint64_t t = (x[3] >> 40) & 1;
		x[3] = ((x[3] << 1) | (x[2] >> 63)) & 0x000001FFFFFFFFFF;
		x[2] = (x[2] << 1) | (x[1] >> 63);
		x[1] = (x[1] << 1) | (x[0] >> 63);
		x[0] = x[0] << 1;
		x[0] ^= t;
		x[1] ^= (t << 10);
	}
	*d = (gf)GFw64be(z[3], z[2], z[1], z[0]);
}

static void
check_eq_gf(const char *msg, const gf *a, const gf *b)
{
	gf x, y;

	gf_norm(&x, a);
	gf_norm(&y, b);
	if (((x.w[0] ^ y.w[0]) | (x.w[1] ^ y.w[1])
		| (x.w[2] ^ y.w[2]) | (x.w[3] ^ y.w[3])) == 0)
	{
		return;
	}
	printf("\nFAIL: %s\n", msg);
	printf("a = %016lX %016lX %016lX %016lX\n",
		a->w[3], a->w[2], a->w[1], a->w[0]);
	printf("b = %016lX %016lX %016lX %016lX\n",
		b->w[3], b->w[2], b->w[1], b->w[0]);
	printf("x = %016lX %016lX %016lX %016lX\n",
		x.w[3], x.w[2], x.w[1], x.w[0]);
	printf("y = %016lX %016lX %016lX %016lX\n",
		y.w[3], y.w[2], y.w[1], y.w[0]);
	exit(EXIT_FAILURE);
}

static void
test_gf(void)
{
	printf("Test gf: ");
	fflush(stdout);

	const char *seed = "test_gf";
	size_t seed_len = strlen(seed);
	for (uint32_t k = 0; k < 1000; k ++) {
		gf a, b, c, d;
		uint64_t za, zb, eab;
		uint8_t tmp[30];

		rndgf(&a, seed, seed_len, 2 * k + 0);
		rndgf(&b, seed, seed_len, 2 * k + 1);
		za = 0;
		zb = 0;
		eab = 0;
		if (k >= 10 && k < 14) {
			switch ((k - 10) & 1) {
			case 0:
				a = GF_ZERO;
				za = (uint64_t)-1;
				break;
			case 1:
				a.w[0] = GF_TEST_RMASK0;
				a.w[1] = GF_TEST_RMASK1;
				a.w[2] = GF_TEST_RMASK2;
				a.w[3] = GF_TEST_RMASK3;
				break;
			}
			switch (((k - 10) >> 1) & 1) {
			case 0:
				b = GF_ZERO;
				zb = (uint64_t)-1;
				break;
			case 1:
				b.w[0] = GF_TEST_RMASK0;
				b.w[1] = GF_TEST_RMASK1;
				b.w[2] = GF_TEST_RMASK2;
				b.w[3] = GF_TEST_RMASK3;
				break;
			}
			eab = -(uint64_t)
				(((k - 10) & 1) == (((k - 10) >> 1) & 1));
		} else if (k == 14) {
			b = a;
			eab = (uint64_t)-1;
		}
		if (gf_equals(&a, &a) != (uint64_t)-1) {
			printf("FAIL: equal 1\n");
			exit(EXIT_FAILURE);
		}
		if (gf_equals(&a, &b) != eab) {
			printf("FAIL: equal 2\n");
			exit(EXIT_FAILURE);
		}
		gf_norm(&c, &a);
		if (gf_equals(&c, &a) != (uint64_t)-1) {
			printf("FAIL: equal 3\n");
			exit(EXIT_FAILURE);
		}
		if (gf_equals(&c, &b) != eab) {
			printf("FAIL: equal 4\n");
			exit(EXIT_FAILURE);
		}
		if (gf_iszero(&a) != za || gf_iszero(&b) != zb) {
			printf("FAIL: iszero 1\n");
			exit(EXIT_FAILURE);
		}
		gf_add(&d, &a, &c);
		if (gf_iszero(&d) != (uint64_t)-1) {
			printf("FAIL: iszero 2\n");
			exit(EXIT_FAILURE);
		}
		gf_mul(&c, &a, &b);
		ref_gf_mul(&d, &a, &b);
		check_eq_gf("mul", &c, &d);
		gf_sqr(&c, &a);
		ref_gf_mul(&d, &a, &a);
		check_eq_gf("sqr", &c, &d);
		gf_sqrt(&c, &b);
		gf_sqr(&d, &c);
		check_eq_gf("sqrt", &b, &d);
		gf_halftrace(&c, &a);
		gf_sqr(&d, &c);
		gf_add(&d, &c, &d);
		if (gf_trace(&a) != 0) {
			gf_add(&d, &d, &GF_ONE);
		}
		check_eq_gf("halftrace", &d, &a);

		gf_inv(&c, &a);
		if (za) {
			check_eq_gf("inv0", &c, &GF_ZERO);
		} else {
			gf_mul(&d, &c, &a);
			check_eq_gf("inv", &d, &GF_ONE);
		}

		gf_div(&c, &a, &b);
		if (zb) {
			check_eq_gf("div0", &c, &GF_ZERO);
		} else {
			gf_mul(&d, &c, &b);
			check_eq_gf("div", &d, &a);
		}

		gf_encode(tmp, &a);
		if (gf_decode(&c, tmp) != (uint64_t)-1) {
			printf("FAIL: decode 1\n");
			exit(EXIT_FAILURE);
		}
		check_eq_gf("enc/dec", &a, &c);
		for (int i = 1; i < 8; i ++) {
			tmp[29] = (tmp[29] & 1) | (1 << i);
			if (gf_decode(&c, tmp) != 0) {
				printf("FAIL: decode 2\n");
				exit(EXIT_FAILURE);
			}
			check_eq_gf("dec2", &c, &GF_ZERO);
		}

		if (k % 100 == 0) {
			printf(".");
			fflush(stdout);
		}
	}

	printf(" done.\n");
	fflush(stdout);
}

int
main(void)
{
	test_gf();
	return 0;
}
