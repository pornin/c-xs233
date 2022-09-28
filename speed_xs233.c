#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "xs233.h"

#if defined __aarch64__
static inline uint64_t
core_cycles(void)
{
	uint64_t d;
	__asm__ __volatile__ ( "dsb sy\nmrs %0, pmccntr_el0" : "=r" (d));
	return d;
}
#else

#include <immintrin.h>

static inline uint64_t
core_cycles(void)
{
	/*
	 * GCC seems not to have the __rdtsc() intrinsic.
	 */
#if defined __GNUC__ && !defined __clang__
	uint32_t hi, lo;

	_mm_lfence();
	__asm__ __volatile__ ("rdtsc" : "=d" (hi), "=a" (lo) : : );
	return ((uint64_t)hi << 32) | (uint64_t)lo;
#else
	_mm_lfence();
	return __rdtsc();
#endif
}
#endif

static int
cmp_u64(const void *p1, const void *p2)
{
	uint64_t v1, v2;

	v1 = *(const uint64_t *)p1;
	v2 = *(const uint64_t *)p2;
	if (v1 < v2) {
		return -1;
	} else if (v1 == v2) {
		return 0;
	} else {
		return 1;
	}
}

static void
init_buf_cycles(uint8_t *buf, size_t len)
{
	for (size_t u = 0; u < len; u ++) {
		buf[u] = (uint8_t)core_cycles();
	}
}

static uint32_t
speed_xsb233_mul(void)
{
	size_t u;
	uint64_t tt[100];
	uint8_t scalar[30];
	xsb233_point P;

#define NUM   ((sizeof tt) / (sizeof tt[0]))

	init_buf_cycles(scalar, sizeof scalar);
	xsb233_mulgen(&P, scalar, sizeof scalar);
	for (u = 0; u < 2 * NUM; u ++) {
		uint64_t begin, end;
		int i;

		begin = core_cycles();
		for (i = 0; i < 40; i ++) {
			xsb233_mul(&P, &P, scalar, sizeof scalar);
			memcpy(scalar, &P, sizeof scalar);
		}
		end = core_cycles();
		if (u >= NUM) {
			tt[u - NUM] = end - begin;
		}
	}
	qsort(tt, NUM, sizeof tt[0], &cmp_u64);
	printf("xsb233_mul:              %9.2f (%.2f .. %.2f)\n",
		(double)tt[NUM / 2] / 40.0,
		(double)tt[NUM / 10] / 40.0,
		(double)tt[(9 * NUM) / 10] / 40.0);
	fflush(stdout);
	return (uint32_t)scalar[0];
}

static uint32_t
speed_xsb233_mul_ladder(void)
{
	size_t u;
	uint64_t tt[100];
	uint8_t scalar[30];
	xsb233_point P;

#define NUM   ((sizeof tt) / (sizeof tt[0]))

	init_buf_cycles(scalar, sizeof scalar);
	xsb233_mulgen(&P, scalar, sizeof scalar);
	for (u = 0; u < 2 * NUM; u ++) {
		uint64_t begin, end;
		int i;

		begin = core_cycles();
		for (i = 0; i < 40; i ++) {
			xsb233_mul_ladder(&P, &P, scalar, sizeof scalar);
			memcpy(scalar, &P, sizeof scalar);
		}
		end = core_cycles();
		if (u >= NUM) {
			tt[u - NUM] = end - begin;
		}
	}
	qsort(tt, NUM, sizeof tt[0], &cmp_u64);
	printf("xsb233_mul_ladder:       %9.2f (%.2f .. %.2f)\n",
		(double)tt[NUM / 2] / 40.0,
		(double)tt[NUM / 10] / 40.0,
		(double)tt[(9 * NUM) / 10] / 40.0);
	fflush(stdout);
	return (uint32_t)scalar[0];
}

static uint32_t
speed_xsb233_mulgen(void)
{
	size_t u;
	uint64_t tt[100];
	uint8_t scalar[30];

#define NUM   ((sizeof tt) / (sizeof tt[0]))

	init_buf_cycles(scalar, sizeof scalar);
	for (u = 0; u < 2 * NUM; u ++) {
		uint64_t begin, end;
		int i;

		begin = core_cycles();
		for (i = 0; i < 40; i ++) {
			xsb233_point P;

			xsb233_mulgen(&P, scalar, sizeof scalar);
			memcpy(scalar, &P, sizeof scalar);
		}
		end = core_cycles();
		if (u >= NUM) {
			tt[u - NUM] = end - begin;
		}
	}
	qsort(tt, NUM, sizeof tt[0], &cmp_u64);
	printf("xsb233_mulgen:           %9.2f (%.2f .. %.2f)\n",
		(double)tt[NUM / 2] / 40.0,
		(double)tt[NUM / 10] / 40.0,
		(double)tt[(9 * NUM) / 10] / 40.0);
	fflush(stdout);
	return (uint32_t)scalar[0];
}

static uint32_t
speed_xsk233_tau_recode(void)
{
	extern void xsk233_tau_recode(int8_t *sd, const void *n);

	size_t u;
	uint64_t tt[100];
	uint8_t scalar[30];

#define NUM   ((sizeof tt) / (sizeof tt[0]))

	init_buf_cycles(scalar, sizeof scalar);
	for (u = 0; u < 2 * NUM; u ++) {
		uint64_t begin, end;
		int i;

		begin = core_cycles();
		for (i = 0; i < 400; i ++) {
			int8_t sd[48];

			xsk233_tau_recode(sd, scalar);
			memcpy(scalar, sd, sizeof scalar);
		}
		end = core_cycles();
		if (u >= NUM) {
			tt[u - NUM] = end - begin;
		}
	}
	qsort(tt, NUM, sizeof tt[0], &cmp_u64);
	printf("xsk233_tau_recode:       %9.2f (%.2f .. %.2f)\n",
		(double)tt[NUM / 2] / 400.0,
		(double)tt[NUM / 10] / 400.0,
		(double)tt[(9 * NUM) / 10] / 400.0);
	fflush(stdout);
	return (uint32_t)scalar[0];
}

static uint32_t
speed_xsk233_mul(void)
{
	size_t u;
	uint64_t tt[100];
	uint8_t scalar[30];
	xsk233_point P;

#define NUM   ((sizeof tt) / (sizeof tt[0]))

	init_buf_cycles(scalar, sizeof scalar);
	xsk233_mulgen(&P, scalar, sizeof scalar);
	for (u = 0; u < 2 * NUM; u ++) {
		uint64_t begin, end;
		int i;

		begin = core_cycles();
		for (i = 0; i < 40; i ++) {
			xsk233_mul(&P, &P, scalar, sizeof scalar);
			memcpy(scalar, &P, sizeof scalar);
		}
		end = core_cycles();
		if (u >= NUM) {
			tt[u - NUM] = end - begin;
		}
	}
	qsort(tt, NUM, sizeof tt[0], &cmp_u64);
	printf("xsk233_mul:              %9.2f (%.2f .. %.2f)\n",
		(double)tt[NUM / 2] / 40.0,
		(double)tt[NUM / 10] / 40.0,
		(double)tt[(9 * NUM) / 10] / 40.0);
	fflush(stdout);
	return (uint32_t)scalar[0];
}

static uint32_t
speed_xsk233_mul_ladder(void)
{
	size_t u;
	uint64_t tt[100];
	uint8_t scalar[30];
	xsk233_point P;

#define NUM   ((sizeof tt) / (sizeof tt[0]))

	init_buf_cycles(scalar, sizeof scalar);
	xsk233_mulgen(&P, scalar, sizeof scalar);
	for (u = 0; u < 2 * NUM; u ++) {
		uint64_t begin, end;
		int i;

		begin = core_cycles();
		for (i = 0; i < 40; i ++) {
			xsk233_mul_ladder(&P, &P, scalar, sizeof scalar);
			memcpy(scalar, &P, sizeof scalar);
		}
		end = core_cycles();
		if (u >= NUM) {
			tt[u - NUM] = end - begin;
		}
	}
	qsort(tt, NUM, sizeof tt[0], &cmp_u64);
	printf("xsk233_mul_ladder:       %9.2f (%.2f .. %.2f)\n",
		(double)tt[NUM / 2] / 40.0,
		(double)tt[NUM / 10] / 40.0,
		(double)tt[(9 * NUM) / 10] / 40.0);
	fflush(stdout);
	return (uint32_t)scalar[0];
}

static uint32_t
speed_xsk233_mul_frob(void)
{
	size_t u;
	uint64_t tt[100];
	uint8_t scalar[30];
	xsk233_point P;

#define NUM   ((sizeof tt) / (sizeof tt[0]))

	init_buf_cycles(scalar, sizeof scalar);
	xsk233_mulgen(&P, scalar, sizeof scalar);
	for (u = 0; u < 2 * NUM; u ++) {
		uint64_t begin, end;
		int i;

		begin = core_cycles();
		for (i = 0; i < 40; i ++) {
			xsk233_mul_frob(&P, &P, scalar, sizeof scalar);
			memcpy(scalar, &P, sizeof scalar);
		}
		end = core_cycles();
		if (u >= NUM) {
			tt[u - NUM] = end - begin;
		}
	}
	qsort(tt, NUM, sizeof tt[0], &cmp_u64);
	printf("xsk233_mul_frob:         %9.2f (%.2f .. %.2f)\n",
		(double)tt[NUM / 2] / 40.0,
		(double)tt[NUM / 10] / 40.0,
		(double)tt[(9 * NUM) / 10] / 40.0);
	fflush(stdout);
	return (uint32_t)scalar[0];
}

static uint32_t
speed_xsk233_mulgen(void)
{
	size_t u;
	uint64_t tt[100];
	uint8_t scalar[30];

#define NUM   ((sizeof tt) / (sizeof tt[0]))

	init_buf_cycles(scalar, sizeof scalar);
	for (u = 0; u < 2 * NUM; u ++) {
		uint64_t begin, end;
		int i;

		begin = core_cycles();
		for (i = 0; i < 40; i ++) {
			xsk233_point P;

			xsk233_mulgen(&P, scalar, sizeof scalar);
			memcpy(scalar, &P, sizeof scalar);
		}
		end = core_cycles();
		if (u >= NUM) {
			tt[u - NUM] = end - begin;
		}
	}
	qsort(tt, NUM, sizeof tt[0], &cmp_u64);
	printf("xsk233_mulgen:           %9.2f (%.2f .. %.2f)\n",
		(double)tt[NUM / 2] / 40.0,
		(double)tt[NUM / 10] / 40.0,
		(double)tt[(9 * NUM) / 10] / 40.0);
	fflush(stdout);
	return (uint32_t)scalar[0];
}

static uint32_t
speed_xsk233_mulgen_frob(void)
{
	size_t u;
	uint64_t tt[100];
	uint8_t scalar[30];

#define NUM   ((sizeof tt) / (sizeof tt[0]))

	init_buf_cycles(scalar, sizeof scalar);
	for (u = 0; u < 2 * NUM; u ++) {
		uint64_t begin, end;
		int i;

		begin = core_cycles();
		for (i = 0; i < 40; i ++) {
			xsk233_point P;

			xsk233_mulgen_frob(&P, scalar, sizeof scalar);
			memcpy(scalar, &P, sizeof scalar);
		}
		end = core_cycles();
		if (u >= NUM) {
			tt[u - NUM] = end - begin;
		}
	}
	qsort(tt, NUM, sizeof tt[0], &cmp_u64);
	printf("xsk233_mulgen_frob:      %9.2f (%.2f .. %.2f)\n",
		(double)tt[NUM / 2] / 40.0,
		(double)tt[NUM / 10] / 40.0,
		(double)tt[(9 * NUM) / 10] / 40.0);
	fflush(stdout);
	return (uint32_t)scalar[0];
}

int
main(void)
{
	uint32_t x = 0;

	x ^= speed_xsb233_mul();
	x ^= speed_xsb233_mul_ladder();
	x ^= speed_xsb233_mulgen();
	printf("\n");
	x ^= speed_xsk233_tau_recode();
	x ^= speed_xsk233_mul();
	x ^= speed_xsk233_mul_ladder();
	x ^= speed_xsk233_mul_frob();
	x ^= speed_xsk233_mulgen();
	x ^= speed_xsk233_mulgen_frob();
	printf("%u\n", (unsigned)x);
	return 0;
}
