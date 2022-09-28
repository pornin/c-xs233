#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gf233.h"

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
init_gf_cycles(gf *x)
{
	x->w[0] = core_cycles() & GF_TEST_RMASK0;
	x->w[1] = core_cycles() & GF_TEST_RMASK1;
	x->w[2] = core_cycles() & GF_TEST_RMASK2;
	x->w[3] = core_cycles() & GF_TEST_RMASK3;
}

static uint32_t
speed_mul(void)
{
	size_t u;
	uint64_t tt[100];
	gf x, y;

#define NUM   ((sizeof tt) / (sizeof tt[0]))

	init_gf_cycles(&x);
	init_gf_cycles(&y);
	for (u = 0; u < 2 * NUM; u ++) {
		uint64_t begin, end;
		int i;

		begin = core_cycles();
		for (i = 0; i < 1000; i ++) {
			gf_mul(&x, &x, &y);
			gf_mul(&y, &x, &y);
			gf_mul(&x, &x, &y);
			gf_mul(&y, &x, &y);
			gf_mul(&x, &x, &y);
			gf_mul(&y, &x, &y);
		}
		end = core_cycles();
		if (u >= NUM) {
			tt[u - NUM] = end - begin;
		}
	}
	qsort(tt, NUM, sizeof tt[0], &cmp_u64);
	printf("gf_mul:                  %9.2f (%.2f .. %.2f)\n",
		(double)tt[NUM / 2] / 6000.0,
		(double)tt[NUM / 10] / 6000.0,
		(double)tt[(9 * NUM) / 10] / 6000.0);
	fflush(stdout);
	return (uint32_t)y.w[0];
}

static uint32_t
speed_sqr(void)
{
	size_t u;
	uint64_t tt[100];
	gf x;

#define NUM   ((sizeof tt) / (sizeof tt[0]))

	init_gf_cycles(&x);
	for (u = 0; u < 2 * NUM; u ++) {
		uint64_t begin, end;
		int i;

		begin = core_cycles();
		for (i = 0; i < 1000; i ++) {
			gf_sqr(&x, &x);
			gf_sqr(&x, &x);
			gf_sqr(&x, &x);
			gf_sqr(&x, &x);
			gf_sqr(&x, &x);
			gf_sqr(&x, &x);
		}
		end = core_cycles();
		if (u >= NUM) {
			tt[u - NUM] = end - begin;
		}
	}
	qsort(tt, NUM, sizeof tt[0], &cmp_u64);
	printf("gf_sqr:                  %9.2f (%.2f .. %.2f)\n",
		(double)tt[NUM / 2] / 6000.0,
		(double)tt[NUM / 10] / 6000.0,
		(double)tt[(9 * NUM) / 10] / 6000.0);
	fflush(stdout);
	return (uint32_t)x.w[0];
}

static uint32_t
speed_sqrt(void)
{
	size_t u;
	uint64_t tt[100];
	gf x;

#define NUM   ((sizeof tt) / (sizeof tt[0]))

	init_gf_cycles(&x);
	for (u = 0; u < 2 * NUM; u ++) {
		uint64_t begin, end;
		int i;

		begin = core_cycles();
		for (i = 0; i < 1000; i ++) {
			gf_sqrt(&x, &x);
			gf_sqrt(&x, &x);
			gf_sqrt(&x, &x);
			gf_sqrt(&x, &x);
			gf_sqrt(&x, &x);
			gf_sqrt(&x, &x);
		}
		end = core_cycles();
		if (u >= NUM) {
			tt[u - NUM] = end - begin;
		}
	}
	qsort(tt, NUM, sizeof tt[0], &cmp_u64);
	printf("gf_sqrt:                 %9.2f (%.2f .. %.2f)\n",
		(double)tt[NUM / 2] / 6000.0,
		(double)tt[NUM / 10] / 6000.0,
		(double)tt[(9 * NUM) / 10] / 6000.0);
	fflush(stdout);
	return (uint32_t)x.w[0];
}

static uint32_t
speed_halftrace(void)
{
	size_t u;
	uint64_t tt[100];
	gf x;

#define NUM   ((sizeof tt) / (sizeof tt[0]))

	init_gf_cycles(&x);
	for (u = 0; u < 2 * NUM; u ++) {
		uint64_t begin, end;
		int i;

		begin = core_cycles();
		for (i = 0; i < 1000; i ++) {
			gf_halftrace(&x, &x);
			gf_halftrace(&x, &x);
			gf_halftrace(&x, &x);
			gf_halftrace(&x, &x);
			gf_halftrace(&x, &x);
			gf_halftrace(&x, &x);
		}
		end = core_cycles();
		if (u >= NUM) {
			tt[u - NUM] = end - begin;
		}
	}
	qsort(tt, NUM, sizeof tt[0], &cmp_u64);
	printf("gf_halftrace:            %9.2f (%.2f .. %.2f)\n",
		(double)tt[NUM / 2] / 6000.0,
		(double)tt[NUM / 10] / 6000.0,
		(double)tt[(9 * NUM) / 10] / 6000.0);
	fflush(stdout);
	return (uint32_t)x.w[0];
}

static uint32_t
speed_div(void)
{
	size_t u;
	uint64_t tt[100];
	gf x, y;

#define NUM   ((sizeof tt) / (sizeof tt[0]))

	init_gf_cycles(&x);
	init_gf_cycles(&y);
	for (u = 0; u < 2 * NUM; u ++) {
		uint64_t begin, end;
		int i;

		begin = core_cycles();
		for (i = 0; i < 200; i ++) {
			gf_div(&x, &x, &y);
			gf_div(&y, &x, &y);
			gf_div(&x, &x, &y);
			gf_div(&y, &x, &y);
			gf_div(&x, &x, &y);
			gf_div(&y, &x, &y);
		}
		end = core_cycles();
		if (u >= NUM) {
			tt[u - NUM] = end - begin;
		}
	}
	qsort(tt, NUM, sizeof tt[0], &cmp_u64);
	printf("gf_div:                  %9.2f (%.2f .. %.2f)\n",
		(double)tt[NUM / 2] / 1200.0,
		(double)tt[NUM / 10] / 1200.0,
		(double)tt[(9 * NUM) / 10] / 1200.0);
	fflush(stdout);
	return (uint32_t)y.w[0];
}

int
main(void)
{
	uint32_t x = 0;
	x ^= speed_mul();
	x ^= speed_sqr();
	x ^= speed_sqrt();
	x ^= speed_halftrace();
	x ^= speed_div();
	printf("%u\n", (unsigned)x);
	return 0;
}
