# The code that uses pclmul requires the compilation target to support
# that opcode and also at least SSE4.1 (AVX or AVX2, when supported,
# imply SSE4.1 support). Some possible option combinations are the
# following:
#
#   -mpclmul -mavx2
#       Enable pclmul and AVX2. This is experimentally good for
#       Skylake cores (and later).
#
#   -mpclmul -msse4.1
#       Enable pclmul and SSE4.1. The output will be able to run on
#       some CPUs that know SSE4.1 but not AVX2, e.g. the Intel Goldmont
#       and Tremont cores (from the Atom/low power line). It seems very
#       slightly slower on Skylake cores.
#
#   -march=native
#       Instruct the compiler to enable opcodes that are supported on the
#       current system (the one where compilation happens). Strangely,
#       on an Intel i5-8259U "Coffee Lake" CPU, with Clang 14.0.0, this
#       yields slightly worse performance than '-mpclmul -mavx2'.
#
# Add '-DGF233_PCLMUL=1' to force use of the pclmul code, bypassing the
# autodetection (of course, this will compile only if the compiler actually
# supports the relevant intrinsic functions). Similarly, '-DGF233_PCLMUL=0'
# disables the pclmul code, even if it would be supported on the target
# architecture.
#
# There are two implementations of the inversion/division in the finite
# field. Normally, the one based on the Itoh-Tsujii algorithm is used;
# if setting '-DGF233_INV_FLT=0', then an alternate implementation is used
# that relies on a binary GCD variant (Brunner-Curiger-Hofstetter). The
# GCD variant is slower, but yields a smaller output binary because it does
# not need the large tables of constants (about 15 kB) used in Itoh-Tsujii.
# This setting has no effect on the non-pclmul code, since the latter
# implements only Itoh-Tsujii.

CC = clang
CFLAGS = -Wall -Wextra -Wundef -Wshadow -O2 -mpclmul -mavx2
LD = clang
LDFLAGS =
LIBS =

.POSIX:

all: test_xs233 speed_xs233 test_gf speed_gf

clean:
	-rm -f test_gf speed_gf test_xs233 speed_xs233 xsb233.o xsk233.o test_gf.o test_xs233.o speed_gf.o speed_xs233.o

test_gf: test_gf.o
	$(LD) $(LDFLAGS) -o test_gf test_gf.o $(LIBS)

speed_gf: speed_gf.o
	$(LD) $(LDFLAGS) -o speed_gf speed_gf.o $(LIBS)

test_xs233: xsb233.o xsk233.o test_xs233.o
	$(LD) $(LDFLAGS) -o test_xs233 xsb233.o xsk233.o test_xs233.o $(LIBS)

speed_xs233: xsb233.o xsk233.o speed_xs233.o
	$(LD) $(LDFLAGS) -o speed_xs233 xsb233.o xsk233.o speed_xs233.o $(LIBS)

xsb233.o: xsb233.c gf233.h gf233_frob116.h gf233_frob58.h xs233.h xs233_common.h
	$(CC) $(CFLAGS) -c -o xsb233.o xsb233.c

xsk233.o: xsk233.c gf233.h gf233_frob116.h gf233_frob58.h xs233.h xs233_common.h
	$(CC) $(CFLAGS) -c -o xsk233.o xsk233.c

test_gf.o: test_gf.c gf233.h gf233_frob116.h gf233_frob58.h
	$(CC) $(CFLAGS) -c -o test_gf.o test_gf.c

speed_gf.o: speed_gf.c gf233.h gf233_frob116.h gf233_frob58.h
	$(CC) $(CFLAGS) -c -o speed_gf.o speed_gf.c

test_xs233.o: test_xs233.c xs233.h
	$(CC) $(CFLAGS) -c -o test_xs233.o test_xs233.c

speed_xs233.o: speed_xs233.c xs233.h
	$(CC) $(CFLAGS) -c -o speed_xs233.o speed_xs233.c
