#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "xs233.h"

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

static size_t
hextobin(void *dst, size_t dst_len, const char *src)
{
	size_t k = 0;
	int z = 0;
	unsigned acc = 0;
	for (;;) {
		int d = *src ++;
		if (d == 0) {
			if (z) {
				printf("invalid hex: odd number of digits\n");
				exit(EXIT_FAILURE);
			}
			return k;
		}
		if (d <= 32 || d == '-' || d == ':') {
			continue;
		}
		if (d >= '0' && d <= '9') {
			d -= '0';
		} else if (d >= 'A' && d <= 'F') {
			d -= 'A' - 10;
		} else if (d >= 'a' && d <= 'f') {
			d -= 'a' - 10;
		} else {
			printf("invalid hex character: %c\n", d);
			exit(EXIT_FAILURE);
		}
		if (z) {
			if (k >= dst_len) {
				printf("hex decode: overflow\n");
				exit(EXIT_FAILURE);
			}
			((uint8_t *)dst)[k ++] = (acc << 4) | d;
		} else {
			acc = d;
		}
		z = !z;
	}
}

#define HEXTOBIN(dst, src)   do { \
		size_t h2b_len = sizeof(dst); \
		if (hextobin(dst, h2b_len, src) != h2b_len) { \
			printf("hex decode: unexpected length\n"); \
			exit(EXIT_FAILURE); \
		} \
	} while (0)

static void
check_eq_buf(const char *msg, const void *b1, const void *b2, size_t len)
{
	if (memcmp(b1, b2, len) == 0) {
		return;
	}
	printf("\nFAIL: %s\n", msg);
	for (size_t u = 0; u < len; u ++) {
		printf("%02x", ((const uint8_t *)b1)[u]);
	}
	printf("\n");
	for (size_t u = 0; u < len; u ++) {
		printf("%02x", ((const uint8_t *)b2)[u]);
	}
	printf("\n");
	exit(EXIT_FAILURE);
}

static void
check_b233_eq(const char *msg, const xsb233_point *p1, const xsb233_point *p2)
{
	uint8_t buf1[30], buf2[30];

	if (xsb233_equals(p1, p2) == 0xFFFFFFFF) {
		return;
	}
	printf("\nFAIL: %s\n", msg);
	xsb233_encode(buf1, p1);
	xsb233_encode(buf2, p2);
	for (size_t u = 0; u < sizeof buf1; u ++) {
		printf("%02x", buf1[u]);
	}
	printf("\n");
	for (size_t u = 0; u < sizeof buf2; u ++) {
		printf("%02x", buf2[u]);
	}
	printf("\n");
	exit(EXIT_FAILURE);
}

static const char *const KAT_XSB233_DECODE_OK[] = {
	"000000000000000000000000000000000000000000000000000000000000",
	"8940bd43927c8b2947b945b3186132156182142eae619a59d2e4485b7500",
	"8fccbab7e7cd1f22ed1e61fe62bbb3482ac9d8d9495ab940d50900673f01",
	"05d37ff3b74b7ce54600866bece0602ff3387f85e4719705a68be2fa5500",
	"d8c9bacab5039f5c9cf60f54852161944c0a03bd4573a1768d2ebe9d8e01",
	"7fb5e0348570597d57c1d2fca09c7f65b31983ab93022740027aa0efba01",
	"ed003b7ab53d132fc0ef1fa36e6535b47290184a3393911009982a819c01",
	"fc78318fc05ebaf2b047f3bc65ffe39b580615a9068496c1977f0bc9cc01",
	"30058d4b2ccd42614fb59b9ef03cf6cd355c140f83fdfe55fe8fb2666501",
	"a44b02f35d234068a67d8498bdcf81098baf54cb49268ce5cb9e2eef5901",
	"0013dffa8bfe74dfb0f54e79ad9b7ce0ce32ade12b58595c6587b67e4401",
	"55b6ad37c6fd5df261e6ac9b2d74536addfb66646bd8e8296a7f6324ce00",
	"47b3b75b1a995e96a672edb8c716839bb589212965e2deb03c7a6ccf6001",
	"82fdb5a8465dd6c59b10dde5432c353edc41577b124f8ef7d0a8f509c500",
	"558f4447507fb52cc08d4fb34f4ad4a2479a3d2246038fd9f0942c9fe600",
	"a24dd168d1d3f75f81a8b8afc31947bdaf166839d8195f9937ff8b581400",
	"bfea156043b9be7dc7c3cd410a0c38f0b2e9a1fd0f414a4fa9ea2cab1001",
	"a97e14af5cf95c8624ba422c57ecec7d898a1e6cc2a10a3e667472cf9e01",
	"5151b9560be3f48042f4adde5d603e0ece23d4d9d7072215b94a78b6ad00",
	"50c08e4f6ec9f4feb32152c19809879983422eb54200eb50792600d73901",
	NULL
};

static const char *const KAT_XSB233_DECODE_BAD[] = {
	"801fb075ad4b6708d7ca674fd29f45118a1e98f34672e6f416db57ba5a01",
	"3c098046550a679c35d1aedb20a4a5d52f4342ba1b79f0a7d2c27a77b600",
	"6771ec2e4a8409b328f591addc14a740502dc293888d761d2a2093b9d501",
	"b18bcf83211455363978c178dcf28c3e8d643d816124067a0c0b0438fd00",
	"7c17bbf7a92f393988465f1d7d1724c8b8510352b077efc82d00115c7400",
	"a260b4bafc91b68f9596cfa114a0813fa57397df479d6a751a204e2f8500",
	"88d7b4a8edd9640f151d9cbd95f77974ae46c79067be72cf9bd421528800",
	"2b19b2b57c4cee44c53c02bb7f2c542000e6f97266a47acd450b3223ff00",
	"97e0a30a19fc45a454d0d4fb703376c80d91751d50b16809e26e2fef3601",
	"646bc8733e38db5f6397f78698f3dc9931eec60344eec49be18fe96d6f01",
	"c01a93fb2fe1e0b9f2e61062a7ec0d5339071c2ce90aebeb10afd79c1600",
	"df7fe9d661debfb6b0de91b204c0dace85eea241d198594b509bb4e49000",
	"3b67cfc6f70edd422c3cc8855e41193563f94074558f5721dcc030b47e01",
	"5d95d619e23e2ffdbb6a780391ea8305e78755b8e80ec94ef71ed3a55001",
	"acdbb3d7959c65c0229a28908222dbc7332b3f643c036f7f2d76bc37bd01",
	"c16c24fdda4895754181cb844438d3a3672b799eeb6d603318c12ca83001",
	"9ef7db1233d0c69ae03b81d3637409487b80d5de2c797b3e9c4ec9bac501",
	"d9d2e06213b9a2873401f030f1800b0af22c41bfa7a5af04a9f5930bfc01",
	"ce36a5a61c4e8559700cbf3f3cf730a4dd663c14badf9b317acf1555af00",
	"1ad612531e5bb81d9c4dccc6725673a28fdf9ea787d5ec3f056c36f65900",
	NULL
};

/*
 * Each group of 6 values encodes points P1 to P6, with:
 *   P3 = P1 + P2
 *   P4 = 2*P1
 *   P5 = 2*P1 + P2 = P4 + P2 = P3 + P1
 *   P6 = 2*P1 + 2*P2 = 2*P3 = P4 + 2*P2 = P5 + P2
 */
static const char *const KAT_XSB233_ADD[] = {
	"02144880b80e145111038b896e3353909604af0b093c53dcb97957c85d01",
	"76d4b7e2c85b59dad633d0b193e4a081640bbcaec951ce756f19bce62c00",
	"9ab06eca75f4afdb3fd838f199ea2efccbedf23632dd2642e64d44c4dc00",
	"eb929217cae07eeabb3a2f73f0529b815c36a0b8dd622e7e112c39728a01",
	"8776ebf4b14066972fdbf06418e2e7e1ef30c7958c7508ee5f12dffd3401",
	"a27faa8d0117ca49b9effa1f4104a42bcb7848014f11c02ff4f4db9a3801",

	"4d58931f05fa526cd4f07de21c71e910cdfd8419a1c3ca455b788c68b201",
	"ca028772e5db4cd219459467d5ca68d27e631e6aa68addffc9834e4ba301",
	"c6f6f201503566437c881d3561f9f8e0029d15169a90852dcb96ef257b01",
	"0c63dc67e9e9087f3b1e857dd769615633eda9849a05409a12347c6c9e00",
	"ac89dd4af3a0208edbcc1e78ebd0eeb0ccd28d56329806f2209ed287a501",
	"8de7a7a4ecc484f8536983c2bb6b6dc483ef06adc103f7116360c07c4301",

	"2cebbc366c6c1f62ee690b370e970c1a718cf016c9d3b03c58485eca5f01",
	"a61388f4a96f1fbe5f19a43ae9a06794fb4120abd03d3116ad434cb7b400",
	"3d033c423aec126e2ccf21eeff0840c413a6054e57a84c5ef4a87811e100",
	"df66d1ff0019d8a0e3cb613caec460c82ef8059039c9a64216a2f9b0e501",
	"48f09a5cefe7b7d56cbf3ed5d19d17d84888b9fcc17d16b8ad3d189c3401",
	"216ae8510928fc21c4bcc9c152fae5be65f2423c202af94482d45e508c00",

	"b1b346bee688b3dcb279315e0f64ad7b71abb413ff26b9cf54accd361300",
	"40af3372e5e28d72bbee0208b9e088407bf35764d218d5c906ff0df98b01",
	"382bc13161d0e039aad71378c0d35f59ee67e969a18c0c1ceba5fc9f2900",
	"8001a637c126818815a11d801ef488777f2d2936a2a1bece94c52abecc00",
	"a569aa7b2165a7144aeb49000c0f1f74e6cbd020d24fe995bad98cc7ad00",
	"a80950085f1b28aaada3753a99f97cd5d5180ffe5d4f7d456efccc9f9600",

	"11e4433fc1402ed10029d41212289ea399d41db9d40a388437ee8a5ebe01",
	"436dee22637a2870237e450541a21828249948585adbd2259a596a22fa01",
	"cd6425529751c0b87b95181b2698460693295760d99713c1abd5cb6d8401",
	"75314bc24edd1afb65bbf2eac38fe8ef1f144fa3930b48de2ea895cb3001",
	"123ce9b917e18a032cc3721f3b1850a16acd4d71bc0f7e6987797645fa01",
	"2b9a20351e0ce435886d22690ddeadf1d2303f38d8737e4471dbc7681c00",

	"e6396f9bcff1106c56e83bd57c68e0c3ce272bcf5f9a6c14e7ffb6e7ca01",
	"c10a69de3fa6394dbbab9128afb2a23fd3cfdbf8f8cf6c456e4710c97300",
	"07ceb73e995462cf4b10979e307d6b8b0d5051a783c6f3ecf4513d7e2e00",
	"ffa05d8cbf717b21e5da33a87c2da01df0f0147861624a93a3f8e62aa601",
	"04fc2cb7f8ee34db20c977ca5e5582581fe515cf4eb8a7ccd7401838a401",
	"f4718c6b783852485210420262b5974e3e635a13a805222d61a072a5d501",

	"7dd5eed5adae63360d282e011ec3596c625a5f91b91dbd85209b550a1000",
	"c91d0d0367efe902b1af07496984acb6c9daa3f973714fd9a5c68741c400",
	"4d4c928d5963c2013f7484eecad069bcc78fb4216009a3b0cc556e59d601",
	"5d902227902048774c922bbaeaa456ec52f4983f56d1d641f96b378ec801",
	"9e4aeed2b6bd766dc5873de2a6c1d93f5a262113a51f0582b460c107af00",
	"51597210f9392667c38a1e4942a7b0b83782288a245336f1606c51924700",

	"3b47769eb490fc588f4251ee38bc8c70423392d6c2381c8803924aa89e01",
	"c33240af7431a1d9cf1bc45d4737f35a9a2b8b8473d6db8459de7b421000",
	"674b3829a9c85ca207e710d762460bef6ddfc92bc47cfde4ddf6206fd901",
	"a37ac566fef860a5a55e1db8599a95867cdf3161726e2f3828d758b0b901",
	"99e73c60a1a9501fcff7b63017cbe8303087a9c12491bfedffa583e48e00",
	"b1c1207cc51c742a7c55afcedcb10fe78951573361e7115e5d1e55522b00",

	"2a852fede5ff207503879030383a20829439f538a7844f61b0e49cee7000",
	"537dea25a66b6995f631b9688117db9470ff4e1032358f9a9020a4f63500",
	"3fa567c96910c09678f9b7609e0aec6f1a40d080d9d2166ad9ef46881701",
	"d7a17d8bd7f84ebde0ddef19beb1165142c4519114e454f1d8017322d200",
	"c551dc502d0915c04459999f5e37ec4a915284371a61396435086ab49201",
	"6444ccb3dc00bdf2103d46e663778b5514ff1d7a848c2df9f3a96cdfe900",

	"77ba1b3758f1ee9efd99a6059224782ba1b92fb00beb3db53d984635d901",
	"fd5029736e2062ad024ccfb5e0e63c6489b537f3587c0c4d8bd47a5a2d00",
	"8036e6f8f53eedf796fa59fae533b743f0c11b955eb8e68c71105bc76600",
	"633c9c5da54168013f4a937aa1d0dc346c579def1a7e54484efa45192301",
	"cfa8569d09b38a25129c61e552e60093d078065901ec007bd422649fdf01",
	"d8a82435a2806938e5f69adbbdd930efc2603e898182d4c1f289127bc300",

	NULL
};

/*
 * Each group of three values is:
 *   point P
 *   scalar n
 *   point n*P
 */
static const char *const KAT_XSB233_MUL[] = {
	"106119152d4352ee7b4f5477bf024d5b03652787ce6d65a1873378cf3c00",
	"e9b6249bbbe83e5acdc45173c92421b785fe8d108c48468f0731fd43a834",
	"dc7990d625cde144b98d41b5bc7745c1f46e07592ce0d7c684b6e8d86500",

	"462375f85cbea6f0c9b6a966e6c6605e0b52c1020263b64deafaf5675801",
	"240944c3a4fb54f9cadd6a182e19715a62e62161632f0ca913540fb58f69",
	"e4a62e179f93d992669c8a763f368b4cd149bac0210ae91c92d66d361001",

	"85f70654a469dcff1802afcf7f65178a9a8d4b9093a076c10fb6a1eb6f01",
	"5c42bdc259deb257843bcae0eaa0d9b96bb43688f2319776f5876a651b03",
	"81ca9c26fa51bdfc603e7d66397e91f1be95ac5467a45eb747b1c21e4700",

	"f514c5c2ff9b864f1ed7b5cd1fff963356656e06f68c29f076a73a661a01",
	"c6f34cd3e3c369e86bd2f3792a3dfab2a51902c36eab85ca0af1779b88f9",
	"107772b42be2ff485d96849e7a0065bb4c6a57759c00d0ec6443dfd53501",

	"ef0469fb216d019185491789b5f2f44b69fd3934473e3cfda7f59d16a501",
	"78def7a7ac41131b0d843339b8af2d0107067d350c75f77fc59709ca6b92",
	"90e1bb7a8f748c09e95224aef5a110adf8c6c058e89167f9d19c68edc901",

	"baff096082f63d76bf21418783b1107529290a640569abccd9178e894900",
	"ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
	"507e493c422fededb0b0128f1d6cf6f295cadfca76c91702eb0260cf4701",

	"baff096082f63d76bf21418783b1107529290a640569abccd9178e894900",
	"d6e0cf03261d0322698a2fe774e913000000000000000000000000000001",
	"bbff096082f63d76bf21418783b1107529290a640569abccd9178e894900",

	"baff096082f63d76bf21418783b1107529290a640569abccd9178e894900",
	"d7e0cf03261d0322698a2fe774e913000000000000000000000000000001",
	"000000000000000000000000000000000000000000000000000000000000",

	"baff096082f63d76bf21418783b1107529290a640569abccd9178e894900",
	"d8e0cf03261d0322698a2fe774e913000000000000000000000000000001",
	"baff096082f63d76bf21418783b1107529290a640569abccd9178e894900",

	"000000000000000000000000000000000000000000000000000000000000",
	"f26efbaa86ce71af4affd3d368abb6ea372d3b426d617313b1916fb48cba",
	"000000000000000000000000000000000000000000000000000000000000",

	NULL
};

static void
test_xsb233_encode(void)
{
	printf("Test xsb233 encode: ");
	fflush(stdout);

	for (int i = 0; KAT_XSB233_DECODE_OK[i] != NULL; i ++) {
		uint8_t buf1[30], buf2[30];
		xsb233_point p;

		HEXTOBIN(buf1, KAT_XSB233_DECODE_OK[i]);
		if (xsb233_decode(&p, buf1) != 0xFFFFFFFF) {
			printf("\ndecode failed\n");
			exit(EXIT_FAILURE);
		}
		if (i == 0) {
			if (xsb233_is_neutral(&p) != 0xFFFFFFFF) {
				printf("\nneutral is not neutral\n");
				exit(EXIT_FAILURE);
			}
		} else {
			if (xsb233_is_neutral(&p) != 0x00000000) {
				printf("\ndecoded to neutral\n");
				exit(EXIT_FAILURE);
			}
		}

		xsb233_encode(buf2, &p);
		check_eq_buf("decode/encode", buf1, buf2, sizeof buf1);

		/*
		 * If we put a non-zero bit in the reserved positions of
		 * the last byte, decoding should fail.
		 */
		for (int j = 1; j < 8; j ++) {
			buf1[29] |= 1 << j;
			if (xsb233_decode(&p, buf1) != 0x00000000) {
				printf("\ndecode should have failed (1)\n");
				exit(EXIT_FAILURE);
			}
			if (xsb233_is_neutral(&p) != 0xFFFFFFFF) {
				printf("\ndid not get neutral (1)\n");
				exit(EXIT_FAILURE);
			}
		}

		printf(".");
		fflush(stdout);
	}

	printf(" ");
	fflush(stdout);

	for (int i = 0; KAT_XSB233_DECODE_BAD[i] != NULL; i ++) {
		uint8_t buf[30];
		xsb233_point p;

		HEXTOBIN(buf, KAT_XSB233_DECODE_BAD[i]);
		if (xsb233_decode(&p, buf) != 0x00000000) {
			printf("\ndecode should have failed (2)\n");
			exit(EXIT_FAILURE);
		}
		if (xsb233_is_neutral(&p) != 0xFFFFFFFF) {
			printf("\ndid not get neutral (2)\n");
			exit(EXIT_FAILURE);
		}

		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

static void
test_xsb233_add(void)
{
	printf("Test xsb233 add: ");
	fflush(stdout);

	for (int i = 0; KAT_XSB233_ADD[i] != NULL; i += 6) {
		uint8_t buf1[30], buf2[30], buf3[30];
		uint8_t buf4[30], buf5[30], buf6[30];
		xsb233_point P1, P2, P3, P4, P5, P6;
		xsb233_point Q, R;

		HEXTOBIN(buf1, KAT_XSB233_ADD[i + 0]);
		HEXTOBIN(buf2, KAT_XSB233_ADD[i + 1]);
		HEXTOBIN(buf3, KAT_XSB233_ADD[i + 2]);
		HEXTOBIN(buf4, KAT_XSB233_ADD[i + 3]);
		HEXTOBIN(buf5, KAT_XSB233_ADD[i + 4]);
		HEXTOBIN(buf6, KAT_XSB233_ADD[i + 5]);
		if (!(xsb233_decode(&P1, buf1)
			&& xsb233_decode(&P2, buf2)
			&& xsb233_decode(&P3, buf3)
			&& xsb233_decode(&P4, buf4)
			&& xsb233_decode(&P5, buf5)
			&& xsb233_decode(&P6, buf6)))
		{
			printf("\nKAT decode failed\n");
			exit(EXIT_FAILURE);
		}
		Q = P1;

		if (xsb233_equals(&P1, &P2) != 0) {
			printf("\npoints should not be equal\n");
			exit(EXIT_FAILURE);
		}

		xsb233_add(&Q, &P1, &P2);
		check_b233_eq("add1", &Q, &P3);
		xsb233_neg(&R, &P1);
		xsb233_add(&Q, &Q, &R);
		check_b233_eq("add2", &Q, &P2);
		xsb233_sub(&Q, &P3, &P1);
		check_b233_eq("add3", &Q, &P2);
		xsb233_double(&Q, &P1);
		check_b233_eq("add4", &Q, &P4);
		xsb233_add(&Q, &P3, &P1);
		xsb233_add(&R, &P4, &P2);
		check_b233_eq("add5", &Q, &P5);
		check_b233_eq("add6", &R, &P5);
		xsb233_double(&Q, &P3);
		xsb233_add(&R, &P5, &P2);
		check_b233_eq("add7", &Q, &P6);
		check_b233_eq("add8", &R, &P6);
		check_b233_eq("add9", &Q, &R);

		for (int j = 0; j <= 10; j ++) {
			Q = P1;
			for (int k = 0; k < j; k ++) {
				xsb233_double(&Q, &Q);
			}
			xsb233_xdouble(&R, &P1, j);
			check_b233_eq("add10", &Q, &R);
		}

		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

static void
test_xsb233_mul(void)
{
	printf("Test xsb233 mul: ");
	fflush(stdout);

	for (int i = 0; KAT_XSB233_MUL[i] != NULL; i += 3) {
		uint8_t buf1[30], buf2[30], buf3[30];
		xsb233_point P, Q, R;

		HEXTOBIN(buf1, KAT_XSB233_MUL[i + 0]);
		HEXTOBIN(buf2, KAT_XSB233_MUL[i + 1]);
		HEXTOBIN(buf3, KAT_XSB233_MUL[i + 2]);
		if (!(xsb233_decode(&P, buf1) && xsb233_decode(&R, buf3))) {
			printf("\nKAT decode failed\n");
			exit(EXIT_FAILURE);
		}
		xsb233_mul(&Q, &P, buf2, sizeof buf2);
		check_b233_eq("mul 1", &Q, &R);
		xsb233_mul_ladder(&Q, &P, buf2, sizeof buf2);
		check_b233_eq("mul_ladder 1", &Q, &R);

		for (int j = 29; j >= 0; j --) {
			buf2[j] = 0;
			xsb233_mul(&Q, &P, buf2, sizeof buf2);
			xsb233_mul(&R, &P, buf2, j);
			check_b233_eq("mul 2", &Q, &R);
			xsb233_mul_ladder(&R, &P, buf2, j);
			check_b233_eq("mul_ladder 2", &Q, &R);
		}

		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

static void
test_xsb233_mulgen(void)
{
	printf("Test xsb233_mulgen: ");
	fflush(stdout);

	for (int i = 0; i < 20; i ++) {
		uint8_t buf[32];
		blake2s_context bc;
		xsb233_point Q, R;

		blake2s_init(&bc, 32);
		blake2s_update(&bc, "test_mulgen", 11);
		buf[0] = (uint8_t)i;
		buf[1] = (uint8_t)(i >> 8);
		blake2s_update(&bc, buf, 2);
		blake2s_final(&bc, buf);
		for (int j = 30; j >= 0; j --) {
			xsb233_mulgen(&Q, buf, j);
			xsb233_mul(&R, &xsb233_generator, buf, j);
			check_b233_eq("mulgen", &Q, &R);
		}

		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

static void
check_k233_eq(const char *msg, const xsk233_point *p1, const xsk233_point *p2)
{
	uint8_t buf1[30], buf2[30];

	if (xsk233_equals(p1, p2) == 0xFFFFFFFF) {
		return;
	}
	printf("\nFAIL: %s\n", msg);
	xsk233_encode(buf1, p1);
	xsk233_encode(buf2, p2);
	for (size_t u = 0; u < sizeof buf1; u ++) {
		printf("%02x", buf1[u]);
	}
	printf("\n");
	for (size_t u = 0; u < sizeof buf2; u ++) {
		printf("%02x", buf2[u]);
	}
	printf("\n");
	exit(EXIT_FAILURE);
}

static const char *const KAT_XSK233_DECODE_OK[] = {
	"000000000000000000000000000000000000000000000000000000000000",
	"478f05b8ff56f97ddcdcdda13b42a773be0c1875bff3e555702e2eb1e900",
	"1123b0cbe6de2d68004ae486854ad503c19853bfd3e6c344f572eb66c001",
	"54df1f3c2160b43bff5440ef33fb35bf1740d9188003cfe82f3bede61d00",
	"9e9230101e312faf1e33cd4fd3dd22b01b0a1be8aa1859b9a1efe8ce9e00",
	"1c9cd9aa57f085717b5c977e9e57917eced2f3cf466b264f656f304ad400",
	"0a8c3050729ca55f86095397b4fc26e4b60483f0e89f0a4e4df71030be01",
	"25039e0e09f5c82541c9e1712474c03392786349162853068150d4549a00",
	"71d29dca39024af194526d9902381d9c6c5cdb2f71ab146be5aed028cc01",
	"e97d39f48852edf98b5a1e712599d277743b68cee5680ebc0f1293f60b01",
	"7b7bab22f163a3aed122145c2cb961efc7b8baae8664faf2ea164b752401",
	"ed874112174ec1a2af4cd3e7c179cff6b2cdd92108883dd800f63d0cbb01",
	"64886d03b7e29e43b2a09c2127ac1c6303051dc8dfe60f41d9b61bd66100",
	"c84a0bb20b4b4940f6d7460f08d9c19c3234a92f96d7938b05ef894f1501",
	"f6c9de32061d6e6747cfdaa36987c7df2a5d96814243af976eafaf447501",
	"161d948d488bd26b3d81beca85f8fc96bb1a36e98e78fcd3c3cb90e3c801",
	"be4339d93b7b5dd5e0dfecaa2e309f3525c8dc7c8090ab618a73c7a95001",
	"9023f2233584b06c3a8e958d8569ec3cac91ca4d7c67516ff74e66d4a200",
	"0d37b500c490e66f1e3646b692579295c0b38dc3a1cc261cdf71b696b000",
	"23b6566c74a3de5cef968835f901b7c218108ec6106096ac8c07afadda00",
	NULL
};

static const char *const KAT_XSK233_DECODE_BAD[] = {
	/* w = 1 yields w^2 + w = 0, which is invalid */
	"010000000000000000000000000000000000000000000000000000000000",

	/* trace(b/(w^2 + w + a)) = 1 */
	"5d7ef866575451e3c90a3ed86c19e435e90a0506ec3e56a1d2070c3b6d01",
	"38d0a8b6e69ec6f70872f4a7a235531a696f4b618f58b6c2a1ea40063400",
	"fbc1529dc5347c086493da1a6505115f8e080feebdaa83fc491774d28400",
	"9262fc6b097a43b8a89a53bbf901ad20d41e29aa754881ed71b74bb12e00",
	"1f3096dcb03ffd92d574ee74d1019adaba91aeaadb5bf69da10e83788701",
	"caa9c228ef23944f86557b325488cca9482d274a9bf024b1bd9038d9ff00",
	"12b4c88a9fe72ae9734f17a8df0a017fd29b549e0adb45336897a8c85e00",
	"649e4366355f8a27fff2cd284eca895881fc1fcf798198a572d3be884201",
	"c0bf607a047e032a878a90601c65f139177bcb8d6c7ac25c100b9a383b00",

	/* w maps to points that cannot be halved */
	"5d3f829d1618989458e400e15aeb484c0a4ba3359b7b1f72a6d9ae442d00",
	"6cdd47f7c7d936f068e8e203c343a075d5306a443b1bb5d1c9df6e11d801",
	"eb86e62958e729cfeae85c32e896469a191fa88a019a6f73138d73a73c00",
	"61cf5b9ee3e3a2529f71d8842cd4cc6bc5d7d12c0defae81739409bb2500",
	"7ff9acb794afb0a1f9511af0ca6314a2502ec3c9702cbe00f39f44677e01",
	"b89fc50d28755305cd98f8bb3882178ac8cf231c58e5c3b8217006c9a201",
	"41b08b4fa81064e52e5969eceff573f083f1f9d27ce1013b7a115a6a6100",
	"952722c25a7128a8fa36d600b1b0583fafe333406dda526dd476e3b50600",
	"763d8e4fe9d109435136f92a6d2f3abb9634c57ff4687ef878717a3f2601",
	"4ab33464767e663f3ed4194f606bed91005f7584195c53dd3c4fbc71d201",

	NULL
};

/*
 * Each group of 6 values encodes points P1 to P6, with:
 *   P3 = P1 + P2
 *   P4 = 2*P1
 *   P5 = 2*P1 + P2 = P4 + P2 = P3 + P1
 *   P6 = 2*P1 + 2*P2 = 2*P3 = P4 + 2*P2 = P5 + P2
 */
static const char *const KAT_XSK233_ADD[] = {
	"ff32132133ff22c3e309ac2a5744ffbe7a745452f9cbc1852591020ef201",
	"285f84c92fe54eba84228362c7ec3b3dc974d3f149d41c53676f0f884c01",
	"b4f7ae8f595620e9993ecec609b220f5a6e55dbcdb72aa193c322a5e5e01",
	"90903e1fb9a714bb872284af141ccb53fe5f36a3b63653e2fdb1f0dd4f00",
	"e6199aa37bb1564c57a0cc3ccd32b921763d724751d98d5917b437022901",
	"a2e81822d3bd47b57ab65618236d6bab4fb9ea267dfd859b409e6292ff00",

	"8f008e388a7dadcfaf7561c973461c0ab20ca1e112bab7d1849ae019ff00",
	"6fddeb568bd884d259b1752fa4075cda41d7858164c48ac2c631ce8b4d01",
	"deb094684ba2425de6145faaf7ff5465b4ec1b23e564d34b8e7464047c01",
	"8342a39c351a804e2a8035ccc9ed27117b7759ac0cb508f6619f48e04b01",
	"50c1897f1cd8e7882e7ac7216a42a6f197221beaf27c6b299b234d073a01",
	"3b288587d25d167c4f61438aca41484608ced390bc5141e7507d17cbbb00",

	"1ab150092668e46283529c2377bfead27461e4ae986039b27bdc892ba000",
	"cf6b27d936d46e2c36bf70f3bf67da166b3a52a292864c369153f84a3f01",
	"fb95dcdcc07046fb561c3c11b9f86364c7defbec037c862d472b72e81500",
	"ee418e647aeb98e694145b3a7858f632f14ec3529375be0ddbd8d11dc501",
	"6d924c79d45f4baede7a60f9c99cd6a99dbbee778f900ac906e8e7c83601",
	"8894f5b516b80f2059100086f3e0e1d300a62cbd5879ee5f871d9d2ef700",

	"bcdd6097170e80beecdf13c9cf01ae15547caf56c5aa8f2fad26864c6100",
	"19e758fc58ed0e700528a1f10f112d4ba8992a1ed60b8a4af76943f22201",
	"fd018b5f0f59dfecab12723018df08705461a3fa2d703cdf7b2696413600",
	"19bdb66fe103ac234b9d3e692df2784bb2ed3b4482485def2bfaf518ac00",
	"3cdd7b97320f9c8f74eb682f4f11a88caa188cdda4e632d3d5748c168700",
	"f5df927b8b1478078e4284b314adfede93412a8a14530a7c115288fd9101",

	"a25a5d30d88913c6d9c510f60b196ae4b347b76bb8b92fe4bc94d9da1900",
	"3998c97a263bd459c126b498c58f0739605c66bd158d1ca9a45abcd30000",
	"dbcef0d1e6e79e54f52193717f670ca85bc7a0ca7d75e8e74642c2454700",
	"49de384e8fe7764bef3a4489af3199472171962016be84759a6c97e3f801",
	"3e2773a25553e35ca382d1c58dc29fa9d148226b980f89827bc3d3887b00",
	"8cf73f76eb32f4dd704119a2527faffc68721fb8e0e94c7bef5202c96100",

	"327879816ab3c276a354ca8e33cc4643b1632496c3df17abb9809a46a701",
	"0af4619c45e23ad644aaf8b5d2a90c8675fd512bb90b492ef88df1e59a00",
	"2f2c5f1ad29893c9ffc669821c5465ebfc8cfedba821da069354ea410401",
	"66f95870eb0c28c92d34dcce95e4b0c12f0c9da13200f93bd4f99986d400",
	"1e14eb7e71e0e52d48e584cfd08ba04f4287ad385240388f1e5758db5900",
	"5701404da4c7e44038c147e4a974471a289bc0cbe8a1780e45969bfb4201",

	"76a96398afbf6be212287eda5671ef8e62116627d9646c0f04081f0f0400",
	"31327483f6edba19761e231d02b53d871c732cd3598f4c8139f716224901",
	"5b85525d41b98e127789fbd15913a3497296d398d558ee28d89df1303401",
	"443d0d33e410953f38ee1e5dfa9666da9bb7bb5a0ed99659d40bd957d600",
	"56de350c2a1c646ae007ba0ff9197bcaf9a664587d9bb2e004489f752201",
	"973d6c27d75380d910dbf25ddf40b0fadc49e1f1a44fa22c4309939d7e01",

	"d7204c14ef665860bef3aaae0fdd50bb6adab7586b76e490c9959a741801",
	"a88e60081f9c9bb3fa7f559e0e02c60b296fb133700bb9e87f688c929700",
	"5ecb998965809d94752288945a8079348f39e37f998d853312b64a77e001",
	"a662c8f83e0e0834a92145a8fa3a3a4031abc5fc73587d4a182f7fff6a00",
	"e2e68ae840d1bfe0646adae98b854f6362c603ad633d83d0fabc8eca1401",
	"3164384a9b38cd09c4d3cad52d58db04dc6127d5fc730ff845156c715c00",

	"d600c7021215c93bbb5074e5ce996ae257078fcdcb79ea614f775d6c0600",
	"5a11626fe4d1dc19afd919c409f6c4fa4f60330ebe6db45e0be25da27401",
	"421381d35463de832e33b077f924e2782e787215418129531d4b83bebb01",
	"1d3e836bc15ea3992f7b6997648b3f4f9eadb2887a38c0852eb89250e401",
	"b2bd23c33575b53b4557d6249ceba44c41efe8de213daef4700c352dab01",
	"819cce21e636b7cc6086aa83cc2e0427921dfcc8364c30a7a077021eae00",

	"2f7dd87f5fd087c349a72d4649f81dbd1f88fc439fc315b0918ac6301501",
	"d3e912c1dfde069bdf37d00ac1932b6fc328c0903d4e15f657407c486500",
	"580bd418c3ee71b271b5024512d4656da8b61ef3e5bd821a9aa133ffe401",
	"00bd8ca49e5be45f9cfa5bf0e22e9e9bb618178e3858b8dbee6659229601",
	"8a3e021f959e3a793e228032fd2330d784ae98b4cc59982675e06c682701",
	"a37b3fadc548f8e54f148a4301e7504a45fd828ef67da67c3851889e2e00",

	NULL
};

/*
 * Each group of three values is:
 *   point P
 *   scalar n
 *   point n*P
 */
static const char *const KAT_XSK233_MUL[] = {
	"65405f2f0a664780b6858a0bc2c458bd3b97eb88df0529c497903eae8601",
	"4e85433e5e2d598f0f644af9cdda327ad9d1eb2f1f00ad73a5d47cf20d5d",
	"87e1fe594929d24b4c5405f5ec03ffedd0d7e30273951ab89a69a06ebb01",

	"8447778f2982c11888b30bbd5a5418a5a5e309519a32872463674e5e2600",
	"6f92c6bc715b2d69f9ee0b55ccf53cbd86c9e5f3af0a5714992401d700df",
	"1a2fe6248157e11eab2141d8755ec880d528460aa756008990aa928fca01",

	"ca8c11ffb1e8196df9f55d9b9c00cc05ff609c19a57d1da59e3d842a7100",
	"81b60819e74f7c43f9d55956cb9af346c791ab634b9e91707c426cd58d8d",
	"7dd17870ceeb3f3f221ef20707b3e2ac5fa6e16139e92ef85b514940e700",

	"88c4737cc61c737d9a3177f4f63804558685b8518a9b4899e822bd3acf00",
	"92b9ee4dd598039761281aac140fd8802b600b0ea98d82c8567d2ee8bfa0",
	"b3e29bd8a8b5b1fb76bcdea37d3fd78ab409d0d982532c7b9f8079a9ba01",

	"cf83e5135b57ded363e44c88cb1d5925a9f18332308d636938357191f500",
	"4c9c704abe0ee618f255596decaf323f9c3c7d6b35d1f0ef9f7f3a963c69",
	"94a8a5040bb6ae2d7bb12a44e677bb6fa843f46c4473a4a8a399391d3a00",

	"08b0df6178ed276829b40449331bc7106dca3475ed6e296305a0017da500",
	"ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff",
	"b8dbd51d8cb04149bcd55196901d77767fca295c87a80a8bb0dfd3abf401",

	"08b0df6178ed276829b40449331bc7106dca3475ed6e296305a0017da500",
	"deab73f1d51afb6ed4bc15b95b9d06000000000000000000000000008000",
	"09b0df6178ed276829b40449331bc7106dca3475ed6e296305a0017da500",

	"08b0df6178ed276829b40449331bc7106dca3475ed6e296305a0017da500",
	"dfab73f1d51afb6ed4bc15b95b9d06000000000000000000000000008000",
	"000000000000000000000000000000000000000000000000000000000000",

	"08b0df6178ed276829b40449331bc7106dca3475ed6e296305a0017da500",
	"e0ab73f1d51afb6ed4bc15b95b9d06000000000000000000000000008000",
	"08b0df6178ed276829b40449331bc7106dca3475ed6e296305a0017da500",

	"000000000000000000000000000000000000000000000000000000000000",
	"ac47ee67e70933727c304bae789bd08aff54ff56d5240f4996ef21976a58",
	"000000000000000000000000000000000000000000000000000000000000",

	NULL
};

static void
test_xsk233_encode(void)
{
	printf("Test xsk233 encode: ");
	fflush(stdout);

	for (int i = 0; KAT_XSK233_DECODE_OK[i] != NULL; i ++) {
		uint8_t buf1[30], buf2[30];
		xsk233_point p;

		HEXTOBIN(buf1, KAT_XSK233_DECODE_OK[i]);
		if (xsk233_decode(&p, buf1) != 0xFFFFFFFF) {
			printf("\ndecode failed\n");
			exit(EXIT_FAILURE);
		}
		if (i == 0) {
			if (xsk233_is_neutral(&p) != 0xFFFFFFFF) {
				printf("\nneutral is not neutral\n");
				exit(EXIT_FAILURE);
			}
		} else {
			if (xsk233_is_neutral(&p) != 0x00000000) {
				printf("\ndecoded to neutral\n");
				exit(EXIT_FAILURE);
			}
		}

		xsk233_encode(buf2, &p);
		check_eq_buf("decode/encode", buf1, buf2, sizeof buf1);

		/*
		 * If we put a non-zero bit in the reserved positions of
		 * the last byte, decoding should fail.
		 */
		for (int j = 1; j < 8; j ++) {
			buf1[29] |= 1 << j;
			if (xsk233_decode(&p, buf1) != 0x00000000) {
				printf("\ndecode should have failed (1)\n");
				exit(EXIT_FAILURE);
			}
			if (xsk233_is_neutral(&p) != 0xFFFFFFFF) {
				printf("\ndid not get neutral (1)\n");
				exit(EXIT_FAILURE);
			}
		}

		printf(".");
		fflush(stdout);
	}

	printf(" ");
	fflush(stdout);

	for (int i = 0; KAT_XSK233_DECODE_BAD[i] != NULL; i ++) {
		uint8_t buf[30];
		xsk233_point p;

		HEXTOBIN(buf, KAT_XSK233_DECODE_BAD[i]);
		if (xsk233_decode(&p, buf) != 0x00000000) {
			printf("\ndecode should have failed (2)\n");
			exit(EXIT_FAILURE);
		}
		if (xsk233_is_neutral(&p) != 0xFFFFFFFF) {
			printf("\ndid not get neutral (2)\n");
			exit(EXIT_FAILURE);
		}

		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

static void
test_xsk233_add(void)
{
	printf("Test xsk233 add: ");
	fflush(stdout);

	for (int i = 0; KAT_XSK233_ADD[i] != NULL; i += 6) {
		uint8_t buf1[30], buf2[30], buf3[30];
		uint8_t buf4[30], buf5[30], buf6[30];
		xsk233_point P1, P2, P3, P4, P5, P6;
		xsk233_point Q, R;

		HEXTOBIN(buf1, KAT_XSK233_ADD[i + 0]);
		HEXTOBIN(buf2, KAT_XSK233_ADD[i + 1]);
		HEXTOBIN(buf3, KAT_XSK233_ADD[i + 2]);
		HEXTOBIN(buf4, KAT_XSK233_ADD[i + 3]);
		HEXTOBIN(buf5, KAT_XSK233_ADD[i + 4]);
		HEXTOBIN(buf6, KAT_XSK233_ADD[i + 5]);
		if (!(xsk233_decode(&P1, buf1)
			&& xsk233_decode(&P2, buf2)
			&& xsk233_decode(&P3, buf3)
			&& xsk233_decode(&P4, buf4)
			&& xsk233_decode(&P5, buf5)
			&& xsk233_decode(&P6, buf6)))
		{
			printf("\nKAT decode failed\n");
			exit(EXIT_FAILURE);
		}
		Q = P1;

		if (xsk233_equals(&P1, &P2) != 0) {
			printf("\npoints should not be equal\n");
			exit(EXIT_FAILURE);
		}

		xsk233_add(&Q, &P1, &P2);
		check_k233_eq("add1", &Q, &P3);
		xsk233_neg(&R, &P1);
		xsk233_add(&Q, &Q, &R);
		check_k233_eq("add2", &Q, &P2);
		xsk233_sub(&Q, &P3, &P1);
		check_k233_eq("add3", &Q, &P2);
		xsk233_double(&Q, &P1);
		check_k233_eq("add4", &Q, &P4);
		xsk233_add(&Q, &P3, &P1);
		xsk233_add(&R, &P4, &P2);
		check_k233_eq("add5", &Q, &P5);
		check_k233_eq("add6", &R, &P5);
		xsk233_double(&Q, &P3);
		xsk233_add(&R, &P5, &P2);
		check_k233_eq("add7", &Q, &P6);
		check_k233_eq("add8", &R, &P6);
		check_k233_eq("add9", &Q, &R);

		for (int j = 0; j <= 10; j ++) {
			Q = P1;
			for (int k = 0; k < j; k ++) {
				xsk233_double(&Q, &Q);
			}
			xsk233_xdouble(&R, &P1, j);
			check_k233_eq("add10", &Q, &R);
		}

		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

static const struct {
	const char *scalar;
	int8_t sd[48];
} KAT_XSK233_TAU[] = {
	{ "c88ec1a2d664c2cf428b2e31d916d0ce4c5c93c84e144ed0848a9141e8b1", {
	   10, -20,  20,  20,  18,  -2,  -6,  -7,  21,   8,   7,  21,
	    0,  -5,  -2,  20, -20,  10,  -3, -14, -14,  -6,  -1,  -4,
	   -6,  -8, -17,  10,  10,  -2,  -6,   6, -13,  10,  12, -18,
	  -20,   0,  14, -16,   4,  -4, -10,   7,  -5,   6,  -2,   0 } },
	{ "e4f84ae53f07e94be2b54a636426b585169fd7494c0ac6c4ed881481b6c9", {
	   21,  -6, -21,  18,   8,  -7,   1,  -2,  20, -16,  -4,  -9,
	    8, -10,   0,  -3,  -9,   0,   9,  -2,   4,  -1,   6,  -3,
	    6,  19,  20,  10,   1,  18,  -2, -16,   8,  18, -18, -14,
	  -10,   4,  11,  -4,  -9,  -3,   6,  -2,   4,  -7,   0,   0 } },
	{ "8fe9c7b11a82304b31566c81bc1e99a47e85806d1614fa70a42c7ccae108", {
	   17, -18, -18,   2,  20,  -6,  -8,  17,  20,  -4,   7,   8,
	    3, -14, -10,  18, -10,  -1,  20,  -8,   6,  19, -14,   0,
	    9,   0,   4,   0,   0, -11,  -8,   5,  20,   2,  -2, -19,
	   10,   8,   4,   3,  15, -12,   4,   1,  18, -20,   0,   0 } },
	{ "dc5bbf4fa58b9d1b18a98d216c7a2442007bb670b4051b92f3216a9ef662", {
	    1, -16,  14, -12,   4,   9,   0, -11, -12,   4,  -6,   3,
	   19, -20,  -2,   7,   6,  11,  -8,  16,   6, -12, -12,  12,
	  -18,  10,   7, -14,  14,  -8,   8,   3, -15,  18,   6,  -2,
	    4,   5,   8, -19,   6, -10, -17, -16, -16,  -8,   0,   0 } },
	{ "8c4fe8d6411731b357d833d0306d4bde77b972ed5ca5d68a39614972ce06", {
	   -3,  -8,   3,  21,   8,  -4, -16,  10,  -5,   6, -10, -18,
	  -12,   6,  10,   9,   8,  -8, -11,   8,  15,   2,  16, -10,
	    1, -17,  -4,   0,  -7,  -5,  -9,  13,  20,  16,  -6,  -7,
	   -9, -14, -10, -11, -18,   4,   7,  -8,   3,   9,   0,   0 } },
	{ "000000000000000000000000000000000000000000000000000000008000", {
	    1,  21,   6,  -1,  -9,  13,   8,  13,   4, -15, -18,  18,
	   -2,   9, -18, -10,   7,   5,   0,  11,  12,   6, -17,   0,
	  -17, -16,  10,   9,   7,  11,  -8,   2,  -7,   3,  17,  -6,
	   -4,  -4,  -9,   0,  14,  12, -10, -13, -20,   6,  -2,   0 } },

	{ NULL, { 0 } }
};

static void
test_xsk233_tau_recode(void)
{
	extern void xsk233_tau_recode(int8_t *sd, const void *n);

	printf("Test xsk233 tau_recode: ");

	for (int i = 0; KAT_XSK233_TAU[i].scalar != NULL; i ++) {
		uint8_t n[30];
		int8_t sd[48];
		const int8_t *sdref;

		HEXTOBIN(n, KAT_XSK233_TAU[i].scalar);
		sdref = KAT_XSK233_TAU[i].sd;
		xsk233_tau_recode(sd, n);
		if (memcmp(sd, sdref, sizeof sd) != 0) {
			printf("FAIL\n");
			printf("got: ");
			for (int j = 0; j < 48; j ++) {
				if (j % 12 == 0) {
					if (j != 0) {
						printf(",\n     ");
					}
				} else {
					printf(", ");
				}
				printf("%3d", sd[j]);
			}
			printf("\n");
			printf("exp: ");
			for (int j = 0; j < 48; j ++) {
				if (j % 12 == 0) {
					if (j != 0) {
						printf(",\n     ");
					}
				} else {
					printf(", ");
				}
				printf("%3d", sdref[j]);
			}
			printf("\n");
			exit(EXIT_FAILURE);
		}

		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

static void
test_xsk233_mul(void)
{
	printf("Test xsk233 mul: ");
	fflush(stdout);

	for (int i = 0; KAT_XSK233_MUL[i] != NULL; i += 3) {
		uint8_t buf1[30], buf2[30], buf3[30];
		xsk233_point P, Q, R;

		HEXTOBIN(buf1, KAT_XSK233_MUL[i + 0]);
		HEXTOBIN(buf2, KAT_XSK233_MUL[i + 1]);
		HEXTOBIN(buf3, KAT_XSK233_MUL[i + 2]);
		if (!(xsk233_decode(&P, buf1) && xsk233_decode(&R, buf3))) {
			printf("\nKAT decode failed\n");
			exit(EXIT_FAILURE);
		}
		xsk233_mul(&Q, &P, buf2, sizeof buf2);
		check_k233_eq("mul 1", &Q, &R);
		xsk233_mul_ladder(&Q, &P, buf2, sizeof buf2);
		check_k233_eq("mul_ladder 1", &Q, &R);
		xsk233_mul_frob(&Q, &P, buf2, sizeof buf2);
		check_k233_eq("mul_frob 1", &Q, &R);

		for (int j = 29; j >= 0; j --) {
			buf2[j] = 0;
			xsk233_mul(&Q, &P, buf2, sizeof buf2);
			xsk233_mul(&R, &P, buf2, j);
			check_k233_eq("mul 2", &Q, &R);
			xsk233_mul_ladder(&R, &P, buf2, j);
			check_k233_eq("mul_ladder 2", &Q, &R);
			xsk233_mul_frob(&R, &P, buf2, j);
			check_k233_eq("mul_frob 2", &Q, &R);
		}

		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

static void
test_xsk233_mulgen(void)
{
	printf("Test xsk233_mulgen: ");
	fflush(stdout);

	for (int i = 0; i < 20; i ++) {
		uint8_t buf[32];
		blake2s_context bc;
		xsk233_point Q, R;

		blake2s_init(&bc, 32);
		blake2s_update(&bc, "test_mulgen", 11);
		buf[0] = (uint8_t)i;
		buf[1] = (uint8_t)(i >> 8);
		blake2s_update(&bc, buf, 2);
		blake2s_final(&bc, buf);
		for (int j = 30; j >= 0; j --) {
			xsk233_mulgen(&Q, buf, j);
			xsk233_mul(&R, &xsk233_generator, buf, j);
			check_k233_eq("mulgen", &Q, &R);
			xsk233_mulgen_frob(&Q, buf, j);
			check_k233_eq("mulgen_frob", &Q, &R);
		}

		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

int
main(void)
{
	test_xsb233_encode();
	test_xsb233_add();
	test_xsb233_mul();
	test_xsb233_mulgen();
	test_xsk233_encode();
	test_xsk233_add();
	test_xsk233_tau_recode();
	test_xsk233_mul();
	test_xsk233_mulgen();
	return 0;
}
