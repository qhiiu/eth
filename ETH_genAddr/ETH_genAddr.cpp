#include <stdint.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <iostream>
#include <vector>
#include <algorithm> 
#include <cstdint>
#include <inttypes.h>
#include <math.h> 
#include <emmintrin.h>

using namespace std;

//============== INT.H - start ================================================================================
//---------- INT.H
#define NB64BLOCK 5
#define NB32BLOCK 10

class Int {
	public:
		Int();
		Int(int64_t i64);
		Int(uint64_t u64);
		Int(Int* a);

		void Add(uint64_t a);
		void Add(Int* a);
		void Add(Int* a, Int* b);
		void AddOne();
		void Sub(uint64_t a);
		void Sub(Int* a);
		void Sub(Int* a, Int* b);
		void Mult(Int* a);
		uint64_t Mult(uint64_t a);
		uint64_t Mult(Int* a, uint64_t b);
		uint64_t IMult(Int* a, int64_t b);
		void Mult(Int* a, Int* b);
		void Neg();

		// Comp 
		bool IsGreaterOrEqual(Int* a);
		bool IsEqual(Int* a);
		bool IsZero();
		bool IsOne();
		bool IsPositive();
		bool IsNegative();
		bool IsEven();

		// Setup field
		static void SetupField(Int* n, Int* R = NULL, Int* R2 = NULL, Int* R3 = NULL, Int* R4 = NULL);

		void ModInv();                             // this <- this^-1 (mod n)
		void MontgomeryMult(Int* a, Int* b);        // this <- a*b*R^-1 (mod n)
		void ModAdd(Int* a);                       // this <- this+a (mod n) [0<a<P]
		void ModAdd(Int* a, Int* b);                // this <- a+b (mod n) [0<a,b<P]
		void ModSub(Int* a);                       // this <- this-a (mod n) [0<a<P]
		void ModSub(Int* a, Int* b);               // this <- a-b (mod n) [0<a,b<P]
		void ModMul(Int* a, Int* b);                // this <- a*b (mod n) 
		void ModNeg();                             // this <- -this (mod n)

		// Specific SecpK1
		void ModMulK1(Int* a, Int* b);
		void ModMulK1(Int* a);
		void ModSquareK1(Int* a);

		// Size
		int GetSize();       // Number of significant 32bit limbs

		// Setter
		void SetInt32(uint32_t value);
		void Set(Int* a);
		void SetBase10(const char* value);
		void SetBase16(const char* value);
		void SetBaseN(int n, const char* charset, const char* value);

		// Getter
		unsigned char GetByte(int n);
		void Get32Bytes(unsigned char* buff);

		// To String
		std::string GetBase10();
		std::string GetBase16();
		std::string GetBaseN(int n, char* charset);

		union {
			uint32_t bits[NB32BLOCK];
			uint64_t bits64[NB64BLOCK];
		};

	private:
		void MatrixVecMul(Int* u, Int* v, int64_t _11, int64_t _12, int64_t _21, int64_t _22, uint64_t* cu, uint64_t* cv);
		void MatrixVecMul(Int* u, Int* v, int64_t _11, int64_t _12, int64_t _21, int64_t _22);
		uint64_t AddCh(Int* a, uint64_t ca, Int* b, uint64_t cb);
		uint64_t AddCh(Int* a, uint64_t ca);
		uint64_t AddC(Int* a);
		void AddAndShift(Int* a, Int* b, uint64_t cH);
		void CLEAR();
		void DivStep62(Int* u, Int* v, int64_t* eta, int* pos, int64_t* uu, int64_t* uv, int64_t* vu, int64_t* vv);
};

// Missing intrinsics
static uint64_t inline _umul128(uint64_t a, uint64_t b, uint64_t* h) {
	uint64_t rhi;
	uint64_t rlo;
	__asm__("mulq  %[b];" :"=d"(rhi), "=a"(rlo) : "1"(a), [b]"rm"(b));
	*h = rhi;
	return rlo;
}

static int64_t inline _mul128(int64_t a, int64_t b, int64_t* h) {
	uint64_t rhi;
	uint64_t rlo;
	__asm__("imulq  %[b];" :"=d"(rhi), "=a"(rlo) : "1"(a), [b]"rm"(b));
	*h = rhi;
	return rlo;
}

static uint64_t inline _udiv128(uint64_t hi, uint64_t lo, uint64_t d, uint64_t* r) {
	uint64_t q;
	uint64_t _r;
	__asm__("divq  %[d];" :"=d"(_r), "=a"(q) : "d"(hi), "a"(lo), [d]"rm"(d));
	*r = _r;
	return q;
}

static uint64_t inline __rdtsc() {
	uint32_t h;
	uint32_t l;
	__asm__("rdtsc;" :"=d"(h), "=a"(l));
	return (uint64_t)h << 32 | (uint64_t)l;
}

#define __shiftright128(a,b,n) ((a)>>(n))|((b)<<(64-(n)))
#define __shiftleft128(a,b,n) ((b)<<(n))|((a)>>(64-(n)))


#define _subborrow_u64(a,b,c,d) __builtin_ia32_sbb_u64(a,b,c,(long long unsigned int*)d);
#define _addcarry_u64(a,b,c,d) __builtin_ia32_addcarryx_u64(a,b,c,(long long unsigned int*)d);
#define _byteswap_uint64 __builtin_bswap64
#define LZC(x) __builtin_clzll(x)
#define TZC(x) __builtin_ctzll(x)

static void inline imm_mul(uint64_t* x, uint64_t y, uint64_t* dst, uint64_t* carryH) {
	unsigned char c = 0;
	uint64_t h, carry;
	dst[0] = _umul128(x[0], y, &h); carry = h;
	c = _addcarry_u64(c, _umul128(x[1], y, &h), carry, dst + 1); carry = h;
	c = _addcarry_u64(c, _umul128(x[2], y, &h), carry, dst + 2); carry = h;
	c = _addcarry_u64(c, _umul128(x[3], y, &h), carry, dst + 3); carry = h;
	c = _addcarry_u64(c, _umul128(x[4], y, &h), carry, dst + 4); carry = h;

	* carryH = carry;
}

static void inline imm_imul(uint64_t* x, uint64_t y, uint64_t* dst, uint64_t* carryH) {
	unsigned char c = 0;
	uint64_t h, carry;
	dst[0] = _umul128(x[0], y, &h); carry = h;
	c = _addcarry_u64(c, _umul128(x[1], y, &h), carry, dst + 1); carry = h;
	c = _addcarry_u64(c, _umul128(x[2], y, &h), carry, dst + 2); carry = h;
	c = _addcarry_u64(c, _umul128(x[3], y, &h), carry, dst + 3); carry = h;

	c = _addcarry_u64(c, _mul128(x[NB64BLOCK - 1], y, (int64_t*)&h), carry, dst + NB64BLOCK - 1); carry = h;
	*carryH = carry;
}

static void inline imm_umul(uint64_t* x, uint64_t y, uint64_t* dst) {
	// Assume that x[NB64BLOCK-1] is 0
	unsigned char c = 0;
	uint64_t h, carry;
	dst[0] = _umul128(x[0], y, &h); carry = h;
	c = _addcarry_u64(c, _umul128(x[1], y, &h), carry, dst + 1); carry = h;
	c = _addcarry_u64(c, _umul128(x[2], y, &h), carry, dst + 2); carry = h;
	c = _addcarry_u64(c, _umul128(x[3], y, &h), carry, dst + 3); carry = h;

	_addcarry_u64(c, 0ULL, carry, dst + (NB64BLOCK - 1));
}

static void inline shiftR(unsigned char n, uint64_t* d) {
	d[0] = __shiftright128(d[0], d[1], n);
	d[1] = __shiftright128(d[1], d[2], n);
	d[2] = __shiftright128(d[2], d[3], n);
	d[3] = __shiftright128(d[3], d[4], n);

	d[NB64BLOCK - 1] = ((int64_t)d[NB64BLOCK - 1]) >> n;
}

static void inline shiftR(unsigned char n, uint64_t* d, uint64_t h) {
	d[0] = __shiftright128(d[0], d[1], n);
	d[1] = __shiftright128(d[1], d[2], n);
	d[2] = __shiftright128(d[2], d[3], n);
	d[3] = __shiftright128(d[3], d[4], n);

	d[NB64BLOCK - 1] = __shiftright128(d[NB64BLOCK - 1], h, n);
}

//----------------- INT.CPP
Int _ONE((uint64_t)1);

Int::Int() {}

Int::Int(Int* a) {	if (a) Set(a);	else CLEAR();  }

Int::Int(int64_t i64) {
	if (i64 < 0) {	memset(bits64, 0xFF, NB64BLOCK * 8); }	else {	CLEAR(); }
	bits64[0] = i64;
}

Int::Int(uint64_t u64) {	CLEAR();	bits64[0] = u64; }
// ------------------------------------------------
void Int::CLEAR() { memset(bits64, 0, NB64BLOCK * 8); }
// ------------------------------------------------
void Int::Set(Int* a) {
	for (int i = 0; i < NB64BLOCK; i++){ 
		bits64[i] = a->bits64[i]; }	
}
// ------------------------------------------------ 
void Int::Add(Int* a) {
	unsigned char c = 0;
	c = _addcarry_u64(c, bits64[0], a->bits64[0], bits64 + 0);
	c = _addcarry_u64(c, bits64[1], a->bits64[1], bits64 + 1);
	c = _addcarry_u64(c, bits64[2], a->bits64[2], bits64 + 2);
	c = _addcarry_u64(c, bits64[3], a->bits64[3], bits64 + 3);
	c = _addcarry_u64(c, bits64[4], a->bits64[4], bits64 + 4);
}
// ------------------------------------------------
void Int::Add(uint64_t a) {
	unsigned char c = 0;
	c = _addcarry_u64(c, bits64[0], a, bits64 + 0);
	c = _addcarry_u64(c, bits64[1], 0, bits64 + 1);
	c = _addcarry_u64(c, bits64[2], 0, bits64 + 2);
	c = _addcarry_u64(c, bits64[3], 0, bits64 + 3);
	c = _addcarry_u64(c, bits64[4], 0, bits64 + 4);
}
// ------------------------------------------------
void Int::AddOne() {
	unsigned char c = 0;
	c = _addcarry_u64(c, bits64[0], 1, bits64 + 0);
	c = _addcarry_u64(c, bits64[1], 0, bits64 + 1);
	c = _addcarry_u64(c, bits64[2], 0, bits64 + 2);
	c = _addcarry_u64(c, bits64[3], 0, bits64 + 3);
	c = _addcarry_u64(c, bits64[4], 0, bits64 + 4);
}
// ------------------------------------------------
void Int::Add(Int* a, Int* b) {
	unsigned char c = 0;
	c = _addcarry_u64(c, b->bits64[0], a->bits64[0], bits64 + 0);
	c = _addcarry_u64(c, b->bits64[1], a->bits64[1], bits64 + 1);
	c = _addcarry_u64(c, b->bits64[2], a->bits64[2], bits64 + 2);
	c = _addcarry_u64(c, b->bits64[3], a->bits64[3], bits64 + 3);
	c = _addcarry_u64(c, b->bits64[4], a->bits64[4], bits64 + 4);
}
// ------------------------------------------------
uint64_t Int::AddCh(Int* a, uint64_t ca, Int* b, uint64_t cb) {
	uint64_t carry;
	unsigned char c = 0;
	c = _addcarry_u64(c, a->bits64[0], b->bits64[0], bits64 + 0);
	c = _addcarry_u64(c, a->bits64[1], b->bits64[1], bits64 + 1);
	c = _addcarry_u64(c, a->bits64[2], b->bits64[2], bits64 + 2);
	c = _addcarry_u64(c, a->bits64[3], b->bits64[3], bits64 + 3);
	c = _addcarry_u64(c, a->bits64[4], b->bits64[4], bits64 + 4);

	_addcarry_u64(c, ca, cb, &carry);
	return carry;
}

uint64_t Int::AddCh(Int* a, uint64_t ca) {
	uint64_t carry;
	unsigned char c = 0;
	c = _addcarry_u64(c, bits64[0], a->bits64[0], bits64 + 0);
	c = _addcarry_u64(c, bits64[1], a->bits64[1], bits64 + 1);
	c = _addcarry_u64(c, bits64[2], a->bits64[2], bits64 + 2);
	c = _addcarry_u64(c, bits64[3], a->bits64[3], bits64 + 3);
	c = _addcarry_u64(c, bits64[4], a->bits64[4], bits64 + 4);

	_addcarry_u64(c, ca, 0, &carry);
	return carry;
}
// ------------------------------------------------
uint64_t Int::AddC(Int* a) {
	unsigned char c = 0;
	c = _addcarry_u64(c, bits64[0], a->bits64[0], bits64 + 0);
	c = _addcarry_u64(c, bits64[1], a->bits64[1], bits64 + 1);
	c = _addcarry_u64(c, bits64[2], a->bits64[2], bits64 + 2);
	c = _addcarry_u64(c, bits64[3], a->bits64[3], bits64 + 3);
	c = _addcarry_u64(c, bits64[4], a->bits64[4], bits64 + 4);

	return c;
}
// ------------------------------------------------
void Int::AddAndShift(Int* a, Int* b, uint64_t cH) {
	unsigned char c = 0;
	c = _addcarry_u64(c, b->bits64[0], a->bits64[0], bits64 + 0);
	c = _addcarry_u64(c, b->bits64[1], a->bits64[1], bits64 + 0);
	c = _addcarry_u64(c, b->bits64[2], a->bits64[2], bits64 + 1);
	c = _addcarry_u64(c, b->bits64[3], a->bits64[3], bits64 + 2);
	c = _addcarry_u64(c, b->bits64[4], a->bits64[4], bits64 + 3);

	bits64[NB64BLOCK - 1] = c + cH;
}
// ------------------------------------------------
void Int::MatrixVecMul(Int* u, Int* v, int64_t _11, int64_t _12, int64_t _21, int64_t _22, uint64_t* cu, uint64_t* cv) {
	Int t1, t2, t3, t4;
	uint64_t c1, c2, c3, c4;
	c1 = t1.IMult(u, _11);
	c2 = t2.IMult(v, _12);
	c3 = t3.IMult(u, _21);
	c4 = t4.IMult(v, _22);
	*cu = u->AddCh(&t1, c1, &t2, c2);
	*cv = v->AddCh(&t3, c3, &t4, c4);
}

void Int::MatrixVecMul(Int* u, Int* v, int64_t _11, int64_t _12, int64_t _21, int64_t _22) {
	Int t1, t2, t3, t4;
	t1.IMult(u, _11);
	t2.IMult(v, _12);
	t3.IMult(u, _21);
	t4.IMult(v, _22);
	u->Add(&t1, &t2);
	v->Add(&t3, &t4);
}

// ------------------------------------------------
bool Int::IsGreaterOrEqual(Int* a) {
	Int p;
	p.Sub(this, a);
	return p.IsPositive();
}
// ------------------------------------------------
bool Int::IsEqual(Int* a) {
	return
		(bits64[4] == a->bits64[4]) &&
		(bits64[3] == a->bits64[3]) &&
		(bits64[2] == a->bits64[2]) &&
		(bits64[1] == a->bits64[1]) &&
		(bits64[0] == a->bits64[0]);
}
// ------------------------------------------------
bool Int::IsOne() {	return IsEqual(&_ONE); }
// ------------------------------------------------
bool Int::IsZero() { return (bits64[4] | bits64[3] | bits64[2] | bits64[1] | bits64[0]) == 0; }
// ------------------------------------------------
void Int::SetInt32(uint32_t value) {	CLEAR();	bits[0] = value;	}
// ------------------------------------------------
unsigned char Int::GetByte(int n) {
	unsigned char* bbPtr = (unsigned char*)bits;
	return bbPtr[n];
}
// ------------------------------------------------
void Int::Get32Bytes(unsigned char* buff) {
	uint64_t* ptr = (uint64_t*)buff;
	ptr[3] = _byteswap_uint64(bits64[0]);
	ptr[2] = _byteswap_uint64(bits64[1]);
	ptr[1] = _byteswap_uint64(bits64[2]);
	ptr[0] = _byteswap_uint64(bits64[3]);
}
// ------------------------------------------------
void Int::Sub(Int* a) {
	unsigned char c = 0;
	c = _subborrow_u64(c, bits64[0], a->bits64[0], bits64 + 0);
	c = _subborrow_u64(c, bits64[1], a->bits64[1], bits64 + 1);
	c = _subborrow_u64(c, bits64[2], a->bits64[2], bits64 + 2);
	c = _subborrow_u64(c, bits64[3], a->bits64[3], bits64 + 3);
	c = _subborrow_u64(c, bits64[4], a->bits64[4], bits64 + 4);
}
// ------------------------------------------------
void Int::Sub(Int* a, Int* b) {
	unsigned char c = 0;
	c = _subborrow_u64(c, a->bits64[0], b->bits64[0], bits64 + 0);
	c = _subborrow_u64(c, a->bits64[1], b->bits64[1], bits64 + 1);
	c = _subborrow_u64(c, a->bits64[2], b->bits64[2], bits64 + 2);
	c = _subborrow_u64(c, a->bits64[3], b->bits64[3], bits64 + 3);
	c = _subborrow_u64(c, a->bits64[4], b->bits64[4], bits64 + 4);
}
// ------------------------------------------------
void Int::Sub(uint64_t a) {
	unsigned char c = 0;
	c = _subborrow_u64(c, bits64[0], a, bits64 + 0);
	c = _subborrow_u64(c, bits64[1], 0, bits64 + 1);
	c = _subborrow_u64(c, bits64[2], 0, bits64 + 2);
	c = _subborrow_u64(c, bits64[3], 0, bits64 + 3);
	c = _subborrow_u64(c, bits64[4], 0, bits64 + 4);
}
// ------------------------------------------------
bool Int::IsPositive() {	return (int64_t)(bits64[NB64BLOCK - 1]) >= 0;  }
// ------------------------------------------------
bool Int::IsNegative() { return (int64_t)(bits64[NB64BLOCK - 1]) < 0;   }
// ------------------------------------------------
bool Int::IsEven() {	return (bits[0] & 0x1) == 0;  }
// ------------------------------------------------
void Int::Neg() {
	unsigned char c = 0;
	c = _subborrow_u64(c, 0, bits64[0], bits64 + 0);
	c = _subborrow_u64(c, 0, bits64[1], bits64 + 1);
	c = _subborrow_u64(c, 0, bits64[2], bits64 + 2);
	c = _subborrow_u64(c, 0, bits64[3], bits64 + 3);
	c = _subborrow_u64(c, 0, bits64[4], bits64 + 4);
}
// ------------------------------------------------
void Int::Mult(Int* a) {	Int b(this);	Mult(a, &b);	}
// ------------------------------------------------
uint64_t Int::Mult(uint64_t a) {
	uint64_t carry;
	imm_mul(bits64, a, bits64, &carry);
	return carry;
}
// ------------------------------------------------
uint64_t Int::IMult(Int* a, int64_t b) {
	uint64_t carry;
	// Make b positive
	if (b < 0LL) {
		unsigned char c = 0;
		c = _subborrow_u64(c, 0, a->bits64[0], bits64 + 0);
		c = _subborrow_u64(c, 0, a->bits64[1], bits64 + 1);
		c = _subborrow_u64(c, 0, a->bits64[2], bits64 + 2);
		c = _subborrow_u64(c, 0, a->bits64[3], bits64 + 3);
		c = _subborrow_u64(c, 0, a->bits64[4], bits64 + 4);

		b = -b;
	} else {	Set(a);	}

	imm_imul(bits64, b, bits64, &carry);
	return carry;
}
// ------------------------------------------------
uint64_t Int::Mult(Int* a, uint64_t b) {
	uint64_t carry;
	imm_mul(a->bits64, b, bits64, &carry);
	return carry;
}
// ------------------------------------------------
void Int::Mult(Int* a, Int* b) {

	unsigned char c = 0;
	uint64_t h;
	uint64_t pr = 0;
	uint64_t carryh = 0;
	uint64_t carryl = 0;

	bits64[0] = _umul128(a->bits64[0], b->bits64[0], &pr);

	for (int i = 1; i < NB64BLOCK; i++) {
		for (int j = 0; j <= i; j++) {
			c = _addcarry_u64(c, _umul128(a->bits64[j], b->bits64[i - j], &h), pr, &pr);
			c = _addcarry_u64(c, carryl, h, &carryl);
			c = _addcarry_u64(c, carryh, 0, &carryh);
		}
		bits64[i] = pr;
		pr = carryl;
		carryl = carryh;
		carryh = 0;
	}
}
// ------------------------------------------------
int Int::GetSize() {
	int i = NB32BLOCK - 1;
	while (i > 0 && bits[i] == 0) i--;
	return i + 1;
}
// ------------------------------------------------
void Int::SetBase10(const char* value) {

	CLEAR();
	Int pw((uint64_t)1);
	Int c;
	int lgth = (int)strlen(value);
	for (int i = lgth - 1; i >= 0; i--) {
		uint32_t id = (uint32_t)(value[i] - '0');
		c.Set(&pw);
		c.Mult(id);
		Add(&c);
		pw.Mult(10);
	}
}

// ------------------------------------------------
void  Int::SetBase16(const char* value) {	SetBaseN(16, "0123456789ABCDEF", value); }
// ------------------------------------------------
std::string Int::GetBase10() {	return GetBaseN(10, "0123456789"); }
// ------------------------------------------------
std::string Int::GetBase16() {	return GetBaseN(16, "0123456789ABCDEF"); }
// ------------------------------------------------
void  Int::SetBaseN(int n, const char* charset, const char* value) {

	CLEAR();

	Int pw((uint64_t)1);
	Int nb((uint64_t)n);
	Int c;

	int lgth = (int)strlen(value);
	for (int i = lgth - 1; i >= 0; i--) {
		const char* p = strchr(charset, toupper(value[i]));
		if (!p) { printf("Invalid charset !!\n");	return;	}
		int id = (int)(p - charset);
		c.SetInt32(id);
		c.Mult(&pw);
		Add(&c);
		pw.Mult(&nb);

	}
}
// ------------------------------------------------
std::string Int::GetBaseN(int n, char* charset) {

	std::string ret;

	Int N(this);
	int isNegative = N.IsNegative();
	if (isNegative) N.Neg();

	// TODO: compute max digit
	unsigned char digits[1024];
	memset(digits, 0, sizeof(digits));

	int digitslen = 1;
	for (int i = 0; i < NB64BLOCK * 8; i++) {
		unsigned int carry = N.GetByte(NB64BLOCK * 8 - i - 1);
		for (int j = 0; j < digitslen; j++) {
			carry += (unsigned int)(digits[j]) << 8;
			digits[j] = (unsigned char)(carry % n);
			carry /= n;
		}
		while (carry > 0) {
			digits[digitslen++] = (unsigned char)(carry % n);
			carry /= n;
		}
	}

	// reverse
	if (isNegative)
		ret.push_back('-');

	for (int i = 0; i < digitslen; i++)
		ret.push_back(charset[digits[digitslen - 1 - i]]);

	if (ret.length() == 0)
		ret.push_back('0');

	return ret;
}
//============== INT.CPP - end ================================================================================

//============== INT->INTMOD.CPP - start ================================================================================
#include <emmintrin.h>

static Int     _P;       // Field characteristic
static Int     _R;       // Montgomery multiplication R
static Int     _R2;      // Montgomery multiplication R2
static Int     _R3;      // Montgomery multiplication R3
static Int     _R4;      // Montgomery multiplication R4
static int32_t  Msize;    // Montgomery mult size
static uint64_t MM64;     // 64bits lsb negative inverse of P
#define MSK62  0x3FFFFFFFFFFFFFFF

extern Int _ONE;

// ------------------------------------------------
void Int::ModAdd(Int* a) {
	Int p;

	Add(a);
	p.Sub(this, &_P);
	if (p.IsPositive())
		Set(&p);
}
// ------------------------------------------------
void Int::ModAdd(Int* a, Int* b) {
	Int p;
	Add(a, b);
	p.Sub(this, &_P);
	if (p.IsPositive()){Set(&p);}
		
}
// ------------------------------------------------
void Int::ModSub(Int* a) {	Sub(a); 	if (IsNegative()){ Add(&_P); }	}
// ------------------------------------------------
void Int::ModSub(Int* a, Int* b) {	Sub(a, b);	if (IsNegative()){ Add(&_P); }   }
// ------------------------------------------------
void Int::ModNeg() {	Neg();	Add(&_P);	}
// ------------------------------------------------
void Int::DivStep62(Int* u, Int* v, int64_t* eta, int* pos, int64_t* uu, int64_t* uv, int64_t* vu, int64_t* vv) {
	int bitCount;
	uint64_t u0 = u->bits64[0];
	uint64_t v0 = v->bits64[0];

#define SWAP(tmp,x,y) tmp = x; x = y; y = tmp;
	uint64_t uh, vh, w, x;
	unsigned char c = 0;

	// Extract 64 MSB of u and v
	// u and v must be positive

	while (*pos >= 1 && (u->bits64[*pos] | v->bits64[*pos]) == 0) (*pos)--;
	if (*pos == 0) {	uh = u->bits64[0];	vh = v->bits64[0];	}
	else {
		uint64_t s = LZC(u->bits64[*pos] | v->bits64[*pos]);
		if (s == 0) {	uh = u->bits64[*pos];	vh = v->bits64[*pos];	}
		else {
			uh = __shiftleft128(u->bits64[*pos - 1], u->bits64[*pos], (uint8_t)s);
			vh = __shiftleft128(v->bits64[*pos - 1], v->bits64[*pos], (uint8_t)s);
		}
	}

	bitCount = 62;
	__m128i _u, _v, _t;

	((int64_t*)&_u)[0] = 1;
	((int64_t*)&_u)[1] = 0;
	((int64_t*)&_v)[0] = 0;
	((int64_t*)&_v)[1] = 1;

	while (true) {
		// Use a sentinel bit to count zeros only up to bitCount
		uint64_t zeros = TZC(v0 | 1ULL << bitCount);
		vh >>= zeros;
		v0 >>= zeros;
		_u = _mm_slli_epi64(_u, (int)zeros);
		bitCount -= (int)zeros;

		if (bitCount <= 0) {	break;	}

		if (vh < uh) {
			SWAP(w, uh, vh);
			SWAP(x, u0, v0);
			SWAP(_t, _u, _v);
		}

		vh -= uh;
		v0 -= u0;
		_v = _mm_sub_epi64(_v, _u);
	}

	*uu = ((int64_t*)&_u)[0];
	*uv = ((int64_t*)&_u)[1];
	*vu = ((int64_t*)&_v)[0];
	*vv = ((int64_t*)&_v)[1];
}

// ------------------------------------------------
void Int::ModInv() {

	Int u(&_P);
	Int v(this);
	Int r((int64_t)0);
	Int s((int64_t)1);

	// Delayed right shift 62bits
	Int r0_P;
	Int s0_P;

	int64_t  eta = -1;
	int64_t uu, uv, vu, vv;
	uint64_t carryS, carryR;
	int pos = NB64BLOCK - 1;
	while (pos >= 1 && (u.bits64[pos] | v.bits64[pos]) == 0) pos--;

	while (!v.IsZero()) { 
		DivStep62(&u, &v, &eta, &pos, &uu, &uv, &vu, &vv);
		// Now update BigInt variables
		MatrixVecMul(&u, &v, uu, uv, vu, vv);

		// Make u,v positive
		// Required only for Pornin's method
		if (u.IsNegative()) {
			u.Neg();
			uu = -uu;
			uv = -uv;
		}
		if (v.IsNegative()) {
			v.Neg();
			vu = -vu;
			vv = -vv;
		}

		MatrixVecMul(&r, &s, uu, uv, vu, vv, &carryR, &carryS);

		// Compute multiple of P to add to s and r to make them multiple of 2^62
		uint64_t r0 = (r.bits64[0] * MM64) & MSK62;
		uint64_t s0 = (s.bits64[0] * MM64) & MSK62;
		r0_P.Mult(&_P, r0);
		s0_P.Mult(&_P, s0);
		carryR = r.AddCh(&r0_P, carryR);
		carryS = s.AddCh(&s0_P, carryS);

		// Right shift all variables by 62bits
		shiftR(62, u.bits64);
		shiftR(62, v.bits64);
		shiftR(62, r.bits64, carryR);
		shiftR(62, s.bits64, carryS);

	}

	if (!u.IsOne()) { CLEAR(); return;	}		// No inverse
	while (r.IsNegative()){ r.Add(&_P); }
	while (r.IsGreaterOrEqual(&_P)){ r.Sub(&_P); }
	Set(&r);
}
// ------------------------------------------------
void Int::ModMul(Int* a, Int* b) {
	Int p;
	p.MontgomeryMult(a, b);
	MontgomeryMult(&_R2, &p);
}
// ------------------------------------------------
void Int::SetupField(Int* n, Int* R, Int* R2, Int* R3, Int* R4) {

	// Size in number of 32bit word
	int nSize = n->GetSize();

	// Last digit inversions (Newton's iteration)
	{
		int64_t x, t;
		x = t = (int64_t)n->bits64[0];
		x = x * (2 - t * x);
		x = x * (2 - t * x);
		x = x * (2 - t * x);
		x = x * (2 - t * x);
		x = x * (2 - t * x);
		MM64 = (uint64_t)(-x);
	}
	_P.Set(n);

	// Size of Montgomery mult (64bits digit)
	Msize = nSize / 2;

	// Compute few power of R
	// R = 2^(64*Msize) mod n
	Int Ri;
	Ri.MontgomeryMult(&_ONE, &_ONE); // Ri = R^-1
	_R.Set(&Ri);                     // R  = R^-1
	_R2.MontgomeryMult(&Ri, &_ONE);  // R2 = R^-2
	_R3.MontgomeryMult(&Ri, &Ri);    // R3 = R^-3
	_R4.MontgomeryMult(&_R3, &_ONE); // R4 = R^-4

	_R.ModInv();                     // R  = R
	_R2.ModInv();                    // R2 = R^2
	_R3.ModInv();                    // R3 = R^3
	_R4.ModInv();                    // R4 = R^4

	if (R){ R->Set(&_R); }
	if (R2){R2->Set(&_R2);}
	if (R3){R3->Set(&_R3);}
	if (R4){R4->Set(&_R4);}
}
// ------------------------------------------------
void Int::MontgomeryMult(Int* a, Int* b) {
	Int pr;
	Int p;
	uint64_t ML;
	uint64_t c;

	// i = 0
	imm_umul(a->bits64, b->bits64[0], pr.bits64);
	ML = pr.bits64[0] * MM64;
	imm_umul(_P.bits64, ML, p.bits64);
	c = pr.AddC(&p);
	memcpy(bits64, pr.bits64 + 1, 8 * (NB64BLOCK - 1));
	bits64[NB64BLOCK - 1] = c;

	for (int i = 1; i < Msize; i++) {
		imm_umul(a->bits64, b->bits64[i], pr.bits64);
		ML = (pr.bits64[0] + bits64[0]) * MM64;
		imm_umul(_P.bits64, ML, p.bits64);
		c = pr.AddC(&p);
		AddAndShift(this, &pr, c);
	}

	p.Sub(this, &_P);
	if (p.IsPositive()){ Set(&p); }
}

// SecpK1 specific section -----------------------------------------------------------------------------
void Int::ModMulK1(Int* a, Int* b) 
{
	unsigned char c;
	uint64_t ah, al;
	uint64_t t[NB64BLOCK];
	uint64_t r512[8];
	r512[5] = 0;
	r512[6] = 0;
	r512[7] = 0;

	// 256*256 multiplier
	imm_umul(a->bits64, b->bits64[0], r512);
	imm_umul(a->bits64, b->bits64[1], t);
	c = _addcarry_u64(0, r512[1], t[0], r512 + 1);
	c = _addcarry_u64(c, r512[2], t[1], r512 + 2);
	c = _addcarry_u64(c, r512[3], t[2], r512 + 3);
	c = _addcarry_u64(c, r512[4], t[3], r512 + 4);
	c = _addcarry_u64(c, r512[5], t[4], r512 + 5);
	imm_umul(a->bits64, b->bits64[2], t);
	c = _addcarry_u64(0, r512[2], t[0], r512 + 2);
	c = _addcarry_u64(c, r512[3], t[1], r512 + 3);
	c = _addcarry_u64(c, r512[4], t[2], r512 + 4);
	c = _addcarry_u64(c, r512[5], t[3], r512 + 5);
	c = _addcarry_u64(c, r512[6], t[4], r512 + 6);
	imm_umul(a->bits64, b->bits64[3], t);
	c = _addcarry_u64(0, r512[3], t[0], r512 + 3);
	c = _addcarry_u64(c, r512[4], t[1], r512 + 4);
	c = _addcarry_u64(c, r512[5], t[2], r512 + 5);
	c = _addcarry_u64(c, r512[6], t[3], r512 + 6);
	c = _addcarry_u64(c, r512[7], t[4], r512 + 7);

	// Reduce from 512 to 320 
	imm_umul(r512 + 4, 0x1000003D1ULL, t);
	c = _addcarry_u64(0, r512[0], t[0], r512 + 0);
	c = _addcarry_u64(c, r512[1], t[1], r512 + 1);
	c = _addcarry_u64(c, r512[2], t[2], r512 + 2);
	c = _addcarry_u64(c, r512[3], t[3], r512 + 3);

	// Reduce from 320 to 256 
	// No overflow possible here t[4]+c<=0x1000003D1ULL
	al = _umul128(t[4] + c, 0x1000003D1ULL, &ah);
	c = _addcarry_u64(0, r512[0], al, bits64 + 0);
	c = _addcarry_u64(c, r512[1], ah, bits64 + 1);
	c = _addcarry_u64(c, r512[2], 0ULL, bits64 + 2);
	c = _addcarry_u64(c, r512[3], 0ULL, bits64 + 3);

	// Probability of carry here or that this>P is very very unlikely
	bits64[4] = 0;
}
// ------------------------------------------------
void Int::ModMulK1(Int* a) 
{
	unsigned char c;
	uint64_t ah, al;
	uint64_t t[NB64BLOCK];
	uint64_t r512[8];
	r512[5] = 0;
	r512[6] = 0;
	r512[7] = 0;

	// 256*256 multiplier
	imm_umul(a->bits64, bits64[0], r512);
	imm_umul(a->bits64, bits64[1], t);
	c = _addcarry_u64(0, r512[1], t[0], r512 + 1);
	c = _addcarry_u64(c, r512[2], t[1], r512 + 2);
	c = _addcarry_u64(c, r512[3], t[2], r512 + 3);
	c = _addcarry_u64(c, r512[4], t[3], r512 + 4);
	c = _addcarry_u64(c, r512[5], t[4], r512 + 5);
	imm_umul(a->bits64, bits64[2], t);
	c = _addcarry_u64(0, r512[2], t[0], r512 + 2);
	c = _addcarry_u64(c, r512[3], t[1], r512 + 3);
	c = _addcarry_u64(c, r512[4], t[2], r512 + 4);
	c = _addcarry_u64(c, r512[5], t[3], r512 + 5);
	c = _addcarry_u64(c, r512[6], t[4], r512 + 6);
	imm_umul(a->bits64, bits64[3], t);
	c = _addcarry_u64(0, r512[3], t[0], r512 + 3);
	c = _addcarry_u64(c, r512[4], t[1], r512 + 4);
	c = _addcarry_u64(c, r512[5], t[2], r512 + 5);
	c = _addcarry_u64(c, r512[6], t[3], r512 + 6);
	c = _addcarry_u64(c, r512[7], t[4], r512 + 7);

	// Reduce from 512 to 320 
	imm_umul(r512 + 4, 0x1000003D1ULL, t);
	c = _addcarry_u64(0, r512[0], t[0], r512 + 0);
	c = _addcarry_u64(c, r512[1], t[1], r512 + 1);
	c = _addcarry_u64(c, r512[2], t[2], r512 + 2);
	c = _addcarry_u64(c, r512[3], t[3], r512 + 3);

	// Reduce from 320 to 256 
	// No overflow possible here t[4]+c<=0x1000003D1ULL
	al = _umul128(t[4] + c, 0x1000003D1ULL, &ah);
	c = _addcarry_u64(0, r512[0], al, bits64 + 0);
	c = _addcarry_u64(c, r512[1], ah, bits64 + 1);
	c = _addcarry_u64(c, r512[2], 0, bits64 + 2);
	c = _addcarry_u64(c, r512[3], 0, bits64 + 3);
	// Probability of carry here or that this>P is very very unlikely
	bits64[4] = 0;
}
// ------------------------------------------------
void Int::ModSquareK1(Int* a) 
{
	unsigned char c;
	uint64_t t[NB64BLOCK], SL, SH, r512[8], t1, t2;

	//k=0
	r512[0] = _umul128(a->bits64[0], a->bits64[0], &t[1]); 

	//k=1
	t[3] = _umul128(a->bits64[0], a->bits64[1], &t[4]);
	c = _addcarry_u64(0, t[3], t[3], &t[3]);
	c = _addcarry_u64(c, t[4], t[4], &t[4]);
	c = _addcarry_u64(c, 0, 0, &t1);
	c = _addcarry_u64(0, t[1], t[3], &t[3]);
	c = _addcarry_u64(c, t[4], 0, &t[4]);
	c = _addcarry_u64(c, t1, 0, &t1);
	r512[1] = t[3];

	//k=2
	t[0] = _umul128(a->bits64[0], a->bits64[2], &t[1]);
	c = _addcarry_u64(0, t[0], t[0], &t[0]);
	c = _addcarry_u64(c, t[1], t[1], &t[1]);
	c = _addcarry_u64(c, 0, 0, &t2);

	SL = _umul128(a->bits64[1], a->bits64[1], &SH);
	c = _addcarry_u64(0, t[0], SL, &t[0]);
	c = _addcarry_u64(c, t[1], SH, &t[1]);
	c = _addcarry_u64(c, t2, 0, &t2);
	c = _addcarry_u64(0, t[0], t[4], &t[0]);
	c = _addcarry_u64(c, t[1], t1, &t[1]);
	c = _addcarry_u64(c, t2, 0, &t2);
	r512[2] = t[0];

	//k=3
	t[3] = _umul128(a->bits64[0], a->bits64[3], &t[4]);
	SL = _umul128(a->bits64[1], a->bits64[2], &SH);

	c = _addcarry_u64(0, t[3], SL, &t[3]);
	c = _addcarry_u64(c, t[4], SH, &t[4]);
	c = _addcarry_u64(c, 0, 0, &t1);
	t1 += t1;
	c = _addcarry_u64(0, t[3], t[3], &t[3]);
	c = _addcarry_u64(c, t[4], t[4], &t[4]);
	c = _addcarry_u64(c, t1, 0, &t1);
	c = _addcarry_u64(0, t[3], t[1], &t[3]);
	c = _addcarry_u64(c, t[4], t2, &t[4]);
	c = _addcarry_u64(c, t1, 0, &t1);
	r512[3] = t[3];

	//k=4
	t[0] = _umul128(a->bits64[1], a->bits64[3], &t[1]);
	c = _addcarry_u64(0, t[0], t[0], &t[0]);
	c = _addcarry_u64(c, t[1], t[1], &t[1]);
	c = _addcarry_u64(c, 0, 0, &t2);

	SL = _umul128(a->bits64[2], a->bits64[2], &SH);
	c = _addcarry_u64(0, t[0], SL, &t[0]);
	c = _addcarry_u64(c, t[1], SH, &t[1]);
	c = _addcarry_u64(c, t2, 0, &t2);
	c = _addcarry_u64(0, t[0], t[4], &t[0]);
	c = _addcarry_u64(c, t[1], t1, &t[1]);
	c = _addcarry_u64(c, t2, 0, &t2);
	r512[4] = t[0];

	//k=5
	t[3] = _umul128(a->bits64[2], a->bits64[3], &t[4]);
	c = _addcarry_u64(0, t[3], t[3], &t[3]);
	c = _addcarry_u64(c, t[4], t[4], &t[4]);
	c = _addcarry_u64(c, 0, 0, &t1);
	c = _addcarry_u64(0, t[3], t[1], &t[3]);
	c = _addcarry_u64(c, t[4], t2, &t[4]);
	c = _addcarry_u64(c, t1, 0, &t1);
	r512[5] = t[3];

	//k=6
	t[0] = _umul128(a->bits64[3], a->bits64[3], &t[1]);
	c = _addcarry_u64(0, t[0], t[4], &t[0]);
	c = _addcarry_u64(c, t[1], t1, &t[1]);
	r512[6] = t[0];

	//k=7
	r512[7] = t[1];

	// Reduce from 512 to 320 
	imm_umul(r512 + 4, 0x1000003D1ULL, t);
	c = _addcarry_u64(0, r512[0], t[0], r512 + 0);
	c = _addcarry_u64(c, r512[1], t[1], r512 + 1);
	c = _addcarry_u64(c, r512[2], t[2], r512 + 2);
	c = _addcarry_u64(c, r512[3], t[3], r512 + 3);

	// Reduce from 320 to 256 
	// No overflow possible here t[4]+c<=0x1000003D1ULL
	SL = _umul128(t[4] + c, 0x1000003D1ULL, &SH);
	c = _addcarry_u64(0, r512[0], SL, bits64 + 0);
	c = _addcarry_u64(c, r512[1], SH, bits64 + 1);
	c = _addcarry_u64(c, r512[2], 0, bits64 + 2);
	c = _addcarry_u64(c, r512[3], 0, bits64 + 3);
	// Probability of carry here or that this>P is very very unlikely
	bits64[4] = 0;
}

//============== INT->INTMOD.CPP - end ================================================================================
//============== POINT - start ================================================================================
class Point
{
	public:
		Point();
		Point(const Point &p);
		~Point();

		void Clear();
		void Reduce();
		std::string toString();

		Int x, y, z;
};

//------------------
Point::Point(){}
Point::~Point(){}
//------------------
Point::Point(const Point &p){
    x.Set((Int *)&p.x);
    y.Set((Int *)&p.y);
    z.Set((Int *)&p.z);
}
//------------------
void Point::Clear(){
    x.SetInt32(0);
    y.SetInt32(0);
    z.SetInt32(0);
}
//------------------
void Point::Reduce(){
    Int i(&z);
    i.ModInv();
    x.ModMul(&x, &i);
    y.ModMul(&y, &i);
    z.SetInt32(1);
}
//------------------
std::string Point::toString()
{
    std::string ret;
    ret  = "\nX=" + x.GetBase16() + "\n";
    ret += "Y=" + y.GetBase16() + "\n";
    ret += "Z=" + z.GetBase16() + "\n";
    return ret;
}
//============== POINT - end ================================================================================
//============== SECP256K - start ================================================================================
class Secp256K1
{
	public:
		Secp256K1();
		~Secp256K1();
		void Init();
		Point ComputePublicKey(Int* privKey);
		// void GetHash160(int type, bool isCompressed, Point& pubKey, unsigned char* hash);

		Point Add2(Point& p1, Point& p2);
		Point AddDirect(Point& p1, Point& p2);
		Point DoubleDirect(Point& p);
		Point G;                 // Generator



		void pubkeyToHash160keccak(Point& pubKey, uint32_t* hash);
		std::string pubkeyToAddr(Point& pubKey);
		std::string hash160keccakToAddr(uint32_t* hash);


	private:
		uint8_t GetByte(std::string& str, int idx);
		Point GTable[256 * 32];     // Generator table
};
//-------------------------------------------------------
Secp256K1::Secp256K1(){}
Secp256K1::~Secp256K1(){}

//-------------------------------------------------------
void Secp256K1::Init() 
{
	// Prime for the finite field
	Int P;
	P.SetBase16("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F");

	// Set up field
	Int::SetupField(&P);

	// Generator point // Điểm cơ sở = 04 79BE667E F9DCBBAC 55A06295 CE870B07 029BFCDB 2DCE28D9 59F2815B 16F81798 483ADA77 26A3C465 5DA4FBFC 0E1108A8 FD17B448 A6855419 9C47D08F FB10D4B8
	G.x.SetBase16("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798");
	G.y.SetBase16("483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8");
	G.z.SetInt32(1);

	// Compute Generator table
	Point N(G);
	for (int i = 0; i < 32; i++) {
		GTable[i * 256] = N;
		N = DoubleDirect(N);
		for (int j = 1; j < 255; j++) {
			GTable[i * 256 + j] = N;
			N = AddDirect(N, GTable[i * 256]);
		}
		GTable[i * 256 + 255] = N; // Dummy point for check function
	}
}
//-------------------------------------------------------

Point Secp256K1::ComputePublicKey(Int* privKey)
{
	int i = 0;
	uint8_t b;
	Point Q;

	Q.Clear(); //Q ( x=0 , y=0 , z=0 )

	// Search first significant byte
	for (i = 0; i < 32; i++) {	
		b = privKey->GetByte(i);
		if (b){ break; }	}
		
	Q = GTable[256 * i + (b - 1)];
	i++;

	for (; i < 32; i++) {	
		b = privKey->GetByte(i);
		if (b){ Q = Add2(Q, GTable[256 * i + (b - 1)]); }	}

	Q.Reduce();

	return Q;
}
//-------------------------------------------------------

uint8_t Secp256K1::GetByte(std::string& str, int idx)
{
	char tmp[3];
	int  val;

	tmp[0] = str.data()[2 * idx];
	tmp[1] = str.data()[2 * idx + 1];
	tmp[2] = 0;

	if (sscanf(tmp, "%X", &val) != 1) {	printf("ParsePublicKeyHex: Error invalid public key specified (unexpected hexadecimal digit)\n");	exit(-1);	}

	return (uint8_t)val;
}
//-------------------------------------------------------
//-------------------------------------------------------

Point Secp256K1::AddDirect(Point& p1, Point& p2)
{
	Int _s, _p, dy, dx;
	Point r;
	r.z.SetInt32(1);

	dy.ModSub(&p2.y, &p1.y);
	dx.ModSub(&p2.x, &p1.x);
	dx.ModInv();
	_s.ModMulK1(&dy, &dx);    // s = (p2.y-p1.y)*inverse(p2.x-p1.x);

	_p.ModSquareK1(&_s);       // _p = pow2(s)

	r.x.ModSub(&_p, &p1.x);
	r.x.ModSub(&p2.x);       // rx = pow2(s) - p1.x - p2.x;

	r.y.ModSub(&p2.x, &r.x);
	r.y.ModMulK1(&_s);
	r.y.ModSub(&p2.y);       // ry = - p2.y - s*(ret.x-p2.x);

	return r;
}
//-------------------------------------------------------

Point Secp256K1::Add2(Point& p1, Point& p2)
{
	// P2.z = 1
	Int u, v, u1, v1, vs2, vs3, us2, a, us2w, vs2v2, vs3u2, _2vs2v2;
	Point r;

	u1.ModMulK1(&p2.y, &p1.z);
	v1.ModMulK1(&p2.x, &p1.z);
	u.ModSub(&u1, &p1.y);
	v.ModSub(&v1, &p1.x);
	us2.ModSquareK1(&u);
	vs2.ModSquareK1(&v);
	vs3.ModMulK1(&vs2, &v);
	us2w.ModMulK1(&us2, &p1.z);
	vs2v2.ModMulK1(&vs2, &p1.x);
	_2vs2v2.ModAdd(&vs2v2, &vs2v2);
	a.ModSub(&us2w, &vs3);
	a.ModSub(&_2vs2v2);

	r.x.ModMulK1(&v, &a);

	vs3u2.ModMulK1(&vs3, &p1.y);
	r.y.ModSub(&vs2v2, &a);
	r.y.ModMulK1(&r.y, &u);
	r.y.ModSub(&vs3u2);

	r.z.ModMulK1(&vs3, &p1.z);

	return r;
}

//-------------------------------------------------------
 
Point Secp256K1::DoubleDirect(Point& p)
{
	Int _s, _p, a;
	Point r;
	r.z.SetInt32(1);

	_s.ModMulK1(&p.x, &p.x);
	_p.ModAdd(&_s, &_s);
	_p.ModAdd(&_s);

	a.ModAdd(&p.y, &p.y);
	a.ModInv();
	_s.ModMulK1(&_p, &a);    // s = (3*pow2(p.x))*inverse(2*p.y);

	_p.ModMulK1(&_s, &_s);
	a.ModAdd(&p.x, &p.x);
	a.ModNeg();
	r.x.ModAdd(&a, &_p);   // rx = pow2(s) + neg(2*p.x);

	a.ModSub(&r.x, &p.x);

	_p.ModMulK1(&a, &_s);
	r.y.ModAdd(&_p, &p.y);
	r.y.ModNeg();           // ry = neg(p.y + s*(ret.x+neg(p.x)));

	return r;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//====================================== keccak160 - start ============================================== 
void keccak160(uint64_t* x, uint64_t* y, uint32_t* hash);


// ---------------------------------------------------------------------------------
// KECCAK/SHA3
// ---------------------------------------------------------------------------------

typedef union {
	uint8_t b[200];
	uint64_t q[25];
	uint32_t d[50];
} _KECCAK_STATE;

const uint64_t _KECCAKF_RNDC[24] = {
	0x0000000000000001ULL, 0x0000000000008082ULL, 0x800000000000808aULL,
	0x8000000080008000ULL, 0x000000000000808bULL, 0x0000000080000001ULL,
	0x8000000080008081ULL, 0x8000000000008009ULL, 0x000000000000008aULL,
	0x0000000000000088ULL, 0x0000000080008009ULL, 0x000000008000000aULL,
	0x000000008000808bULL, 0x800000000000008bULL, 0x8000000000008089ULL,
	0x8000000000008003ULL, 0x8000000000008002ULL, 0x8000000000000080ULL,
	0x000000000000800aULL, 0x800000008000000aULL, 0x8000000080008081ULL,
	0x8000000000008080ULL, 0x0000000080000001ULL, 0x8000000080008008ULL
};

#define ROTL64(a,b) (((a) << (b)) | ((a) >> (64 - b)))

#define bswap32(v) (((v) >> 24) | (((v) >> 8) & 0xff00) | (((v) << 8) & 0xff0000) | ((v) << 24))

void keccak160(uint64_t* x, uint64_t* y, uint32_t* hash)
{
	_KECCAK_STATE e;
	memset(&e, 0, 200);
	uint32_t* X = (uint32_t*)x;
	uint32_t* Y = (uint32_t*)y;
	e.d[0] = bswap32(X[7]);
	e.d[1] = bswap32(X[6]);
	e.d[2] = bswap32(X[5]);
	e.d[3] = bswap32(X[4]);
	e.d[4] = bswap32(X[3]);
	e.d[5] = bswap32(X[2]);
	e.d[6] = bswap32(X[1]);
	e.d[7] = bswap32(X[0]);
	e.d[8] = bswap32(Y[7]);
	e.d[9] = bswap32(Y[6]);
	e.d[10] = bswap32(Y[5]);
	e.d[11] = bswap32(Y[4]);
	e.d[12] = bswap32(Y[3]);
	e.d[13] = bswap32(Y[2]);
	e.d[14] = bswap32(Y[1]);
	e.d[15] = bswap32(Y[0]);
	uint64_t* s = e.q; 
	e.d[16] ^= 0x01;
	e.d[33] ^= 0x80000000;
	int i;
	uint64_t v, w, t[5], u[5];
	for (i = 0; i < 24; i++) {
		/* theta: c = a[0,i] ^ a[1,i] ^ .. a[4,i] */
		t[0] = s[0] ^ s[5] ^ s[10] ^ s[15] ^ s[20];
		t[1] = s[1] ^ s[6] ^ s[11] ^ s[16] ^ s[21];
		t[2] = s[2] ^ s[7] ^ s[12] ^ s[17] ^ s[22];
		t[3] = s[3] ^ s[8] ^ s[13] ^ s[18] ^ s[23];
		t[4] = s[4] ^ s[9] ^ s[14] ^ s[19] ^ s[24];
		/* theta: d[i] = c[i+4] ^ ROTL64(c[i+1],1) */
		u[0] = t[4] ^ ROTL64(t[1], 1);
		u[1] = t[0] ^ ROTL64(t[2], 1);
		u[2] = t[1] ^ ROTL64(t[3], 1);
		u[3] = t[2] ^ ROTL64(t[4], 1);
		u[4] = t[3] ^ ROTL64(t[0], 1);
		/* theta: a[0,i], a[1,i], .. a[4,i] ^= d[i] */
		s[0] ^= u[0]; s[5] ^= u[0]; s[10] ^= u[0]; s[15] ^= u[0]; s[20] ^= u[0];
		s[1] ^= u[1]; s[6] ^= u[1]; s[11] ^= u[1]; s[16] ^= u[1]; s[21] ^= u[1];
		s[2] ^= u[2]; s[7] ^= u[2]; s[12] ^= u[2]; s[17] ^= u[2]; s[22] ^= u[2];
		s[3] ^= u[3]; s[8] ^= u[3]; s[13] ^= u[3]; s[18] ^= u[3]; s[23] ^= u[3];
		s[4] ^= u[4]; s[9] ^= u[4]; s[14] ^= u[4]; s[19] ^= u[4]; s[24] ^= u[4];
		/* rho pi: b[..] = ROTL64(a[..], ..) */
		v = s[1];
		s[1] = ROTL64(s[6], 44);
		s[6] = ROTL64(s[9], 20);
		s[9] = ROTL64(s[22], 61);
		s[22] = ROTL64(s[14], 39);
		s[14] = ROTL64(s[20], 18);
		s[20] = ROTL64(s[2], 62);
		s[2] = ROTL64(s[12], 43);
		s[12] = ROTL64(s[13], 25);
		s[13] = ROTL64(s[19], 8);
		s[19] = ROTL64(s[23], 56);
		s[23] = ROTL64(s[15], 41);
		s[15] = ROTL64(s[4], 27);
		s[4] = ROTL64(s[24], 14);
		s[24] = ROTL64(s[21], 2);
		s[21] = ROTL64(s[8], 55);
		s[8] = ROTL64(s[16], 45);
		s[16] = ROTL64(s[5], 36);
		s[5] = ROTL64(s[3], 28);
		s[3] = ROTL64(s[18], 21);
		s[18] = ROTL64(s[17], 15);
		s[17] = ROTL64(s[11], 10);
		s[11] = ROTL64(s[7], 6);
		s[7] = ROTL64(s[10], 3);
		s[10] = ROTL64(v, 1);
		/* chi: a[i,j] ^= ~b[i,j+1] & b[i,j+2] */
		v = s[0]; w = s[1]; s[0] ^= (~w) & s[2]; s[1] ^= (~s[2]) & s[3]; s[2] ^= (~s[3]) & s[4]; s[3] ^= (~s[4]) & v; s[4] ^= (~v) & w;
		v = s[5]; w = s[6]; s[5] ^= (~w) & s[7]; s[6] ^= (~s[7]) & s[8]; s[7] ^= (~s[8]) & s[9]; s[8] ^= (~s[9]) & v; s[9] ^= (~v) & w;
		v = s[10]; w = s[11]; s[10] ^= (~w) & s[12]; s[11] ^= (~s[12]) & s[13]; s[12] ^= (~s[13]) & s[14]; s[13] ^= (~s[14]) & v; s[14] ^= (~v) & w;
		v = s[15]; w = s[16]; s[15] ^= (~w) & s[17]; s[16] ^= (~s[17]) & s[18]; s[17] ^= (~s[18]) & s[19]; s[18] ^= (~s[19]) & v; s[19] ^= (~v) & w;
		v = s[20]; w = s[21]; s[20] ^= (~w) & s[22]; s[21] ^= (~s[22]) & s[23]; s[22] ^= (~s[23]) & s[24]; s[23] ^= (~s[24]) & v; s[24] ^= (~v) & w;
		/* iota: a[0,0] ^= round constant */
		s[0] ^= _KECCAKF_RNDC[i];
	}
	hash[0] = e.d[3];
	hash[1] = e.d[4];
	hash[2] = e.d[5];
	hash[3] = e.d[6];
	hash[4] = e.d[7];

}



//====================================== keccak160 - end ============================================== 
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

void Secp256K1::pubkeyToHash160keccak(Point& pubKey, uint32_t* hash)
{
	keccak160(pubKey.x.bits64, pubKey.y.bits64, (uint32_t*)hash);
}
//----------------------------------------------------------------------------------

std::string Secp256K1::hash160keccakToAddr(uint32_t* hash) //hiiu - here - codenow
{
	char tmp[3];
	std::string ret;

	ret.append("0x");
	for (int i = 0; i < 20; i++) {
		sprintf(tmp, "%02x", ((uint8_t*)hash)[i]);
		ret.append(tmp);
	}

	return ret;
}
//----------------------------------------------------------------------------------

std::string Secp256K1::pubkeyToAddr(Point& pubKey) // hiiu - here - codenow
{
	uint32_t hash[5];
	char tmp[3];
	std::string ret;

	keccak160(pubKey.x.bits64, pubKey.y.bits64, hash);

	ret.append("0x");
	for (int i = 0; i < 20; i++) {
		sprintf(tmp, "%02x", ((uint8_t*)hash)[i]);
		ret.append(tmp);
	}

	return ret;
}
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

//==============================================================================================

//============== HIIU::BITCOIN - start =============================================================================
#include <cassert>

class ETH
{
public:
	Point privToPubkey(char* p_hex);
	void privToHash160keccak(char* p_hex,uint32_t* hash);
	std::string privToAddr(char* p_hex);

	void pubkeyToHash160keccak(Point& pubKey, uint32_t* hash);
	std::string pubkeyToAddr(Point& pubKey);

	std::string hash160keccakToAddr(uint32_t* hash);

	void decodeAddrToHash160keccak(char* p_hex);
};
//----------------------------------------------------------------------------------

Point ETH::privToPubkey(char* p_hex)
{
	char* priv_hex = p_hex;

	Int priv; 
	priv.SetBase16(priv_hex); 
	 
	Secp256K1*	secp = new Secp256K1();   
	secp->Init();	

	Point pubKey = secp->ComputePublicKey(&priv); //

	delete secp;

	return pubKey;
}
//----------------------------------------------------------------------------------

void ETH::privToHash160keccak(char* p_hex,uint32_t* hash)
{
	char* priv_hex = p_hex;

	Int priv; 
	priv.SetBase16(priv_hex); 
	 
	Secp256K1*	secp = new Secp256K1();   
	secp->Init();	

	Point pubKey = secp->ComputePublicKey(&priv); //
	keccak160(pubKey.x.bits64, pubKey.y.bits64, (uint32_t*)hash);

	delete secp;
}
//----------------------------------------------------------------------------------

std::string ETH::privToAddr(char* p_hex)
{
	char* priv_hex = p_hex;

	Int priv; 
	priv.SetBase16(priv_hex); 
	 
	Secp256K1*	secp = new Secp256K1();   
	secp->Init();	
	
	//priv -> pubkey
	Point pubKey = secp->ComputePublicKey(&priv); //

	delete secp;

	uint32_t hash[5];
	char tmp[3];
	std::string addr;

	//pubkey -> hash160keccak
	keccak160(pubKey.x.bits64, pubKey.y.bits64, hash);

	//hash160keccak -> addr 
	addr.append("0x");
	for (int i = 0; i < 20; i++) {
		sprintf(tmp, "%02x", ((uint8_t*)hash)[i]);
		addr.append(tmp);
	}

	return addr;
}
//----------------------------------------------------------------------------------
void ETH::pubkeyToHash160keccak(Point& pubKey, uint32_t* hash)
{
	keccak160(pubKey.x.bits64, pubKey.y.bits64, (uint32_t*)hash);
}
//----------------------------------------------------------------------------------
std::string ETH::pubkeyToAddr(Point& pubKey)
{
	uint32_t hash[5];
	char tmp[3];
	std::string addr;

	keccak160(pubKey.x.bits64, pubKey.y.bits64, hash);

	addr.append("0x");
	for (int i = 0; i < 20; i++) {
		sprintf(tmp, "%02x", ((uint8_t*)hash)[i]);
		addr.append(tmp);
	}

	return addr;
}

//----------------------------------------------------------------------------------
std::string ETH::hash160keccakToAddr(uint32_t* hash)
{
	char tmp[3];
	std::string addr;

	addr.append("0x");
	for (int i = 0; i < 20; i++) {
		sprintf(tmp, "%02x", ((uint8_t*)hash)[i]);
		addr.append(tmp);
	}

	return addr;
}
//----------------------------------------------------------------------------------
void ETH::decodeAddrToHash160keccak(char* addr)
{
	string address = addr;

	std::vector<unsigned char> hashORxpoint;
	address.erase(0, 2); // 
	for (int i = 0; i < 40; i += 2) {
		uint8_t c = 0;
		for (size_t j = 0; j < 2; j++) {
			uint32_t c0 = (uint32_t)address[i + j];
			uint8_t c2 = (uint8_t)strtol((char*)&c0, NULL, 16);
			if (j == 0) 
				c2 = c2 << 4; 
			c |= c2;
		}
		hashORxpoint.push_back(c);
	}
	
	uint32_t hash160Keccak[5];
	for (size_t i = 0; i < hashORxpoint.size(); i++) {
		((uint8_t*)hash160Keccak)[i] = hashORxpoint.at(i);
	}

	for (int i = 0; i < 5; i++){	printf("\n --- hash160Keccak[] : %d ", hash160Keccak[i]); }
}
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

//============== HIIU::BITCOIN - end =============================================================================

int main()
{

	std::cout<<std::endl<<" -------------------------------------------- ";

	ETH eth;
	eth.decodeAddrToHash160keccak("0x73B1D0D7EEa322bAc6cfE7B9329daf79ca1Ac76A");

	// std::cout << "\nAddr : " << addr;

	printf("\n\n\n\n\n");

	//---------------------------------------- 
	//priv --> hash160 
	uint32_t h[5];
	eth.privToHash160keccak("5a8425eac0e9479f3f48cd9cd5d10b0ebba6b22bcd49d2caf678622ae362eb40",h);
	for (int i = 0; i < 5; i++){	printf("\n --- h[] : %d ", h[i]); }


	printf("\n\n\n\n\n");

	//---------------------------------------- 

	printf("\n\n\n\n\n");


	//-- addr -> hash160[5] -------------------------------------- 



	printf("\n\n\n\n\n");
	//---------------------------------------- 

    return 0;
}
