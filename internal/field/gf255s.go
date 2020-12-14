package field

// This file implements computations in the fields of integers
// modulo p = 2^255 - 3957.

// =======================================================================
// Field GF255s: integers modulo p = 2^255 - 3957
type GF255s [4]uint64

const mq255s uint64 = 3957

// Field element of value 0.
var GF255s_ZERO = GF255s{0, 0, 0, 0}

// Field element of value 1.
var GF255s_ONE = GF255s{1, 0, 0, 0}

// Field element of value 2.
var GF255s_TWO = GF255s{2, 0, 0, 0}

// Field element of value 1/2^508 (used internally for inversions).
var GF255s_INVT508 = GF255s{
	0xC7E0DEC400D7BDB6, 0xCCABD4771F6FB10F,
	0x940F23A06B74BE6E, 0x1C45852F33548365}

// d <- a
func (d *GF255s) Set(a *GF255s) *GF255s {
	copy(d[:], a[:])
	return d
}

// d <- a + b
func (d *GF255s) Add(a, b *GF255s) *GF255s {
	gf_add((*[4]uint64)(d), (*[4]uint64)(a), (*[4]uint64)(b), mq255s)
	return d
}

// d <- a - b
func (d *GF255s) Sub(a, b *GF255s) *GF255s {
	gf_sub((*[4]uint64)(d), (*[4]uint64)(a), (*[4]uint64)(b), mq255s)
	return d
}

// d <- -a
func (d *GF255s) Neg(a *GF255s) *GF255s {
	gf_neg((*[4]uint64)(d), (*[4]uint64)(a), mq255s)
	return d
}

// If ctl == 1:  d <- a
// If ctl == 0:  d <- b
// ctl MUST be 0 or 1
func (d *GF255s) Select(a, b *GF255s, ctl uint64) *GF255s {
	gf_select((*[4]uint64)(d), (*[4]uint64)(a), (*[4]uint64)(b), ctl)
	return d
}

// d <- d OR (a AND mask)
// mask value should be 0 or 0xFFFFFFFFFFFFFFFF
func (d *GF255s) CondOrFrom(a *GF255s, mask uint64) *GF255s {
	d[0] |= mask & a[0]
	d[1] |= mask & a[1]
	d[2] |= mask & a[2]
	d[3] |= mask & a[3]
	return d
}

// If ctl == 1:  d <- -a
// If ctl == 0:  d <- a
// ctl MUST be 0 or 1
func (d *GF255s) CondNeg(a *GF255s, ctl uint64) *GF255s {
	gf_condneg((*[4]uint64)(d), (*[4]uint64)(a), mq255s, ctl)
	return d
}

// d <- a*b
func (d *GF255s) Mul(a, b *GF255s) *GF255s {
	gf_mul((*[4]uint64)(d), (*[4]uint64)(a), (*[4]uint64)(b), mq255s)
	return d
}

// d <- a^2
func (d *GF255s) Sqr(a *GF255s) *GF255s {
	gf_sqr((*[4]uint64)(d), (*[4]uint64)(a), mq255s)
	return d
}

// d <- a^(2^n) for any n >= 0
// This is constant-time with regard to a and d, but not to n.
func (d *GF255s) SqrX(a *GF255s, n uint) *GF255s {
	gf_sqr_x((*[4]uint64)(d), (*[4]uint64)(a), n, mq255s)
	return d
}

// d <- a/2
func (d *GF255s) Half(a *GF255s) *GF255s {
	gf_half((*[4]uint64)(d), (*[4]uint64)(a), mq255s)
	return d
}

// d <- a*2^n
// n must be between 1 and 15 (inclusive)
func (d *GF255s) Lsh(a *GF255s, n uint) *GF255s {
	gf_lsh((*[4]uint64)(d), (*[4]uint64)(a), n, mq255s)
	return d
}

// d <- 1/a   (if a == 0, this sets d to 0)
func (d *GF255s) Inv(a *GF255s) *GF255s {
	gf_inv_scaled((*[4]uint64)(d), (*[4]uint64)(a), mq255s)
	gf_mul((*[4]uint64)(d), (*[4]uint64)(d), (*[4]uint64)(&GF255s_INVT508), mq255s)
	return d
}

// Returns 1 if d == 0, or 0 otherwise.
func (d *GF255s) IsZero() uint64 {
	return gf_iszero((*[4]uint64)(d), mq255s)
}

// Returns 1 if d == a, or 0 otherwise.
func (d *GF255s) Eq(a *GF255s) uint64 {
	return gf_eq((*[4]uint64)(d), (*[4]uint64)(a), mq255s)
}

// Legendre symbol computation; return value:
//    0  if d == 0
//    1  if d != 0 and is a quadratic residue
//   -1  if d != 0 and is a not a quadratic residue
// Value is returned as uint64, i.e. 0xFFFFFFFFFFFFFFFF for non-squares.
func (d *GF255s) Legendre() uint64 {
	return gf_legendre((*[4]uint64)(d), mq255s)
}

// Square root computation. If the source value (a) is a quadratic
// residue, then this function sets this object (d) to a square root
// of a, and returns 1; otherwise, it sets d to zero and returns 0.
// When the input is a square, this function returns the square root
// whose least significant bit (as an integer in the 0..p-1 range) is zero.
func (d *GF255s) Sqrt(a *GF255s) uint64 {
	// Since p = 3 mod 4, we get a candidate square root by raising
	// the input to the power (p+1)/4. We use an addition chain with
	// 252 squarings and 12 extra multiplications.
	var x, x2, y [4]uint64

	// x2 <- a^3
	gf_sqr(&x2, (*[4]uint64)(a), mq255s)
	gf_mul(&x2, &x2, (*[4]uint64)(a), mq255s)

	// x <- a^(2^3-1)
	gf_sqr(&x, &x2, mq255s)
	gf_mul(&x, &x, (*[4]uint64)(a), mq255s)

	// x <- a^(2^9-1)
	gf_sqr_x(&y, &x, 3, mq255s)
	gf_mul(&y, &y, &x, mq255s)
	gf_sqr_x(&y, &y, 3, mq255s)
	gf_mul(&x, &y, &x, mq255s)

	// x <- a^(2^27-1)
	gf_sqr_x(&y, &x, 9, mq255s)
	gf_mul(&y, &y, &x, mq255s)
	gf_sqr_x(&y, &y, 9, mq255s)
	gf_mul(&x, &y, &x, mq255s)

	// x <- a^(2^81-1)
	gf_sqr_x(&y, &x, 27, mq255s)
	gf_mul(&y, &y, &x, mq255s)
	gf_sqr_x(&y, &y, 27, mq255s)
	gf_mul(&x, &y, &x, mq255s)

	// x <- a^(2^243-1)
	gf_sqr_x(&y, &x, 81, mq255s)
	gf_mul(&y, &y, &x, mq255s)
	gf_sqr_x(&y, &y, 81, mq255s)
	gf_mul(&x, &y, &x, mq255s)

	// x <- a^(2^253 - 1024 + 35)
	gf_sqr_x(&x, &x, 5, mq255s)
	gf_mul(&x, &x, (*[4]uint64)(a), mq255s)
	gf_sqr_x(&x, &x, 5, mq255s)
	gf_mul(&x, &x, &x2, mq255s)

	// Verify the result. If not square, then set the result to 0.
	gf_sqr(&y, &x, mq255s)
	qr := gf_eq(&y, (*[4]uint64)(a), mq255s)
	for i := 0; i < 4; i++ {
		x[i] &= -qr
	}

	// Normalize the result, and negate the value if the least
	// significant bit is 1.
	gf_norm(&x, &x, mq255s)
	gf_condneg(&x, &x, mq255s, x[0]&1)

	// Return the result.
	copy(d[:], x[:])
	return qr
}

// Encode element into exactly 32 bytes. The encoding is appended to the
// provided slice, and the resulting slice is returned. The extension is
// done in place if the provided slice has enough capacity.
func (d *GF255s) Encode(dst []byte) []byte {
	return gf_encode(dst, (*[4]uint64)(d), mq255s)
}

// Decode element from 32 bytes. If the source is invalid (out of range),
// then the decoded value is zero, and 0 is returned; otherwise, 1 is
// returned.
func (d *GF255s) Decode(src []byte) uint64 {
	return gf_decode((*[4]uint64)(d), src, mq255s)
}

// Decode element from bytes. The input bytes are interpreted as an
// integer (unsigned little-endian convention) which is reduced modulo
// the field modulus. By definition, this process cannot fail.
func (d *GF255s) DecodeReduce(src []byte) *GF255s {
	gf_decodeReduce((*[4]uint64)(d), src, mq255s)
	return d
}
