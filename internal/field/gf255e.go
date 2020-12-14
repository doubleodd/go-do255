package field

// This file implements computations in the field of integers
// modulo p = 2^255 - 18651.

// =======================================================================
// Field GF255e: integers modulo p = 2^255 - 18651
type GF255e [4]uint64

const mq255e uint64 = 18651

// Field element of value 0.
var GF255e_ZERO = GF255e{0, 0, 0, 0}

// Field element of value 1.
var GF255e_ONE = GF255e{1, 0, 0, 0}

// Field element of value 3.
var GF255e_THREE = GF255e{3, 0, 0, 0}

// Field element of value 7.
var GF255e_SEVEN = GF255e{7, 0, 0, 0}

// Field element of value 8.
var GF255e_EIGHT = GF255e{8, 0, 0, 0}

// Field element of value 27.
var GF255e_TWENTYSEVEN = GF255e{27, 0, 0, 0}

// Field element of value 176.
var GF255e_HUNDREDSEVENTYSIX = GF255e{176, 0, 0, 0}

// Field element of value 308.
var GF255e_THREEHUNDREDEIGHT = GF255e{308, 0, 0, 0}

// Field element of value 343.
var GF255e_THREEHUNDREDFORTYTHREE = GF255e{343, 0, 0, 0}

// Field element of value 1/2^508 (used internally for inversions).
var GF255e_INVT508 = GF255e{
	0xD40D5B5D2BE1CF5D, 0x3B3573987282DD51,
	0x3ECCB22800EED6AE, 0x44F35C558E8FAC0B}

// d <- a
func (d *GF255e) Set(a *GF255e) *GF255e {
	copy(d[:], a[:])
	return d
}

// d <- a + b
func (d *GF255e) Add(a, b *GF255e) *GF255e {
	gf_add((*[4]uint64)(d), (*[4]uint64)(a), (*[4]uint64)(b), mq255e)
	return d
}

// d <- a - b
func (d *GF255e) Sub(a, b *GF255e) *GF255e {
	gf_sub((*[4]uint64)(d), (*[4]uint64)(a), (*[4]uint64)(b), mq255e)
	return d
}

// d <- -a
func (d *GF255e) Neg(a *GF255e) *GF255e {
	gf_neg((*[4]uint64)(d), (*[4]uint64)(a), mq255e)
	return d
}

// If ctl == 1:  d <- a
// If ctl == 0:  d <- b
// ctl MUST be 0 or 1
func (d *GF255e) Select(a, b *GF255e, ctl uint64) *GF255e {
	gf_select((*[4]uint64)(d), (*[4]uint64)(a), (*[4]uint64)(b), ctl)
	return d
}

// d <- d OR (a AND mask)
// mask value should be 0 or 0xFFFFFFFFFFFFFFFF
func (d *GF255e) CondOrFrom(a *GF255e, mask uint64) *GF255e {
	d[0] |= mask & a[0]
	d[1] |= mask & a[1]
	d[2] |= mask & a[2]
	d[3] |= mask & a[3]
	return d
}

// If ctl == 1:  d <- -a
// If ctl == 0:  d <- a
// ctl MUST be 0 or 1
func (d *GF255e) CondNeg(a *GF255e, ctl uint64) *GF255e {
	gf_condneg((*[4]uint64)(d), (*[4]uint64)(a), mq255e, ctl)
	return d
}

// d <- a*b
func (d *GF255e) Mul(a, b *GF255e) *GF255e {
	gf_mul((*[4]uint64)(d), (*[4]uint64)(a), (*[4]uint64)(b), mq255e)
	return d
}

// d <- a^2
func (d *GF255e) Sqr(a *GF255e) *GF255e {
	gf_sqr((*[4]uint64)(d), (*[4]uint64)(a), mq255e)
	return d
}

// d <- a^(2^n) for any n >= 0
// This is constant-time with regard to a and d, but not to n.
func (d *GF255e) SqrX(a *GF255e, n uint) *GF255e {
	gf_sqr_x((*[4]uint64)(d), (*[4]uint64)(a), n, mq255e)
	return d
}

// d <- a/2
func (d *GF255e) Half(a *GF255e) *GF255e {
	gf_half((*[4]uint64)(d), (*[4]uint64)(a), mq255e)
	return d
}

// d <- a*2^n
// n must be between 1 and 15 (inclusive)
func (d *GF255e) Lsh(a *GF255e, n uint) *GF255e {
	gf_lsh((*[4]uint64)(d), (*[4]uint64)(a), n, mq255e)
	return d
}

// d <- 1/a   (if a == 0, this sets d to 0)
func (d *GF255e) Inv(a *GF255e) *GF255e {
	gf_inv_scaled((*[4]uint64)(d), (*[4]uint64)(a), mq255e)
	gf_mul((*[4]uint64)(d), (*[4]uint64)(d), (*[4]uint64)(&GF255e_INVT508), mq255e)
	return d
}

// Returns 1 if d == 0, or 0 otherwise.
func (d *GF255e) IsZero() uint64 {
	return gf_iszero((*[4]uint64)(d), mq255e)
}

// Returns 1 if d == a, or 0 otherwise.
func (d *GF255e) Eq(a *GF255e) uint64 {
	return gf_eq((*[4]uint64)(d), (*[4]uint64)(a), mq255e)
}

// Legendre symbol computation; return value:
//    0  if d == 0
//    1  if d != 0 and is a quadratic residue
//   -1  if d != 0 and is a not a quadratic residue
// Value is returned as uint64, i.e. 0xFFFFFFFFFFFFFFFF for non-squares.
func (d *GF255e) Legendre() uint64 {
	return gf_legendre((*[4]uint64)(d), mq255e)
}

// Square root computation. If the source value (a) is a quadratic
// residue, then this function sets this object (d) to a square root
// of a, and returns 1; otherwise, it sets d to zero and returns 0.
// When the input is a square, this function returns the square root
// whose least significant bit (as an integer in the 0..p-1 range) is zero.
func (d *GF255e) Sqrt(a *GF255e) uint64 {
	// Since p = 5 mod 8, we use Atkin's algorithm:
	//   b <- (2*a)^((p-5)/8)
	//   c <- 2*a*b^2
	//   return a*b*(c - 1)
	var b, c, e, x, x2, x96, y [4]uint64

	// e <- 2*a
	gf_lsh(&e, (*[4]uint64)(a), 1, mq255e)

	// Raise e to the power (p-5)/8. We use an addition chain with
	// 251 squarings and 13 extra multiplications:
	//   (p-5)/8 = (2^240-1)*2^12 + (2^2-1)*2^9 + (2^3-1)*2^5 + 2^2

	// x2 <- e^3
	gf_sqr(&x2, &e, mq255e)
	gf_mul(&x2, &x2, &e, mq255e)

	// x <- e^(2^4-1)
	gf_sqr_x(&x, &x2, 2, mq255e)
	gf_mul(&x, &x, &x2, mq255e)

	// x <- e^(2^8-1)
	gf_sqr_x(&y, &x, 4, mq255e)
	gf_mul(&x, &y, &x, mq255e)

	// x <- e^(2^16-1)
	gf_sqr_x(&y, &x, 8, mq255e)
	gf_mul(&x, &y, &x, mq255e)

	// x <- e^(2^48-1)
	gf_sqr_x(&y, &x, 16, mq255e)
	gf_mul(&y, &y, &x, mq255e)
	gf_sqr_x(&y, &y, 16, mq255e)
	gf_mul(&x, &y, &x, mq255e)

	// x96 <- e^(2^96-1)
	gf_sqr_x(&y, &x, 48, mq255e)
	gf_mul(&x96, &y, &x, mq255e)

	// x <- e^(2^240-1)
	gf_sqr_x(&y, &x96, 96, mq255e)
	gf_mul(&y, &y, &x96, mq255e)
	gf_sqr_x(&y, &y, 48, mq255e)
	gf_mul(&x, &y, &x, mq255e)

	// x <- e^((p-5)/8)
	gf_sqr_x(&x, &x, 3, mq255e)
	gf_mul(&x, &x, &x2, mq255e)
	gf_sqr_x(&x, &x, 2, mq255e)
	gf_mul(&x, &x, &e, mq255e)
	gf_sqr_x(&x, &x, 2, mq255e)
	gf_mul(&x, &x, &x2, mq255e)
	gf_sqr_x(&x, &x, 3, mq255e)
	gf_mul(&x, &x, &e, mq255e)
	gf_sqr_x(&b, &x, 2, mq255e)

	// We now have b = (2*a)^((p-5)/8).

	// c <- 2*a*b^2
	gf_sqr(&c, &b, mq255e)
	gf_mul(&c, &c, &e, mq255e)

	// x <- a*b*(c - 1)
	gf_sub(&x, &c, (*[4]uint64)(&GF255e_ONE), mq255e)
	gf_mul(&x, &x, (*[4]uint64)(a), mq255e)
	gf_mul(&x, &x, &b, mq255e)

	// Verify the result. If not square, then set the result to 0.
	gf_sqr(&y, &x, mq255e)
	qr := gf_eq(&y, (*[4]uint64)(a), mq255e)
	for i := 0; i < 4; i++ {
		x[i] &= -qr
	}

	// Normalize the result, and negate the value if the least
	// significant bit is 1.
	gf_norm(&x, &x, mq255e)
	gf_condneg(&x, &x, mq255e, x[0]&1)

	// Return the result.
	copy(d[:], x[:])
	return qr
}

// Encode element into exactly 32 bytes. The encoding is appended to the
// provided slice, and the resulting slice is returned. The extension is
// done in place if the provided slice has enough capacity.
func (d *GF255e) Encode(dst []byte) []byte {
	return gf_encode(dst, (*[4]uint64)(d), mq255e)
}

// Decode element from 32 bytes. If the source is invalid (out of range),
// then the decoded value is zero, and 0 is returned; otherwise, 1 is
// returned.
func (d *GF255e) Decode(src []byte) uint64 {
	return gf_decode((*[4]uint64)(d), src, mq255e)
}

// Decode element from bytes. The input bytes are interpreted as an
// integer (unsigned little-endian convention) which is reduced modulo
// the field modulus. By definition, this process cannot fail.
func (d *GF255e) DecodeReduce(src []byte) *GF255e {
	gf_decodeReduce((*[4]uint64)(d), src, mq255e)
	return d
}
