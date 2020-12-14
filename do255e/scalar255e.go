package do255e

import (
	"github.com/doubleodd/go-do255/internal/scalar"
	"math/bits"
)

// This file defines types and functions for do255e scalars, i.e.
// integers modulo r = 2^254 - 131528281291764213006042413802501683931
//
// Unless explicitly documented, all functions here are constant-time.

// Do255eScalar is the type for an integer modulo the prime order of the
// do255e group. Default value is zero.
type Do255eScalar [4]uint64

// Decode a scalar from exactly 32 bytes. Returned value is:
//   1   scalar properly decoded, value is not zero
//   0   scalar properly decoded, value is zero
//  -1   source bytes were not a valid scalar encoding
// If the decoding fails, then the scalar value is forced to zero.
func (s *Do255eScalar) Decode(src []byte) int {
	return scalar.Decode((*[4]uint64)(s), src, &do255eOrder)
}

// Decode a scalar from some bytes. All provided bytes are read and
// interpreted as an integer in unsigned little endian convention, which
// is reduced modulo the curve subgroup order. This process cannot fail.
func (s *Do255eScalar) DecodeReduce(src []byte) {
	scalar.DecodeReduce((*[4]uint64)(s), src, do255eModrReduce384Partial)
}

// Encode a scalar into exactly 32 bytes. The bytes are appended to the
// provided slice; the new slice is returned. The extension is done in
// place if the provided slice has enough capacity.
func (s *Do255eScalar) Encode(dst []byte) []byte {
	return scalar.Encode(dst, (*[4]uint64)(s), do255eModrReduce256)
}

// Encode a scalar into exactly 32 bytes.
func (s *Do255eScalar) Bytes() [32]byte {
	return scalar.ToBytes((*[4]uint64)(s), do255eModrReduce256)
}

// Compare a scalar with zero. Returned value is 1 if the scalar is zero,
// 0 otherwise.
func (s *Do255eScalar) IsZero() int {
	var t [4]uint64
	do255eModrReduce256(&t, (*[4]uint64)(s))
	z := t[0] | t[1] | t[2] | t[3]
	return int(1 - ((z | -z) >> 63))
}

// Compare two scalars together. Returned value is 1 if the scalars are
// equal to each other, 0 otherwise.
func (s *Do255eScalar) Equal(a *Do255eScalar) int {
	var t Do255eScalar
	t.Sub(s, a)
	return t.IsZero()
}

// Scalar addition: s is set to a + b (mod r).
// A pointer to s is returned.
func (s *Do255eScalar) Add(a, b *Do255eScalar) *Do255eScalar {
	scalar.Add((*[4]uint64)(s), (*[4]uint64)(a), (*[4]uint64)(b), do255eModrReduce256Partial)
	return s
}

// Scalar subtraction: s is set to a - b (mod r).
// A pointer to s is returned.
func (s *Do255eScalar) Sub(a, b *Do255eScalar) *Do255eScalar {
	scalar.Sub((*[4]uint64)(s), (*[4]uint64)(a), (*[4]uint64)(b), do255eModrReduce256Partial, &do255eOrder)
	return s
}

// Scalar negation: s is set to -a (mod r).
// A pointer to s is returned.
func (s *Do255eScalar) Neg(a *Do255eScalar) *Do255eScalar {
	var t = [4]uint64{0, 0, 0, 0}
	scalar.Sub((*[4]uint64)(s), &t, (*[4]uint64)(a), do255eModrReduce256Partial, &do255eOrder)
	return s
}

// Scalar multiplication: s is set to a*b (mod r).
// A pointer to s is returned.
func (s *Do255eScalar) Mul(a, b *Do255eScalar) *Do255eScalar {
	scalar.Mul((*[4]uint64)(s), (*[4]uint64)(a), (*[4]uint64)(b), do255eModrReduce384Partial)
	return s
}

// Group order is r = 2^254 - r0, with:
//    r0 = 131528281291764213006042413802501683931
const do255e_r0_lo uint64 = 0xE0AD37518B27BADB
const do255e_r0_hi uint64 = 0x62F36CF0ABF873AC

var do255eOrder = [4]uint64{
	0x1F52C8AE74D84525,
	0x9D0C930F54078C53,
	0xFFFFFFFFFFFFFFFF,
	0x3FFFFFFFFFFFFFFF,
}

// Given input 'a' (up to 2^258-1), perform a partial reduction modulo r;
// output (into 'd') fits on 255 bits and is (much) lower than 2*r. The
// high bits of 'a' are provided as extra parameter ah.
func do255eModrReduce256PartialWithExtra(d, a *[4]uint64, ah uint64) {
	// Truncate to 254 bits and get extra bits into ah.
	ah = (ah << 2) | (a[3] >> 62)

	// Compute ah*r0 into u0:u1:u2.
	u1, u0 := bits.Mul64(ah, do255e_r0_lo)
	u2, lo := bits.Mul64(ah, do255e_r0_hi)
	var cc uint64
	u1, cc = bits.Add64(u1, lo, 0)
	u2 += cc

	// 2^254 = r0 mod r
	d[0], cc = bits.Add64(a[0], u0, 0)
	d[1], cc = bits.Add64(a[1], u1, cc)
	d[2], cc = bits.Add64(a[2], u2, cc)
	d[3] = (a[3] & 0x3FFFFFFFFFFFFFFF) + cc
}

// Partial reduction ensures that the output fits on 255 bits and is
// less than 2*r.
func do255eModrReduce256Partial(d, a *[4]uint64) {
	do255eModrReduce256PartialWithExtra(d, a, 0)
}

// Given a partially reduced input 'a' (less than 2*r), finish reduction
// (conditional subtraction of r).
func do255eModrReduce256Finish(d, a *[4]uint64) {
	// Subtracting r is equivalent to adding r0, and subtracting
	// 2^254.
	var t [4]uint64
	var cc uint64
	t[0], cc = bits.Add64(a[0], do255e_r0_lo, 0)
	t[1], cc = bits.Add64(a[1], do255e_r0_hi, cc)
	t[2], cc = bits.Add64(a[2], 0, cc)
	t[3], cc = bits.Add64(a[3], 0, cc)
	t[3] -= 0x4000000000000000

	// Since the result fits on 255 bits, the top bit is a sign bit,
	// which we use to decide whether we use t[] or a[] as result.
	m := -(t[3] >> 63)
	for i := 0; i < 4; i++ {
		d[i] = t[i] ^ (m & (a[i] ^ t[i]))
	}
}

// Perform full reduction of a scalar.
func do255eModrReduce256(d, a *[4]uint64) {
	do255eModrReduce256Partial(d, a)
	do255eModrReduce256Finish(d, d)
}

// Given a 384-bit input 'a', perform a partial reduction modulo r;
// output fits on 255 bits and is less than 2*r.
func do255eModrReduce384Partial(d *[4]uint64, a *[6]uint64) {
	// Multiply the high third (a4:a5) by r0 into tw.
	var t1, t2 [2]uint64
	var tw [4]uint64
	t1[0] = do255e_r0_lo
	t1[1] = do255e_r0_hi
	t2[0] = a[4]
	t2[1] = a[5]
	scalar.Mul128x128(&tw, &t1, &t2)

	// Compute 4*tw and add to the low part of 'a'.
	// Since 4*r0 =~ 2^128.63, the result fits on 258 bits.
	var th uint64
	th = tw[3] >> 62
	tw[3] = (tw[3] << 2) | (tw[2] >> 62)
	tw[2] = (tw[2] << 2) | (tw[1] >> 62)
	tw[1] = (tw[1] << 2) | (tw[0] >> 62)
	tw[0] = tw[0] << 2
	var cc uint64
	tw[0], cc = bits.Add64(tw[0], a[0], 0)
	tw[1], cc = bits.Add64(tw[1], a[1], cc)
	tw[2], cc = bits.Add64(tw[2], a[2], cc)
	tw[3], cc = bits.Add64(tw[3], a[3], cc)
	th += cc

	// Perform partial reduction.
	do255eModrReduce256PartialWithExtra(d, &tw, th)
}

// Given an input k (fully reduced, in 0..r-1) and e < 2^127-2,
// compute d = round(k*e / r)
func do255eMulDivrRounded(d *[2]uint64, k *[4]uint64, e *[2]uint64) {
	// z <- k*e
	var z [6]uint64
	scalar.Mul256x128(&z, k, e)

	// z <- z + (r-1)/2
	var cc uint64
	z[0], cc = bits.Add64(z[0], 0x8FA964573A6C2292, 0)
	z[1], cc = bits.Add64(z[1], 0xCE864987AA03C629, cc)
	z[2], cc = bits.Add64(z[2], 0xFFFFFFFFFFFFFFFF, cc)
	z[3], cc = bits.Add64(z[3], 0x1FFFFFFFFFFFFFFF, cc)
	z[4], cc = bits.Add64(z[4], 0, cc)
	z[5] += cc

	// y <- floor(z / 2^254) + 1
	var y [2]uint64
	y[0] = (z[3] >> 62) | (z[4] << 2)
	y[1] = (z[4] >> 62) | (z[5] << 2)
	y[0], cc = bits.Add64(y[0], 1, 0)
	y[1] += cc

	// t <- y*r0
	var r0 [2]uint64
	r0[0] = do255e_r0_lo
	r0[1] = do255e_r0_hi
	var t [4]uint64
	scalar.Mul128x128(&t, &y, &r0)

	// t <- t + z0
	// We are only interested in the high limb.
	z[3] &= 0x3FFFFFFFFFFFFFFF
	_, cc = bits.Add64(z[0], t[0], 0)
	_, cc = bits.Add64(z[1], t[1], cc)
	_, cc = bits.Add64(z[2], t[2], cc)
	th, _ := bits.Add64(z[3], t[3], cc)

	// The high limb is in th and it is lower than 2^63. If it
	// is lower than 2^62, then y is too large and we must
	// decrement it; otherwise, we keep it unchanged.
	d[0], cc = bits.Sub64(y[0], 1-(th>>62), 0)
	d[1] = y[1] - cc
}

// Constants for scalar splitting.
var do255e_u = [2]uint64{0x2ACCF9DEC93F6111, 0x1A509F7A53C2C6E6}
var do255e_v = [2]uint64{0x0B7A31305466F77E, 0x7D440C6AFFBB3A93}

// Split scalar k (256 bits) into k0 and k1 (128 bits each, signed),
// such that k = k0 + k1*mu mod r, where mu is a square root of -1
// modulo r.
func (k *Do255eScalar) SplitMu(k0, k1 *[2]uint64) {
	// Ensure that k is fully reduced modulo r.
	var t [4]uint64
	do255eModrReduce256(&t, (*[4]uint64)(k))

	// c = round(k*v / r)
	// d = round(k*u / r)
	var c, d [2]uint64
	do255eMulDivrRounded(&c, &t, &do255e_v)
	do255eMulDivrRounded(&d, &t, &do255e_u)

	// k0 = k - d*u - c*v
	var y [2]uint64
	var cc uint64
	scalar.Mul128x128trunc(&y, &d, &do255e_u)
	t[0], cc = bits.Sub64(t[0], y[0], 0)
	t[1], _ = bits.Sub64(t[1], y[1], cc)
	scalar.Mul128x128trunc(&y, &c, &do255e_v)
	k0[0], cc = bits.Sub64(t[0], y[0], 0)
	k0[1], _ = bits.Sub64(t[1], y[1], cc)

	// k1 = d*v - c*u
	scalar.Mul128x128trunc(k1, &d, &do255e_v)
	scalar.Mul128x128trunc(&y, &c, &do255e_u)
	k1[0], cc = bits.Sub64(k1[0], y[0], 0)
	k1[1], _ = bits.Sub64(k1[1], y[1], cc)
}

// Recode a scalar with 5-bit Booth encoding. Output is a sequence of
// small integers in the -15..+16 range such that:
//   a = \sum_{i=0}^{51} d[i]*2^(5*i)
// Top digit d[51] is nonnegative. If the input value is less than 2^255,
// then the top digit can only be 0 or 1.
// Each output digit is encoded in a byte as sign+mantissa: the low 5 bits
// of the byte are the absolute value of the digit (in the 0..16 range),
// and the high bit of the byte is set to 1 for a negative digit, 0 otherwise.
// When the digit is 0, the function may encode it as -0 (0x80) or +0 (0x00)
// (the top digit d[51] cannot be -0, only +0).
func (a *Do255eScalar) recode5(d *[52]byte) {
	scalar.Recode5(d, (*[4]uint64)(a))
}
