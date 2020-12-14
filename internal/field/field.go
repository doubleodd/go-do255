package field

import (
	"encoding/binary"
	"math/bits"
)

// This file implements computations on some finite fields of integers
// modulo 2^255 - mq (for some small integers mq). This implementation
// is portable (no assembly) but should be decently efficient on 64-bit
// architectures. It is safe (constant-time) as long as 64-bit operations
// (especially 64x64->128 multiplication, using math/bits.Mul64()) are
// constant-time, which should be true on most modern systems.

// =======================================================================
// Internal functions
// =======================================================================

// Unless otherwise stated, all functions below accept source and destination
// operands to be the same objects. Parameter order is destination first
// (similar to mathematical notation: "d = a + b").
// The 'mq' parameter is the small integer such that modulus is p = 2^255 - mq.
// For all fields supported by this module, mq < 2^15 = 32767.
//
// Storage format: an array of four 64-bit unsigned integers, which encode
// the value in base 2^64 (little-endian order: first limb is least
// significant). Values are not necessarily reduced on output; all functions
// accept inputs in the whole 0..2^256-1 range.

// Internal function for field addition.
// Parameters:
//   d    destination
//   a    first operand
//   b    second operand
//   mq   modulus definition parameter
func gf_add(d, a, b *[4]uint64, mq uint64) {
	// First pass: sum over 256 bits + carry
	var cc uint64 = 0
	for i := 0; i < 4; i++ {
		d[i], cc = bits.Add64(a[i], b[i], cc)
	}

	// Second pass: if there is a carry, subtract 2*p = 2^256 - 2*mq;
	// i.e. we add 2*mq.
	d[0], cc = bits.Add64(d[0], (mq<<1)&-cc, 0)
	for i := 1; i < 4; i++ {
		d[i], cc = bits.Add64(d[i], 0, cc)
	}

	// If there is an extra carry, then this means that the initial
	// sum was at least 2^257 - 2*mq, in which case the current low
	// limb is necessarily lower than 2*mq, and adding 2*mq again
	// won't trigger an extra carry.
	d[0] += (mq << 1) & -cc
}

// Internal function for field subtraction.
// Parameters:
//   d    destination
//   a    first operand
//   b    second operand
//   mq   modulus definition parameter
func gf_sub(d, a, b *[4]uint64, mq uint64) {
	// First pass: difference over 256 bits + borrow
	var cc uint64 = 0
	for i := 0; i < 4; i++ {
		d[i], cc = bits.Sub64(a[i], b[i], cc)
	}

	// Second pass: if there is a borrow, add 2*p = 2^256 - 2*mq;
	// i.e. we subtract 2*mq.
	d[0], cc = bits.Sub64(d[0], (mq<<1)&-cc, 0)
	for i := 1; i < 4; i++ {
		d[i], cc = bits.Sub64(d[i], 0, cc)
	}

	// If there is an extra borrow, then this means that the
	// subtraction of 2*mq above triggered a borrow, and the first
	// limb is at least 2^64 - 2*mq; we can subtract 2*mq again without
	// triggering another borrow.
	d[0] -= (mq << 1) & -cc
}

// Internal function for field negation.
// Parameters:
//   d    destination
//   a    operand
//   mq   modulus definition parameter
func gf_neg(d, a *[4]uint64, mq uint64) {
	// First pass: compute 2*p - a over 256 bits.
	var cc uint64
	d[0], cc = bits.Sub64(-(mq << 1), a[0], 0)
	for i := 1; i < 4; i++ {
		d[i], cc = bits.Sub64(0xFFFFFFFFFFFFFFFF, a[i], cc)
	}

	// Second pass: if there is a borrow, add back p = 2^255 - mq.
	var e uint64 = -cc
	d[0], cc = bits.Add64(d[0], e&-mq, 0)
	for i := 1; i < 3; i++ {
		d[i], cc = bits.Add64(d[i], e, cc)
	}
	d[3], _ = bits.Add64(d[3], e>>1, cc)
}

// Internal function for constant-time selection. Output d is set to
// the value of a if ctl == 1, or to the value of b if ctl == 0.
// ctl MUST be 0 or 1.
// Parameters:
//   d     destination
//   a     first source
//   b     second source
//   ctl   1 to use the first source, 0 for the second source
// ctl MUST be 0 or 1
func gf_select(d, a, b *[4]uint64, ctl uint64) {
	ma := -ctl
	mb := ^ma
	for i := 0; i < 4; i++ {
		d[i] = (a[i] & ma) | (b[i] & mb)
	}
}

// Conditional negation: if ctl == 1, then d is set to -a; otherwise,
// if ctl == 0, then d is set to a. ctl MUST be 0 or 1.
//   d     destination
//   a     operand
//   mq    modulus definition parameter
//   ctl   control parameter
func gf_condneg(d, a *[4]uint64, mq uint64, ctl uint64) {
	var t [4]uint64
	gf_neg(&t, a, mq)
	gf_select(d, &t, a, ctl)
}

// Internal function for multiplication.
// Parameters:
//   d    destination
//   a    first operand
//   b    second operand
//   mq   modulus definition parameter
func gf_mul(d, a, b *[4]uint64, mq uint64) {
	var t [8]uint64
	var hi, lo, cc uint64

	// Step 1: multiply the two operands as plain integers, 512-bit
	// result goes to t[]. We have 16 products a[i]*b[j] to compute
	// and add at the right place; sequence below tries to do them
	// in an order that minimizes carry propagation steps.

	// a0*b0, a1*b1, a2*b2, a3*b3
	t[1], t[0] = bits.Mul64(a[0], b[0])
	t[3], t[2] = bits.Mul64(a[1], b[1])
	t[5], t[4] = bits.Mul64(a[2], b[2])
	t[7], t[6] = bits.Mul64(a[3], b[3])

	// a0*b1, a0*b3, a2*b3
	hi, lo = bits.Mul64(a[0], b[1])
	t[1], cc = bits.Add64(t[1], lo, 0)
	t[2], cc = bits.Add64(t[2], hi, cc)
	hi, lo = bits.Mul64(a[0], b[3])
	t[3], cc = bits.Add64(t[3], lo, cc)
	t[4], cc = bits.Add64(t[4], hi, cc)
	hi, lo = bits.Mul64(a[2], b[3])
	t[5], cc = bits.Add64(t[5], lo, cc)
	t[6], cc = bits.Add64(t[6], hi, cc)
	t[7] += cc

	// a1*b0, a3*b0, a3*b2
	hi, lo = bits.Mul64(a[1], b[0])
	t[1], cc = bits.Add64(t[1], lo, 0)
	t[2], cc = bits.Add64(t[2], hi, cc)
	hi, lo = bits.Mul64(a[3], b[0])
	t[3], cc = bits.Add64(t[3], lo, cc)
	t[4], cc = bits.Add64(t[4], hi, cc)
	hi, lo = bits.Mul64(a[3], b[2])
	t[5], cc = bits.Add64(t[5], lo, cc)
	t[6], cc = bits.Add64(t[6], hi, cc)
	t[7] += cc

	// a0*b2, a1*b3
	hi, lo = bits.Mul64(a[0], b[2])
	t[2], cc = bits.Add64(t[2], lo, 0)
	t[3], cc = bits.Add64(t[3], hi, cc)
	hi, lo = bits.Mul64(a[1], b[3])
	t[4], cc = bits.Add64(t[4], lo, cc)
	t[5], cc = bits.Add64(t[5], hi, cc)
	t[6], cc = bits.Add64(t[6], 0, cc)
	t[7] += cc

	// a2*b0, a3*b1
	hi, lo = bits.Mul64(a[2], b[0])
	t[2], cc = bits.Add64(t[2], lo, 0)
	t[3], cc = bits.Add64(t[3], hi, cc)
	hi, lo = bits.Mul64(a[3], b[1])
	t[4], cc = bits.Add64(t[4], lo, cc)
	t[5], cc = bits.Add64(t[5], hi, cc)
	t[6], cc = bits.Add64(t[6], 0, cc)
	t[7] += cc

	// a1*b2, a2*b1
	var x0, x1, x2 uint64
	x1, x0 = bits.Mul64(a[1], b[2])
	hi, lo = bits.Mul64(a[2], b[1])
	x0, cc = bits.Add64(x0, lo, 0)
	x1, x2 = bits.Add64(x1, hi, cc)
	t[3], cc = bits.Add64(t[3], x0, 0)
	t[4], cc = bits.Add64(t[4], x1, cc)
	t[5], cc = bits.Add64(t[5], x2, cc)
	t[6], cc = bits.Add64(t[6], 0, cc)
	t[7] += cc

	// Step 2: fold upper half into lower half, multiplied by 2*mq.
	// Each high word (t[4..7]) is multipied by 2*mq, yielding a
	// low half (64 bits, added into the low words t[0..3]) and a
	// high half (h0..h3, value at most 2*mq-1 < 2^16).

	var h0, h1, h2, h3 uint64
	h0, lo = bits.Mul64(t[4], mq<<1)
	t[0], cc = bits.Add64(t[0], lo, 0)
	h1, lo = bits.Mul64(t[5], mq<<1)
	t[1], cc = bits.Add64(t[1], lo, cc)
	h2, lo = bits.Mul64(t[6], mq<<1)
	t[2], cc = bits.Add64(t[2], lo, cc)
	h3, lo = bits.Mul64(t[7], mq<<1)
	t[3], cc = bits.Add64(t[3], lo, cc)
	h3 += cc

	// We must still add the upper words h0..h3 into the result, at
	// their proper place. h3 is to be folded again; we also include
	// bit 255 into h3 so that this step triggers no further carry.
	// Note that (2*h3+1)*mq <= 2*mq^2 < 2^31, hence we can do that
	// multiplication with the basic operator instead of Mul64().
	// Since this step produces the final output words, we can write
	// them into the destination directly.

	h3 = (h3 << 1) | (t[3] >> 63)
	t[3] &= 0x7FFFFFFFFFFFFFFF
	d[0], cc = bits.Add64(t[0], h3*mq, 0)
	d[1], cc = bits.Add64(t[1], h0, cc)
	d[2], cc = bits.Add64(t[2], h1, cc)
	d[3], cc = bits.Add64(t[3], h2, cc)
}

// Internal function for squaring.
// Parameters:
//   d    destination
//   a    operand
//   mq   modulus definition parameter
func gf_sqr(d, a *[4]uint64, mq uint64) {
	var t [8]uint64
	var hi, lo, cc uint64

	// Step 1: square the operand as a plain integer, 512-bit
	// result goes to t[]. Sequence below tries to do them
	// in an order that minimizes carry propagation steps.

	// First the non-square products:
	//   a0*a1, a0*a2, a0*a3, a1*a2, a1*a3, a2*a3
	// This partial sum is necessarily lower than 2^448, so there
	// is no carry to spill into t[7].
	t[2], t[1] = bits.Mul64(a[0], a[1])
	t[4], t[3] = bits.Mul64(a[0], a[3])
	t[6], t[5] = bits.Mul64(a[2], a[3])
	hi, lo = bits.Mul64(a[0], a[2])
	t[2], cc = bits.Add64(t[2], lo, 0)
	t[3], cc = bits.Add64(t[3], hi, cc)
	hi, lo = bits.Mul64(a[1], a[3])
	t[4], cc = bits.Add64(t[4], lo, cc)
	t[5], cc = bits.Add64(t[5], hi, cc)
	t[6] += cc
	hi, lo = bits.Mul64(a[1], a[2])
	t[3], cc = bits.Add64(t[3], lo, 0)
	t[4], cc = bits.Add64(t[4], hi, cc)
	t[5], cc = bits.Add64(t[5], 0, cc)
	t[6] += cc

	// Double the current sum.
	t[7] = t[6] >> 63
	t[6] = (t[6] << 1) | (t[5] >> 63)
	t[5] = (t[5] << 1) | (t[4] >> 63)
	t[4] = (t[4] << 1) | (t[3] >> 63)
	t[3] = (t[3] << 1) | (t[2] >> 63)
	t[2] = (t[2] << 1) | (t[1] >> 63)
	t[1] = t[1] << 1

	// Add the squares: a0*a0, a1*a1, a2*a2, a3*a3
	hi, t[0] = bits.Mul64(a[0], a[0])
	t[1], cc = bits.Add64(t[1], hi, 0)
	hi, lo = bits.Mul64(a[1], a[1])
	t[2], cc = bits.Add64(t[2], lo, cc)
	t[3], cc = bits.Add64(t[3], hi, cc)
	hi, lo = bits.Mul64(a[2], a[2])
	t[4], cc = bits.Add64(t[4], lo, cc)
	t[5], cc = bits.Add64(t[5], hi, cc)
	hi, lo = bits.Mul64(a[3], a[3])
	t[6], cc = bits.Add64(t[6], lo, cc)
	t[7], _ = bits.Add64(t[7], hi, cc)

	// Step 2: we now have the 512-bit result in t[0..7]. We apply
	// reduction modulo p. This is the same code as in gf_mul();
	// see the comments in that function.

	var h0, h1, h2, h3 uint64
	h0, lo = bits.Mul64(t[4], mq<<1)
	t[0], cc = bits.Add64(t[0], lo, 0)
	h1, lo = bits.Mul64(t[5], mq<<1)
	t[1], cc = bits.Add64(t[1], lo, cc)
	h2, lo = bits.Mul64(t[6], mq<<1)
	t[2], cc = bits.Add64(t[2], lo, cc)
	h3, lo = bits.Mul64(t[7], mq<<1)
	t[3], cc = bits.Add64(t[3], lo, cc)
	h3 += cc

	h3 = (h3 << 1) | (t[3] >> 63)
	t[3] &= 0x7FFFFFFFFFFFFFFF
	d[0], cc = bits.Add64(t[0], h3*mq, 0)
	d[1], cc = bits.Add64(t[1], h0, cc)
	d[2], cc = bits.Add64(t[2], h1, cc)
	d[3], cc = bits.Add64(t[3], h2, cc)
}

// Internal multiplication of multiple squarings: d = a^(2^n)
// Parameters:
//   d    destination
//   a    operand
//   n    number of squarings to perform
//   mq   modulus definition parameter
func gf_sqr_x(d, a *[4]uint64, n uint, mq uint64) {
	if n == 0 {
		copy(d[:], a[:])
		return
	}
	gf_sqr(d, a, mq)
	for n -= 1; n != 0; n-- {
		gf_sqr(d, d, mq)
	}
}

// Internal function for halving (division by 2).
// Parameters:
//   d    destination
//   a    operand
//   mq   modulus definition parameter
func gf_half(d, a *[4]uint64, mq uint64) {
	// We right shift, and add (p+1)/2 = 2^254 - ((mq-1)/2) conditionally
	// on the least significant bit of the source.
	var e uint64 = -(a[0] & 1)
	var cc uint64
	d[0], cc = bits.Add64((a[0]>>1)|(a[1]<<63), e&-((mq-1)>>1), 0)
	for i := 1; i < 3; i++ {
		d[i], cc = bits.Add64((a[i]>>1)|(a[i+1]<<63), e, cc)
	}
	d[3], _ = bits.Add64(a[3]>>1, e>>2, cc)
}

// Internal function for left-shifting by some bits.
// Parameters:
//   d    destination
//   a    operand
//   n    shift count (at least 1, at most 15).
//   mq   modulus definition parameter
func gf_lsh(d, a *[4]uint64, n uint, mq uint64) {
	// First pass: left shift, extra bits in g.
	var g uint64 = a[0] >> (64 - n)
	d[0] = a[0] << n
	for i := 1; i < 4; i++ {
		w := a[i]
		d[i] = (w << n) | g
		g = w >> (64 - n)
	}

	// Second pass: reduction of extra bits (with the top bit of the
	// value).
	g = (g << 1) | (d[3] >> 63)
	var cc uint64
	d[0], cc = bits.Add64(d[0], g*mq, 0)
	for i := 1; i < 3; i++ {
		d[i], cc = bits.Add64(d[i], 0, cc)
	}
	d[3] = (d[3] & 0x7FFFFFFFFFFFFFFF) + cc
}

// Internal function for normalization. This function ensures that the
// output is in the 0..p-1 range. It is meant to be called prior to
// encoding, or for comparisons.
//   d    destination
//   a    operand
//   mq   modulus definition parameter
func gf_norm(d, a *[4]uint64, mq uint64) {
	// Fold the top bit to ensure a value of at most 2^255 + mq-1.
	var cc uint64
	d[0], cc = bits.Add64(a[0], mq&-(a[3]>>63), 0)
	for i := 1; i < 3; i++ {
		d[i], cc = bits.Add64(a[i], 0, cc)
	}
	d[3] = (a[3] & 0x7FFFFFFFFFFFFFFF) + cc

	// Subtract p.
	d[0], cc = bits.Sub64(d[0], -mq, 0)
	for i := 1; i < 3; i++ {
		d[i], cc = bits.Sub64(d[i], 0xFFFFFFFFFFFFFFFF, cc)
	}
	d[3], cc = bits.Sub64(d[3], 0x7FFFFFFFFFFFFFFF, cc)

	// If there is a borrow, add p back.
	var e uint64 = -cc
	d[0], cc = bits.Add64(d[0], e&-mq, 0)
	for i := 1; i < 3; i++ {
		d[i], cc = bits.Add64(d[i], e, cc)
	}
	d[3], cc = bits.Add64(d[3], e>>1, cc)
}

// Internal function for comparing a value with zero. This function
// returns 1 if the value is equal to 0 modulo p; otherwise, it returns 0.
//   a    operand
//   mq   modulus definition parameter
func gf_iszero(a *[4]uint64, mq uint64) uint64 {
	// There are three possible representations for zero: 0, p and 2*p.
	t0 := a[0]
	t1 := a[0] + mq
	t2 := a[0] + (mq << 1)
	for i := 1; i < 3; i++ {
		t0 |= a[i]
		t1 |= ^a[i]
		t2 |= ^a[i]
	}
	t0 |= a[3]
	t1 |= a[3] ^ 0x7FFFFFFFFFFFFFFF
	t2 |= ^a[3]
	return 1 - (((t0 | -t0) & (t1 | -t1) & (t2 | -t2)) >> 63)
}

// Internal function for comparing two values. This function returns 1
// the values are equal modulo p, 0 otherwise.
//   a    first operand
//   b    second operand
//   mq   modulus definition parameter
func gf_eq(a, b *[4]uint64, mq uint64) uint64 {
	var t [4]uint64
	gf_sub(&t, a, b, mq)
	return gf_iszero(&t, mq)
}

// Internal function for encoding a field element into 32 bytes. The
// encoded element is appended to the specified slice; the new slice
// (with the appended data) is returned.
func gf_encode(b []byte, a *[4]uint64, mq uint64) []byte {
	len1 := len(b)
	len2 := len1 + 32
	var b2 []byte
	if cap(b) >= len2 {
		b2 = b[:len2]
	} else {
		b2 = make([]byte, len2)
		copy(b2, b)
	}
	dst := b2[len1:]
	var t [4]uint64
	gf_norm(&t, a, mq)
	for i := 0; i < 4; i++ {
		binary.LittleEndian.PutUint64(dst[8*i:], t[i])
	}
	return b2
}

// Internal function for decoding a field element from 32 bytes. If the
// source is not in the valid range (0..p-1), then the destination is
// set to all zeros, and 0 is returned; otherwise, 1 is returned.
func gf_decode(d *[4]uint64, src []byte, mq uint64) uint64 {
	for i := 0; i < 4; i++ {
		d[i] = binary.LittleEndian.Uint64(src[8*i:])
	}
	// Compare with the modulus. If there is a borrow (cc == 1),
	// then the value is correct; otherwise (cc == 0) it is out of
	// range and shall be cleared.
	_, cc := bits.Sub64(d[0], -mq, 0)
	_, cc = bits.Sub64(d[1], 0xFFFFFFFFFFFFFFFF, cc)
	_, cc = bits.Sub64(d[2], 0xFFFFFFFFFFFFFFFF, cc)
	_, cc = bits.Sub64(d[3], 0x7FFFFFFFFFFFFFFF, cc)
	for i := 0; i < 4; i++ {
		d[i] &= -cc
	}
	return cc
}

// Internal function for decoding a field element from bytes, with
// reduction. An arbitrary number of input bytes can be used. This
// process cannot fail.
func gf_decodeReduce(d *[4]uint64, src []byte, mq uint64) {
	var t [8]uint64

	// Initialize the low half of t with the rightmost bytes; we use
	// j bytes such that len(src)-j is a multiple of 32.
	n := len(src)
	j := n & 31
	if j == 0 && n != 0 {
		j = 32
	}
	n -= j
	var buf [32]byte
	copy(buf[:], src[n:])
	for i := 0; i < 4; i++ {
		t[i] = binary.LittleEndian.Uint64(buf[8*i:])
	}

	// For all remaining chunks of 32 bytes (right-to-left order),
	// shift the current value, add the next chunk, and reduce.
	for n > 0 {
		n -= 32
		copy(t[4:], t[:4])
		for i := 0; i < 4; i++ {
			t[i] = binary.LittleEndian.Uint64(src[n+8*i:])
		}

		// Fold upper half into lower half, multiplied by 2*mq.
		// Each high word (t[4..7]) is multipied by 2*mq,
		// yielding a low half (64 bits, added into the low
		// words t[0..3]) and a high half (h0..h3, value at most
		// 2*mq-1 < 2^16).
		var h0, h1, h2, h3 uint64
		var lo, cc uint64
		h0, lo = bits.Mul64(t[4], mq<<1)
		t[0], cc = bits.Add64(t[0], lo, 0)
		h1, lo = bits.Mul64(t[5], mq<<1)
		t[1], cc = bits.Add64(t[1], lo, cc)
		h2, lo = bits.Mul64(t[6], mq<<1)
		t[2], cc = bits.Add64(t[2], lo, cc)
		h3, lo = bits.Mul64(t[7], mq<<1)
		t[3], cc = bits.Add64(t[3], lo, cc)
		h3 += cc

		// We must still add the upper words h0..h3 into the
		// result, at their proper place. h3 is to be folded
		// again; we also include bit 255 into h3 so that this
		// step triggers no further carry. Note that
		// (2*h3+1)*mq <= 2*mq^2 < 2^31, hence we can do that
		// multiplication with the basic operator instead of
		// Mul64(). Since this step produces the final output
		// words, we can write them into the destination
		// directly.
		h3 = (h3 << 1) | (t[3] >> 63)
		t[3] &= 0x7FFFFFFFFFFFFFFF
		t[0], cc = bits.Add64(t[0], h3*mq, 0)
		t[1], cc = bits.Add64(t[1], h0, cc)
		t[2], cc = bits.Add64(t[2], h1, cc)
		t[3], cc = bits.Add64(t[3], h2, cc)
	}

	// Copy the result.
	copy(d[:], t[:4])
}

// =======================================================================
// Internal support code for inversion and Legendre symbol.
//
// Both inversion and Legendre symbol could be more easily implemented
// with exponentiations:
//    1/y = y^(p-2) mod p
//    Legendre(y) = y^((p-1)/2) mod p
// However, fully optimized implementations will prefer to use the
// algorithms employed below, since they are faster (even for 64-bit
// architectures with efficient 64x64->128 multiplications). The
// binary GCD algorithm is described here:
//    https://eprint.iacr.org/2020/972
// The adaptation to Legendre symbol is straightforward, and has been
// described here:
//    https://research.nccgroup.com/2020/09/28/faster-modular-inversion-and-legendre-symbol-and-an-x25519-speed-record/

// Count leading zeros in a 64-bit value.
// Output is in the 0..64 range.
// We do not use bits.LeadingZeros64() because the default implementation
// is not constant-time (on some architectures, there is a constant-time
// efficient opcode for that, e.g. LZCNT on recent x86, but there still are
// many systems without such facilities).
func countLeadingZeros(x uint64) uint64 {
	var r, c uint64
	r = 0
	c = -(((x >> 32) - 1) >> 63)
	r += c & 32
	x ^= c & (x ^ (x << 32))
	c = -(((x >> 48) - 1) >> 63)
	r += c & 16
	x ^= c & (x ^ (x << 16))
	c = -(((x >> 56) - 1) >> 63)
	r += c & 8
	x ^= c & (x ^ (x << 8))
	c = -(((x >> 60) - 1) >> 63)
	r += c & 4
	x ^= c & (x ^ (x << 4))
	c = -(((x >> 62) - 1) >> 63)
	r += c & 2
	x ^= c & (x ^ (x << 2))
	c = -(((x >> 63) - 1) >> 63)
	r += c & 1
	x ^= c & (x ^ (x << 1))
	r += 1 - ((x | -x) >> 63)
	return r
}

// Compute d <- (a*f+b*g)/2^31. Factors f and g are provided as uint64, but
// they really are signed integers in the -2^31..+2^31 range. Values
// a and b are plain 255-bit nonnegative integers. The division is
// assumed to be exact, and the signed result is assumed to fit in
// 256 bits (with its sign bit). If the result is negative, then it is
// negated. Returned value is 1 if that final negation had to be applied,
// 0 otherwise.
func gf_lin_div31_abs(d, a, b *[4]uint64, f, g uint64) uint64 {
	// If f < 0, replace f with -f, but keep the sign in sf.
	// Same treatment for g.
	sf := f >> 63
	f = (f ^ -sf) + sf
	sg := g >> 63
	g = (g ^ -sg) + sg

	// Apply signs sf and sg to a and b, respectively.
	var ta, tb [4]uint64
	var cc uint64
	ta[0], cc = bits.Add64(a[0]^-sf, sf, 0)
	for i := 1; i < 4; i++ {
		ta[i], cc = bits.Add64(a[i]^-sf, 0, cc)
	}
	tb[0], cc = bits.Add64(b[0]^-sg, sg, 0)
	for i := 1; i < 4; i++ {
		tb[i], cc = bits.Add64(b[i]^-sg, 0, cc)
	}

	// Compute a*f+b*g into d, with extra word in t.
	// Note: f and g are at most 2^31 here.
	z1, z0 := bits.Mul64(ta[0], f)
	hi, lo := bits.Mul64(tb[0], g)
	d[0], cc = bits.Add64(z0, lo, 0)
	t, _ := bits.Add64(z1, hi, cc)
	for i := 1; i < 4; i++ {
		z1, z0 = bits.Mul64(ta[i], f)
		hi, lo = bits.Mul64(tb[i], g)
		z0, cc = bits.Add64(z0, lo, 0)
		z1, _ = bits.Add64(z1, hi, cc)
		d[i], cc = bits.Add64(z0, t, 0)
		t = z1 + cc
	}

	// If a < 0, then the result is overestimated by 2^256*f; similarly
	// for the case b < 0. We adjust the result here.
	t -= -(ta[3] >> 63) & f
	t -= -(tb[3] >> 63) & g

	// Do the division by 2^31 (right-shift, since the division is
	// assumed to be exact).
	for i := 0; i < 3; i++ {
		d[i] = (d[i] >> 31) | (d[i+1] << 33)
	}
	d[3] = (d[3] >> 31) | (t << 33)

	// If the result is negative, negate it.
	t >>= 63
	d[0], cc = bits.Add64(d[0]^-t, t, 0)
	for i := 1; i < 4; i++ {
		d[i], cc = bits.Add64(d[i]^-t, 0, cc)
	}
	return t
}

// Compute u*f+v*g (mod p). Parameters f and g are provided with an
// unsigned type, but they really are signed integers in the
// -2^62..+2^62 range.
func gf_lin(d, u, v *[4]uint64, f, g uint64, mq uint64) {
	// If f < 0, replace f with -f, but keep the sign in sf.
	// Same treatment for g.
	sf := f >> 63
	f = (f ^ -sf) + sf
	sg := g >> 63
	g = (g ^ -sg) + sg

	// Apply signs sf and sg to u and v (in the field).
	var tu, tv [4]uint64
	gf_condneg(&tu, u, mq, sf)
	gf_condneg(&tv, v, mq, sg)

	// Compute u*f+v*g into d, with extra word in t.
	// Since |f| <= 2^62 and |g| <= 2^62, the 64-bit addition cannot
	// overflow.
	z1, z0 := bits.Mul64(tu[0], f)
	hi, lo := bits.Mul64(tv[0], g)
	var cc uint64
	d[0], cc = bits.Add64(z0, lo, 0)
	t, _ := bits.Add64(z1, hi, cc)
	for i := 1; i < 4; i++ {
		z1, z0 = bits.Mul64(tu[i], f)
		hi, lo = bits.Mul64(tv[i], g)
		z0, cc = bits.Add64(z0, lo, 0)
		z1, _ = bits.Add64(z1, hi, cc)
		d[i], cc = bits.Add64(z0, t, 0)
		t = z1 + cc
	}

	// Upper word can be up to 63 bits. We apply reduction.
	t = (t << 1) | (d[3] >> 63)
	d[3] &= 0x7FFFFFFFFFFFFFFF
	z1, z0 = bits.Mul64(t, mq)
	d[0], cc = bits.Add64(d[0], z0, 0)
	d[1], cc = bits.Add64(d[1], z1, cc)
	for i := 2; i < 4; i++ {
		d[i], cc = bits.Add64(d[i], 0, cc)
	}
}

// Internal function for inversion in the field: it computes 2^508/y.
// If the input is y = 0, then 0 is returned. This function is
// constant-time. This uses the optimized binary GCD from:
//    https://eprint.iacr.org/2020/972
// For a complete inversion, the caller must still divide the result
// by 2^508, which is normally done with a multiplication with the
// precomputed constant 1/2^508 mod p.
// Parameters:
//   d    destination
//   y    operand
//   mq   modulus definition parameter
func gf_inv_scaled(d, y *[4]uint64, mq uint64) {
	// Binary GCD starts with:
	//    a <- y
	//    b <- p
	//    u <- 1
	//    v <- 0
	// Then, at each step:
	//    if a is even, then:
	//        a <- a/2, u <- u/2 mod p
	//    else:
	//        if a < b:
	//            (a, u, b, v) <- (b, v, a, u)
	//        a <- (a-b)/2, u <- (u-v)/2 mod p
	// When a reaches 0, it stays there; at that point, b contains
	// the GCD of y and p (normally 1, since p is prime) and v
	// is the inverse of y modulo p.
	//
	// In the optimized version, we group iterations by chunks of 31.
	// In each chunk, we work over approximations of a and b, with
	// only the top 33 bits and bottom 31 bits. Updates to the actual
	// a, b, u and v are grouped: each chunk of 31 iterations computes
	// "update factors" which are then applied en masse.

	var a, b, u, v [4]uint64

	gf_norm(&a, y, mq)
	b[0] = -mq
	u[0] = 1
	v[0] = 0
	for i := 1; i < 3; i++ {
		b[i] = 0xFFFFFFFFFFFFFFFF
		u[i] = 0
		v[i] = 0
	}
	b[3] = 0x7FFFFFFFFFFFFFFF
	u[3] = 0
	v[3] = 0

	// First do 15*31 = 465 iterations.
	for i := 0; i < 15; i++ {
		// Extract approximations of a and b over 64 bits:
		//  - If len(a) <= 64 and len(b) <= 64, then we just use
		//    their values (low limb of each).
		//  - Otherwise, with n = max(len(a), len(b)), we use:
		//       (a mod 2^31) + 2^31*floor(a / 2^(n-33))
		//       (b mod 2^31) + 2^31*floor(b / 2^(n-33))
		// We first locate the top two words to use (i.e. we
		// skip limbs which are zero for both values).
		m3 := a[3] | b[3]
		m2 := a[2] | b[2]
		m1 := a[1] | b[1]
		tnz3 := -((m3 | -m3) >> 63)
		tnz2 := -((m2 | -m2) >> 63) & ^tnz3
		tnz1 := -((m1 | -m1) >> 63) & ^tnz3 & ^tnz2
		tnzm := (m3 & tnz3) | (m2 & tnz2) | (m1 & tnz1)
		tnza := (a[3] & tnz3) | (a[2] & tnz2) | (a[1] & tnz1)
		tnzb := (b[3] & tnz3) | (b[2] & tnz2) | (b[1] & tnz1)
		snza := (a[2] & tnz3) | (a[1] & tnz2) | (a[0] & tnz1)
		snzb := (b[2] & tnz3) | (b[1] & tnz2) | (b[0] & tnz1)

		// If len(a) <= 64 and len(b) <= 64, then all tnz* and
		// snz* are 0. Otherwise:
		//   tnzm != 0, length yields value of n
		//   tnza contains top limb of a, snza the second limb
		//   tnzb contains top limb of b, snzb the second limb
		//
		// Shifting is delicate here; on some architectures,
		// whether any shift count is greater than 31 or not may
		// leak through timing side-channels. Moreover, it is
		// possible that tnzm == 0. Thus:
		//   - we extract s = number of leading zeros in tnzm
		//     (s <= 0 <= 64);
		//   - if s <= 31, then we keep tnza unchanged; otherwise,
		//     we bring in 32 bits from snza; idem for tnzb/snzb;
		//   - we ensure s is in 0..31.
		// If tnzm == 0, then we end up with s == 0, which is not
		// a problem.
		s := countLeadingZeros(tnzm)
		sm := -(s >> 5)
		tnza ^= sm & (tnza ^ ((tnza << 32) | (snza >> 32)))
		tnzb ^= sm & (tnzb ^ ((tnzb << 32) | (snzb >> 32)))
		s &= 31
		tnza <<= s
		tnzb <<= s

		// At this point, if len(a) <= 64 and len(b) <= 64, then
		// all tnz* are zero. Otherwise, one of tnz1, tnz2 or tnz3
		// is 0xFFFFFFFFFFFFFFFF, and we need the top 33 bits of
		// tnza and tnzb.
		tnza |= a[0] & ^(tnz1 | tnz2 | tnz3)
		tnzb |= b[0] & ^(tnz1 | tnz2 | tnz3)
		xa := (a[0] & 0x000000007FFFFFFF) | (tnza & 0xFFFFFFFF80000000)
		xb := (b[0] & 0x000000007FFFFFFF) | (tnzb & 0xFFFFFFFF80000000)

		// We now have our approximations in xa and xb. We run
		// the chunk of 31 iterations, keeping track of updates
		// in the "update factors" fg0 and fg1.
		var fg0 uint64 = 1
		var fg1 uint64 = uint64(1) << 32
		for j := 0; j < 31; j++ {
			a_odd := -(xa & 1)
			_, cc := bits.Sub64(xa, xb, 0)
			swap := a_odd & -cc
			t := swap & (xa ^ xb)
			xa ^= t
			xb ^= t
			t = swap & (fg0 ^ fg1)
			fg0 ^= t
			fg1 ^= t
			xa -= a_odd & xb
			fg0 -= a_odd & fg1
			xa >>= 1
			fg1 <<= 1
		}

		// Extract individual update factors from the packed
		// representations fg0 and fg1.
		fg0 += 0x7FFFFFFF7FFFFFFF
		fg1 += 0x7FFFFFFF7FFFFFFF
		f0 := (fg0 & 0xFFFFFFFF) - 0x7FFFFFFF
		g0 := (fg0 >> 32) - 0x7FFFFFFF
		f1 := (fg1 & 0xFFFFFFFF) - 0x7FFFFFFF
		g1 := (fg1 >> 32) - 0x7FFFFFFF

		// Update a and b. Corresponding update factors are
		// conditionally negated if the update found a negative
		// output.
		var na, nb, nu, nv [4]uint64
		nega := gf_lin_div31_abs(&na, &a, &b, f0, g0)
		negb := gf_lin_div31_abs(&nb, &a, &b, f1, g1)
		f0 = (f0 ^ -nega) + nega
		g0 = (g0 ^ -nega) + nega
		f1 = (f1 ^ -negb) + negb
		g1 = (g1 ^ -negb) + negb
		gf_lin(&nu, &u, &v, f0, g0, mq)
		gf_lin(&nv, &u, &v, f1, g1, mq)
		a = na
		b = nb
		u = nu
		v = nv
	}

	// If y is invertible, then len(a) + len(b) <= 45 at this
	// point, so we can avoid the extraction, and 43 iterations are
	// sufficient to always arrive at b = 1. Since the values are
	// exact (no approximation), we can run all 43 iterations in one go.
	// However, we cannot pack update factors as we did previously,
	// since they will no longer fit in 32 bits each.
	xa := a[0]
	xb := b[0]
	var f0 uint64 = 1
	var g0 uint64 = 0
	var f1 uint64 = 0
	var g1 uint64 = 1
	for j := 0; j < 43; j++ {
		a_odd := -(xa & 1)
		_, cc := bits.Sub64(xa, xb, 0)
		swap := a_odd & -cc
		t := swap & (xa ^ xb)
		xa ^= t
		xb ^= t
		t = swap & (f0 ^ f1)
		f0 ^= t
		f1 ^= t
		t = swap & (g0 ^ g1)
		g0 ^= t
		g1 ^= t
		xa -= a_odd & xb
		f0 -= a_odd & f1
		g0 -= a_odd & g1
		xa >>= 1
		f1 <<= 1
		g1 <<= 1
	}
	gf_lin(&v, &u, &v, f1, g1, mq)

	// If the original value was zero, then v contains zero at this
	// point. Otherwise, it contains 2^508/y mod p. Either way, v
	// has the correct result.
	copy(d[:], v[:])
}

// Internal function for Legendre symbol. Return value is:
//    0  if y == 0
//    1  if y != 0 and is a quadratic residue
//   -1  if y != 0 and is not a quadratic residue
// Value is returned as an uint64, so -1 is 0xFFFFFFFFFFFFFFFF.
// Parameters:
//   y    operand
//   mq   modulus definition parameter
func gf_legendre(y *[4]uint64, mq uint64) uint64 {
	// This follows the same steps as the binary GCD used in
	// gf_inv_scaled(), with the following differences:
	//   - We do not keep track of the u and v values.
	//   - In the inner loop, we update the current symbol value
	//     as we apply the operations on a and b.
	//   - The last two iterations of the inner loop must access up
	//     to three bottom bits of b, so we compute the updated a and b
	//     (low bits only) at that point to get the correct bits.
	//
	// The running symbol value is held in the least significant bit
	// of 'ls' (other bits should be ignored). This is, in fact, the
	// Kronecker symbol (extension of the Jacobi symbol, which
	// itself extends the Legendre symbol). The algorithm relies on
	// the following well-known properties:
	//
	//   If x = y mod n, then (x|n) = (y|n), as long as either n > 0,
	//   or x and y have the same sign.
	//
	//   If x and y are not both negative, then (x|y) = (y|x), unless
	//   both x = 3 mod 4 and y = 3 mod 4, in which case (x|y) = -(y|x).
	//
	//   (2|n) = 1 if n = 1 or 7 mod 8, or -1 if n = 3 or 5 mod 8.
	//
	// We use these properties to keep track of symbol updates
	// through the binary GCD operations. The crucial observation is
	// that while it may happen that a or b becomes negative at some
	// point, it never happens that they are both negative.
	// Therefore, when replacing a with a-b, either b > 0, or a and
	// a-b have the same sign, so the symbol is conserved unchanged.
	// Also, when a and b are swapped, they are not both negative,
	// and we can update the symbol by looking at the two low bits
	// only.

	var a, b [4]uint64
	var ls uint64 = 0

	gf_norm(&a, y, mq)
	b[0] = -mq
	for i := 1; i < 3; i++ {
		b[i] = 0xFFFFFFFFFFFFFFFF
	}
	b[3] = 0x7FFFFFFFFFFFFFFF

	// First do 15*31 = 465 iterations.
	for i := 0; i < 15; i++ {
		// Extract approximations of a and b over 64 bits.
		m3 := a[3] | b[3]
		m2 := a[2] | b[2]
		m1 := a[1] | b[1]
		tnz3 := -((m3 | -m3) >> 63)
		tnz2 := -((m2 | -m2) >> 63) & ^tnz3
		tnz1 := -((m1 | -m1) >> 63) & ^tnz3 & ^tnz2
		tnzm := (m3 & tnz3) | (m2 & tnz2) | (m1 & tnz1)
		tnza := (a[3] & tnz3) | (a[2] & tnz2) | (a[1] & tnz1)
		tnzb := (b[3] & tnz3) | (b[2] & tnz2) | (b[1] & tnz1)
		snza := (a[2] & tnz3) | (a[1] & tnz2) | (a[0] & tnz1)
		snzb := (b[2] & tnz3) | (b[1] & tnz2) | (b[0] & tnz1)

		s := countLeadingZeros(tnzm)
		sm := -(s >> 5)
		tnza ^= sm & (tnza ^ ((tnza << 32) | (snza >> 32)))
		tnzb ^= sm & (tnzb ^ ((tnzb << 32) | (snzb >> 32)))
		s &= 31
		tnza <<= s
		tnzb <<= s

		tnza |= a[0] & ^(tnz1 | tnz2 | tnz3)
		tnzb |= b[0] & ^(tnz1 | tnz2 | tnz3)
		xa := (a[0] & 0x000000007FFFFFFF) | (tnza & 0xFFFFFFFF80000000)
		xb := (b[0] & 0x000000007FFFFFFF) | (tnzb & 0xFFFFFFFF80000000)

		// Run 29 iterations.
		var fg0 uint64 = 1
		var fg1 uint64 = uint64(1) << 32
		for j := 0; j < 29; j++ {
			a_odd := -(xa & 1)
			_, cc := bits.Sub64(xa, xb, 0)
			swap := a_odd & -cc
			ls += swap & ((xa & xb) >> 1)
			t := swap & (xa ^ xb)
			xa ^= t
			xb ^= t
			t = swap & (fg0 ^ fg1)
			fg0 ^= t
			fg1 ^= t
			xa -= a_odd & xb
			fg0 -= a_odd & fg1
			xa >>= 1
			fg1 <<= 1
			ls += (xb + 2) >> 2
		}

		// For the last two iterations, we need to recompute the
		// updated a and b (low words only) to get their low bits.
		f0 := ((fg0 + 0x7FFFFFFF7FFFFFFF) & 0xFFFFFFFF) - 0x7FFFFFFF
		g0 := ((fg0 + 0x7FFFFFFF7FFFFFFF) >> 32) - 0x7FFFFFFF
		f1 := ((fg1 + 0x7FFFFFFF7FFFFFFF) & 0xFFFFFFFF) - 0x7FFFFFFF
		g1 := ((fg1 + 0x7FFFFFFF7FFFFFFF) >> 32) - 0x7FFFFFFF
		a0 := (a[0]*f0 + b[0]*g0) >> 29
		b0 := (a[0]*f1 + b[0]*g1) >> 29
		for j := 0; j < 2; j++ {
			a_odd := -(xa & 1)
			_, cc := bits.Sub64(xa, xb, 0)
			swap := a_odd & -cc
			ls += swap & ((a0 & b0) >> 1)
			t := swap & (xa ^ xb)
			xa ^= t
			xb ^= t
			t = swap & (fg0 ^ fg1)
			fg0 ^= t
			fg1 ^= t
			t = swap & (a0 ^ b0)
			a0 ^= t
			b0 ^= t
			xa -= a_odd & xb
			fg0 -= a_odd & fg1
			a0 -= a_odd & b0
			xa >>= 1
			fg1 <<= 1
			a0 >>= 1
			ls += (b0 + 2) >> 2
		}

		// Extract individual update factors from the packed
		// representations fg0 and fg1.
		f0 = ((fg0 + 0x7FFFFFFF7FFFFFFF) & 0xFFFFFFFF) - 0x7FFFFFFF
		g0 = ((fg0 + 0x7FFFFFFF7FFFFFFF) >> 32) - 0x7FFFFFFF
		f1 = ((fg1 + 0x7FFFFFFF7FFFFFFF) & 0xFFFFFFFF) - 0x7FFFFFFF
		g1 = ((fg1 + 0x7FFFFFFF7FFFFFFF) >> 32) - 0x7FFFFFFF

		// Update a and b.
		var na, nb [4]uint64
		nega := gf_lin_div31_abs(&na, &a, &b, f0, g0)
		gf_lin_div31_abs(&nb, &a, &b, f1, g1)
		ls += nega & (nb[0] >> 1)
		a = na
		b = nb
	}

	// If y is invertible, then len(a) + len(b) <= 45 at this
	// point, so we can avoid the extraction, and 43 iterations are
	// sufficient. We do not need to compute update factors. Since
	// values are exact, we always have the proper bit values. We
	// can stop at 43 iterations, because at the last iteration,
	// b = 1 and a = 0 or 1; in either case, no modification of the
	// Legendre symbol will occur.
	xa := a[0]
	xb := b[0]
	for j := 0; j < 43; j++ {
		a_odd := -(xa & 1)
		_, cc := bits.Sub64(xa, xb, 0)
		swap := a_odd & -cc
		ls += swap & ((xa & xb) >> 1)
		t := swap & (xa ^ xb)
		xa ^= t
		xb ^= t
		xa -= a_odd & xb
		xa >>= 1
		ls += (xb + 2) >> 2
	}

	// If y != 0, then the low bit of ls contains its QR status
	// (0 = square, 1 = non-square). If y == 0, then we replace
	// the output value with 0, as per the API.
	return (1 - ((ls & 1) << 1)) & (gf_iszero(y, mq) - 1)
}
