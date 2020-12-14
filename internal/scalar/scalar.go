package scalar

import (
	"encoding/binary"
	"math/bits"
)

// This file contains some helper functions that are used by implementations
// of scalars for both do255e and do255s. A scalar is an integer modulo
// the prime r = 2^254 + r0, where |r0| < 2^127. Some generic functions
// operating on big integers (of a given size) are also provided.
//
// As a general rule, operations on scalars are not critical for performance
// and can tolerate suboptimal implementations; but they must be strictly
// constant-time because some scalar values are secret (e.g. private keys).

// Extend a slice for appending n bytes. The two returned values are the
// new extended slice (no extra allocation if the original slice was large
// enough), and the sub-slice where data should be written.
// (Inspired by https://github.com/gtank/ristretto255 )
func prepareAppend(b []byte, n int) (head, tail []byte) {
	len1 := len(b)   // current length
	len2 := len1 + n // new length after extension
	if cap(b) >= len2 {
		head = b[:len2]
	} else {
		head = make([]byte, len2)
		copy(head, b)
	}
	tail = head[len1:]
	return
}

// 128x128->128 multiplication.
func Mul128x128trunc(d, a, b *[2]uint64) {
	t1, t0 := bits.Mul64(a[0], b[0])
	t1 += a[0]*b[1] + a[1]*b[0]
	d[0] = t0
	d[1] = t1
}

// 128x128->256 multiplication.
func Mul128x128(d *[4]uint64, a, b *[2]uint64) {
	var lo, hi, cc uint64
	d[1], d[0] = bits.Mul64(a[0], b[0])
	d[3], d[2] = bits.Mul64(a[1], b[1])
	hi, lo = bits.Mul64(a[0], b[1])
	d[1], cc = bits.Add64(d[1], lo, 0)
	d[2], cc = bits.Add64(d[2], hi, cc)
	d[3] += cc
	hi, lo = bits.Mul64(a[1], b[0])
	d[1], cc = bits.Add64(d[1], lo, 0)
	d[2], cc = bits.Add64(d[2], hi, cc)
	d[3] += cc
}

// 256x128->384 multiplication.
func Mul256x128(d *[6]uint64, a *[4]uint64, b *[2]uint64) {
	var c0, c1 [2]uint64
	var t0, t1 [4]uint64
	c0[0] = a[0]
	c0[1] = a[1]
	Mul128x128(&t0, &c0, b)
	c1[0] = a[2]
	c1[1] = a[3]
	Mul128x128(&t1, &c1, b)
	var cc uint64
	d[0] = t0[0]
	d[1] = t0[1]
	d[2], cc = bits.Add64(t0[2], t1[0], 0)
	d[3], cc = bits.Add64(t0[3], t1[1], cc)
	d[4], cc = bits.Add64(0, t1[2], cc)
	d[5] = t1[3] + cc
}

// 256x256->512 multiplication.
func Mul256x256(d *[8]uint64, a *[4]uint64, b *[4]uint64) {
	var c0, c1 [2]uint64
	var t0, t1 [6]uint64
	c0[0] = b[0]
	c0[1] = b[1]
	Mul256x128(&t0, a, &c0)
	c1[0] = b[2]
	c1[1] = b[3]
	Mul256x128(&t1, a, &c1)
	var cc uint64
	d[0] = t0[0]
	d[1] = t0[1]
	d[2], cc = bits.Add64(t0[2], t1[0], 0)
	d[3], cc = bits.Add64(t0[3], t1[1], cc)
	d[4], cc = bits.Add64(t0[4], t1[2], cc)
	d[5], cc = bits.Add64(t0[5], t1[3], cc)
	d[6], cc = bits.Add64(0, t1[4], cc)
	d[7] = t1[5] + cc
}

// Decode a scalar value from bytes. Modulus r is provided. Returned
// value:
//   1   decode successful, value is in range and non-zero
//   0   decode successful, value is zero
//  -1   decode failed, value is out of range.
// On error, output value (in d[]) is forced to zero.
func Decode(d *[4]uint64, src []byte, r *[4]uint64) int {
	// Decode in little-endian.
	for i := 0; i < 4; i++ {
		d[i] = binary.LittleEndian.Uint64(src[8*i:])
	}

	// Check whether all bytes were zero.
	zz := d[0] | d[1] | d[2] | d[3]
	zz = 1 - ((zz | -zz) >> 63)

	// Compare value with r; if not lower (borrow is zero), then
	// this is invalid.
	var cc uint64 = 0
	for i := 0; i < 4; i++ {
		_, cc = bits.Sub64(d[i], r[i], cc)
	}
	for i := 0; i < 4; i++ {
		d[i] &= -cc
	}

	// If input was valid, then cc == 1; otherwise, cc == 0. If
	// input was zero, then cc == 1 (it was valid) and zz == 1;
	// otherwise, zz == 0.
	return int(int64(((cc << 1) - zz) - 1))
}

// Type for a scalar reduction function: input is a 256-bit integer, output
// is normalized into the 0..r-1 range.
type Reduce256 func(*[4]uint64, *[4]uint64)

// Type for a scalar reduction function: input is a 384-bit integer, output
// fits on 256 bits (but is not necessarily normalized to 0..r-1).
type Reduce384 func(*[4]uint64, *[6]uint64)

// Encode a scalar into exactly 32 bytes. The scalar is reduced by
// invoking the provided reduction function. The bytes are appended
// to the provided slice. The extension is done in place if the
// provided slice has enough capacity. The new slice is returned.
func Encode(b []byte, s *[4]uint64, rf Reduce256) []byte {
	b2, dst := prepareAppend(b, 32)
	var t [4]uint64
	rf(&t, s)
	for i := 0; i < 4; i++ {
		binary.LittleEndian.PutUint64(dst[8*i:], t[i])
	}
	return b2
}

// Encode a scalar into exactly 32 bytes. The scalar is reduced by
// invoking the provided reduction function.
func ToBytes(s *[4]uint64, rf Reduce256) [32]byte {
	var dst [32]byte
	Encode(dst[:0], s, rf)
	return dst
}

// Decode a scalar from bytes; the bytes are interpreted with unsigned
// little-endian convention into a big integer, which is reduced modulo
// the curve subgroup order r. All bytes from the input slice are used.
// If the input slice is empty, then the obtained value is 0. The
// reduction is applied with the provided function (rf) for reduction
// 384->256.
func DecodeReduce(d *[4]uint64, src []byte, rf Reduce384) {
	n := len(src)

	// Set output to 0.
	for i := 0; i < 4; i++ {
		d[i] = 0
	}

	// Special case: empty slice.
	if n == 0 {
		return
	}

	// Fill the scalar with the last chunk. We put as many bytes as
	// we can in it, provided that the remaining number of bytes (j)
	// is a multiple of 32.
	var j int
	if n >= 32 {
		j = n - (n & 15) - 16
		if j == (n - 16) {
			j = n - 32
		}
	} else {
		j = 0
	}
	for i := 0; i < (n - j); i++ {
		d[i>>3] |= uint64(src[j+i]) << uint((i&7)<<3)
	}

	// For all remaining chunks of 16 bytes, multiply the current
	// value by 2^128 (left shift), add the new chunk, and do a
	// reduction round.
	for j > 0 {
		j -= 16
		var t [6]uint64
		t[0] = binary.LittleEndian.Uint64(src[j:])
		t[1] = binary.LittleEndian.Uint64(src[j+8:])
		copy(t[2:], d[:])
		rf(d, &t)
	}
}

// Scalar addition; partial reduction function is provided (rf). The
// reduction function must ensure that the result fits on 255 bits.
func Add(d, a, b *[4]uint64, rf Reduce256) {
	var t1, t2 [4]uint64
	rf(&t1, a)
	rf(&t2, b)
	var cc uint64 = 0
	for i := 0; i < 4; i++ {
		d[i], cc = bits.Add64(t1[i], t2[i], cc)
	}
	// No output carry is possible, since both inputs were reduced
	// to less than 2^255.
}

// Scalar subtraction; partial reduction function (rf) and order (r)
// are provided. The reductin function must ensure that the results
// is less than 2*r.
func Sub(d, a, b *[4]uint64, rf Reduce256, r *[4]uint64) {
	// Reduce second operand to less than 2*r.
	var t2 [4]uint64
	rf(&t2, b)

	// Perform subtraction.
	var cc uint64 = 0
	for i := 0; i < 4; i++ {
		d[i], cc = bits.Sub64(a[i], t2[i], cc)
	}

	// If there is an output borrow, then we must add 2*r. Since
	// the second input was reduced to less than 2*r, adding 2*r
	// once is enough. Moreover, r < 2^255, so 2*r fits on 256 bits.
	var r2 [4]uint64
	r2[0] = -cc & (r[0] << 1)
	r2[1] = -cc & ((r[1] << 1) | (r[0] >> 63))
	r2[2] = -cc & ((r[2] << 1) | (r[1] >> 63))
	r2[3] = -cc & ((r[3] << 1) | (r[2] >> 63))
	cc = 0
	for i := 0; i < 4; i++ {
		d[i], cc = bits.Add64(d[i], r2[i], cc)
	}
}

// Scalar multiplication; partial reduction function (rf, for 384->256)
// is prodived.
func Mul(d, a, b *[4]uint64, rf Reduce384) {
	var t6 [6]uint64
	var t8 [8]uint64
	Mul256x256(&t8, a, b)
	copy(t6[:], t8[2:])
	rf(d, &t6)
	copy(t6[:], t8[0:2])
	copy(t6[2:], d[:])
	rf(d, &t6)
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
func Recode5(d *[52]byte, a *[4]uint64) {
	acc := a[0]
	acc_len := 64
	j := 1
	var cc uint = 0
	for i := 0; i < 51; i++ {
		var b uint
		if acc_len < 5 {
			next := a[j]
			j++
			b = uint(acc|(next<<uint(acc_len))) & 31
			acc = next >> uint(5-acc_len)
			acc_len = 59 + acc_len
		} else {
			b = uint(acc) & 31
			acc >>= 5
			acc_len -= 5
		}
		b += cc
		m := (16 - b) >> 8
		b ^= m & (b ^ (160 - b))
		cc = m & 1
		d[i] = byte(b)
	}
	d[51] = byte(uint(acc) + cc)
}

// Recode a small _unsigned_ 128-bit integer with 5-bit Booth encoding.
// Output is 26 digits, the top digit is necessarily nonnegative (and
// cannot be -0).
func Recode5Small(d *[26]byte, k *[2]uint64) {
	// First 12 digits from the low limb.
	var db uint64 = 0
	t := k[0]
	for i := 0; i < 12; i++ {
		b := (t & 0x1F) + db
		m := (16 - b) >> 8
		b ^= m & (b ^ (160 - b))
		db = m & 1
		d[i] = byte(b)
		t >>= 5
	}

	// Get more bits from the high limb for the next 12 digits.
	t |= k[1] << 4
	for i := 12; i < 24; i++ {
		b := (t & 0x1F) + db
		m := (16 - b) >> 8
		b ^= m & (b ^ (160 - b))
		db = m & 1
		d[i] = byte(b)
		t >>= 5
	}

	// Last two digits.
	t = k[1] >> 56
	b := (t & 0x1F) + db
	m := (16 - b) >> 8
	b ^= m & (b ^ (160 - b))
	db = m & 1
	d[24] = byte(b)
	t >>= 5
	d[25] = byte(t + db)
}

// Recode a small _signed_ 128-bit integer with 5-bit Booth encoding.
// If the source value is negative, then what is recoded is its absolute
// value; the source sign is returned (1 for negative, 0 for zero or
// positive). Output is 26 digits, the top digit is necessarily
// nonnegative (and cannot be -0).
func Recode5SmallSigned(d *[26]byte, k *[2]uint64) uint64 {
	// Compute abs(k) (in x0:x1) and record its sign (in sk).
	sk := k[1] >> 63
	x0, cc := bits.Add64(k[0]^-sk, sk, 0)
	x1 := (k[1] ^ -sk) + cc

	// First 12 digits from the low limb.
	var db uint64 = 0
	t := x0
	for i := 0; i < 12; i++ {
		b := (t & 0x1F) + db
		m := (16 - b) >> 8
		b ^= m & (b ^ (160 - b))
		db = m & 1
		d[i] = byte(b)
		t >>= 5
	}

	// Get more bits from the high limb for the next 12 digits.
	t |= x1 << 4
	for i := 12; i < 24; i++ {
		b := (t & 0x1F) + db
		m := (16 - b) >> 8
		b ^= m & (b ^ (160 - b))
		db = m & 1
		d[i] = byte(b)
		t >>= 5
	}

	// Last two digits.
	t = x1 >> 56
	b := (t & 0x1F) + db
	m := (16 - b) >> 8
	b ^= m & (b ^ (160 - b))
	db = m & 1
	d[24] = byte(b)
	t >>= 5
	d[25] = byte(t + db)

	// Return original sign.
	return sk
}
