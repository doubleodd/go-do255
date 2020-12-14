package do255s

import (
	"github.com/doubleodd/go-do255/internal/scalar"
	"math/bits"
)

// This file defines types and functions for do255e scalars, i.e.
// integers modulo r = 2^254 + 56904135270672826811114353017034461895
//
// Unless explicitly documented, all functions here are constant-time.

// Do255sScalar is the type for an integer modulo the prime order of the
// do255s group. Default value is zero.
type Do255sScalar [4]uint64

// Decode a scalar from exactly 32 bytes. Returned value is:
//   1   scalar properly decoded, value is not zero
//   0   scalar properly decoded, value is zero
//  -1   source bytes were not a valid scalar encoding
// If the decoding fails, then the scalar value is forced to zero.
func (s *Do255sScalar) Decode(src []byte) int {
	return scalar.Decode((*[4]uint64)(s), src, &do255sOrder)
}

// Decode a scalar from some bytes. All provided bytes are read and
// interpreted as an integer in unsigned little endian convention, which
// is reduced modulo the curve subgroup order. This process cannot fail.
func (s *Do255sScalar) DecodeReduce(src []byte) {
	scalar.DecodeReduce((*[4]uint64)(s), src, do255sModrReduce384Partial)
}

// Encode a scalar into exactly 32 bytes. The bytes are appended to the
// provided slice; the new slice is returned. The extension is done in
// place if the provided slice has enough capacity.
func (s *Do255sScalar) Encode(dst []byte) []byte {
	return scalar.Encode(dst, (*[4]uint64)(s), do255sModrReduce256)
}

// Encode a scalar into exactly 32 bytes; the newly allocated slice
// is returned.
func (s *Do255sScalar) Bytes() [32]byte {
	return scalar.ToBytes((*[4]uint64)(s), do255sModrReduce256)
}

// Compare a scalar with zero. Returned value is 1 if the scalar is zero,
// 0 otherwise.
func (s *Do255sScalar) IsZero() int {
	var t [4]uint64
	do255sModrReduce256(&t, (*[4]uint64)(s))
	z := t[0] | t[1] | t[2] | t[3]
	return int(1 - ((z | -z) >> 63))
}

// Compare two scalars together. Returned value is 1 if the scalars are
// equal to each other, 0 otherwise.
func (s *Do255sScalar) Equal(a *Do255sScalar) int {
	var t Do255sScalar
	t.Sub(s, a)
	return t.IsZero()
}

// Scalar addition: s is set to a + b (mod r).
// A pointer to s is returned.
func (s *Do255sScalar) Add(a, b *Do255sScalar) *Do255sScalar {
	scalar.Add((*[4]uint64)(s), (*[4]uint64)(a), (*[4]uint64)(b), do255sModrReduce256Partial)
	return s
}

// Scalar subtraction: s is set to a - b (mod r).
// A pointer to s is returned.
func (s *Do255sScalar) Sub(a, b *Do255sScalar) *Do255sScalar {
	scalar.Sub((*[4]uint64)(s), (*[4]uint64)(a), (*[4]uint64)(b), do255sModrReduce256Partial, &do255sOrder)
	return s
}

// Scalar negation: s is set to -a (mod r).
// A pointer to s is returned.
func (s *Do255sScalar) Neg(a *Do255sScalar) *Do255sScalar {
	var t = [4]uint64{0, 0, 0, 0}
	scalar.Sub((*[4]uint64)(s), &t, (*[4]uint64)(a), do255sModrReduce256Partial, &do255sOrder)
	return s
}

// Scalar multiplication: s is set to a*b (mod r).
// A pointer to s is returned.
func (s *Do255sScalar) Mul(a, b *Do255sScalar) *Do255sScalar {
	scalar.Mul((*[4]uint64)(s), (*[4]uint64)(a), (*[4]uint64)(b), do255sModrReduce384Partial)
	return s
}

// Group order is r = 2^254 + r0, with:
//    r0 = 56904135270672826811114353017034461895
const do255s_r0_lo uint64 = 0xDCF2AC65396152C7
const do255s_r0_hi uint64 = 0x2ACF567A912B7F03
const do255s_r_top uint64 = 0x4000000000000000

const do255s_r0x4_lo uint64 = 0x73CAB194E5854B1C
const do255s_r0x4_hi uint64 = 0xAB3D59EA44ADFC0F

var do255sOrder = [4]uint64{do255s_r0_lo, do255s_r0_hi, 0, do255s_r_top}

var do255s_r2 = [8]uint64{
	0xA31F34E2739216B1, 0x86A297C9835B5211,
	0x95DCE66BF04303AD, 0x8728B04D2F0F9E3C,
	0xEE7956329CB0A963, 0x1567AB3D4895BF81,
	0x0000000000000000, 0x1000000000000000,
}

// Given input 'a' (up to 2^286-1), perform a partial reduction modulo r;
// output (into 'd') fits on 255 bits and is lower than 2*r. The
// high bits of 'a' are provided as extra parameter ah.
func do255sModrReduce256PartialWithExtra(d, a *[4]uint64, ah uint64) {
	// Truncate to 254 bits and get extra bits into ah.
	ah = (ah << 2) | (a[3] >> 62)

	// Compute ah*r0 into u0:u1:u2.
	u1, u0 := bits.Mul64(ah, do255s_r0_lo)
	u2, lo := bits.Mul64(ah, do255s_r0_hi)
	var cc uint64
	u1, cc = bits.Add64(u1, lo, 0)
	u2 += cc

	// 2^254 = -r0 mod r
	d[0], cc = bits.Sub64(a[0], u0, 0)
	d[1], cc = bits.Sub64(a[1], u1, cc)
	d[2], cc = bits.Sub64(a[2], u2, cc)
	d[3], cc = bits.Sub64(a[3]&0x3FFFFFFFFFFFFFFF, 0, cc)

	// If we got a borrow, then we must add back r. Since ah*r0 < 2^192,
	// the result will be nonnegative, but less than r.
	m := -cc
	d[0], cc = bits.Add64(d[0], m&do255sOrder[0], 0)
	d[1], cc = bits.Add64(d[1], m&do255sOrder[1], cc)
	d[2], cc = bits.Add64(d[2], m&do255sOrder[2], cc)
	d[3] = d[3] + (m & do255sOrder[3]) + cc
}

// Partial reduction ensures that the output is fits on 255 bits and is
// less than 2*r.
func do255sModrReduce256Partial(d, a *[4]uint64) {
	do255sModrReduce256PartialWithExtra(d, a, 0)
}

// Given a partially reduced input 'a' (less than 2*r), finish reduction
// (conditional subtraction of r).
func do255sModrReduce256Finish(d, a *[4]uint64) {
	// Try to subtract r.
	var t [4]uint64
	var cc uint64
	t[0], cc = bits.Sub64(a[0], do255s_r0_lo, 0)
	t[1], cc = bits.Sub64(a[1], do255s_r0_hi, cc)
	t[2], cc = bits.Sub64(a[2], 0, cc)
	t[3], cc = bits.Sub64(a[3], do255s_r_top, cc)

	// If the result is nonnegative, then keep it; otherwise, use the
	// original value.
	m := -cc
	for i := 0; i < 4; i++ {
		d[i] = t[i] ^ (m & (a[i] ^ t[i]))
	}
}

// Perform full reduction of a scalar.
func do255sModrReduce256(d, a *[4]uint64) {
	do255sModrReduce256Partial(d, a)
	do255sModrReduce256Finish(d, d)
}

// Given a 384-bit input 'a', perform a partial reduction modulo r;
// output fits on 255 bits and is less than 2*r.
func do255sModrReduce384Partial(d *[4]uint64, a *[6]uint64) {
	// Multiply the high third (a4:a5) by 4*r0 into tw.
	var t1, t2 [2]uint64
	var tw [4]uint64
	t1[0] = do255s_r0x4_lo
	t1[1] = do255s_r0x4_hi
	t2[0] = a[4]
	t2[1] = a[5]
	scalar.Mul128x128(&tw, &t1, &t2)

	// Subtract 4*r0*ah from the low part of 'a', then
	// add back 4*r. Since 4*r0 =~ 2^127.42, the result may be
	// slightly above 2^257, but will fit on 258 bits.
	var cc uint64
	tw[0], cc = bits.Sub64(a[0], tw[0], 0)
	tw[1], cc = bits.Sub64(a[1], tw[1], cc)
	tw[2], cc = bits.Sub64(a[2], tw[2], cc)
	tw[3], cc = bits.Sub64(a[3], tw[3], cc)
	tw4 := -cc
	tw[0], cc = bits.Add64(tw[0], do255s_r0x4_lo, 0)
	tw[1], cc = bits.Add64(tw[1], do255s_r0x4_hi, cc)
	tw[2], cc = bits.Add64(tw[2], 0, cc)
	tw[3], cc = bits.Add64(tw[3], 0, cc)
	tw4 += (cc + 1)

	// Perform partial reduction.
	do255sModrReduce256PartialWithExtra(d, &tw, tw4)
}

// Compare two 512-bit nonnegative integers; return value is true if
// a < b, false otherwise.
// THIS IS NOT CONSTANT-TIME.
func cmp512ltVartime(a, b *[8]uint64) bool {
	for i := 7; i >= 0; i-- {
		if a[i] < b[i] {
			return true
		}
		if a[i] > b[i] {
			return false
		}
	}
	return false
}

// Get the bit length of a 512-bit signed integer.
// THIS IS NOT CONSTANT-TIME.
func bitlength512Vartime(a *[8]uint64) int {
	m := -(a[7] >> 63)
	for i := 7; i >= 0; i-- {
		aw := a[i] ^ m
		if aw != 0 {
			return (i << 6) + 64 - bits.LeadingZeros64(aw)
		}
	}
	return 0
}

// Add a*2^s to d.
// THIS IS NOT CONSTANT-TIME.
func addLshift192Vartime(d, a *[3]uint64, s int) {
	if s < 64 {
		if s == 0 {
			var cc uint64
			d[0], cc = bits.Add64(d[0], a[0], 0)
			d[1], cc = bits.Add64(d[1], a[1], cc)
			d[2], _ = bits.Add64(d[2], a[2], cc)
		} else {
			var t [3]uint64
			t[0] = a[0] << uint(s)
			t[1] = (a[1] << uint(s)) | (a[0] >> uint(64-s))
			t[2] = (a[2] << uint(s)) | (a[1] >> uint(64-s))
			var cc uint64
			d[0], cc = bits.Add64(d[0], t[0], 0)
			d[1], cc = bits.Add64(d[1], t[1], cc)
			d[2], _ = bits.Add64(d[2], t[2], cc)
		}
	} else {
		// This case is rare.
		if s >= 192 {
			return
		}
		if (s & 63) == 0 {
			if s == 64 {
				var cc uint64
				d[1], cc = bits.Add64(d[1], a[0], 0)
				d[2], _ = bits.Add64(d[2], a[1], cc)
			} else { // s == 128
				d[2] += a[0]
			}
			return
		} else {
			if s < 128 {
				a0 := a[0] << uint(s-64)
				a1 := (a[1] << uint(s-64)) | (a[0] >> uint(128-s))
				var cc uint64
				d[1], cc = bits.Add64(d[1], a0, 0)
				d[2], _ = bits.Add64(d[2], a1, cc)
			} else {
				d[2] += a[0] << uint(s-128)
			}
		}
	}
}

// Subtract a*2^s from d.
// THIS IS NOT CONSTANT-TIME.
func subLshift192Vartime(d, a *[3]uint64, s int) {
	if s < 64 {
		if s == 0 {
			var cc uint64
			d[0], cc = bits.Sub64(d[0], a[0], 0)
			d[1], cc = bits.Sub64(d[1], a[1], cc)
			d[2], _ = bits.Sub64(d[2], a[2], cc)
		} else {
			var t [3]uint64
			t[0] = a[0] << uint(s)
			t[1] = (a[1] << uint(s)) | (a[0] >> uint(64-s))
			t[2] = (a[2] << uint(s)) | (a[1] >> uint(64-s))
			var cc uint64
			d[0], cc = bits.Sub64(d[0], t[0], 0)
			d[1], cc = bits.Sub64(d[1], t[1], cc)
			d[2], _ = bits.Sub64(d[2], t[2], cc)
		}
	} else {
		// This case is rare.
		if s >= 192 {
			return
		}
		if (s & 63) == 0 {
			if s == 64 {
				var cc uint64
				d[1], cc = bits.Sub64(d[1], a[0], 0)
				d[2], _ = bits.Sub64(d[2], a[1], cc)
			} else { // s == 128
				d[2] -= a[0]
			}
			return
		} else {
			if s < 128 {
				a0 := a[0] << uint(s-64)
				a1 := (a[1] << uint(s-64)) | (a[0] >> uint(128-s))
				var cc uint64
				d[1], cc = bits.Sub64(d[1], a0, 0)
				d[2], _ = bits.Sub64(d[2], a1, cc)
			} else {
				d[2] -= a[0] << uint(s-128)
			}
		}
	}
}

// Add a*2^s to d.
// THIS IS NOT CONSTANT-TIME.
func addLshift512Vartime(d, a *[8]uint64, s int) {
	if s < 64 {
		if s == 0 {
			var cc uint64 = 0
			for i := 0; i < 8; i++ {
				d[i], cc = bits.Add64(d[i], a[i], cc)
			}
		} else {
			var t [8]uint64
			t[0] = a[0] << uint(s)
			for i := 1; i < 8; i++ {
				t[i] = (a[i] << uint(s)) | (a[i-1] >> uint(64-s))
			}
			var cc uint64 = 0
			for i := 0; i < 8; i++ {
				d[i], cc = bits.Add64(d[i], t[i], cc)
			}
		}
	} else {
		// This case is rare.
		if s >= 512 {
			return
		}
		j := s >> 6
		s &= 63
		var t [8]uint64
		if s == 0 {
			for i := j; i < 8; i++ {
				t[i] = a[i-j]
			}
		} else {
			t[j] = a[0] << uint(s)
			for i := j + 1; i < 8; i++ {
				t[i] = (a[i-j] << uint(s)) | (a[i-j-1] >> uint(64-s))
			}
		}
		var cc uint64 = 0
		for i := j; i < 8; i++ {
			d[i], cc = bits.Add64(d[i], t[i], cc)
		}
	}
}

// Subtract a*2^s from d.
// THIS IS NOT CONSTANT-TIME.
func subLshift512Vartime(d, a *[8]uint64, s int) {
	if s < 64 {
		if s == 0 {
			var cc uint64 = 0
			for i := 0; i < 8; i++ {
				d[i], cc = bits.Sub64(d[i], a[i], cc)
			}
		} else {
			var t [8]uint64
			t[0] = a[0] << uint(s)
			for i := 1; i < 8; i++ {
				t[i] = (a[i] << uint(s)) | (a[i-1] >> uint(64-s))
			}
			var cc uint64 = 0
			for i := 0; i < 8; i++ {
				d[i], cc = bits.Sub64(d[i], t[i], cc)
			}
		}
	} else {
		// This case is rare.
		if s >= 512 {
			return
		}
		j := s >> 6
		s &= 63
		var t [8]uint64
		if s == 0 {
			for i := j; i < 8; i++ {
				t[i] = a[i-j]
			}
		} else {
			t[j] = a[0] << uint(s)
			for i := j + 1; i < 8; i++ {
				t[i] = (a[i-j] << uint(s)) | (a[i-j-1] >> uint(64-s))
			}
		}
		var cc uint64 = 0
		for i := j; i < 8; i++ {
			d[i], cc = bits.Sub64(d[i], t[i], cc)
		}
	}
}

// Split an input scalar k into shorter integers c0 and c1 such that
// k = c0/c1 mod r, with |c0| < 2^128 and |c1| < 2^128. The absolute
// values of c0 and c1 are returned in c0 and c1, respectively, while
// the signs of c0 and c1 are the returned values (true for negative,
// false for nonnegative).
// THIS IS NOT CONSTANT-TIME.
func (k *Do255sScalar) ReduceBasisVartime(c0, c1 *[2]uint64) (negc0 bool, negc1 bool) {
	// Algorithm is described in: https://eprint.iacr.org/2020/454
	//
	// Init:
	//   reduce k modulo r
	//   u = [r, 0]
	//   v = [k, 1]
	//   nu = r^2
	//   nv = k^2 + 1
	//   sp = r*k
	//
	// Loop:
	//   - if nu < nv then:
	//        (u, v) <- (v, u)
	//        (nu, nv) <- (nv, nu)
	//   - if bitlength(nv) <= 255 then:
	//        return (v0, v1)
	//   - s <- max(0, bitlength(sp) - bitlength(nv))
	//   - if sp > 0 then:
	//        u <- u - lshift(v, s)
	//        nu <- nu + lshift(nv, 2*s) - lshift(sp, s+1)
	//        sp <- sp - lshift(nv, s)
	//     else:
	//        u <- u + lshift(v, s)
	//        nu <- nu + lshift(nv, 2*s) + lshift(sp, s+1)
	//        sp <- sp + lshift(nv, s)

	// Reduce k mod r.
	var vk [4]uint64
	do255sModrReduce256(&vk, (*[4]uint64)(k))

	// u <- (r, 0)
	var u0, u1 [3]uint64
	copy(u0[:], do255sOrder[:])

	// v <- (k, 1)
	var v0, v1 [3]uint64
	copy(v0[:], vk[:])
	v1[0] = 1

	// nu <- r^2
	var nu [8]uint64
	copy(nu[:], do255s_r2[:])

	// nv <- k^2 + 1
	var nv [8]uint64
	scalar.Mul256x256(&nv, &vk, &vk)
	var cc uint64 = 1
	for i := 0; i < 8; i++ {
		nv[i], cc = bits.Add64(nv[i], 0, cc)
	}

	// sp <- r*k
	var sp [8]uint64
	scalar.Mul256x256(&sp, &vk, &do255sOrder)

	// Main algorithm loop.
	for {
		// If nu < nv, then swap(u,v) and swap(nu,nv)
		if cmp512ltVartime(&nu, &nv) {
			for i := 0; i < 3; i++ {
				tt := u0[i]
				u0[i] = v0[i]
				v0[i] = tt
				tt = u1[i]
				u1[i] = v1[i]
				v1[i] = tt
			}
			for i := 0; i < 8; i++ {
				tt := nu[i]
				nu[i] = nv[i]
				nv[i] = tt
			}
		}

		// Get bit length of nv; if it is no more than 255,
		// then v is small enough to be returned. We know that
		// we can get ||v||^2 down to about 1.075*r, i.e. only
		// slightly above 2^254, so this condition can always
		// be reached.
		bl_nv := bitlength512Vartime(&nv)
		if bl_nv <= 255 {
			if v0[2] == 0 {
				c0[0] = v0[0]
				c0[1] = v0[1]
				negc0 = false
			} else {
				var cc uint64
				c0[0], cc = bits.Sub64(0, v0[0], 0)
				c0[1], _ = bits.Sub64(0, v0[1], cc)
				negc0 = true
			}
			if v1[2] == 0 {
				c1[0] = v1[0]
				c1[1] = v1[1]
				negc1 = false
			} else {
				var cc uint64
				c1[0], cc = bits.Sub64(0, v1[0], 0)
				c1[1], _ = bits.Sub64(0, v1[1], cc)
				negc1 = true
			}
			return
		}

		// Get shift count:
		//    s = min(0, bitlength(sp) - bitlength(nv))
		bl_sp := bitlength512Vartime(&sp)
		s := bl_sp - bl_nv
		s &= ^(s >> 31)

		// Subtract or add, depending on the sign of sp.
		if (sp[7] >> 63) == 0 {
			subLshift192Vartime(&u0, &v0, s)
			subLshift192Vartime(&u1, &v1, s)
			addLshift512Vartime(&nu, &nv, 2*s)
			subLshift512Vartime(&nu, &sp, s+1)
			subLshift512Vartime(&sp, &nv, s)
		} else {
			addLshift192Vartime(&u0, &v0, s)
			addLshift192Vartime(&u1, &v1, s)
			addLshift512Vartime(&nu, &nv, 2*s)
			addLshift512Vartime(&nu, &sp, s+1)
			addLshift512Vartime(&sp, &nv, s)
		}
	}
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
func (a *Do255sScalar) recode5(d *[52]byte) {
	scalar.Recode5(d, (*[4]uint64)(a))
}
