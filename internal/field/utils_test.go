package field

import (
	"crypto/sha512"
	"encoding/binary"
	"fmt"
	"math/big"
)

// =====================================================================
// Custom PRNG (based on SHA-512) for reproducible tests.

type prng struct {
	buf [64]byte
	ptr int
}

// Initialize the PRNG with an explicit seed.
func (p *prng) init(seed string) {
	hv := sha512.Sum512([]byte(seed))
	copy(p.buf[:], hv[:])
	p.ptr = 0
}

// Fill the provided slice with pseudorandom bytes from the PRNG.
func (p *prng) generate(d []byte) {
	n := len(d)
	for n > 0 {
		c := 32 - p.ptr
		if c == 0 {
			hv := sha512.Sum512(p.buf[:])
			copy(p.buf[:], hv[:])
			p.ptr = 0
			c = 32
		}
		if c > n {
			c = n
		}
		copy(d, p.buf[p.ptr:p.ptr+c])
		d = d[c:]
		n -= c
		p.ptr += c
	}
}

// Generate a random 256-bit integer from the PRNG.
func (p *prng) mk256(d *[4]uint64) {
	var bb [32]byte
	p.generate(bb[:])
	for i := 0; i < 4; i++ {
		d[i] = binary.LittleEndian.Uint64(bb[8*i:])
	}
}

// Make a new random field element from the PRNG.
func (p *prng) mkgf(d *[4]uint64) {
	var t [4]uint64
	p.mk256(&t)
	copy(d[:], t[:])
}

// Create a new big integer by reducing the provided 256-bit integer a[]
// modulo m.
func int256ToBigMod(a *[4]uint64, m *big.Int) big.Int {
	var x, y big.Int
	for i := 3; i >= 0; i-- {
		y.SetUint64(a[i])
		x.Lsh(&x, 64).Add(&x, &y)
	}
	for x.Cmp(m) >= 0 {
		x.Sub(&x, m)
	}
	return x
}

// Get the string representation of a 256-bit integer (hexadecimal, with
// '0x' prefix).
func int256ToString(a *[4]uint64) string {
	return fmt.Sprintf("0x%016X%016X%016X%016X", a[3], a[2], a[1], a[0])
}

// Convert a field element to a string which can be copy-pasted into Sage
// for easier verifications (provided that 'K' was defined in Sage as the
// relevant field).
func gfToString(a *[4]uint64) string {
	var t [4]uint64
	copy(t[:], a[:])
	return "K(" + int256ToString(&t) + ")"
}

// Convert an (internal) field element representation to a big integer
// modulo the provided integer p.
func gfToBig(a *[4]uint64, p *big.Int) big.Int {
	var t [4]uint64
	copy(t[:], a[:])
	return int256ToBigMod(&t, p)
}

// Decode a sequence of bytes into a big integer, with unsigned little-endian
// convention.
func decodeToBigLE(src []byte) big.Int {
	n := len(src)
	tt := make([]byte, n)
	for i := 0; i < n; i++ {
		tt[i] = src[n-1-i]
	}
	var x big.Int
	x.SetBytes(tt)
	return x
}
