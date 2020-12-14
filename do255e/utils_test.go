package do255e

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
