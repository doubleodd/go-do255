package field

import (
	"fmt"
	"math/big"
	"testing"
)

// Tests for GF255s (integers modulo p = 2^255 - 3957).

// =====================================================================

func TestGF255sAdd(t *testing.T) {
	var rng prng
	rng.init("test add GF255s")
	var a, b, c GF255s
	var p, mq big.Int
	mq.SetUint64(3957)
	p.SetUint64(1).Lsh(&p, 255).Sub(&p, &mq)
	for i := 0; i < 100000; i++ {
		if i < 30000 {
			for j := 0; j < 4; j++ {
				a[j] = 0xFFFFFFFFFFFFFFFF
				b[j] = 0xFFFFFFFFFFFFFFFF
			}
			a[3] = 0xFFFFFFFFFFFFFFFF - uint64(i)
		} else {
			rng.mkgf((*[4]uint64)(&a))
			rng.mkgf((*[4]uint64)(&b))
		}
		c.Add(&a, &b)

		za := gfToBig((*[4]uint64)(&a), &p)
		zb := gfToBig((*[4]uint64)(&b), &p)
		zc := gfToBig((*[4]uint64)(&c), &p)
		var zd big.Int
		zd.Add(&za, &zb)
		if zd.Cmp(&p) >= 0 {
			zd.Sub(&zd, &p)
		}
		if zc.Cmp(&zd) != 0 {
			t.Fatalf("ERR add:\na = %s\nb = %s\nc = %s\n", gfToString((*[4]uint64)(&a)), gfToString((*[4]uint64)(&b)), gfToString((*[4]uint64)(&c)))
		}
	}
}

func TestGF255sSub(t *testing.T) {
	var rng prng
	rng.init("test sub GF255s")
	var a, b, c GF255s
	var p, mq big.Int
	mq.SetUint64(3957)
	p.SetUint64(1).Lsh(&p, 255).Sub(&p, &mq)
	for i := 0; i < 100000; i++ {
		if i < 30000 {
			for j := 0; j < 4; j++ {
				a[j] = 0xFFFFFFFFFFFFFFFF
				b[j] = 0xFFFFFFFFFFFFFFFF
			}
			a[3] = 0xFFFFFFFFFFFFFFFF - uint64(i)
		} else {
			rng.mkgf((*[4]uint64)(&a))
			rng.mkgf((*[4]uint64)(&b))
		}

		c.Sub(&a, &b)
		za := gfToBig((*[4]uint64)(&a), &p)
		zb := gfToBig((*[4]uint64)(&b), &p)
		zc := gfToBig((*[4]uint64)(&c), &p)
		var zd big.Int
		zd.Sub(&za, &zb)
		if zd.Sign() < 0 {
			zd.Add(&zd, &p)
		}
		if zc.Cmp(&zd) != 0 {
			t.Fatalf("ERR sub:\na = %s\nb = %s\nc = %s\n", gfToString((*[4]uint64)(&a)), gfToString((*[4]uint64)(&b)), gfToString((*[4]uint64)(&c)))
		}

		c.Sub(&b, &a)
		zc = gfToBig((*[4]uint64)(&c), &p)
		zd.Sub(&zb, &za)
		if zd.Sign() < 0 {
			zd.Add(&zd, &p)
		}
		if zc.Cmp(&zd) != 0 {
			t.Fatalf("ERR sub:\na = %s\nb = %s\nc = %s\n", gfToString((*[4]uint64)(&b)), gfToString((*[4]uint64)(&a)), gfToString((*[4]uint64)(&c)))
		}
	}
}

func TestGF255sNeg(t *testing.T) {
	var rng prng
	rng.init("test neg GF255s")
	var a, c GF255s
	var p, mq big.Int
	mq.SetUint64(3957)
	p.SetUint64(1).Lsh(&p, 255).Sub(&p, &mq)
	for i := 0; i < 100000; i++ {
		if i < 30000 {
			for j := 0; j < 3; j++ {
				a[j] = 0xFFFFFFFFFFFFFFFF
			}
			a[3] = 0xFFFFFFFFFFFFFFFF - uint64(i)
		} else if i < 60000 {
			a[0] = uint64(i - 30000)
			for j := 1; j < 4; j++ {
				a[j] = 0
			}
		} else {
			rng.mkgf((*[4]uint64)(&a))
		}

		c.Neg(&a)
		za := gfToBig((*[4]uint64)(&a), &p)
		zc := gfToBig((*[4]uint64)(&c), &p)
		var zd big.Int
		zd.Neg(&za)
		if zd.Sign() < 0 {
			zd.Add(&zd, &p)
		}
		if zc.Cmp(&zd) != 0 {
			t.Fatalf("ERR neg:\na = %s\nc = %s\n", gfToString((*[4]uint64)(&a)), gfToString((*[4]uint64)(&c)))
		}
	}
}

func TestGF255sDecode(t *testing.T) {
	var rng prng
	rng.init("test decode GF255s")
	var p, mq big.Int
	mq.SetUint64(3957)
	p.SetUint64(1).Lsh(&p, 255).Sub(&p, &mq)
	var bb [100]byte
	for i := 0; i < 1000; i++ {
		rng.generate(bb[:])
		for j := 0; j <= len(bb); j++ {
			var a GF255s
			a.DecodeReduce(bb[:j])
			za := gfToBig((*[4]uint64)(&a), &p)
			zb := decodeToBigLE(bb[:j])
			zb.Mod(&zb, &p)
			if za.Cmp(&zb) != 0 {
				sbb := "0x"
				for k := j - 1; k >= 0; k-- {
					sbb += fmt.Sprintf("%02X", bb[k])
				}
				t.Fatalf("ERR decode:\nsrc = %s\na = %s\nb = %s\n", sbb, gfToString((*[4]uint64)(&a)), "0x"+zb.Text(16))
			}
		}
	}
}

func TestGF255sMul(t *testing.T) {
	var rng prng
	rng.init("test mul GF255s")
	var a, b, c GF255s
	var p, mq big.Int
	mq.SetUint64(3957)
	p.SetUint64(1).Lsh(&p, 255).Sub(&p, &mq)
	for i := 0; i < 100000; i++ {
		if i < 30000 {
			for j := 0; j < 4; j++ {
				a[j] = 0xFFFFFFFFFFFFFFFF
				b[j] = 0xFFFFFFFFFFFFFFFF
			}
			a[3] = 0xFFFFFFFFFFFFFFFF - uint64(i)
		} else {
			rng.mkgf((*[4]uint64)(&a))
			rng.mkgf((*[4]uint64)(&b))
		}

		c.Mul(&a, &b)
		za := gfToBig((*[4]uint64)(&a), &p)
		zb := gfToBig((*[4]uint64)(&b), &p)
		zc := gfToBig((*[4]uint64)(&c), &p)
		var zd big.Int
		zd.Mul(&za, &zb).Mod(&zd, &p)
		if zc.Cmp(&zd) != 0 {
			t.Fatalf("ERR mul:\na = %s\nb = %s\nc = %s\n", gfToString((*[4]uint64)(&a)), gfToString((*[4]uint64)(&b)), gfToString((*[4]uint64)(&c)))
		}

		c.Mul(&b, &a)
		zc = gfToBig((*[4]uint64)(&c), &p)
		if zc.Cmp(&zd) != 0 {
			t.Fatalf("ERR mul:\na = %s\nb = %s\nc = %s\n", gfToString((*[4]uint64)(&b)), gfToString((*[4]uint64)(&a)), gfToString((*[4]uint64)(&c)))
		}

		c.Sqr(&a)
		zc = gfToBig((*[4]uint64)(&c), &p)
		zd.Mul(&za, &za).Mod(&zd, &p)
		if zc.Cmp(&zd) != 0 {
			t.Fatalf("ERR sqr:\na = %s\nc = %s\n", gfToString((*[4]uint64)(&a)), gfToString((*[4]uint64)(&c)))
		}
	}
}

func TestGF255sShift(t *testing.T) {
	var rng prng
	rng.init("test shift GF255s")
	var a, c GF255s
	var p, mq big.Int
	mq.SetUint64(3957)
	p.SetUint64(1).Lsh(&p, 255).Sub(&p, &mq)
	for i := 0; i < 100000; i++ {
		if i < 30000 {
			for j := 0; j < 3; j++ {
				a[j] = 0xFFFFFFFFFFFFFFFF
			}
			a[3] = 0xFFFFFFFFFFFFFFFF - uint64(i)
		} else if i < 60000 {
			a[0] = uint64(i - 30000)
			for j := 1; j < 4; j++ {
				a[j] = 0
			}
		} else {
			rng.mkgf((*[4]uint64)(&a))
		}

		c.Half(&a)
		za := gfToBig((*[4]uint64)(&a), &p)
		zc := gfToBig((*[4]uint64)(&c), &p)
		var zd big.Int
		zd.Set(&za)
		if zd.Bit(0) != 0 {
			zd.Add(&zd, &p)
		}
		zd.Rsh(&zd, 1)
		if zc.Cmp(&zd) != 0 {
			t.Fatalf("ERR half:\na = %s\nc = %s\n", gfToString((*[4]uint64)(&a)), gfToString((*[4]uint64)(&c)))
		}

		for n := 1; n <= 15; n++ {
			c.Lsh(&a, uint(n))
			zc := gfToBig((*[4]uint64)(&c), &p)
			zd.Lsh(&za, uint(n)).Mod(&zd, &p)
			if zc.Cmp(&zd) != 0 {
				t.Fatalf("ERR lsh:\na = %s\nn = %d\nc = %s\n", gfToString((*[4]uint64)(&a)), n, gfToString((*[4]uint64)(&c)))
			}
		}
	}
}

func TestGF255sInv(t *testing.T) {
	var rng prng
	rng.init("test inv GF255s")
	var a, c, d GF255s
	for i := 0; i < 100000; i++ {
		if i < 100 {
			a[0] = -uint64(3957 + i)
			for j := 1; j < 3; j++ {
				a[j] = 0xFFFFFFFFFFFFFFFF
			}
			a[3] = 0x7FFFFFFFFFFFFFFF
		} else {
			rng.mkgf((*[4]uint64)(&a))
		}

		c.Inv(&a)
		d.Mul(&a, &c)
		var good uint64
		if i == 0 {
			good = d.IsZero()
		} else {
			good = d.Eq(&GF255s_ONE)
		}
		if good != 1 {
			t.Fatalf("ERR inv:\na = %s\nc = %s\nd = %s\n", gfToString((*[4]uint64)(&a)), gfToString((*[4]uint64)(&c)), gfToString((*[4]uint64)(&d)))
		}
	}
}

func TestGF255sLegendre(t *testing.T) {
	var rng prng
	rng.init("test Legendre GF255s")
	var a GF255s
	for i := 0; i < 100000; i++ {
		// Since 2^255 - 3957 = 3 mod 8, -2 is a square but -1 is not.
		if i < 100 {
			a[0] = -uint64(3957 + 2*i*i)
			for j := 1; j < 3; j++ {
				a[j] = 0xFFFFFFFFFFFFFFFF
			}
			a[3] = 0x7FFFFFFFFFFFFFFF
		} else {
			rng.mkgf((*[4]uint64)(&a))
			a.Sqr(&a)
		}

		ls := a.Legendre()
		if i == 0 {
			if ls != 0 {
				t.Fatalf("ERR Legendre: returned %d for 0\n", ls)
			}
			continue
		}
		if ls != 1 {
			t.Fatalf("ERR Legendre:\na = %s\nls = %d\n", gfToString((*[4]uint64)(&a)), ls)
		}
		a.Neg(&a)
		ls = a.Legendre()
		if ls != 0xFFFFFFFFFFFFFFFF {
			t.Fatalf("ERR Legendre:\na = %s\nls = %d\n", gfToString((*[4]uint64)(&a)), ls)
		}
	}
}

func TestGF255sSqrt(t *testing.T) {
	var rng prng
	rng.init("test Sqrt GF255s")
	var a, c, d GF255s
	for i := 0; i < 100000; i++ {
		// Since 2^255 - 3957 = 3 mod 8, -2 is a square but -1 is not.
		if i < 100 {
			a[0] = -uint64(3957 + 2*i*i)
			for j := 1; j < 3; j++ {
				a[j] = 0xFFFFFFFFFFFFFFFF
			}
			a[3] = 0x7FFFFFFFFFFFFFFF
		} else {
			rng.mkgf((*[4]uint64)(&a))
			a.Sqr(&a)
		}

		if c.Sqrt(&a) != 1 {
			t.Fatalf("ERR sqrt: failed on QR\na = %s\n", gfToString((*[4]uint64)(&a)))
		}
		if i == 0 {
			continue
		}
		d.Sqr(&c)
		if d.Eq(&a) != 1 {
			t.Fatalf("ERR sqrt: not a square root:\na = %s\nc = %s\n", gfToString((*[4]uint64)(&a)), gfToString((*[4]uint64)(&c)))
		}
		var ec [32]byte
		c.Encode(ec[:])
		if (ec[0] & 1) != 0 {
			t.Fatalf("ERR sqrt: wrong square root sign\na = %s\nc = %s\n", gfToString((*[4]uint64)(&a)), gfToString((*[4]uint64)(&c)))
		}

		a.Neg(&a)
		if c.Sqrt(&a) != 0 {
			t.Fatalf("ERR sqrt: non-QR was not detected\na = %s\nc = %s\n", gfToString((*[4]uint64)(&a)), gfToString((*[4]uint64)(&c)))
		}
		if c.IsZero() != 1 {
			t.Fatalf("ERR sqrt: non-QR was not normalized to zero\n")
		}
	}
}

func BenchmarkGF255sSqrt(b *testing.B) {
	var x GF255s
	x[0] = 0x875F15FEE9E2BCD3
	x[1] = 0xAF045BDB0DB8D55B
	x[2] = 0xFCECE0C9FAEA5119
	x[3] = 0x5CE56163CA711773
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		x.Sqrt(&x)
	}
}

func BenchmarkGF255sInv(b *testing.B) {
	var x GF255s
	x[0] = 0x875F15FEE9E2BCD3
	x[1] = 0xAF045BDB0DB8D55B
	x[2] = 0xFCECE0C9FAEA5119
	x[3] = 0x5CE56163CA711773
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		x.Inv(&x)
	}
}

func BenchmarkGF255sLegendre(b *testing.B) {
	var x GF255s
	x[0] = 0x875F15FEE9E2BCD3
	x[1] = 0xAF045BDB0DB8D55B
	x[2] = 0xFCECE0C9FAEA5119
	x[3] = 0x5CE56163CA711773
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		x.Legendre()
	}
}
