package do255s

import (
	"fmt"
	"math/big"
	"testing"
)

func TestDo255sScalar(t *testing.T) {
	var rng prng
	rng.init("test scalars do255s")
	var r, r0, tt big.Int
	tt.SetUint64(0xDCF2AC65396152C7)
	r0.SetUint64(0x2ACF567A912B7F03).Lsh(&r0, 64).Add(&r0, &tt)
	r.SetUint64(1).Lsh(&r, 254).Add(&r, &r0)
	r_enc := r.Bytes()
	for j := 0; (j + j) < len(r_enc); j++ {
		t := r_enc[j]
		r_enc[j] = r_enc[len(r_enc)-1-j]
		r_enc[len(r_enc)-1-j] = t
	}

	var a, b, c Do255sScalar
	for i := 0; i < 1000; i++ {
		var bb, bb2 [70]byte
		rng.generate(bb[:])
		a.DecodeReduce(bb[:])
		za := int256ToBigMod((*[4]uint64)(&a), &r)
		for j := 0; j < len(bb); j++ {
			bb2[len(bb)-1-j] = bb[j]
		}
		var zt big.Int
		zt.SetBytes(bb2[:]).Mod(&zt, &r)
		if za.Cmp(&zt) != 0 {
			s := "0x"
			for j := 0; j < len(bb2); j++ {
				s += fmt.Sprintf("%02X", bb2[j])
			}
			t.Fatalf("Reduction failed:\nsrc = %s\ndst = %s\n", s, int256ToString((*[4]uint64)(&a)))
		}

		rng.mk256((*[4]uint64)(&a))
		rng.mk256((*[4]uint64)(&b))
		za = int256ToBigMod((*[4]uint64)(&a), &r)
		zb := int256ToBigMod((*[4]uint64)(&b), &r)

		c.Add(&a, &b)
		zc := int256ToBigMod((*[4]uint64)(&c), &r)
		zt.Add(&za, &zb).Mod(&zt, &r)
		if zc.Cmp(&zt) != 0 {
			t.Fatalf("Wrong add:\na = %s\nb = %s\nc = %s\n", int256ToString((*[4]uint64)(&a)), int256ToString((*[4]uint64)(&b)), int256ToString((*[4]uint64)(&c)))
		}

		c.Sub(&a, &b)
		zc = int256ToBigMod((*[4]uint64)(&c), &r)
		zt.Sub(&za, &zb).Mod(&zt, &r)
		if zc.Cmp(&zt) != 0 {
			t.Fatalf("Wrong sub:\na = %s\nb = %s\nc = %s\n", int256ToString((*[4]uint64)(&a)), int256ToString((*[4]uint64)(&b)), int256ToString((*[4]uint64)(&c)))
		}

		c.Mul(&a, &b)
		zc = int256ToBigMod((*[4]uint64)(&c), &r)
		zt.Mul(&za, &zb).Mod(&zt, &r)
		if zc.Cmp(&zt) != 0 {
			t.Fatalf("Wrong mul:\na = %s\nb = %s\nc = %s\n", int256ToString((*[4]uint64)(&a)), int256ToString((*[4]uint64)(&b)), int256ToString((*[4]uint64)(&c)))
		}

		var aa [32]byte
		copy(aa[:], r_enc[:])
		w := int(aa[0]) + (int(aa[1]) << 8)
		w += i - 500
		aa[0] = byte(w)
		aa[1] = byte(w >> 8)
		ds := a.Decode(aa[:])
		if i < 500 {
			if ds != 1 {
				s := "0x"
				for j := 0; j < len(aa); j++ {
					s += fmt.Sprintf("%02X", aa[31-j])
				}
				t.Fatalf("In-range scalar decoding failed:\nsrc = %s\nout = %d\n", s, ds)
			}
			zb.SetUint64(uint64(500 - i))
			zt.Sub(&r, &zb)
		} else {
			if ds != -1 {
				s := "0x"
				for j := 0; j < len(aa); j++ {
					s += fmt.Sprintf("%02X", aa[31-j])
				}
				t.Fatalf("Out of range scalar was decoded:\nsrc = %s\ndst = %s\nout = %d\n", s, int256ToString((*[4]uint64)(&a)), ds)
			}
			zt.SetUint64(0)
		}
		za = int256ToBigMod((*[4]uint64)(&a), &r)
		if za.Cmp(&zt) != 0 {
			s := "0x"
			for j := 0; j < len(aa); j++ {
				s += fmt.Sprintf("%02X", aa[31-j])
			}
			t.Fatalf("Wrong decoded value:\nsrc = %s\ndst = %s\n", s, int256ToString((*[4]uint64)(&a)))
		}
	}
}

func TestDo255sLagrange(t *testing.T) {
	var rng prng
	rng.init("test Lagrange do255s")
	for i := 0; i < 10000; i++ {
		var bb [32]byte
		rng.generate(bb[:])
		var k Do255sScalar
		k.DecodeReduce(bb[:])
		var c0, c1 [2]uint64
		negc0, negc1 := k.ReduceBasisVartime(&c0, &c1)
		var d0, d1 Do255sScalar
		copy(d0[:], c0[:])
		copy(d1[:], c1[:])
		if negc0 {
			d0.Neg(&d0)
		}
		if negc1 {
			d1.Neg(&d1)
		}
		var tt Do255sScalar
		tt.Mul(&k, &d1)
		if tt.Equal(&d0) != 1 {
			t.Fatalf("Lattice basis reduction failed:\nk = %s\nc0 = %s\nc1 = %s\n", int256ToString((*[4]uint64)(&k)), int256ToString((*[4]uint64)(&d0)), int256ToString((*[4]uint64)(&d1)))
		}
	}
}
