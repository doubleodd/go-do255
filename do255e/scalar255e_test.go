package do255e

import (
	"fmt"
	"math/big"
	"testing"
)

func TestDo255eScalar(t *testing.T) {
	var rng prng
	rng.init("test scalars do255e")
	var r, r0, tt big.Int
	tt.SetUint64(0xE0AD37518B27BADB)
	r0.SetUint64(0x62F36CF0ABF873AC).Lsh(&r0, 64).Add(&r0, &tt)
	r.SetUint64(1).Lsh(&r, 254).Sub(&r, &r0)
	r_enc := r.Bytes()
	for j := 0; (j + j) < len(r_enc); j++ {
		t := r_enc[j]
		r_enc[j] = r_enc[len(r_enc)-1-j]
		r_enc[len(r_enc)-1-j] = t
	}

	var a, b, c Do255eScalar
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
