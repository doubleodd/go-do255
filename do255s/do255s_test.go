package do255s

import (
	"bytes"
	"encoding/hex"
	"testing"
)

func TestDo255sDecode(t *testing.T) {
	for i := 0; i < len(kat_DO255S_DECODE_OK); i++ {
		bb, _ := hex.DecodeString(kat_DO255S_DECODE_OK[i])
		rc := Do255sCheckPoint(bb)
		var P Do255sPoint
		rd := P.Decode(bb)
		if i == 0 {
			if rc != 0 || rd != 0 {
				t.Fatalf("Neutral not decoded: %d, %d\n", rc, rd)
			}
			if P.IsNeutral() == 0 {
				t.Fatalf("Neutral not decoded as neutral\n")
			}
		} else {
			if rc != 1 || rd != 1 {
				t.Fatalf("Valid non-neutral not decoded: %d, %d\n", rc, rd)
			}
		}
		e2 := P.Bytes()
		if !bytes.Equal(bb, e2[:]) {
			t.Fatalf("Point not reencoded properly:\nsrc = %s\ndst = %s\n", hex.EncodeToString(bb), hex.EncodeToString(e2[:]))
		}
	}

	bzz := do255sNeutral.Bytes()

	for i := 0; i < len(kat_DO255S_DECODE_BAD); i++ {
		bb, _ := hex.DecodeString(kat_DO255S_DECODE_BAD[i])
		rc := Do255sCheckPoint(bb)
		if rc != -1 {
			t.Fatalf("Invalid point reported as decodable (1):\nsrc = %s\n", hex.EncodeToString(bb))
		}
		var P Do255sPoint
		P.Generator()
		rd := P.Decode(bb)
		if rd != -1 {
			t.Fatalf("Invalid point reported as decodable (2):\nsrc = %s\n", hex.EncodeToString(bb))
		}
		e2 := P.Bytes()
		if !bytes.Equal(e2[:], bzz[:]) {
			t.Fatalf("Invalid point not normalized to neutral:\ndst = %s\n", hex.EncodeToString(e2[:]))
		}
	}
}

func TestDo255sMapBytes(t *testing.T) {
	for i := 0; i < len(kat_DO255S_POINT_MAP); i += 2 {
		bb1, _ := hex.DecodeString(kat_DO255S_POINT_MAP[i])
		bb2, _ := hex.DecodeString(kat_DO255S_POINT_MAP[i+1])
		var P Do255sPoint
		P.MapBytes(bb1)
		bb3 := P.Bytes()
		if !bytes.Equal(bb2, bb3[:]) {
			t.Fatalf("Mapping failed:\nexp = %s\ngot = %s\n", hex.EncodeToString(bb2), hex.EncodeToString(bb3[:]))
		}
	}
}

func TestDo255sPointAdd(t *testing.T) {
	for i := 0; i < len(kat_DO255S_POINT_ADD); i += 6 {
		var P1, P2, P3, P4, P5, P6 Do255sPoint
		bb1, _ := hex.DecodeString(kat_DO255S_POINT_ADD[i])
		bb2, _ := hex.DecodeString(kat_DO255S_POINT_ADD[i+1])
		bb3, _ := hex.DecodeString(kat_DO255S_POINT_ADD[i+2])
		bb4, _ := hex.DecodeString(kat_DO255S_POINT_ADD[i+3])
		bb5, _ := hex.DecodeString(kat_DO255S_POINT_ADD[i+4])
		bb6, _ := hex.DecodeString(kat_DO255S_POINT_ADD[i+5])
		P1.Decode(bb1)
		P2.Decode(bb2)
		P3.Decode(bb3)
		P4.Decode(bb4)
		P5.Decode(bb5)
		P6.Decode(bb6)

		if P1.Equal(&P2) != 0 || P1.Equal(&P3) != 0 || P1.Equal(&P4) != 0 || P1.Equal(&P5) != 0 || P1.Equal(&P6) != 0 || P2.Equal(&P3) != 0 || P2.Equal(&P4) != 0 || P2.Equal(&P5) != 0 || P2.Equal(&P6) != 0 || P3.Equal(&P4) != 0 || P3.Equal(&P5) != 0 || P3.Equal(&P6) != 0 || P4.Equal(&P5) != 0 || P4.Equal(&P6) != 0 || P5.Equal(&P6) != 0 {
			t.Fatalf("Equal() malfunction (1)\n")
		}
		if !(P1.Equal(&P1) == 1 && P2.Equal(&P2) == 1 && P3.Equal(&P3) == 1 && P4.Equal(&P4) == 1 && P5.Equal(&P5) == 1 && P6.Equal(&P6) == 1) {
			t.Fatalf("Equal() malfunction (2)\n")
		}
		if P1.IsNeutral() != 0 || P2.IsNeutral() != 0 || P3.IsNeutral() != 0 || P4.IsNeutral() != 0 || P5.IsNeutral() != 0 || P6.IsNeutral() != 0 || do255sNeutral.IsNeutral() == 0 {
			t.Fatalf("IsNeutral() malfunction\n")
		}

		// P3 = P1 + P2
		var Q3 Do255sPoint
		Q3.Generator()
		Q3.Add(&P1, &P2)
		if Q3.Equal(&P3) == 0 {
			t.Fatalf("Addition failed (1):\nexp = %s\ngot = %s\n", hex.EncodeToString(P3.Encode(nil)), hex.EncodeToString(Q3.Encode(nil)))
		}
		Q3.Generator()
		Q3.Add(&P2, &P1)
		if Q3.Equal(&P3) == 0 {
			t.Fatalf("Addition failed (2):\nexp = %s\ngot = %s\n", hex.EncodeToString(P3.Encode(nil)), hex.EncodeToString(Q3.Encode(nil)))
		}

		// P4 = 2*P1
		var Q4 Do255sPoint
		Q4.Generator()
		Q4.Add(&P1, &P1)
		if Q4.Equal(&P4) == 0 {
			t.Fatalf("Addition failed (3):\nexp = %s\ngot = %s\n", hex.EncodeToString(P4.Encode(nil)), hex.EncodeToString(Q4.Encode(nil)))
		}
		Q4.Generator()
		Q4.Double(&P1)
		if Q4.Equal(&P4) == 0 {
			t.Fatalf("Addition failed (4):\nexp = %s\ngot = %s\n", hex.EncodeToString(P4.Encode(nil)), hex.EncodeToString(Q4.Encode(nil)))
		}

		// P5 = P4 + P2 = P3 + P1
		var Q5 Do255sPoint
		Q5.Generator()
		Q5.Add(&Q4, &P2)
		if Q5.Equal(&P5) == 0 {
			t.Fatalf("Addition failed (5):\nexp = %s\ngot = %s\n", hex.EncodeToString(P5.Encode(nil)), hex.EncodeToString(Q5.Encode(nil)))
		}
		Q5.Generator()
		Q5.Add(&Q3, &P1)
		if Q5.Equal(&P5) == 0 {
			t.Fatalf("Addition failed (6):\nexp = %s\ngot = %s\n", hex.EncodeToString(P5.Encode(nil)), hex.EncodeToString(Q5.Encode(nil)))
		}

		// P6 = P5 + P2 = P4 + 2*P2
		var Q6 Do255sPoint
		Q6.Generator()
		Q6.Add(&Q5, &P2)
		if Q6.Equal(&P6) == 0 {
			t.Fatalf("Addition failed (7):\nexp = %s\ngot = %s\n", hex.EncodeToString(P6.Encode(nil)), hex.EncodeToString(Q6.Encode(nil)))
		}
		Q6.Generator()
		Q6.Double(&P2)
		Q6.Add(&Q6, &Q4)
		if Q6.Equal(&P6) == 0 {
			t.Fatalf("Addition failed (8):\nexp = %s\ngot = %s\n", hex.EncodeToString(P6.Encode(nil)), hex.EncodeToString(Q6.Encode(nil)))
		}
		Q6.Generator()
		Q6.Add(&P2, &Q4)
		Q6.Add(&Q6, &P2)
		if Q6.Equal(&P6) == 0 {
			t.Fatalf("Addition failed (9):\nexp = %s\ngot = %s\n", hex.EncodeToString(P6.Encode(nil)), hex.EncodeToString(Q6.Encode(nil)))
		}

		// Adding the neutral should not change the point.
		var Q7 Do255sPoint
		Q7.Generator()
		Q7.Add(&Q6, &do255sNeutral)
		if Q7.Equal(&Q6) == 0 {
			t.Fatalf("Addition of neutral failed (1):\nexp = %s\ngot = %s\n", hex.EncodeToString(Q6.Encode(nil)), hex.EncodeToString(Q7.Encode(nil)))
		}
		Q7.Generator()
		Q7.Add(&do255sNeutral, &Q6)
		if Q7.Equal(&Q6) == 0 {
			t.Fatalf("Addition of neutral failed (2):\nexp = %s\ngot = %s\n", hex.EncodeToString(Q6.Encode(nil)), hex.EncodeToString(Q7.Encode(nil)))
		}

		// Testing negation.
		var Q8, Q9 Do255sPoint
		Q8.Generator()
		Q9.Generator()
		Q8.Neg(&Q6)
		Q9.Add(&Q8, &Q6)
		if Q9.IsNeutral() == 0 {
			t.Fatalf("Addition of negation failed:\ngot = %s\n", hex.EncodeToString(Q9.Encode(nil)))
		}

		// Testing sequences of doublings.
		var Q10, Q11 Do255sPoint
		for j := 0; j < 10; j++ {
			Q10.Generator()
			Q10.DoubleX(&Q6, uint(j))
			Q11.Set(&Q6)
			for k := 0; k < j; k++ {
				Q11.Double(&Q11)
			}
			if Q10.Equal(&Q11) == 0 {
				t.Fatalf("Successive doublings failed (n = %d):\nexp = %s\ngot = %s\n", j, hex.EncodeToString(Q11.Encode(nil)), hex.EncodeToString(Q10.Encode(nil)))
			}
		}
	}
}

func TestDo255sPointMul(t *testing.T) {
	for i := 0; i < len(kat_DO255S_POINT_MUL); i += 3 {
		var P1, P2, P3 Do255sPoint
		var n Do255sScalar
		bb1, _ := hex.DecodeString(kat_DO255S_POINT_MUL[i])
		bb2, _ := hex.DecodeString(kat_DO255S_POINT_MUL[i+1])
		bb3, _ := hex.DecodeString(kat_DO255S_POINT_MUL[i+2])
		P1.Decode(bb1)
		n.DecodeReduce(bb2)
		P2.Decode(bb3)
		P3.Mul(&P1, &n)
		if P2.Equal(&P3) == 0 {
			t.Fatalf("Wrong point multiplication result:\nexp = %s\ngot = %s\n", hex.EncodeToString(P2.Encode(nil)), hex.EncodeToString(P3.Encode(nil)))
		}
	}

	var rng prng
	rng.init("test mulgen do255s")
	for i := 0; i < 1000; i++ {
		var n Do255sScalar
		if i == 0 {
			n[0] = 0
			n[1] = 0
			n[2] = 0
			n[3] = 0
		} else {
			rng.mk256((*[4]uint64)(&n))
		}
		var P1, P2 Do255sPoint
		P1.MulGen(&n)
		P2.Generator().Mul(&P2, &n)
		if P1.Equal(&P2) == 0 {
			t.Fatalf("Wrong point mulgen result:\nexp = %s\ngot = %s\n", hex.EncodeToString(P1.Encode(nil)), hex.EncodeToString(P2.Encode(nil)))
		}
	}

	bb, _ := hex.DecodeString(kat_DO255S_MC_POINT_MUL[0])
	var P Do255sPoint
	P.Decode(bb)
	for i := 1; i < len(kat_DO255S_MC_POINT_MUL); i++ {
		for j := 0; j < 1000; j++ {
			var n Do255sScalar
			n.DecodeReduce(bb)
			P.Mul(&P, &n)
			P.Encode(bb[:0])
		}
		str := hex.EncodeToString(bb)
		exp := kat_DO255S_MC_POINT_MUL[i]
		if str != exp {
			t.Fatalf("Wrong MC mul result:\nexp = %s\ngot = %s\n", exp, str)
		}
	}
}

func TestDo255sVerifyHelper(t *testing.T) {
	var rng prng
	rng.init("test verify helper do255s")
	for i := 0; i < 1000; i++ {
		var n, k0, k1 Do255sScalar

		rng.mk256((*[4]uint64)(&n))
		rng.mk256((*[4]uint64)(&k0))
		rng.mk256((*[4]uint64)(&k1))
		if i < 16 {
			// Special cases.

			// P == neutral
			if (i & 1) == 0 {
				n[0] = 0
				n[1] = 0
				n[2] = 0
				n[3] = 0
			}

			// k0 == 0
			if (i & 2) == 0 {
				k0[0] = 0
				k0[1] = 0
				k0[2] = 0
				k0[3] = 0
			}

			// k1 == 0
			if (i & 4) == 0 {
				k1[0] = 0
				k1[1] = 0
				k1[2] = 0
				k1[3] = 0
			}

			// R == neutral
			if (i & 8) == 0 {
				k0.Mul(&n, &k1)
				k0.Neg(&k0)
			}
		}
		var P, R, M Do255sPoint
		P.MulGen(&n)
		R.MulGen(&k0)
		M.Mul(&P, &k1)
		R.Add(&R, &M)
		if P.VerifyHelper(&k0, &k1, R.Encode(nil)) != 1 {
			t.Fatalf("Verify helper (ct) failed\n")
		}
		if !P.VerifyHelperVartime(&k0, &k1, R.Encode(nil)) {
			t.Fatalf("Verify helper (vartime) failed\n")
		}
		M.Generator()
		R.Add(&R, &M)
		if P.VerifyHelper(&k0, &k1, R.Encode(nil)) != 0 {
			t.Fatalf("Verify helper (ct) should have failed\n")
		}
		if P.VerifyHelperVartime(&k0, &k1, R.Encode(nil)) {
			t.Fatalf("Verify helper (vartime) should have failed\n")
		}
	}
}

func BenchmarkMul255s(b *testing.B) {
	var P Do255sPoint
	bb, _ := hex.DecodeString("c591e486e9f0f87589ba486b51629c358cddf4bd8c50db5c53d5631a3f785074")
	P.Decode(bb)
	var s Do255sScalar
	bb, _ = hex.DecodeString("81d3066b23528ff3728a200ae411d815ba612b45c2556fa1000cb665de692a73")
	s.DecodeReduce(bb)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		P.Mul(&P, &s)
	}
}

func BenchmarkMulGen255s(b *testing.B) {
	var P Do255sPoint
	var s Do255sScalar
	bb, _ := hex.DecodeString("81d3066b23528ff3728a200ae411d815ba612b45c2556fa1000cb665de692a73")
	s.DecodeReduce(bb)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		P.MulGen(&s)
	}
}

func BenchmarkVerify255s(b *testing.B) {
	encR := do255sGenerator.Bytes()
	var rng prng
	rng.init("bench verify do255s")
	var k0, k1 Do255sScalar
	rng.mk256((*[4]uint64)(&k0))
	rng.mk256((*[4]uint64)(&k1))
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		do255sGenerator.VerifyHelper(&k0, &k1, encR[:])
	}
}

func BenchmarkVerifyVartime255s(b *testing.B) {
	encR := do255sGenerator.Bytes()
	var rng prng
	rng.init("bench verify do255s")
	var k0, k1 [100]Do255sScalar
	for i := 0; i < len(k0); i++ {
		rng.mk256((*[4]uint64)(&k0[i]))
		rng.mk256((*[4]uint64)(&k1[i]))
	}
	b.ResetTimer()
	j := 0
	for i := 0; i < b.N; i++ {
		do255sGenerator.VerifyHelperVartime(&k0[j], &k1[j], encR[:])
		j++
		if j >= len(k0) {
			j = 0
		}
	}
}

func BenchmarkLagrange255s(b *testing.B) {
	var rng prng
	rng.init("bench Lagrange do255s")
	var k [100]Do255sScalar
	for i := 0; i < len(k); i++ {
		rng.mk256((*[4]uint64)(&k[i]))
	}
	b.ResetTimer()
	j := 0
	for i := 0; i < b.N; i++ {
		var c0, c1 [2]uint64
		k[j].ReduceBasisVartime(&c0, &c1)
		j++
		if j >= len(k) {
			j = 0
		}
	}
}

// Strings that decode to valid points.
var kat_DO255S_DECODE_OK = []string{
	"0000000000000000000000000000000000000000000000000000000000000000",
	"c591e486e9f0f87589ba486b51629c358cddf4bd8c50db5c53d5631a3f785074",
	"3280b9a23784102fb56338df91bc417efa9d04a4718a36a55552d0e55a8d6a25",
	"5357ab7306780b1b2fac2f107ff1effe3f3c39ad3c9fa4e77702b88059d52f25",
	"ed12d7281a7638d58234df32ec1c4985d1ae099c5f7ffe066df6693c5b242f04",
	"252f89df45487c2d67b066bd2bb5ecad1e8f6970a63ce456fe2ba430ea665a0d",
	"5cd741e897f34f6769f49e6d192e091b021394b7c29fcd0f92211dc473b22a23",
	"736f644265f9d4ab1d3520c004ef6f46c6ac48eb36ca216c0a3f35a4eb55d347",
	"b58cc302de47ce281bd63edfd6464ff56b0742d1b99c540b63e7427d823c277c",
	"54f9b832bc6ba95f300235f0593d38564b9041df76b27c76c4bbb170c5892404",
	"ded481afe5d5ee743be61c06f6306dc6885bb1241022d7feed4f7654eafb5c50",
	"468e63bd5efe7dd6b88855804c465c94a929c6051bb411dc5259e02ff6d3e072",
	"b3e8e5972e218178596832d80b946f89da310666fb3751160ea7fc8a6b2de163",
	"03b720f4d6b9e76485eeb8fc82c183679ee25d0d69d4511c91d742ea60391a74",
	"83ca018c52abee2354aed77f55d344dce48c9a30408e9531904999a7b28bec0c",
	"f57f1503088132b922c7888b742e82565e4826e01968c035fde90d11b456353c",
	"b32fdcf6a19cc515f6b819c431ba5d2f3b80a111ea3dc4adff749d802060d640",
	"e38a8462fb16cd79cebe97c582294e9019a801b21c28578c036d35dae7b6e11a",
	"096ef5f368322cce5a902cc9e316007352519a52c37c2ddf7e59b7f1697f0f00",
	"b89e8f053f295e983f03bddd2761faf59eba3c1f69071697a388da00b1457e3f",
	"6982b80b863d40f177f0c48402cb246ef85644354eb84003f84d9ffe4754592a",
}

// Strings that do not decode to valid points.
var kat_DO255S_DECODE_BAD = []string{
	// These values cannot be decoded (w is out of range).
	"3a92541c88f5ca5620985ec23400c1149969685ff403f9ff8fa63349a9f518a9",
	"7dcebc2ba5e16650eef4f69d4fb444786858d42576e2fcbb867938df7451e4c6",
	"f60268231297120bced5af37a3513caa0e2b5733a30cd04c3034084e5e0f5acc",
	"10f225adbaa42e20723ee481a584269c0b573603b65565ffa8f39d8c96237683",
	"02f8239691757d594850446ffb6248b889af05fb9580789fdb0d827faef2bbd9",
	"8cfcce45e6979c13b7a20129f49bf857dac14fb6dd6311237111c278747860af",
	"07321e9fa5affc3dae33da1af30ee0eb05db0fbbccee8bc990fb2cf0c44ec2d2",
	"65f3c2f7ee8c16f053c047d10b185507cea824a35f9b38ed9db54f3bd2616c83",
	"ff57953b5b8497d3266be9e27b05e974a8dec70452ddc421df115dbc0670b7a0",
	"a0741124d8c6cd1aa30528026bed28b9ee00096d57aa12048f0f7b67dd6e0c93",
	"d0745a09c13838989c83ec2ad350fdb94ea9517d20527354067ad1f43b37fef1",
	"36c2472d2a704844e64024dfcb1077fa715bf394994915143157f90337a3bdb8",
	"830271462ef4b2430012cb6b089a2d087d2b982e91ecb194d4997e643cad0487",
	"9632aa4a38a6dddd15d8f7877ba46ccfbe6b80f7254c008a973cd0f3ade8eaf5",
	"f3db83e51c406741d20710fe4aa5f04e06c0ed4794b21d6931946007c68fcbb9",
	"dc8c9ca8b34e132d1ac61694d8e457f296cd21f3d75f640ede46267f3b7c80bf",
	"922c234934c93ae72825d2167430d5939305fea8019d63c014527e7f2727a2c2",
	"ea6fbe2ec7177b20ad1b7bab4c61bdff817e9afb9dc9b7f7ce8ff4df2ae721cf",
	"4b0c025553ed615139140e7852682cd87fda491eec9ffafa5e515da9e71a6ad6",
	"dff93a41bef66c4725e63b20634a758c61a0fff9d92e5874fbed4d44cb317fc5",

	// These values cannot be decoded (w matches no point).
	"01e244ae27c9bf0406559c04ab7b25587f56ff0b79968ed05cee35ac7e527e31",
	"0904f4a4622446cd7f453e13968f22aaf1841fe17461e0998753e8a6c3f1db1b",
	"abf5bb06f321a68875a285c8ae8a1f10352c566a9a347c8758b385b05151c147",
	"e0102e5cb6cf5e1f4b9d027d147d7cfc161ae6441f33f55a0acc2737d3117229",
	"2d81ecbfec945c26df0d4cbe401f55fd550bea5ca5facf26f592e18f3e398957",
	"7d781c87ecd3ae16b13d72ced0bf9fad2ab00014e1c5035f95cb32bb8d295e60",
	"3cbba04f2ec5f32c2fc68183de992759d0b6600983ea5e836c190bbaddef1b5b",
	"f53d4f896d6cb6c0ac1f36939df609419c90c27b0533cee21d0893173be4a111",
	"b2ec464f0ce61e3015e68e2f8c8ca55da784b605694a740391f0b6d405756b63",
	"e8e297612bb58431e08168f6b27739461cc5d8e5b8e02edc7dcf0209440bdb62",
	"28c94ceb6c5b83b802cfbd4ae3e04455061269a78b09299c2767c90ed2aea714",
	"5dc5a2ee12aaaf093fa2c1fa845fa45431ad630760424ce24c4e15951fb69645",
	"e1434027d7502a74417c663183df9a570746c873e02d0ec2cbf7a33d5024fa33",
	"4120309222b712ffe3390b263c018eb3957154337578affef09b8dc918c57e14",
	"f045fcd44e2c4582bc929cac99b8fe62013a48b7a15a8c7102b283eb18ec8068",
	"548982ec9b6a9d396bf36c3d7df1fb24975cdb750b7c316f218bbc388aebc732",
	"3309ba8cad6523b7984f19342e5f63b83b83964a798a1388ecbf2673aa15fb40",
	"60db10a949d2e9ed81a4fe4b6b0d6e4454ac420f22e7dfcaf4503aefd7113542",
	"b93e2061042a1e65640bc8632e880a565e459ee4ac4755705cf9f0769cf58c57",
	"80493859a49bf1ebc8a3e7195157b2b98a0ed135d77959ef8a94a64b709f7628",
}

// Mapping of bytes to points.
var kat_DO255S_POINT_MAP = []string{
	// Map-to-curve test vectors for Do255s
	// Each group of two values is: input bytes, mapped point
	"0100000000000000000000000000000000000000000000000000000000000000",
	"0000000000000000000000000000000000000000000000000000000000000000",

	"ef8ddec9c3906319fb74f84036829ba34a1f51f7700b78c81bc91389b8adcb29",
	"8c09ce78bc83500190e7059310f5e8637170769d844233bc557e5e2497dcf820",

	"09aa856c6506ed9848570423dc8fbf4e916fe7f36ee895109893895a73415f76",
	"336f8f5872a29e0defd157585246037e0e349e6ef3bd1c990efdeda52d9d181a",

	"2615809de9b07f3c5cf0170ceccdab6843798f628743daab41b235ebd93f0916",
	"f801553749eaab1c17ab414754bdf0214ae534a6254cc3f74767b5c86732b472",

	"0cc3f41eef7051a5bd1b32ff447618375a6efef1adebb9331abf4bd7b172d672",
	"17a43022c6b25554d98355bb937ec9d3179c999c2082c96da20d9cdd44666b00",

	"d706eec3dc6f112f1fa8ee531dc6af144c03a7ef261c117f7d3ead67f4f6a835",
	"702cf2796ae47b5b0892c8eac02567c28fab77f50ae42a01167899aa08ac0646",

	"927fadf4f0357ee03d8f52503692c404c9b3ac05443990619852098ed1ffe741",
	"a6815504de2232cc30e4f68b7329e5824339d47429fcd9ac374c9f312a239108",

	"3ed2a27bb19a5d0e4015e310ea2c7d5be9363d2e3a699c50ff1073d89d984b3e",
	"d30a3eadf2a120cb196da94463c216e0a5b5c239892c4c21caa1995a2c7cb86e",

	"184cffce57049542aab9a816dc806b4419a6e98a538685c627bc1cb86acbebad",
	"74ea8109db1b8c0829bd1fcae8ccdbae2df1683909992c2bbeff1b9338c2c47a",

	"299039d72adb0f9afb61e59c3d7f0fdb102cf152e3cac224922ca5be9087b12e",
	"350b451b72b66f4d0995d468cc51716a36aaa45f4280c6e378de3d277e2e5422",

	"27411f5d3c65de7d6c3a34bfe03985c3b9ab86d67a33af1851f025c4dc2e45c8",
	"8272b971dbdad55b02c5ba6bbffd02c41a9979a7a6373514c10d4a861a81404e",

	"27caf691b8003ef948ffbbf0520bf129d47024e0576ac3f28e73911b75afed33",
	"c88444c61ca57d298a2c40e59a26d2fe7261e12ea8e04a2b19e5f10b5c61a26f",

	"aa09eba29f536707455c34328d74b417f5100f3cb499cb360fc7cc368153d4a1",
	"7a70f2f0076852e9eb0331333021ab2d41eab2dd72ec2e143f25a77796a4f647",

	"8c8f7c49e620854b260b5f5ad5e3d221580b24aaac46377b13024f7d3198d6f2",
	"f406e61169965432bdc314b8dfc35e72ad886c76c72ccdcca8abd227709c2c35",

	"6fddae9d444eaf137aa2831b1713c93dbade201ea2a811a608ae0cfbcacffb6a",
	"a829e7440a42b80665eea59eb1cf5480e348be004a8b4b4b07e493122b0e5513",

	"0193677a09da4f1d819ef8418d15e0cae8eb9a89aafb96f13a262a3e6d950f11",
	"dd027f5ce90e7afa06c890fad453f76005e9cdf827342e91ca7042b3e722120b",

	"1d947f967ddc8fe741eb3533f7f25f78c541f53d30fcb7d6a5b5d8f269020030",
	"8c12d16626cde555cb922dccaa8ff8feed931df1bcce5f4f85a240ae2dd26821",

	"32049bcf13de6ac721a02a4bdc0737abc4c3c0f28b9a6123ccb33a9e7647d136",
	"70b18ef2e9c0a33a8c6a79f6823e405cf8451c293b5999eb7a7c1a67cdfba460",

	"69cf1efafc674d1a41ba28a16b58759bb42d30bc77e8ec5abe19eae43f150085",
	"1234cd73cffb6350e8b816d349d79064ae1d5c01869a08d11384070b53027f75",

	"3453d8c0de0389d588ae9791b6105f17c7407fabd3e2a665194e591600392a64",
	"4fe4ff767a01552863bc22363b2431f19c3321e2bb0323a2c5bca6dc94be560d",

	"530195eb8bf148eecdd97ecebf7580ebc685a2a7a1a7aa54bfa27cf87a203ec4",
	"c26b2c32991126301254bed67ca20674dd8f42151a7b80898eab4b2a03aa712b",

	"5f915e83780fdf6cc484b2d890c2a5f053f2eb845722882c074a271ce3bde224",
	"f81dac9fa4565539ab631a282c86d2cd49b9f18fd18357cf6ec77c93330c6a7b",

	"c412fee4f10c9b8950077e31438f9bfdc644a28cc9f3cfd83ef480a82569a41f",
	"178e85745e49482e6ca1fa345ec910f1e2f8724ea03c5093fa47b6491d54c174",

	"3121bd886339c8004216c982f690b50f4e633d60c9988b669d6d2c1fccab4ddf",
	"30e95dca620183bb105b09b13b8f9e2ed3ad5856824975e9728019e18d57fc1b",

	"14b087ed68b14798d2d1bc4ea65ff28305862b94ca6601d9bb8f807043974e67",
	"38f3860c889afc82974d2aa0bf33afdba9bc4b40cbfe59e331ee94e4dba9167e",

	"b26fcbc480eb4d531fcdc50bed5d292bfca679ffc358655e0135b169a351ddd6",
	"92f06ee33ee6c62b24b8546b9e235bf21dc924e9e2709ce1e912271c3156a237",

	"39f63d9002b06a64fbd6e22e541f6374e6731fc831c1f1c2e46c5391c959b8c6",
	"32c7847e59e8e0c39d5bd9efe2cd656b2d7654c733316dc7e5fa4dd1fccbd314",

	"7b7c740f9fffb909e45fb3d5ef8b705f009590f7cd55a96047a30407f2553e5e",
	"bb9eef67a061a575c3b71a3e29993eb8c523fe9ae11a707a30c048f0fe2d8007",

	"0932c9c6acde82775d0bb254731bb096833362990f68528c9b14f9d076e93c02",
	"906894b1b036386b4ddf2a8a58d45f308c9321109fc12a34341e50f0e80df216",

	"ba4ad4152fd792db8ea1575c87494d23036ee200ceef8d8c1bc99f6df38ab514",
	"dd9ab7a2c8778a7cd89769a71246551ed70f7606fdecd717565184740c1d4f03",

	"5896d16b48209ca651ba8f686688d091072805b245f136871b43bb8167fca960",
	"ed209723b4d1d06ffec847ae89708474dc4fc133fe793614347fab61e67b2c74",

	"ce5fdf706ea9c16e47d4cb9fab4ba6173858acfefec7096669aab4c40b928346",
	"39fbb562fefcbeaa069d6ba0c975e99bc36ea3f017c24b3bcf36fa6ab84e0c17",

	"a5f64a58f56dfb212f63cbde27a246e96f05dcc60687a8eb51bb84200e3be5e9",
	"ff8094c66be12b8008b3360fd7ca828798f497797f879c4b72fb1575e4b85b69",

	"70b2451d5465130a6cf7998d3b752dc8ccd7547598149435e77573a86bff1ea5",
	"9969658848ac05f24f3387ef5f6ffb7e4aa3339171537143a23984205341cb3a",

	"b87e45cd4369dba7cd7b877a20325b027ba87c90d23939844c0c9a194ee39d3a",
	"2e6b5b348687087b30c13bb9ed59fd9e324490ea09b66d1a55a7601a8fb25b77",

	"163eec28be72a61876479a7e5838e438c5d053ac09e4a9b6c4f5e91ece2a9244",
	"06aa16b19bf14ef4607f31f7a782ec38da885f912b24fa19ffc834a16861d805",

	"9e3b7632ae68362ed9c716bd1e52c040315531a4f75c67ce3e3c8253d36c98e5",
	"f42ec2b036a422009273ca5d4edaa6e1a910cf288f6e1435856769feddc99564",

	"f1ec32a37b6a82cbe322e049687af65ab32383e5b4ff1a40b40876c4c9876ad8",
	"ce2661546a9ffe2c7bcfb6782290870b68544b2624a52031c179be82d530497b",

	"8f6b4c01d3090814d79008367f0aec97e27a6aecf751eb2cd678b744018649ab",
	"46bf5c2ad04470fd09f209af2f65f1c8063bf58ae55fbb06758d8c150a78c105",

	"395c4cd8ca231c8a20169c456908ce8b88905e740f389837b33e693ce6d9e1f1",
	"227fcab45ba482a85ac0b4f30b0e1c4fbca733aae08022a017a814cb20601c2d",
}

// Point addition tests. Each group of 6 strings is P1..P6, with:
//   P3 = P1 + P2
//   P4 = 2*P1
//   P5 = P4 + P2 = P3 + P1
//   P6 = P5 + P2 = P4 + 2*P2
var kat_DO255S_POINT_ADD = []string{
	"a640cf808233d9ba7d04df6eb7f79bfca282e16788a06268458a23a67c786915",
	"6227d95a323ef8dcec8f963ee75b0c3610aaa85eb141eaf0c9351c3bb7614357",
	"3e662bf0ac88715ece7ae60455d8e3cd62ad6015eb15f3e9bfb1affa13ff0e2a",
	"d96d674f7963b382007d7969307150c4cc95b1063e4b073a0183a7519ab67752",
	"61e387b94adf41923549c7178bb6b8c924423eb96d6d0fee123203400dd0a44c",
	"2b8e33dd62d6f3572ea7c53770089717c25bac5279bf2f43c0fe2f1b0932f27f",

	"1a3a1178a2bd338bdf7ca49fa59ecb9a652bc96c6a48577aa881b51adf4d282e",
	"6040ce129c4c483fb92cf164a63a6a2d6e4848999391d4b94147fb8d9fe59551",
	"4dad0c9677b1fd43be0e79d2e7b6d90e1a5c55b925bbdfb50bd52352bbcfe11a",
	"dd11bd55b0b4079954655b36ffb20a96f1f4f163904daffe14b8e0efd750727a",
	"a250018185e4123367020d813e5b728b55c73c5e1ea7ac2f0116d54ba701ee7d",
	"32b7903197d16ab2fde54a282910aae877b858918c5627e9c4126e4e5dc8c479",

	"4ddc18a7d00a0bdda5a875971e6c731f2a2959130cdec37ee4913e54fd674b1a",
	"13e378c4c73a554ce2b28b9b26cd0f5dcb8fdd73cddbc71194d98da886e2fd69",
	"4b1d837735e6f1785f16ae1112989d538fb11ed2a1d26de666016d65be5a0934",
	"2d324d8c027cb2618369443b1dac74b2fcc99bb5d1cf24082b33fec22d3c906f",
	"97dfecbad133037a8b2848e3c8e29beb398fe331f8ac9c1e5646ef847c379a36",
	"1a5a97a862136f2f36777005d9e82947faa03b8cbc76520ddcd833c0f0300d74",

	"980578f78bf511d15c4bef18686a0ead48cd74fd15631e3300949e64d507865a",
	"e305b704759244b0fddbd82cfb5496b72a7eb8a76943fd73467df6073662d948",
	"7bc92979187737d63345c79420ba1251a011872357e023aac8539b9c0ec81331",
	"3e35c89f8742a9bc007e1fcbf5e736267b8f4fb3b1145b6c397662233e05f145",
	"8931d362632b2c3c19b14bd425d2c91d356573bf544c628e586b909eb849306b",
	"9a5bc0169afdad2c9e22478ce8b0b10dd881b4278e0b19561d42529c9c721432",

	"ef4cf925122ab9574e5c8e96da58adca49faadca63fc5049f8e61c89988b3d27",
	"14d1512620692e704b03839103075cc452aca3ebef3f9ab4dcff7bd0d62ef723",
	"dd6e8f2a38be87cbefece4adbc8fd3dfb0eb901da1adda39318f04333eb49a6f",
	"6af14986a4b3abe3caf085e3e7c373737cafbaffb8bb7f828d6c830806c2125f",
	"4d1c54d59ab47e71bc0e02354f1826bf9090d5d727badbbcb5f4a68f75a9dc13",
	"bf4ca4c88648f0a8bc6d608219da5c32e3c23f0bb0e0c65df40eba1c0f4b4116",

	"5a866e0e5b955a299d93417d16532948b06579c37b22b829f6751fcb1f222f5e",
	"c81e212420d832cae7575262188d2303e7a245c22725ad79d10ceda84ecf6c4a",
	"a5fb129ef16915a67975e42c6eaf9c766a4943d8983b104dfb2095e28daebe43",
	"a4aa6b8bbb5e54d4b59dd4695ba9369ab207981a5bfeb9e10d67630bcae67a3a",
	"44e8e363febe36514664a9f57713a10d5aa31b2acc13a93b7a411d386712146c",
	"e78f9b616bd6438c894cf37740e574163487b266fec78a7b962fc69372ddcb0d",

	"f4e77d81028ad7ea08fbff56163fc8b3b58956f09e0fef5c1dbe58e8b630c277",
	"83478fef1acff3b8e642530c6002e0c4093b39ec127fa5f568c88bd8a2e0fd46",
	"c358d1adb7f48a0c3baf35dc1c12734ae63132071380f26b0237bd0b6527b17a",
	"059a8aba8d5c99995ec568c9b799d4025e973088fa08665b6e844d32d754e027",
	"308e175d7a39cc88a1881c4476501c097cef861b65786c0e890b6fb6cee7890a",
	"9e7ea170a7c12bed18c6ce91455f14e90a310d18a60de655aec9106555683301",

	"d94b17ae0d156cd7701980163adfe43f53a52e4b493c24231aaacd1994ee1e2f",
	"217a0af174b9c48e16c18183323db11e7e152ba85ea64681fd5f3caa489c1470",
	"b538611b042a80842504533c0fc968eed008509773f9b97fa8b93416e8a96925",
	"480817661be3ea8a9f6cbf5213c7ef7f1007b207ec016ca141139f503212c005",
	"1928cebc87cfab8f727ee04c4f296604001f93e9110fd963f4e162f925c97517",
	"4950d17b374f18b9db841928058f58eff7f2d41208cd7066f8a14dab091bfd3d",

	"3451b6465903cd506a60b5d19e828f8fba9a30ae1e4a9a75a58dcced87db0a58",
	"87449490955acf02069a0677d351f2afcb5cc5e8aa2bbb2cd8b2fd53bac7a138",
	"4909536e0e2ce88ab529803022ba132a12ec1e140b59b1bdefc01273b6c00405",
	"7268449992fc72695f3d36908cfb70d9789940f03eaf1491d7f5fbdd93367158",
	"3b930f298efa8eb77007d1817c042bf1d439e5a1a6f9969160f5a17e9115fa29",
	"3a7e6b21aff84f0e95a627f9915a9b6c4393eef8e06c2184d6fd346006e6ce42",

	"01f3d5e31341701444bf1fff853e10ef39cd5e9da454e136190aa0c38370e46b",
	"c20930495669cf2ce86a3870a7341e881dd763b919cd17d9c6dcea5f9b1cef3f",
	"70c85608f178aaf809689cbe0e6be6344bf4c8ebfd04cedd822b809c34606172",
	"05e6d3585f6bc0e5bb6dbc12b6a3544cc16657738407674743a5c7f4d4356212",
	"3e1a83ab6c56128c9c789d65263b757072ab98731245430b979fa5a50479a917",
	"5962a17a505e76348a07a06fea32ddd380cd4efa4fca24bec30b18110281fc43",

	"9ba12f4e0eafb933c308fd5da5f4d3d55cdad1c197d848165d1bdf4a8d01b04b",
	"c48ac3ce9e271636bde728cfbf88bce1c7d4df91068aaabdf124c377f52b8a7c",
	"183a4140f451b54a43df9a6147a4846cb63cf829bb0e87547e082ad49f92f16f",
	"23e790504acaa7e2d5927dfe2cb98eee7ae47c1817c30334d5a5bd3ef4db1d2f",
	"1dd985dc467a5187d2ba6687597153085460c1dee9cff40705db34434c1d4f30",
	"0de054052572e7696e766ab3851e3c4fabf04d9c3fd1e4872f7f5736c1b34426",

	"4c3ffae19afb8620c6e867983a06e3cd3cc238323d93fba4142b5c801faae25d",
	"86a0d63d0cb41621a1f1c56030799ae34a57932e92d5081c3733f9496116101c",
	"a54b96b61371fcf4c3a2fa68c9bcfa804ff28b70d64dee42dd64cefa28849c47",
	"f0c949a2897ff7a9009f7d0f9b6464a82a64f8aaf16a452eefdff25e69503a71",
	"909236068a12858179fe0dbbdb67466b861469c835b3ac1b57eb9229dd083b60",
	"b18c21a34114f914cd8d155e7e6089b3c8b72e262aa83fd388c420a78a056c5b",

	"7c7e0d9f029cebc1efba013595a28b06f8ebcb5286a0d4611d336dfe4fbef73d",
	"ae068c18b240bc106b43e7aec262e441c700296284f83677f2141a6d5d7f6a67",
	"7bbc38a70a3064d8f52d8b6a0bbfa49578a595b0af49c8ddeac09fbdfd3f544d",
	"ea98385f84356007e7ffba03702e12a89086926936506fde777cb197d166bb10",
	"c6099de75d6711fed77a5344e8c7ee5355b608043ce9ee4c0c35e5f1e52b8543",
	"1fdf7f1fa1faf2041d58e757f65f251b03eccc093ab8fb206cbb70e0732b9527",

	"c062787a188868f598801bb2b45244d3a64987bc05aed3b2ba0238cbf1902a0e",
	"6bfaf1582a3518bf52c297cd0419aa4a197972c8f40dd9862008bf596f49d841",
	"8aafb17e7cc0ba1054dd06ffaa81803e20429b2dc9ed1d6594df8f5081add408",
	"9540226944ce86fbca7ec3d0041288c93f8de3ce714742dc9e16a43a631b8650",
	"c16ff9ff1e6bec20ac97addf4fb0027a7dd04841384c4a37a12e03c3083e3179",
	"a4bb87c6c9b78269929d0f24f428f118ec605400a5c26a22ac22b0b15f0efc0c",

	"8b50fbc3ef43eaacfd453ca05c20b74e22ab95db9507f9c9109d8c37af407215",
	"c877e5939c705a57fbaf5f75506a6172bb4a4436b23cf83e80c4a9edf242022f",
	"c8982f4dff39bc51419625d17bf0cce215f520e1a33d8d3b73cc21b1b2bbc658",
	"1e9077eee79ef046056ce2e3057d720f1870c59696f116cd54971152e639876e",
	"13bedf646bef836bf86f4aa03f08b3cd2f2532b45111053e92aba56986accc2f",
	"47cc22449fdc2af3e88904cd7714670505773b3d5483c4d332d25496555bb978",

	"8c691786576a60074ba1156ca5a4d0d7a71133b0c4cb7bbb8c87efbf0811b45e",
	"189d3e7a5ae4e805a4a1b9f89e410a70aad1927e1bce117fc7e365a89206822e",
	"03a508232dc80744e628fc8eca32b7a0abb85aaa7f58e49da8e3843a4ae58a42",
	"a1d95266f2390bcec21ac2ed4eba6ea5a506c9bbdc89c543dc894bfeb4457a77",
	"30a601da257fa562472a2f785f6a9518cf0fdf4de683ca2651ef12df2058bf32",
	"bb652bd7f4e631c461ad7c251a5e638c49b76af05daee1626e55cb7e12f0d146",

	"24156a730debf2b2a0f5efe84a6cfdd1d07ab61b6a16ea56048c99294d254721",
	"0f2722057dbcd9cd3124a52166ec5eb73e01bedadb16dce73860ddfbda1ddc6b",
	"2ff1a64d2ec074fe9d060d11015f8fd835dcc885ed266968c7baf437dbbc3977",
	"88c9a3bf21647f769fabe5440170f43f34e720785ef17a52bb487873894a3130",
	"7c5d1b2632051bc4304d2ad4353f484ad2500c52c2bf5d98e222337b699f6e6e",
	"df6fb6079fae6e04be70a8b0f1ee5ecb89392f14c564b698b96d99f4bbfabb0c",

	"48dff31893d9a18b9d140f545b6eaf96a03e752cdf7c099f4cd4e85dc3190845",
	"4fcb08eca6160c48057efcab290f3232d1f5f50b1b7f9be529f71ec5e40dfa09",
	"0a314c221238d38e37aa1600d319a88c5dda964809ab43d07ba4033f670ad46f",
	"d54657b2b5ba22cdd1c9ebe846ae4a53d6484feb641578dd26bc01c6e0d94409",
	"01816d3e0361f5f0bf007d7bded050b1f923881ed6be118542ee45ba9d95ae07",
	"089170ef5470057a02db477cb9e993eb00563459dd68e5ee4b1eda7c9be09038",

	"e4276557bb7938928a671a1d740a6eb174f46aae97da5b27a71ea3321bf57267",
	"664554c073fb6ba6c95721802949ee9dfb22adee718b5dd1fc1205a14aaacc7b",
	"a28a97509a17d7b05b73c57eb334257b05ef220f009fa9e7ea5d1b1df617486c",
	"2dda9bd368aa4eb61737bf91e9913f0f50d6aa10e441695c81f96306a48e506e",
	"7b128440f290897d78b553f29df73ab4831d4de1f617ef7d9759177e3bb35e17",
	"5af72fe3775742e97db400783f033da4a0620cb8305ffd8570788decc3035f65",

	"c9d8d4d274e82877e09ed0c0b7c32c46779bce29f22b7ca3852bacf297ec2c46",
	"a56593047ac3ccca17975e80655a44f32cebe8272c8d0415e70d65a8b3d69762",
	"ea59dc471f4de1f0937b20d2748c324b8d0e1b5f2592b44c769f44362545b852",
	"68ddfa06d3116251c20b464046e072811beb0351a2107b475a4515605bff9f19",
	"df74d55747c447d7ffa55545207ab3a26654590ad631174171548734572e8215",
	"89c435bd0da90d2007629c58db2f5220226363dc2c324e8242241a14c7f7de10",
}

// Point multiplication tests. Each group of 3 strings is: P, n, n*P
// Scalar n is not reduced modulo r.
var kat_DO255S_POINT_MUL = []string{
	"5dc536d2997635a8a3bbff585cebccb8b05f4689535d464a687daf26b0a85e50",
	"ee5f3f239a06a05bb113cd4f568e4f85ec3a6fd7f394b0d35c30af7ef25cb088",
	"294d9711c65e9d92d6a461390c74af1e6736d9eb1beab7956c95c5f93f7b0016",

	"22cdffbbb8f87841aab2ec1afa2c5fb5cb8f67e286d0fdcaf6124951a21a597a",
	"24c84823084e64146d09b64a58ed5f39256b568834c4215687ac642436a5d17e",
	"1ea8dba51fd35925488a671414cbb276c1e9bae3731967ca63cf53cae6428b0d",

	"2ca39783b2d55263d0b1c8f0cdce65f35536959e0127d8522b9624535d750e7c",
	"0fe8dd3dacef7df60e4ba71de3d35e5e6aeb6fbf7c960f72626fd3aae43c8b88",
	"ba32589f0ec3475f0ec62059e65a370df09ddcca269c0548d14e66967b95f664",

	"b514f908b623f7858ab324cdef9eea187b094192cf3c5c55e0a657cc13faeb14",
	"84d3e49037e45ecde991d0742d2f82c3ba9e2a0f0902fde2596369d9b0e5dd94",
	"8715245a20000eb3af5970f87237c4bba588af8e703143089b6b9f86fa1df51b",

	"73156faf9fbd65fd4af9ad2fec458760cb13ae044adf35bc4212ea0698e46d4e",
	"d350128f36b85d365060770fa450cb4a37a95842b5b84de39528c2a66b5425e6",
	"89daddb7da6f0379478aa027e963d140c6a3696b791140fee1e6b50ded5c1c56",

	"e74f7f5f0deaebbecac8ec772db37e5d495ba594e4f97a42e36ae544d1b6b748",
	"8ee094aed7a63b40981dddc40ec2e8ac9e1b8d3d3e40bdb4fdf6e35ad2afb262",
	"0a5fd0843da111a93b585a877690136cc610a0042d90937cbd5c88a7b9bfd76d",

	"ea9fab33e79dc4c47f83d95eb0f8cc0d486c78d1cb512603cd83aa7f7ecda23a",
	"a8e4a6496f60ea4cc879bfce649ffde1651c6258e6fdb6109dc431faf1f1676a",
	"7fd8e7464f704e0e4476ff6f188ee8b89ba8687d2dcf2097e30926289859c45c",

	"1bcb940392743b88e7bd5162aba895bd215e6726ae560d7c50162ccf60d61229",
	"078d80335d06b2f52557ba6382a32d6531bef339e7b6fd96881a332800f4d85a",
	"270537a047f33d4b82c4d9f12c3c3e1a3ee981dbaa78c7715f5874add68b9f7f",

	"fed66edc33f4bab3130c952ad4593bdf572b2b05ee0c422f4f7681d14549a55b",
	"cce339e0bfd3f9e30f14a44217e1c0bafdcb03120031a9eb234db4b25a6b359e",
	"433eed46a714dee7c7ec6b28fadec8d4ec67d9e0064fd12cc82e11a0b23d815a",

	"95603f9756b59436ba21c374e2ec38d179d719de24f9500bad1b4d06a0f7bd55",
	"0f0653922e02f8b71c7deb00d61dd6bb88ac94d158cac6f6018fef6cc932eb34",
	"87e70397eb175f392cd150dbf05ac8660ef8edc5ce9063fdc9a8d3002c90013d",

	"e01ca2fc36b41a801ac05c3d9f07befd1c08e25ff5a7b5f43956fb49a69b8e76",
	"dda5d59e8d25b5f6819d70bc3121120f2ab43cf61a676adeaa2237d16cbd35f7",
	"b55d2dffb68ee646b2b65bcd7164ee73412d8b2d81c7c53bf0b0f3fdc5367a14",

	"f989d9044394cb23198341262a09b701546acd9a096eebb4ca6c6fce9e790207",
	"e41ea9fdf13d5740295d30b7546a15bf74cc6fa08329839e11d1e902a1e8e338",
	"98a8193f4ac5ef3025992819a146e8d506e9c0896180cc89c9dc9f92f81b3f7e",

	"634e2ec626c1b57aa596122b196cd6cc4bb22ad2ade86703d9a45d24c8f19c1f",
	"55519d76b70e2a4fcf5b7f2161dfbdcb2ec0366f33f614343999bf5a43f07d4c",
	"92b4f42339ae4c0d290c4ffc543cbc13df85de7de8307126a595e6a32ea83349",

	"6046ca0c0c42306a8b600a0ead94364d22249b4a422d363c9858d512f6fc8b6a",
	"7e88fd3c2ff9fe0675e6b0e660c20a6ce0c111f096fe492a411edd9ffeb6cd48",
	"d317b8695e7797f2148695e2235efd7b764cf700d9010e28ccba1383d615864e",

	"dbd2df5aee4e9c100cb43fa8882df37d26bfb1cdfd732f9023cd9b626ac01a5b",
	"6e1aea679e9989af043f8b8ec0501fd62fe5f912d281f162f5702fceb659ba71",
	"a0a0c612850f2f761d47732b12f1b870218f77f01cbd9ceb9e96239e939b1a51",

	"a4b5c4f695e5a7da5081f2d00037220184cbe49cc69f3f0b6a6f9196c55ab924",
	"3f4d9db8fb2b46ea67329dc8110333dc2cc7c4d71a00a2fdda4ccd32fa204e23",
	"354446a1a2e23ece38fab88a85c926bdd4e338107386aba1f76192fc9565c178",

	"12c6501fba2cbbda12625bb61c197132a81fa722ef3a1401b5fca5104fa65e40",
	"1a688a922425c286ef6ff140ae11ca95cc175412233ce0a3abb22500c52c5bf2",
	"01598c40b7500e3c6f406d70932e5e2e116cf1eb93d0dc4f6fda6fa1ac5ad77c",

	"d00c41d984b324d821942f265bcfddef1b92d63d702470f50cd60015c52c012f",
	"2c1c60d6ad54a3a0e7b2f1b2dcf220ff51780217d6a6f1a88e677ef676b6580d",
	"99f8e913e519e0c94dea397dd62b2f80394cb0e3fd8e4c0eb66e040199fee156",

	"9c4401b972d38bc839bcc503950ce8d5cbdd1bd7b389124a21798bfb94dcf659",
	"198c7dabc8da7215fa19efb2d633f75790b02ef2e9b5d15117da2dac07aa8e7e",
	"a969d581f759f9551941c0f6ea801258ee7590ba0dce346d200ebc745a79677a",

	"137c9e0678b34e923b0e96e03d42181812105228f541b2aa2ca047eabf0fb96c",
	"a5f3e423c54bfb63b6db2908dadccb141ec9c8da8f20f635bd7542fe4bc311ba",
	"b9263ce912675c3602b6941f735610bd1ea46749b9f8f86af48adb203d631b6d",
}

var kat_DO255S_MC_POINT_MUL = []string{
	"a90da7dccc5a50916e7ce76e8292ddb9722a8bc210fb1a71258c50ea17e15722",

	"ae754a06df209494ece70bcffb24f0f4ddcd89c2674ddf9965d7ee14f69e4978",
	"fc34479bf9839114bc53daf081a155d68187d55809e684015a53f0f0cfdc9c0e",
	"c8c05caece52b2c29d6dc9f0bdd904ba42bcefe3f401098675a8ddb219ba3663",
	"88a4e1c958f1a24c016787c1caff2f6879b82380d5b02d92a29e2068bfa66328",
	"d50512d8b9580075fe99f77f78710900cc0bc21ec95cbe67e5de9e42e8c0f53c",
	"6da965ad8b472a62411bee2ebfc334868e3d6d4e7d82356f63cfbaaa86be3727",
	"85342d60cf733489500a66b5f130b5ed05e4fd35be61f177c4b8a30dddc8677f",
	"67d7674f7febe3e613f958010747a92d5724434f992b7a4446ef9c099b78b028",
	"eb0f2cb73d44fc21f8af9bb30502e44c8a07c3578feee67f0774d6541bf0ec13",
	"5ed693571fa186f9370e4eb4cea05ce4c28c7bba60d11f8c7ed0726d14a73f56",
}
