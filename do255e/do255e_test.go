package do255e

import (
	"bytes"
	"encoding/hex"
	"testing"
)

func TestDo255eDecode(t *testing.T) {
	for i := 0; i < len(kat_DO255E_DECODE_OK); i++ {
		bb, _ := hex.DecodeString(kat_DO255E_DECODE_OK[i])
		rc := Do255eCheckPoint(bb)
		var P Do255ePoint
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

	bzz := do255eNeutral.Bytes()

	for i := 0; i < len(kat_DO255E_DECODE_BAD); i++ {
		bb, _ := hex.DecodeString(kat_DO255E_DECODE_BAD[i])
		rc := Do255eCheckPoint(bb)
		if rc != -1 {
			t.Fatalf("Invalid point reported as decodable (1):\nsrc = %s\n", hex.EncodeToString(bb))
		}
		var P Do255ePoint
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

func TestDo255eMapBytes(t *testing.T) {
	for i := 0; i < len(kat_DO255E_POINT_MAP); i += 2 {
		bb1, _ := hex.DecodeString(kat_DO255E_POINT_MAP[i])
		bb2, _ := hex.DecodeString(kat_DO255E_POINT_MAP[i+1])
		var P Do255ePoint
		P.MapBytes(bb1)
		bb3 := P.Bytes()
		if !bytes.Equal(bb2, bb3[:]) {
			t.Fatalf("Mapping failed:\nexp = %s\ngot = %s\n", hex.EncodeToString(bb2), hex.EncodeToString(bb3[:]))
		}
	}
}

func TestDo255ePointAdd(t *testing.T) {
	for i := 0; i < len(kat_DO255E_POINT_ADD); i += 6 {
		var P1, P2, P3, P4, P5, P6 Do255ePoint
		bb1, _ := hex.DecodeString(kat_DO255E_POINT_ADD[i])
		bb2, _ := hex.DecodeString(kat_DO255E_POINT_ADD[i+1])
		bb3, _ := hex.DecodeString(kat_DO255E_POINT_ADD[i+2])
		bb4, _ := hex.DecodeString(kat_DO255E_POINT_ADD[i+3])
		bb5, _ := hex.DecodeString(kat_DO255E_POINT_ADD[i+4])
		bb6, _ := hex.DecodeString(kat_DO255E_POINT_ADD[i+5])
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
		if P1.IsNeutral() != 0 || P2.IsNeutral() != 0 || P3.IsNeutral() != 0 || P4.IsNeutral() != 0 || P5.IsNeutral() != 0 || P6.IsNeutral() != 0 || do255eNeutral.IsNeutral() == 0 {
			t.Fatalf("IsNeutral() malfunction\n")
		}

		// P3 = P1 + P2
		var Q3 Do255ePoint
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
		var Q4 Do255ePoint
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
		var Q5 Do255ePoint
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
		var Q6 Do255ePoint
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
		var Q7 Do255ePoint
		Q7.Generator()
		Q7.Add(&Q6, &do255eNeutral)
		if Q7.Equal(&Q6) == 0 {
			t.Fatalf("Addition of neutral failed (1):\nexp = %s\ngot = %s\n", hex.EncodeToString(Q6.Encode(nil)), hex.EncodeToString(Q7.Encode(nil)))
		}
		Q7.Generator()
		Q7.Add(&do255eNeutral, &Q6)
		if Q7.Equal(&Q6) == 0 {
			t.Fatalf("Addition of neutral failed (2):\nexp = %s\ngot = %s\n", hex.EncodeToString(Q6.Encode(nil)), hex.EncodeToString(Q7.Encode(nil)))
		}

		// Testing negation.
		var Q8, Q9 Do255ePoint
		Q8.Generator()
		Q9.Generator()
		Q8.Neg(&Q6)
		Q9.Add(&Q8, &Q6)
		if Q9.IsNeutral() == 0 {
			t.Fatalf("Addition of negation failed:\ngot = %s\n", hex.EncodeToString(Q9.Encode(nil)))
		}

		// Testing sequences of doublings.
		var Q10, Q11 Do255ePoint
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

func TestDo255ePointMul(t *testing.T) {
	for i := 0; i < len(kat_DO255E_POINT_MUL); i += 3 {
		var P1, P2, P3 Do255ePoint
		var n Do255eScalar
		bb1, _ := hex.DecodeString(kat_DO255E_POINT_MUL[i])
		bb2, _ := hex.DecodeString(kat_DO255E_POINT_MUL[i+1])
		bb3, _ := hex.DecodeString(kat_DO255E_POINT_MUL[i+2])
		P1.Decode(bb1)
		n.DecodeReduce(bb2)
		P2.Decode(bb3)
		P3.Mul(&P1, &n)
		if P2.Equal(&P3) == 0 {
			t.Fatalf("Wrong point multiplication result:\nexp = %s\ngot = %s\n", hex.EncodeToString(P2.Encode(nil)), hex.EncodeToString(P3.Encode(nil)))
		}
	}

	var rng prng
	rng.init("test mulgen do255e")
	for i := 0; i < 1000; i++ {
		var n Do255eScalar
		if i == 0 {
			n[0] = 0
			n[1] = 0
			n[2] = 0
			n[3] = 0
		} else {
			rng.mk256((*[4]uint64)(&n))
		}
		var P1, P2 Do255ePoint
		P1.MulGen(&n)
		P2.Generator().Mul(&P2, &n)
		if P1.Equal(&P2) == 0 {
			t.Fatalf("Wrong point mulgen result:\nexp = %s\ngot = %s\n", hex.EncodeToString(P1.Encode(nil)), hex.EncodeToString(P2.Encode(nil)))
		}
	}

	bb, _ := hex.DecodeString(kat_DO255E_MC_POINT_MUL[0])
	var P Do255ePoint
	P.Decode(bb)
	for i := 1; i < len(kat_DO255E_MC_POINT_MUL); i++ {
		for j := 0; j < 1000; j++ {
			var n Do255eScalar
			n.DecodeReduce(bb)
			P.Mul(&P, &n)
			P.Encode(bb[:0])
		}
		str := hex.EncodeToString(bb)
		exp := kat_DO255E_MC_POINT_MUL[i]
		if str != exp {
			t.Fatalf("Wrong MC mul result:\nexp = %s\ngot = %s\n", exp, str)
		}
	}
}

func TestDo255eVerifyHelper(t *testing.T) {
	var rng prng
	rng.init("test verify helper do255e")
	for i := 0; i < 1000; i++ {
		var n, k0, k1 Do255eScalar

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
		var P, R, M Do255ePoint
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

func BenchmarkMul255e(b *testing.B) {
	var P Do255ePoint
	bb, _ := hex.DecodeString("5863d164f563f2fe6b720174c440d41666397b4aaf6d2243c57f5ac77b6f5361")
	P.Decode(bb)
	var s Do255eScalar
	bb, _ = hex.DecodeString("81d3066b23528ff3728a200ae411d815ba612b45c2556fa1000cb665de692a73")
	s.DecodeReduce(bb)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		P.Mul(&P, &s)
	}
}

func BenchmarkMulGen255e(b *testing.B) {
	var P Do255ePoint
	var s Do255eScalar
	bb, _ := hex.DecodeString("81d3066b23528ff3728a200ae411d815ba612b45c2556fa1000cb665de692a73")
	s.DecodeReduce(bb)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		P.MulGen(&s)
	}
}

func BenchmarkVerify255e(b *testing.B) {
	encR := do255eGenerator.Bytes()
	var rng prng
	rng.init("bench verify do255e")
	var k0, k1 Do255eScalar
	rng.mk256((*[4]uint64)(&k0))
	rng.mk256((*[4]uint64)(&k1))
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		do255eGenerator.VerifyHelper(&k0, &k1, encR[:])
	}
}

func BenchmarkVerifyVartime255e(b *testing.B) {
	encR := do255eGenerator.Bytes()
	var rng prng
	rng.init("bench verify do255e")
	var k0, k1 [100]Do255eScalar
	for i := 0; i < 100; i++ {
		rng.mk256((*[4]uint64)(&k0[i]))
		rng.mk256((*[4]uint64)(&k1[i]))
	}
	b.ResetTimer()
	j := 0
	for i := 0; i < b.N; i++ {
		do255eGenerator.VerifyHelperVartime(&k0[j], &k1[j], encR[:])
		j++
		if j >= len(k0) {
			j = 0
		}
	}
}

// Strings that decode to valid points.
var kat_DO255E_DECODE_OK = []string{
	"0000000000000000000000000000000000000000000000000000000000000000",
	"5863d164f563f2fe6b720174c440d41666397b4aaf6d2243c57f5ac77b6f5361",
	"81d3066b23528ff3728a200ae411d815ba612b45c2556fa1000cb665de692a73",
	"eb01802d010a2e2dd6af94ee08ce72f89b9a7996f197174e3aac3086d45c5c6b",
	"8b38c5596c5f334db20143d106813e8b47573452098af9dca9fb35d3a0f9c77d",
	"17f96c91779fa4e8de9efe799dc513bb3830a9295b28f7ec78d1f651eeeeb861",
	"901818436d2422c5c4b9120506922ca9209c422b3465acf4b499e416679c076a",
	"946b683865a5d8daf4c643bf378ba87175e333ea53f5c28803e7da1e8e37c87f",
	"061958de1ff3560bfc352b6b1fa4e4169ff443f540f4a3356274efe9a504ee37",
	"add42bb95400b3d2a07f6e096bc1a762578083c38c67d6a053648973e950cc62",
	"c081d24b00b5e41ae485411795241b46799d4c8e6c9c92c372e220021628423d",
	"e9e9b4b160a1c391b0abbc321bdfb81168960249c785e812910bc7e110017065",
	"b1463c47e5d9015a6df289dd3c72b2f8f92b619683561aa4b518ea2c9346a06d",
	"e66400cdde69bd0f4740c18dc734b5a9c8a73a731a9fc2a2d2bc2d82607e8a64",
	"ed31ae129c1d51524ac74b0e880f9c4e899c9d01246fd1e14e6647c6b07e4f4f",
	"402e84ac250247b758aea5bcf350342d351f4c1d136be05bbaf3b67536e87f00",
	"04f27a42e40833d894f41830a2bb69bb704e0b964a2ba6241b916bc5b2f9e442",
	"d611039fb328c085ba83e560a1f4b6fb302a2a19e27afcf79b9eded0cb91761a",
	"c370faf1b9aa1ea5768ae36bd1712d7673cad20416d5325b5135b3acde3ad25a",
	"a0e880d8fdd28fa5d90fed2c78454c7c1a7ee7254e4166a626e9ddbe8aac9d54",
	"3a9418d0fdc91071c347682c78cc113b6c3c5a8e6790817e54ad798e958e1f31",
}

// Strings that do not decode to valid points.
var kat_DO255E_DECODE_BAD = []string{
	// These values cannot be decoded (w is out of range).
	"08761e618ba81e1cb145e406e8c6d63b115d2ca460a585a85275ffed14ec1ea9",
	"6743c991104ac09dd95f81607645aa69892ac9c79e186fb666c201b465f40e94",
	"435f64012646f31d93edc5f20a4bd03c51b30cd314f18ac23ad74779e2190eb1",
	"59a13eaefc3e1043f748adef252abf956ca58287be56a2a8990fbe77c6790689",
	"d9001e28d9468785a32fa452a99a5fa617e3cc25dac80c7b4dc668883ee19eef",
	"4b80e321fc3b29f3dee7ee0450eb309bbc3f628e0b0e799e1161a3c0f0ff8aa6",
	"632d62209619b3a388d00001fcb25d627f830517814267f3cff381c6d83b1ef2",
	"7ad205df3300a6350b03c8337c89550564fd8c91779e9d0342dbed07346e3a80",
	"a9747cda67e5ad2bfe85253f094c105856c4ca94de7f58ec1ca64353eb4eafa1",
	"0049a4befd5678f03d3b23a986dc3350fee4cf41140cbaf78db10919ecd288a9",
	"45c00fec5f4ee698f8fb70fa6fbb102bd2707878cdc22dbec2a073d7670752ad",
	"f32961c81d13acd47d56e10abd4cf21086ee9ba34cea03384bcdd19b8336cc8d",
	"120106c537b05cd7ec4a24b83e99c399f212776c0e3f43ed923a7cf8370316b1",
	"7160d3059ab1eef1cb3f33d4296e0ca05594f39f59b11e32e1e61ba478ce05ed",
	"d700eecb2c254be73f7d19aae82c883fea2362225769b4cdeba60f63cccc7aed",
	"22c91e4ae9566785b1fa88bcf6d57011300c63af5f4c4f16a6366d57db0973e1",
	"b5ec84834affc9d9026d2051da7b4584fda756e39f1b69df8a36a8b35902e3ee",
	"bdf310230866f6e8ed2497956f148e01f7d712ba46255de4e6aa2eaec6d25adb",
	"0b03203dd0d310b4f203ff7216bb0b666e91c2680ed06aef96c78f7a4e72a5f2",
	"b7dadcce34e857c65acf1e2c3f36376aeb906cf48f37559582d1b6067ea0d5a8",

	// These values cannot be decoded (w matches no point).
	"7d901cb1f92cd676a64241f368a415e48c6a8a662566b25a309d8f1823eab363",
	"e4f8ba84d29399e935a8f8dfee14f1ffb5b5bc20c4889ab391c0de6464a5ca46",
	"e3603e6b73b33cb33c9cea3119d0e223520ca31dee0c9238a10f7036b0a62f5a",
	"d0c7fb74ee0c7315b689cf84efff6c4407d86a1620138840512090c26ad58972",
	"630a1f53aa31c5809408d113a3537d2d63815c221378469b309ac4544cd33877",
	"bda70de7c625827b184dc76747e1f52f218259720f8790082c409e66fc771a2e",
	"25128aecc25469119568d04c274e1cb523bf5007fce1487957db96e14954f11d",
	"20c393a3b0f36c5461a5edc04b574105d2d9c8ec615626810b838cc2b77f9a43",
	"808fab937bda301503242d8d6f8da0809b1e9c9d958b46f1c9b2cde0ba365d29",
	"bcc21cf263eb2840bb6e7a974c91ecbfdb6be13f9cf37c749652086abf4c7526",
	"be0e14ae5cec6ec93c83274c83f84dfd9296243da9209151f5e020af6b29103b",
	"ad82e5d3f63b9507e8cb5faa79052259da755a7f1437fc7cb07a5398330b6147",
	"e7fa1c035f07d020e6b91d6b65f08f6385e3b39cd449c7f67729b1e058103038",
	"f855695bf1c0fc5618437e1c10d30f40ed0b9b823d9424d291dc3bf3f0203953",
	"1d41b6b8376bbd4316623b800fa253d674feafe3b2a51bc709b7b3619fd5f025",
	"af249c269234142d7d87ac6a23e2c65ae26d18e5493ae15e3469ee46c73b9125",
	"3044fa93f5d5a1acbb4bbdd050aabc1b15d6c95fe4c0ab74cf7cb0cb8e503d3b",
	"b88a6ffe9e90bd7546a22ea54bef4cde2e5e6cb36bd004fef175bc4567ce3120",
	"e28e882d49e4bd2dc3e3301ffc2e8e234810d067ed7783e3f84515985a65f812",
	"f326a01c98b91363347b6e6816366be025f1437a680b99993f28c0e41f01b323",
}

// Mapping of bytes to points.
var kat_DO255E_POINT_MAP = []string{
	// Map-to-curve test vectors for Do255e
	// Each group of two values is: input bytes, mapped point
	"0000000000000000000000000000000000000000000000000000000000000000",
	"0000000000000000000000000000000000000000000000000000000000000000",

	"7c77393af403852b621bbae660adeb83939dd0a5f50a986438f1cee5aa076f2b",
	"f22b8c41ef41d627ccb0315b63a344b4bbd0c468126d48c5fb464c2a10234e25",

	"2a24450b202dbd95efdf0926740f914b548057b52840dd6954d3361c135ababe",
	"28a3522313142c44e6495452ece93e213b4a794f9238f9bcc80389dead84d955",

	"82a075600b10b133622709204280900fc0386d666676fee3708b129838c6bed8",
	"bdcbf7989c2d5e9f266c766bdeb8290a0e0d0f23faa5d21906d290f541ca8847",

	"e21baddbe4292c21084e50113560b5b6192a34c5a349d55eafc368fd62ff904c",
	"dd3cf979869404d984555aa660a60cc9fbea3c8c7efb6bd67bc366f63acedd38",

	"e6dd4f7067fb4c0cd3a4dcd4e6c48ee0bbfdb4745fbdd548cd16405bac258181",
	"07a37c1c788e206e713df4318daaee2997427f3016849fd09319ccb2584b1117",

	"3a396c8cab2a781765b7ee28485320ee8e662ffad75b6e45ae16eca1cfc9d62b",
	"d64c989d65d58f1d5f3a0d214305246fe91a2573877125a01ff92c2bc465d76f",

	"296e39ee957395e29a2053c7d7222555030500f88b9135cd7d1b48a95e654a93",
	"4927392fa7d750f00559acec28f480ac4937b7d881c59f63170292435cc6ed69",

	"c5562800cc863d83f735c7ac9fbe1c36f2e69a3b5af9ffb8dc4a7e094a771455",
	"5597a65c7391679627df68674303b900e0625a499a822251f271f1a51ee80e5d",

	"b3f09c45abf1c3b8e1198a7572df6dfb4f62df6cdba0f76ec55b225e0fe77ce5",
	"5bd3b7683a0cbbf0b9f5dba9893d74bff57db1ec1ce56071cb14d00d60e0cd58",

	"0a19a23eb556a74203b62a06ed663afacf691ee13c62b384148e56e6c2b2b37b",
	"0e630f3953e0132cfd8db0decdaaa33060948d8e3b36ada436b5078b451ed361",

	"c28a168b21c3a1107533c998cecd45597b03c2be70554013521918256f5baf15",
	"9f780a60fe7435024e455049ca6c4e0154d4848e6172c3fbf1396249b2cc206e",

	"552ab32fc9e4d4226334a31b08c05e5e5c0ddaff8b5f291bc1762a9f8d920ffc",
	"be5de20456aa8904c98ae2dfcab82fd47ab492c59aa3a9af26a9bbf52683205d",

	"cab9ab03b89c099894e759a01fc0faececae6eaa7a0f67cc166bb5a4c422634f",
	"1896d0514ed6f9795175e22b34d19248b117dc684cc3ffeed96dfb28eaeb633d",

	"471e3c38e2550fddc47246ce8583d0a6a2316c7aa207861a357f6e8995c5d3bf",
	"8f6085b268cef2d676ee879f608abd613eb593716d9b82f6e659f95c443cf013",

	"b4d07e65ae3a7b808e619350b7f39a118cf779cd2fdd691668d6cc0d6684cd01",
	"f3ac72e9b70cfb74f1f53d90abe98123fbc8cf19af9d1c04eb4906455f5b3c26",

	"c081b2e772e34361a99001ce607c75271a615588da226a0095ef3f6fca15ea06",
	"7ae24441ae9477d82e4d8bbfd815c120d6def83ab992f6a13e30b47b59587377",

	"e5b8e5d2aa21225f8bf2a219952233b328fb975360c8c7adb1e56335abc6ad41",
	"ae3840621877c76f4572ba238ccc4411a68fe4596d8634b5925bade40357274b",

	"79e264d13a72e2c262000f89429828ed66040cca0a1e7c949270221e12dbb80c",
	"71628b93d86ba6b8a7e72daa6c44ff6df82a80f91c3cab3e6124594cf284396c",

	"16788253ca890ad995049642ce00bf04c1d8083d25d974f7780ba963034d7141",
	"2d77a65d15122d977149e7012bd70899514545f42e9d3585251dd39fb90f222b",

	"0d80f3763d03cf4746439badfea8f28bae79220d636b374b31c401088e83ef8b",
	"1a3e56d52ca4d670dad6de6a34db373cbcc3f0b763638006dedac85882570127",

	"7f0004fc41e569b73041691b7974c51a9b1b04ea8f76fa1caf628e34bfb1f4f6",
	"3c6545119f989679e9e363c322e50828861d9a9e25ff8a21a26aa5bbfb7ec172",

	"88b2470104159787e3db1e024c3db4eba7bb795f347270b29d5ec7f02aad49c8",
	"f6e875b6f81dda89c73b79621ec942756b383db0a510fa931354c65c1d50871a",

	"0262670011cb98108e1c1e776eb7cefdda204c66f5e6bd1664229e61ed08817b",
	"2f307bbc3a3181dc432e869ac40468c587e1f469748892add24bcdbbc9971046",

	"828a4cca2cc9cccf7164aaf20770536c3174921d0b55c650b10faa588df4f2ac",
	"dbfc23a3b67028ff389609d6af5177c29da61b34452c7ccf8d6e6562fdaf2b2a",

	"8d7bd596075453d76c4562bbab98589c8e21f5623d60189335362329d15f9265",
	"f0097dc3b33994fcdaaeb6bf0629a3c6da6adf03c497ffd30768d8065cf07303",

	"a1770a277066fac683ac2c2df85e21e657f12c2e818e54b80e3bb1ce2c8d5191",
	"8e5207791c09593db77edef73767aa2ad9108ddcff5ea0fa0a2e7d711d9b5f67",

	"8bfaf9e763343eecde43e3a22ea33626a4a016d23148c25fe3ff48386f3d8836",
	"0e7847eae28288c147064370beaed29895bfdc22e5b1e9e89e0c7d6069889d21",

	"e7a2946b955bb7be27adcf4c59867c24abd018d46e57eaef53a0d531d1def7c7",
	"66d6c79a8ef0651c10228c80634254d5e948ea301fba87b76c72e9c1a51b6e41",

	"2583f37d0878c0cbf4bb6478a804c9f2f7b779c9c22104a4ebb320a3b501c707",
	"5c90c5abc2bc789c82bca2aa7fa247084a9514d4da8c7db077ed27aa7b105e4d",

	"8381b36d9e53d498d93446895736efd43a9d932201925d7bb851269347736e10",
	"1966062ec5e0de2601e9615a08d2eca562c4b55c44ecaa7abb487cc6b5b2955d",

	"a9f08891e93c5750784f818e885917e05e60877c8f6a80f2f756f68a7c19ea22",
	"c9b1b5412d6741e59aecf34a9b4f3c967a75ea47cdb31953d7a88a793a808525",

	"5672edfc5480a3fe3ecea6d3421efa43bb23f74a549df865e3f31b3ba083be72",
	"9755fe7a49517013d6a0a9c300ee4e9d363df8426a96acd7b1b4245282f95a3f",

	"8ef25822dea237140212c8fa43a6ac5b4c369f032c020b06a9d4dfb337612ada",
	"5d2f36e163362a8ef2050f7a56e5ebcc5cf47cb1841ee89097f0864eb942996e",

	"114606d85fdca97ef2717c6819212ce8823e1ba09157e96f1e10239c90ea731a",
	"085e7474203f9f8d10069edb5cd0eb53d9bcda19402a239ddeec703c9e92f534",

	"2fa6413288fd64a5390a01cc2bde69181113f2b3d4ac87b024577e0e8036f9c5",
	"0e9de1be3279625b8013c4a1646aedf1b0be75bedd4e73f049816dce6651126e",

	"dcdf5d997984a90b61be027b9d0aa5a9e6d4f93c22ce72162fe47c3ff6d92294",
	"a8ba440f9460ab0667b50119ea613300fb941ddc82b291ebefaa5e73a7bdff00",

	"041a24caf2b427154320009bcdb12f8a859dc02f4100cd9f731d5ab995d25c16",
	"1e560340a15e33c02edf628e2b4bf3e33f3ec50c469d88600b31e2aeae67a12e",

	"f1617c38e005c0d8452edcdab243e7da129484e87970c462a1582773e5166874",
	"43886629b6814a7f449b12077afd9e1990d9b50f4488fc6c8e755d669575d220",

	"f6c657965b8a999310799f4c55f9e352f9f78602a4a37c2e937aa117a153f2fe",
	"025834988b1c7719797b8020e7bba4a3b449c7ee3a74119c396c066db8d2992c",
}

// Point addition tests. Each group of 6 strings is P1..P6, with:
//   P3 = P1 + P2
//   P4 = 2*P1
//   P5 = P4 + P2 = P3 + P1
//   P6 = P5 + P2 = P4 + 2*P2
var kat_DO255E_POINT_ADD = []string{
	"8bf6108baa2a18e584a9d9ddb84d60d1fb7152c0b070dc32211454aa63182028",
	"71cb77fb7b314b5917f984007e4b9e8da84e30df23b6361f3d5748859b3c3d46",
	"05efef11a57d4bddc59507920f67c702744485e17bca3960b9c75a6bd19ac535",
	"32a7f632e99ef0e67af738cbbb35cf6e3191b614840ff5071a67f0cd92a66169",
	"d3e338bd7a3d585b1d77795426a8ddd08c08b92a981087887d29e4643599316c",
	"9d9a190314a3a2f0640faafe0d63821380be7ce4b41bf2cab0bc0f3c1978196f",

	"7827dbeea7d4909ce79a9cca5c1683196353d9e15ce378e3b503c6f8fb937156",
	"a73ad61994300390f82e1a1bc730cfab1b47bb15e9feacf6831bfbc075af2554",
	"abc43d03dc36abdec998f1076a57cf673e9ab6ea064e3787d652fcbc84267f6d",
	"cbb97881af8dd4c407f4bfed6e61f02c74e87a16bcd8046336e324bc8c1cbe3f",
	"6553437175dd8e34e6a977f4abac327624245d7237d57133af319b46c76baf71",
	"b890a1ce8535b4b7191b666aaa8740f29ec858e6b417b287ba786b0eef881f1f",

	"16da527e26afc209d365df1cd4ea3613f2d40272d51d03e16414f662b4645748",
	"ab7791b15b2a06b5aed5ae8ff815db1e6a62495c30b8e4021a0d3296b93c4c41",
	"c439ec6d77de43c6ec146c843bf098b668200abf23831f4bab5aa1c4738e7c19",
	"2f82160afb2b3c6b394ae6861bab9e3077c810c7f136213d6ba94cdfdddc1d42",
	"fc69d563fa3c7b5b00710d7e5b940520e8990474c7bb60bc07d823ef5ae3837f",
	"a900860305b2398ac56d3a01c7823eaa80551b55f47e70eee3b9186be1da5662",

	"949f530930d463015efbd3aac6cee9e11e3522785f889ee5c8163a35fa06a73d",
	"b35f9e085825f3d62ff28700db54be2eb4efe30d4ca65ce16ea1b76b6c44fd2f",
	"06cc9347b401d20480fe0e2e85dcfe5b72bcefe30f5eb1c7cd900f94c084141b",
	"339e4e5c897420df32f25ebe5838b581e02f811685ce80f57c4406def28cd911",
	"efbf2eea194160d60259e400f546028f1de1ce57e41c24070e37855748d9b92a",
	"ef9fe3a03923d1a53a3f72f8f92bb2f538096e5414bacd9ef1b82baf4650c717",

	"522645d388ba218cbf28ab87b1ae6358dd1c609c46e1e174aff0d44c38855b56",
	"32cede21959628e69241141709677523a162fd0f9e7872cb8e4981e618293f69",
	"a813de7a5dea2df59bd5db75322fe62aada06cd0ca6e4fab52c364c3f9a02f35",
	"a1bb202cf9cc7237c8a8529e7512b3fc931e4b83feb4dde7ffaebaef4f5b4d46",
	"013a508f242ad1e205eb9c9437d49a80a355fbc5ecde40d14ed2e37918a2c651",
	"0984e413fd315d3ceedfc33b25d1119c4e5ac3153f481fd6304cfa05b69ca112",

	"0375359b6d6614120e445bfc974b888858fcc20469a1362c386f530bc76f8f07",
	"fb26cc0184231ee374cc0fb319bbebe07727e954fcd5316c9d752b2d5360733b",
	"4928c88e81caf3722037798d20bfbfac13ca3f89faf1fb0275ae334328a0d532",
	"d3caa6985f1d538c973a4c5bb65961f4c6b210fbad79f3f4a7a2875c591bee18",
	"022f0fe0313a6f07246906d42508aaeac38d7b6fc2b15faee197a757ef89a90f",
	"5b1e42fcdae23b9c3c3b6f76b75a151ac634eccc003afaa8fe28951c282de648",

	"4440458c78397f91bb5b8d55e62fa2ca76de17c15d9451605b9f05a482fdbd38",
	"e6ddb801f2dc0b11ad0c0f22216bee8bd32af5a706168e1bb6f57d8bbf492312",
	"5764a6d58f7b80a2f7ef3a5cdb4c15eff9fd2636e25a99db632658f5342fb925",
	"186e3ebd608eeff2980ae50b4fc7704890d83830d7f503d2d3c5f0b0521a4e28",
	"7f21526d9a0e2f70cef68007fc7f424429792c7a10025cfff6761242c0aa6954",
	"18ebdb169f70b5a97846b5534f2d92c5200a9362dc9680a53c80ff29bdef485d",

	"019c01e44f1b3533557ff98e13b51e9f821f672e410db1aa94ea0009e0474a74",
	"804bbe8b9847e5f2b79bbe85704f6db18e67c56d6d9f8e47169d7b923a27d07d",
	"e31925dcb78adc7ff9931dfb10e472096819de19e9f9ecbd9b564fdff1b2bd5e",
	"780068e4858c04e2bf6f5eb4ffa0a6b244440314f926ca3b64502590dfab442e",
	"781bbbb46597f4f8df52549bbcad7121f058ea9b9af2606fbb39f177c5734d69",
	"330839d3cdbc637542af3cb3d883a8c0df814f99f4f65a627ab55bd5c0ca9418",

	"97127912865b7a61972c2a26529afc606792395dad31ccc63e886b3baffdcd07",
	"cd613277de9313f158807e1b5a7aca4d4f88ad8acac8d562db7dfb7406d74970",
	"126b737a30166bad6b2b7067914969b83d299f94e2a0eac859dbe4fa622e6b16",
	"31dd812b5d40f895f411d7a212995dff94c4906a8904149f23293e7d3ca47e2c",
	"1ffab5efdfeb44526bf8a88c22efd5e79288c80f9271f835faa214d64060380b",
	"0a80970ae1707caadd76acd4fc3c1dc9f1f7f9db97542c3836842e27b72c6b1d",

	"d875a8b4ed944713a54541f8be96f6b5d240cb08d02ef13059ca883be716d53c",
	"795671ed1a8167aca295b1ce4f876839488a0515142b18b806671ad8820dd25a",
	"444fb1cd52b2e78f456848213770c024cc63c26dee0c95d780f393f6dd0da65a",
	"3105e2ecf53d159145db28f117415229aaf4f78249cffda30f6ae112b25de554",
	"0a612d272a741a0975064b575252da0f46d8e05f31e0d87610e7e51f376a5509",
	"06b816c9d80bf4dda119218f10d979c88147cf66dbcffc8494fde5e0d4df6314",

	"6d1136b1d187c6375e0330127b89308972cd648ca8d4c0b3b3084faef6a0ae7b",
	"4361d9f6b265d5759ba21f16b7ec7a69408273380e275540b6af742a1708c259",
	"ae8a77a19adf2dbc82c5c04d92c902d04f4a85ee4c80f64053a8843285df800a",
	"6327e44a9b9796328256748870f3f073f73746085f16a97fe053ec34e98aff10",
	"5d7fb07b9d2f53ec351bbd163b0b4aeb0685f1901b249d4d7cef4fc0c956e25d",
	"1884a71c2496a3c9f5898c31cec62c941b7a1bafb584cf8a325ca7823fd85847",

	"91cc0ccdd5d8d2b910458051630886f0673132a6b76a2c23eb645a89f99e3717",
	"9f0a91c033ac72c60d572c1f344ef9fc1e9c0172a5a7e281e6639fd7e4cedd38",
	"0e5fbc3bf39d05a7e147541e1ca171e7a15675b2373b0db070dcbaf80d5dc253",
	"1e3a2d1bb41c39867070d81067434ead00a750ba8627d0d480410e3184731c02",
	"fb58beb3611f9355ce14d9087fc2259928f58b2fd462f1d3836455810985a118",
	"24be6c7b725b3af82bdb75d71e97f105cff5cc5da34be848eefe48c9342ffd31",

	"cf817a9bd2d22183b03272ff3276e42c4634051f1315987dd9e0ae5a445ac617",
	"8ab17277cfed468c8c652b27856d24793590690b40be91a8468dc03abd7ad61c",
	"cd824cc69fb91c47a9ab94058e208eb1733629d3dc3f625435753b1fdbc47563",
	"f18abf4199b14cdb8ad2b326da607592604fd59c9704560ab188ca0dcf8f4370",
	"b244624df14cc53ff5a01bcf25e873c5c0c4dbeee33fbf55d84283f142cbf27b",
	"9c547dad9a911895a33189e6885849baba64ca7c63d9c377e6728a86d901f339",

	"5330654d131b9def2c388564f5f78b059fbb99ba7c6f06741972dcea83e75d67",
	"f0783bcf87b5924f38ab75810a8de6f588f2b222680f07b6f75988f38b1a2d2b",
	"68bdb7f049133af828eab84184a632c9aa0c0ed441584bd18547ff03ba7d6a47",
	"dc32f7696728eb88eba440d2f543b25d4ac9e8433dcab897671a9663d0e1ea08",
	"9051f8a29c4072a28075b1e4e4b44c50d2a502f9c25f05af85609311da7e5d21",
	"9ce241d02fd4ad88677bd5bb887e1d05631e672cdab4f49d390fe38d1b419202",

	"23aa1c0155bb5e9759ed2a0532db1f8e51f89dab6c1e75cc258445d53c227223",
	"e8a502c19b44324b448bf78d402e915000b607adb9b438af5fa869a687774c15",
	"e8e5546ebba258b43bc7f28b8bcb52f9e07ea7971271ee7fcb0ca5d43b824861",
	"2833cad71970a23547cd2cb00914420cef26dc28a1b950fc2d773bec352ba05c",
	"f82abc383a3aa6f4465fb4cafdbe8614a824337a86aa5ba53cefce8820772336",
	"68136ffba26af05621be27b015755821fad002b0f55e586786ebe358592cd702",

	"03148f09d7a2648c00abd3da21dcae3df8473d1b612bfb9a401ffc049b0c8354",
	"80c087a47761b790653e1d86a518115163df2857950391e8af9ca97c5a280f60",
	"3de574edf3a9f2fcbf23390226235a10a1459a298b62adabde41c2cd11396e61",
	"73c30df91adc16faf853ed42b81113caf584d6ae4eea7fa6e9ce6c951e5d130b",
	"f96be16653d7c87ca872603fd1b2fb61419c72beb843a6073d6dd0553b82cc0f",
	"f700f1d273e0d817fd8280ac57d1e502e99889713e0ebbbcbe59bf0b150fab56",

	"4e58f499625b53cc971e69a28145f03295b903fd9f3f84178fc6b8b634334414",
	"f6b39e655d6778231927cb055855d748253bb1b19b4d428096d36f66df9c4617",
	"581e8de1b87d2db07b44deca726b42fe92f124c35b401f3a1eabcd60a95de93a",
	"b5d527e23c48a98f62013d300ab15b9109d7017dbd372fc46aeca2d98009ba6e",
	"9e84aa553f5e341d425ae035a78c4f596675d1a1b6f7a778600bd6198e960e13",
	"635e4679ba93e1ff5aabe350cbfd9f5632c6175d603d3c010d90b4d62cb0bb4b",

	"8de5c537d165cb2d1e19f6935e6b0b4cd0f683434a6dbf7d7d23d3a868b26b27",
	"7db3664ddea75e39c231e003d3a2986e9030c02e796c14886a9152daa135572f",
	"7bf14a8eb8d4f292863b40d7a6fb24517cfd8412f10cd7006d0a8b5270609431",
	"df19b33b430b5b3cce3c3db41bf4726b3256c9e1ff532ad752cc9e24aa08370b",
	"4cab0b529c42814c8d479ab3f658482114ca5a91e0684741985ba945fe8c687d",
	"a281c63147dd8fb29faddfc2bc7c6f37ee37e6feec7c8acf46b7d7608267140f",

	"c024f64a2411676d42122c4b79854ef4e0dc00a22c547227e1557f236971dd5d",
	"50c8657ccabec0d17de623a2c41d34cf88cf9796fc7e76df820612b70c0cdf03",
	"fa7cdfc900600c4f1549a40deaa7bb692621e127ad23ba096175eaa1bcf44730",
	"093b873c1c294550da42e04c77521198d5e8b94e7e6c4bd6b980083f2324c067",
	"677c68c92caf8986d442a2cc3160f3906262c5e7a3467369fbaa1feba4e34033",
	"76bc5021fecd8c27bf08cf9a32ca5a7495b40d97d07d9eda73e7e221c7aaa50f",

	"9c354f93143666937ebed8402584e0a60e9ba754bc4e8ef7a7ca1dc5d38d0c7e",
	"1046ef3436f1cae1294503bd0c107a170206b9eadc25f2c6c7fdb27df5b31b74",
	"bc76ec59ae200a641878917894207f2d59f4045fec7bfdbb5a120f07b289b879",
	"30e61ad45b8c89cbec642f83cb6d7ac0b44f8687396cb39f747a084b57111f02",
	"57e17f471ef8b2d3568d59a1330f0ccbbcc15de23f43d2e6222ab7bfc58dc934",
	"7897bb18009af0f6c0bc83662dc8df9a00e26e030904fdfff13047e9648b3206",
}

// Point multiplication tests. Each group of 3 strings is: P, n, n*P
// Scalar n is not reduced modulo r.
var kat_DO255E_POINT_MUL = []string{
	"ef599ab0984eba6f1aee6c55415e9c30d1cb856b58492b645cec774b8f5d6e1c",
	"fef20582903e6fac0f8c4c4f09598d4feca23851250200204b3988dc73d3713e",
	"3fd5579d954e26abcf2cc8f5b452c3bc060eb012ecbfb81c8a061bd1ab54cc4d",

	"735b7c942443c90187b57707ef3e70cd6ebee4119e3e27a22eb27e19ba2dea23",
	"4606e325b06211bc19d2ea57e5413f7b7c3c9a2b85787af4164289f80239197d",
	"a7c405102137937c2589cdac9c73ec8afb939ef5e4508d159855d4c4cbfdc10e",

	"1d0128c75bc84f26393d0e8ed60ac65c49823aff69a03a57058f2a9f99d9fc7d",
	"4ff84fe04d1b5e9348138dbc179125c6b8dc8bfe67ffdb759cd317b794e00697",
	"40637e11cd493c67ce9664b91c1587f9199a1bffb62182c1de099c270f78da63",

	"659490c6898a71e9f7ddd11fc3d58a26f3673e4b8ffc73c62384d032393d522f",
	"beb7868f12ab5d17d83db56959f126cdde4d38d27a09dff17d64e6e917ea01b7",
	"6cdf044d38e8ee606e77be414255c4885c584d231f921bb9c04f88e5c4f88226",

	"3bc044686531f9a9475083a43fccb2186cd4b55c7a404aae1528d3912f65c028",
	"376c73d83abb5c746ae5cc34f86d2e782c109f3a1fed8fd0ad65c988e0fbca9a",
	"4f1aac7047637ac387c71bc2024c471e51639a9d9cc9e4aa046dfdd341322360",

	"62ab48031c59fa983709d9f8026c1caabc5d9576574e00e6a6a98d951b8f353b",
	"26e99a6e445c253168a96528b46aac08d94132a1780cb03c6646bc40a5aba447",
	"da26085b3f34ac6927ed9762b2ce5875e974282d8e8288b704c7149643d9a735",

	"6d3e16376b59f7ffc4a348e990e9ad2b66a0b04679b9cdd210bd59073ff55411",
	"3122b512408cb45009b26c266e2e868cc68c6bf5c91b1e292b56643eec5cb7b5",
	"800d80f6386bf5f5fcb51e7471ab05f15c486dbcc3c074c69ab0627c4845d56f",

	"b74aad34274fc6c328845899e1504400f82f632c0b79e5a72a91970a3a34f06f",
	"2e0ab2403e9cf3b909fc4f2c74715244c2a2d48f4b54ce158977bb1a5a44b58c",
	"fbb794ea32f07eece120ec7c516a136ae96f1bfe19a75393cf9c8476c7a4de2e",

	"dae2b85e5b6c6550cefb5b478cace23f7192d05d3bc4902f3fd7251f0e8fba6a",
	"a6d2fbd77623a7b52683e491d2a6d3e9c4c7bc78a96c95060faa12792c9f1a21",
	"54f67774b511a4301e18ceddfbc5d66690aaed090296a3f2bdd38fa59764e640",

	"138b69ce3c02c42384505d9e2e7c7ee9e2a45f8d825df084be48ced0b1b6fb4f",
	"c3849f4207ee75cba24a641a6ab9ebc7e613854bccbc385be35f27ea84123803",
	"9cb32db76b1073f0e6338288070f3090ceb1e59ac24a716840aaf7f408217a1e",

	"ac85820a0e9659376f241adbbab2da3fe797983096e7ee475a9ba2e49e74964d",
	"87c16396c42af29948cd1594f2396e59c45778ed36c739da0166c9cca2de2bd6",
	"bec4a7126f927e98807316f948fed476c2cd6ce131a8e8d1521e55cb1bf34f22",

	"51b4a3e99927cbc8bbb8584abafc0d1e915d2e5f0bf93617522593d5029ab935",
	"011916d017f7eeef2c95493d7b040f2fb5fa093fe0538294c854e74d861e5771",
	"38e5eff913b154c4b536de092f1734d59fdf0ae0f8e81af9fadd307b17259f59",

	"899c66330df8d48a7da39517f0babe7298fd18f79165f65218297258d2094879",
	"44c5603c739facc3a392c1277b39044e27c027eb2499e44796f38245e1f5affd",
	"4c6ec4dc4882256af0e362f849c3ebf5f79ff79bb79626d8516b301efbae247a",

	"e342ded1ed0f07137d7b2992e7fa4ea8e08bf45d0811c5e1d613e47938d5d21f",
	"d51a875780f9335ffb5a7aba312312388ba951bbe3980b8c79e390c0dcabb8c4",
	"dd72b3efd123911595d40909d5f27ed6b5b40e7677855f834366ad545f1fb34c",

	"9db6e81ed4fcf00e7b4a13509a976b25eb674b44ce6d47645cea7aa7f77b8e27",
	"03560cf6003a1a5bf34148c081a44ad345e6c3639c355b8ec63e396c77546fbc",
	"78f32da15f9ebe01b8d028eb3f4879204e1a141059d8d85b755947aa3b369073",

	"d8e11751381ad5cdccd6815ee6bdac4e386113dd24d2e17a48626a4280fcc94b",
	"a7b8933f86798c4bbba9589c100033ea2d0954802d1cc4d68f9282ea84bffb0e",
	"c961bcd739f530db91cbc83becc997c82e716dfe9ddf30c2946dd7eb57f62728",

	"6df18f24e7e180d23dcbac1f785023b8a11e55a2068438a745dc4301db1db301",
	"c843a3be8090195d566690801760a39626b12af2d398db96d24701197814cc7b",
	"46ed7c43da0dee7a80226e6e67f408f38a526c22e7593ad6e5933e4700af4954",

	"91d114745fa07c9ea408c091b0590ad2c921b83acba95270c948955c86619c55",
	"ff0f743577a47f9f73832ebf0e9b0d6dff89b208fc4b60e1ba18316af401b61d",
	"004397186d5503ed6eec0d6204c4462a12d590d2226d504c83add12bac59c414",

	"847610b28dfc4d72d52fdfc4e343e4508e65356319505868b643cb368775f07d",
	"ca761a610323e709e3014c5ec22272ce78eb4d5fafa0bd01b5d61aea7f3ffdc3",
	"ac114c9a50893d4b333cb05487bbf92b3c16228aacacbf51a41ea786b2f0ed25",

	"e1e1e70147d942f1ff27e837032a9bf34350209f4b562c6753e76039558ae728",
	"dc6b2d0b39ea749c36fed438e2e6f2d54f87f318a0ac80e8beb68c72f707945b",
	"0b1e250899cac9ea52e2837fbe1541a9a8d8ad5a29e7dca67cf1e9aa2b846f5e",
}

// Monte-Carlo test vectors for Do255e
// Starting point P[0] = (2**120)*G.
// Then, for i = 1...10000 (inclusive): P[i] = int(P[i-1].w)*P[i-1]
// Below are P[0], P[1000], P[2000],... P[10000]
var kat_DO255E_MC_POINT_MUL = []string{
	"fae3f207045cdc8e1a73b743fbfd81b9d23bec1392603a524e4e53677ac8646d",

	"7c13b0849e3e5e26756e471cd2e8dd06c4d96e24c2edb547d39c9813579d0670",
	"bb140aaa8912e7379464635b220241b7eecd1cb21ac9ae63b9d6f2f2c076945a",
	"8a73932369250c2ced3c2c123558f150cf7b67320f18ea3e0a5271cb12dfca51",
	"a15c7e3bc340143bfa773b01530a034dfb0e9182d2bf5df120fbb8442e12ad7c",
	"4a2f85396ee4ef8b18c6fb5ac7ba7cdf4abb63e0084a675aa23ace64424c2f78",
	"11b2b72661412f1825b3c045906858345932542e04aaf916f5d0d65978883722",
	"6374912bf9223322bc6bab34b11440c779f1696f8d5031a4a09d51b6b4ce001d",
	"b5cd6da3d2465fa2ec920e09be4ff44857d80f62949ffa4bcdad91b60bfc8224",
	"f187c40e8f29bd7e825be662377d39f50e5d34e22a27ea910f5c8ac311ce1f71",
	"91bafd2379679082436e7c7ef2459f06e967225cea6ea37491b1ac901fab6f68",
}
