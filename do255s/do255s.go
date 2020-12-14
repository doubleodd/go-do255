package do255s

import (
	gf "github.com/doubleodd/go-do255/internal/field"
	"github.com/doubleodd/go-do255/internal/scalar"
)

// This file implements operations on curve points for do255s
// (specifically, on elements of the prime order group defined over
// do255s). It also provides function for computing with scalar
// values (integers modulo the group order).
//
// API: a point is represented in memory by a Do255sPoint structure.
// Contents of such a structure are opaque. These structures are
// mutable; the various functions such as setAdd() modify the point
// on which they are called. It is always acceptable to also use the
// destination structure as one of the operands. All such functions
// return a pointer to the structure on which they were called, so
// that calls may be syntactically chained.
//
// A point can be encoded to, and decoded from, a sequence of 32 bytes.
// Encoding is unique and verified. Decoding an invalid sequence of
// bytes yields an error. The group neutral element can be encoded (as
// a sequence of 32 bytes of value 0x00). Whether the neutral element
// is acceptable or not when decoding depends on usage context; the
// decoding functions returns a specific flag value when the point is
// the neutral, that the caller may test explicitly.
//
// A scalar is similarly represented in memory by a Do255sScalar
// mutable structure. Scalars can be encoded to, and decoded from,
// 32 bytes. Encoding is canonical: a given scalar value can be encoded
// in a unique way, and only one sequence of bytes may decode to that
// scalar.
//
// Unless explcitly documented, all functions here are constant-time.

// Do255sPoint is the type for a do255s point.
//
// Default value for a point structure is not valid. The NewDo255sPoint()
// function makes sure to return only initialized structures. If allocating
// a point structure manually, make sure to properly set it to a valid point
// before using it as source.
type Do255sPoint struct {
	// Internally, we use fractional (x,u) coordinates, which have
	// complete and efficient formulas.
	x, z, u, t gf.GF255s
}

// Preallocated neutral point. Do not modify.
var do255sNeutral = Do255sPoint{
	x: gf.GF255s{0, 0, 0, 0},
	z: gf.GF255s{1, 0, 0, 0},
	u: gf.GF255s{0, 0, 0, 0},
	t: gf.GF255s{1, 0, 0, 0},
}

// Preallocated conventional generator point. Do not modify.
var do255sGenerator = Do255sPoint{
	x: gf.GF255s{
		0x4803AC7D33B156B1, 0x3EF832265840B591,
		0x213759ECCB010B9D, 0x39BD72651783FB6D},
	z: gf.GF255s{1, 0, 0, 0},
	u: gf.GF255s{3, 0, 0, 0},
	t: gf.GF255s{1, 0, 0, 0},
}

// Create a new point. The point is set to the group neutral element (N).
func NewDo255sPoint() *Do255sPoint {
	P := new(Do255sPoint)
	*P = do255sNeutral
	return P
}

// Set the point P to the neutral element (N).
// A pointer to this structure is returned.
func (P *Do255sPoint) Neutral() *Do255sPoint {
	*P = do255sNeutral
	return P
}

// Set the point P to the conventional generator (G).
// A pointer to this structure is returned.
func (P *Do255sPoint) Generator() *Do255sPoint {
	*P = do255sGenerator
	return P
}

// Encode a point into exactly 32 bytes. The bytes are appended to the
// provided slice; the new slice is returned. The extension is done in
// place if the provided slice has enough capacity.
func (P *Do255sPoint) Encode(dst []byte) []byte {
	// Encoded value is the representation of w = 1/u. If u == 0,
	// we encode 0.
	var w gf.GF255s
	w.Inv(&P.u).Mul(&w, &P.t)
	return w.Encode(dst)
}

// Encode the square of the w coordinate of a point into exactly 32 bytes.
// The bytes are appended to the provided slice; the new slice is returned.
// The extension is done in place if the provided slice has enough capacity.
// This function is meant to support ECDH.
func (P *Do255sPoint) EncodeSquaredW(dst []byte) []byte {
	// Encoded value is the representation of w^2 = 1/u^2. If u == 0,
	// we encode 0.
	var w gf.GF255s
	w.Inv(&P.u).Mul(&w, &P.t).Sqr(&w)
	return w.Encode(dst)
}

// Encode a point into exactly 32 bytes.
func (P *Do255sPoint) Bytes() [32]byte {
	var d [32]byte
	P.Encode(d[:0])
	return d
}

// Test whether a given chunk of 32 bytes is a valid representation of
// do255s group element. This is faster than actually decoding it.
// Returned value is:
//    1   valid encoding of a non-neutral group element
//    0   valid encoding of the neutral point N
//   -1   invalid encoding
func Do255sCheckPoint(src []byte) int {
	var d gf.GF255s
	r := d.Decode(src)
	zz := d.IsZero()

	// d <- sqrt((w^2 - a)^2 - 4*b)   (with a = -1 and b = 1/2 for do255s)
	// Note: d can never be 0, because b is not a quadratic residue.
	d.Sqr(&d)
	d.Add(&d, &gf.GF255s_ONE)
	d.Sqr(&d)
	d.Sub(&d, &gf.GF255s_TWO)

	// Encoding is a valid non-neutral element if and only if
	// d is a quadratic residue. Since d != 0, this can only be 1 or -1.
	qr := d.Legendre()

	// If r == 0 then input is not a valid field element encoding
	// and we must return -1.
	// If r == 1:
	//    If zz == 1, then the point is the neutral. In that case,
	//    we have qr == -1 because a^2-4*b is not a quadratic residue.
	//    If zz == 0, then qr == 1 for valid points, -1 for invalid.
	r = (r - 1) + (-r & (qr + zz))
	return int(int64(r))
}

// Decode a point from exactly 32 bytes. Returned value is 1 if the
// point could be successfully decoded into a non-neutral group element,
// 0 if it could be successfully decoded as the neutral point N, or -1
// if it could not be decoded. If the decoding was not successful, then
// the destination structure is set to the neutral N.
//
// This function is constant-time with regard to the decoded value
// and also with regard to the validity status (timing-based side
// channels do not leak whether the value was found to be a valid point).
//
// Returned value is:
//    1   valid encoding of a non-neutral group element
//    0   valid encoding of the neutral point N
//   -1   invalid encoding
func (P *Do255sPoint) Decode(src []byte) int {
	var w, x, d gf.GF255s

	// Decode field element and test for zero. We want zz == 1 if
	// and only if the field element decoding worked AND the value
	// is zero.
	r := w.Decode(src)
	zz := r & w.IsZero()

	// x <- w^2 - a   (with a = -1 for do255s)
	x.Sqr(&w)
	x.Add(&x, &gf.GF255s_ONE)

	// d <- sqrt((w^2 - a)^2 - 4*b)   (with b = 1/2 for do255s)
	d.Sqr(&x)
	d.Sub(&d, &gf.GF255s_TWO)
	r &= d.Sqrt(&d)

	// x <- ((w^2 - a) + d)/2
	// We cannot get x == 0 here.
	x.Add(&x, &d)
	x.Half(&x)

	// If x is a square, then we must use the other solution, i.e.
	// ((w^2 - a) - d)/2, which we obtain by subtracting d.
	qr := x.Legendre()
	d.Select(&d, &gf.GF255s_ZERO, (qr+1)>>1)
	x.Sub(&x, &d)

	// If decoding failed, or if the value was 0, then r == 0 at
	// this point. In that case, we want to set the point to the
	// neutral (0:1:0:1). Otherwise, we set it to (x:1:1:w)
	// (because u = 1/w).
	P.x.Select(&x, &gf.GF255s_ZERO, r)
	P.z.Set(&gf.GF255s_ONE)
	P.u.Select(&gf.GF255s_ONE, &gf.GF255s_ZERO, r)
	P.t.Select(&w, &gf.GF255s_ONE, r)

	// If the point was the neutral, then r == 0 and zz == 1.
	// Otherwise, zz == 0, and r == 0 or 1, depending on point
	// validity.
	return int(int64((zz - 1) & ((r << 1) - 1)))
}

// Test whether a point is the neutral element N.
// Returned value is 1 for the neutral, 0 otherwise.
func (P *Do255sPoint) IsNeutral() int {
	return int(P.u.IsZero())
}

// Test whether this structure (P) represents the same point as the
// provided other structure (Q).
// Returned value is true if both points are the same, false otherwise.
func (P *Do255sPoint) Equal(Q *Do255sPoint) int {
	var t1, t2 gf.GF255s
	t1.Mul(&P.u, &Q.t)
	t2.Mul(&P.t, &Q.u)
	return int(t1.Eq(&t2))
}

// Copy a point structure into another.
// A pointer to this structure is returned.
func (P *Do255sPoint) Set(Q *Do255sPoint) *Do255sPoint {
	P.x.Set(&Q.x)
	P.z.Set(&Q.z)
	P.u.Set(&Q.u)
	P.t.Set(&Q.t)
	return P
}

// If ctl == 1, then copy point Q1 into P.
// If ctl == 0, then copy point Q2 into P.
// ctl MUST be 0 or 1. This is a constant-time selection primitive.
func (P *Do255sPoint) Select(P1, P2 *Do255sPoint, ctl int) {
	P.x.Select(&P1.x, &P2.x, uint64(ctl))
	P.z.Select(&P1.z, &P2.z, uint64(ctl))
	P.u.Select(&P1.u, &P2.u, uint64(ctl))
	P.t.Select(&P1.t, &P2.t, uint64(ctl))
}

// Set this point to the sum of the two provided points.
// A pointer to this structure (P) is returned.
func (P *Do255sPoint) Add(P1, P2 *Do255sPoint) *Do255sPoint {
	var t1, t2, t3, t4, t5, t6, t7, t8, t9, t10 gf.GF255s

	// t1 <- X1*X2
	// t2 <- Z1*Z2
	// t3 <- U1*U2
	// t4 <- T1*T2
	t1.Mul(&P1.x, &P2.x)
	t2.Mul(&P1.z, &P2.z)
	t3.Mul(&P1.u, &P2.u)
	t4.Mul(&P1.t, &P2.t)

	// t5 <- (X1 + Z1)*(X2 + Z2) - t1 - t2 = X1*Z2 + X2*Z1
	t5.Add(&P1.x, &P1.z)
	t8.Add(&P2.x, &P2.z)
	t5.Mul(&t5, &t8)
	t5.Sub(&t5, &t1)
	t5.Sub(&t5, &t2)

	// t6 <- (U1 + T1)*(U2 + T2) - t3 - t4 = U1*T2 + U2*T1
	t6.Add(&P1.u, &P1.t)
	t9.Add(&P2.u, &P2.t)
	t6.Mul(&t6, &t9)
	t6.Sub(&t6, &t3)
	t6.Sub(&t6, &t4)

	// t7 <- t1 + b*t2  (with b = 1/2)
	t7.Half(&t2)
	t7.Add(&t1, &t7)

	// t8 <- t4*t7
	t8.Mul(&t4, &t7)

	// t9 <- t3*(2*b*t5 + a*t7) = t3*(t5 - t7)  (since a = -1 and b = 1/2)
	t9.Sub(&t5, &t7)
	t9.Mul(&t3, &t9)

	// t10 <- (t4 + alpha*t3)*(t5 + t7)  (with alpha = 1/2)
	t5.Add(&t5, &t7)
	t3.Half(&t3)
	t3.Add(&t3, &t4)
	t10.Mul(&t3, &t5)

	// U3 <- -t6*(t1 - b*t2) = t6*((1/2)*t2 - t1)
	t2.Half(&t2)
	t1.Sub(&t2, &t1)
	P.u.Mul(&t6, &t1)

	// Z3 <- t8 - t9
	P.z.Sub(&t8, &t9)

	// T3 <- t8 + t9
	P.t.Add(&t8, &t9)

	// X3 <- b*(t10 - t8 + beta*t9)  (with b = 1/2 and beta = -3/2)
	t10.Sub(&t10, &t8)
	t10.Sub(&t10, &t9)
	t9.Half(&t9)
	t10.Sub(&t10, &t9)
	P.x.Half(&t10)

	return P
}

// Set this point to the difference of the two provided points (P1 - P2).
// A pointer to this structure (P) is returned.
func (P *Do255sPoint) Sub(P1, P2 *Do255sPoint) *Do255sPoint {
	var P2n Do255sPoint
	P2n.x = P2.x
	P2n.z = P2.z
	P2n.u.Neg(&P2.u)
	P2n.t = P2.t
	return P.Add(P1, &P2n)
}

// Internal type for a point in affine (x, u) coordinates. This is used
// to speed up some computations but we do not make it public so that
// the API remains simple.
type do255sPointAffine struct {
	x, u gf.GF255s
}

// Set this point to the sum of the two provided points, the second of
// which being in affine (x, u) coordinates.
// A pointer to this structure (P) is returned.
func (P *Do255sPoint) addMixed(P1 *Do255sPoint, P2 *do255sPointAffine) *Do255sPoint {
	var t1, t3, t5, t6, t7, t8, t9, t10 gf.GF255s

	// t1 <- X1*X2
	// t2 <- Z1*Z2 = Z1
	// t3 <- U1*U2
	// t4 <- T1*T2 = T1
	t1.Mul(&P1.x, &P2.x)
	t3.Mul(&P1.u, &P2.u)

	// t5 <- X1*Z2 + X2*Z1 = X1 + X2*Z1
	t5.Mul(&P1.z, &P2.x)
	t5.Add(&t5, &P1.x)

	// t6 <- U1*T2 + U2*T1 = U1 + U2*T1
	t6.Mul(&P1.t, &P2.u)
	t6.Add(&t6, &P1.u)

	// t7 <- t1 + b*t2  (with b = 1/2 and t2 = Z1)
	t7.Half(&P1.z)
	t7.Add(&t7, &t1)

	// t8 <- t4*t7 = T1*t7
	t8.Mul(&P1.t, &t7)

	// t9 <- t3*(2*b*t5 + a*t7) = t3*(t5 - t7)  (since a = -1 and b = 1/2)
	t9.Sub(&t5, &t7)
	t9.Mul(&t3, &t9)

	// t10 <- (t4 + alpha*t3)*(t5 + t7)  (with alpha = 1/2 and t4 = T1)
	t5.Add(&t5, &t7)
	t3.Half(&t3)
	t3.Add(&t3, &P1.t)
	t10.Mul(&t3, &t5)

	// U3 <- -t6*(t1 - b*t2) = t6*((1/2)*Z1 - t1)
	t5.Half(&P1.z)
	t1.Sub(&t5, &t1)
	P.u.Mul(&t6, &t1)

	// Z3 <- t8 - t9
	P.z.Sub(&t8, &t9)

	// T3 <- t8 + t9
	P.t.Add(&t8, &t9)

	// X3 <- b*(t10 - t8 + beta*t9)  (with b = 1/2 and beta = -3/2)
	t10.Sub(&t10, &t8)
	t10.Sub(&t10, &t9)
	t9.Half(&t9)
	t10.Sub(&t10, &t9)
	P.x.Half(&t10)

	return P
}

// Set this point to the difference of the two provided points, the second of
// which being in affine (x, u) coordinates.
// A pointer to this structure (P) is returned.
func (P *Do255sPoint) subMixed(P1 *Do255sPoint, P2 *do255sPointAffine) *Do255sPoint {
	var P2n do255sPointAffine
	P2n.x = P2.x
	P2n.u.Neg(&P2.u)
	return P.addMixed(P1, &P2n)
}

// Set this point (P) to the double of the provided point Q.
// A pointer to this structure (P) is returned.
func (P *Do255sPoint) Double(Q *Do255sPoint) *Do255sPoint {
	// 3M+6S formulas

	// X' = (a^2-4*b)*X*Z
	// Z' = X^2 + a*X*Z + b*Z^2
	var t1, t2, xp, zp gf.GF255s
	t1.Sqr(&Q.x)
	t2.Sqr(&Q.z)
	xp.Add(&Q.x, &Q.z).Sqr(&xp).Sub(&t1, &xp).Add(&xp, &t2).Half(&xp)
	zp.Half(&t2).Add(&zp, &t1).Add(&zp, &xp)

	// t1 = X^2
	// t2 = Z^2

	// X'' = 4*b*X'*Z'
	// Z'' = X'^2 - 2*a*X'*Z' + (a^2-4*b)*Z'^2
	var t3, t4, t5 gf.GF255s
	t3.Sqr(&xp)
	t4.Sqr(&zp)
	t5.Add(&xp, &zp).Sqr(&t5).Sub(&t5, &t3).Sub(&t5, &t4)
	P.x.Set(&t5)
	P.z.Add(&t3, &t5).Sub(&P.z, &t4)

	// t3 = X'^2
	// t4 = -(a^2-4*b)*Z'^2

	// U'' = 2*(a^2-4*b)*(X^2 - b*Z^2)*Z'*U
	// T'' = (X'^2 - (a^2-4*b)*Z'^2)*T
	zp.Mul(&zp, &Q.u)
	t1.Lsh(&t1, 1).Sub(&t2, &t1)
	P.u.Mul(&zp, &t1)
	t3.Add(&t3, &t4)
	P.t.Mul(&Q.t, &t3)

	return P
}

// Set this point (P) to (2^n)*Q (i.e. perform n successive doublings).
// This function is constant-time with regard to the point values, but
// not to the number of doublings (n); computation time is proportional
// to n.
// A pointer to this structure (P) is returned.
func (P *Do255sPoint) DoubleX(Q *Do255sPoint, n uint) *Do255sPoint {
	var tX, tW, tZ, t1, t2 gf.GF255s

	if n == 0 {
		P.Set(Q)
		return P
	}
	if n == 1 {
		P.Double(Q)
		return P
	}

	// First half-doubling, combined with conversion from fractional
	// (x,u) to Jacobian (x,w); output in E[r](-2*a,a^2-4*b).
	//   X' = Z^2*T^4
	//   W' = Z*T^2 - (2*X + a*Z)*U^2
	//   Z' = Z*U*T
	// Note that a = -1 for curve do255s.
	// Cost: 4M+2S
	t2.Sqr(&Q.u)       // t2 <- U^2
	t1.Mul(&Q.z, &Q.t) // t1 <- Z*T
	tW.Lsh(&Q.x, 1)    // tW <- 2*X
	tW.Sub(&tW, &Q.z)  // tW <- 2*X + a*Z
	tX.Mul(&t1, &Q.t)  // tX <- Z*T^2
	tW.Mul(&tW, &t2)   // tW <- (2*X + a*Z)*U^2
	tZ.Mul(&t1, &Q.u)  // tZ <- Z*U*T
	tW.Sub(&tX, &tW)   // tW <- Z*T^2 - (2*X + a*Z)*U^2
	tX.Sqr(&tX)        // tX <- Z^2*T^4

	// Second half-doubling, with output in Jacobian (x, w)
	//   X = 16*b*Z'^4
	//   W = 2*X' - 2*a*Z'^2 - W'2
	//   Z = 2*W'*Z'
	// With a = -1 and b = 1/2 for curve do255s.
	// Cost: 4S
	t1.Sqr(&tZ)      // t1 <- Z'^2
	tZ.Add(&tW, &tZ) // tZ <- W' + Z'
	t2.Sqr(&tW)      // t2 <- W'^2
	tW.Add(&tX, &t1) // tW <- X' - a*Z'^2
	tZ.Sqr(&tZ)      // tZ <- (W' + Z')^2
	tW.Lsh(&tW, 1)   // tW <- 2*X' - 2*a*Z'^2
	tX.Sqr(&t1)      // tX <- Z'^4
	tW.Sub(&tW, &t2) // tW <- 2*X' - 2*a*Z'^2 - W'^2
	t2.Add(&t2, &t1) // t2 <- W'^2 + Z'^2
	tX.Lsh(&tX, 3)   // tX <- 16*b*Z'^4
	tZ.Sub(&tZ, &t2) // tZ <- 2*W'*Z'

	// n-2 inner doublings in Jacobian (x, w) coordinates.
	for n -= 2; n > 0; n-- {
		// Formulas:
		//   X' = 16*b*(W*Z)^4
		//   W' = -(W^4 - (a^2-4*b)*Z^4)
		//   Z' = 2*W*Z*(2*X + a*Z^2 - W^2)
		// With a = -1 and b = 1/2 for curve do255s; in particular:
		//   W^4 - (a^2-4*b)*Z^4 = W^4 + Z^4
		//                       = (W^2 - a*Z^2)^2 - 2*(W*Z)^2
		// Cost: 2M+4S
		t1.Mul(&tW, &tZ) // t1 <- W*Z
		tW.Add(&tW, &tZ) // tW <- W + Z
		tZ.Lsh(&tX, 1)   // tZ <- 2*X
		t2.Sqr(&t1)      // t2 <- (W*Z)^2
		t1.Lsh(&t1, 1)   // t1 <- 2*W*Z
		tW.Sqr(&tW)      // tW <- (W + Z)^2
		t2.Lsh(&t2, 1)   // t2 <- 2*(W*Z)^2
		tW.Sub(&tW, &t1) // tW <- W^2 - a*Z^2
		tZ.Sub(&tZ, &tW) // tZ <- 2*X + a*Z^2 - W^2
		tW.Sqr(&tW)      // tW <- (W^2 - a*Z^2)^2
		tZ.Mul(&tZ, &t1) // tZ <- 2*W*Z*(2*X + a*Z^2 - W^2)
		tW.Sub(&t2, &tW) // tW <- -(W^4 - (a^2-4*b)*Z^4)
		tX.Sqr(&t2)      // tX <- 4*(W*Z)^4
		tX.Lsh(&tX, 1)   // tX <- 16*b*(W*Z)^4
	}

	// Final doubling with conversion back into fractional coordinates.
	//   X' = 4*b*(W*Z)^2
	//   Z' = (2*X + a*Z^2 - W^2)^2
	//   U' = 2*W*Z*(2*X + a*Z^2 - W^2)
	//   T' = -(W^4 - (a^2-4*b)*Z^4)
	// Cost: 2M+4S
	t1.Mul(&tW, &tZ)   // t1 <- W*Z
	tW.Add(&tW, &tZ)   // tW <- W + Z
	tZ.Lsh(&tX, 1)     // tZ <- 2*X
	t2.Sqr(&t1)        // t2 <- (W*Z)^2
	t1.Lsh(&t1, 1)     // t1 <- 2*W*Z
	tW.Sqr(&tW)        // tW <- (W + Z)^2
	P.x.Lsh(&t2, 1)    // X' <- 4*b*(W*Z)^2
	tW.Sub(&tW, &t1)   // tW <- W^2 - a*Z^2
	tZ.Sub(&tZ, &tW)   // tZ <- 2*X + a*Z^2 - W^2
	tW.Sqr(&tW)        // tW <- W^4 - (a^2-4*b)*Z^4 + 2*(W*Z)^2
	P.u.Mul(&t1, &tZ)  // U' <- 2*W*Z*(2*X + a*Z^2 - W^2)
	P.z.Sqr(&tZ)       // Z' <- (2*X + a*Z^2 - W^2)^2
	P.t.Sub(&P.x, &tW) // T' <- -(W^4 - (a^2-4*b)*Z^4)

	return P
}

// Set P to the opposite of point Q.
// A pointer to this structure (P) is returned.
func (P *Do255sPoint) Neg(Q *Do255sPoint) *Do255sPoint {
	P.x.Set(&Q.x)
	P.z.Set(&Q.z)
	P.u.Neg(&Q.u)
	P.t.Set(&Q.t)
	return P
}

// Multiply a point Q by a given scalar n.
// A pointer to this structure (P) is returned.
func (P *Do255sPoint) Mul(Q *Do255sPoint, n *Do255sScalar) *Do255sPoint {
	// Recode input scalar in 52 signed digits in the -15..+16 range.
	var sd [52]byte
	n.recode5(&sd)

	// Initialize a window with i*Q for i in 1..16
	var win [16]Do255sPoint
	win[0] = *Q
	win[1].Double(Q)
	for i := 3; i <= 15; i += 2 {
		win[i-1].Add(&win[i-2], Q)
		win[i].Double(&win[((i+1)>>1)-1])
	}

	// Use the top digit to initialize the accumulator point M.
	P.do255sLookupWindow(&win, uint(sd[51]))

	// Process other digits from top to bottom.
	for i := 50; i >= 0; i-- {
		P.DoubleX(P, 5)
		var M Do255sPoint
		M.do255sLookupWindow(&win, uint(sd[i]&0x1F))
		M.u.CondNeg(&M.u, uint64(sd[i]>>7))
		P.Add(P, &M)
	}

	return P
}

// Multiply the conventional generator by a given scalar n. This is
// functionally equivalent (but faster) to P.Generator().Mul(&P, n).
// A pointer to this structure (P) is returned.
func (P *Do255sPoint) MulGen(n *Do255sScalar) *Do255sPoint {
	// Recode input scalar into 5-bit Booth encoding.
	var sd [52]byte
	n.recode5(&sd)

	// Lookup initial accumulator by using the top digit (which is
	// guaranteed nonnegative).
	var Ma do255sPointAffine
	Ma.do255sLookupWindowAffine(&do255sWin_G195_xu, uint(sd[51]))
	P.x = Ma.x
	P.z = gf.GF255s_ONE
	P.u = Ma.u
	P.t = gf.GF255s_ONE

	// Add points corresponding to top digits of the three other
	// quarter-scalars.
	Ma.do255sLookupWindowAffine(&do255sWin_G_xu, uint(sd[12]&0x1F))
	Ma.u.CondNeg(&Ma.u, uint64(sd[12]>>7))
	P.addMixed(P, &Ma)
	Ma.do255sLookupWindowAffine(&do255sWin_G65_xu, uint(sd[25]&0x1F))
	Ma.u.CondNeg(&Ma.u, uint64(sd[25]>>7))
	P.addMixed(P, &Ma)
	Ma.do255sLookupWindowAffine(&do255sWin_G130_xu, uint(sd[38]&0x1F))
	Ma.u.CondNeg(&Ma.u, uint64(sd[38]>>7))
	P.addMixed(P, &Ma)

	// Process all other digits from high to low. We process the
	// four quarter-scalars in parallel.
	for i := 11; i >= 0; i-- {
		P.DoubleX(P, 5)

		Ma.do255sLookupWindowAffine(&do255sWin_G_xu, uint(sd[i]&0x1F))
		Ma.u.CondNeg(&Ma.u, uint64(sd[i]>>7))
		P.addMixed(P, &Ma)
		Ma.do255sLookupWindowAffine(&do255sWin_G65_xu, uint(sd[i+13]&0x1F))
		Ma.u.CondNeg(&Ma.u, uint64(sd[i+13]>>7))
		P.addMixed(P, &Ma)
		Ma.do255sLookupWindowAffine(&do255sWin_G130_xu, uint(sd[i+26]&0x1F))
		Ma.u.CondNeg(&Ma.u, uint64(sd[i+26]>>7))
		P.addMixed(P, &Ma)
		Ma.do255sLookupWindowAffine(&do255sWin_G195_xu, uint(sd[i+39]&0x1F))
		Ma.u.CondNeg(&Ma.u, uint64(sd[i+39]>>7))
		P.addMixed(P, &Ma)
	}

	return P
}

// Constant-time lookup of a point in a window. Provided window has
// 16 elements. Input offset ('index') is in the 0..16 range. This
// function sets P to a copy of win[index - 1] if index != 0, or to
// the neutral if index == 0.
func (P *Do255sPoint) do255sLookupWindow(win *[16]Do255sPoint, index uint) {
	// Initialize P to all-zeros.
	P.x = gf.GF255s_ZERO
	P.z = gf.GF255s_ZERO
	P.u = gf.GF255s_ZERO
	P.t = gf.GF255s_ZERO

	// Lookup all values.
	for i := 0; i < 16; i++ {
		m := int64(index) - int64(i+1)
		mm := ^uint64((m | -m) >> 63)
		P.x.CondOrFrom(&win[i].x, mm)
		P.z.CondOrFrom(&win[i].z, mm)
		P.u.CondOrFrom(&win[i].u, mm)
		P.t.CondOrFrom(&win[i].t, mm)
	}

	// Set Z and T to 1 if the index is zero.
	mz := uint64((int64(index) - 1) >> 63)
	P.z.CondOrFrom(&gf.GF255s_ONE, mz)
	P.t.CondOrFrom(&gf.GF255s_ONE, mz)
}

// Constant-time lookup of a point in a window. This is similar to
// do255sLookupWindow(), except that this function works on points in
// affine (x, u) coordinates.
func (P *do255sPointAffine) do255sLookupWindowAffine(win *[16]do255sPointAffine, index uint) {
	// Initialize P to all-zeros (which is the valid representation
	// of the neutral element).
	P.x = gf.GF255s_ZERO
	P.u = gf.GF255s_ZERO

	// Lookup all values.
	for i := 0; i < 16; i++ {
		m := int64(index) - int64(i+1)
		mm := ^uint64((m | -m) >> 63)
		P.x.CondOrFrom(&win[i].x, mm)
		P.u.CondOrFrom(&win[i].u, mm)
	}
}

// Add to point P a point from a window, given an encoded index.
// THIS IS NOT CONSTANT-TIME.
func (P *Do255sPoint) addFromWindowVartime(win *[16]Do255sPoint, ej byte) {
	j := int(ej & 0x1F)
	if j != 0 {
		if ej < 0x80 {
			P.Add(P, &win[j-1])
		} else {
			P.Sub(P, &win[j-1])
		}
	}
}

// Add to point P a point from a window (affine), given an encoded index.
// THIS IS NOT CONSTANT-TIME.
func (P *Do255sPoint) addFromWindowAffineVartime(win *[16]do255sPointAffine, ej byte) {
	j := int(ej & 0x1F)
	if j != 0 {
		if ej < 0x80 {
			P.addMixed(P, &win[j-1])
		} else {
			P.subMixed(P, &win[j-1])
		}
	}
}

// Check whether k0*G + k1*P (with G being the conventional generator)
// yields a point which would encode to the specified sequence of bytes
// encR. This function is meant to support signature verification.
// IT IS NOT CONSTANT-TIME; thus, it should be used only on public
// elements (which is normally the case when verifying signatures).
// Returned value is true on match, false otherwise.
func (P *Do255sPoint) VerifyHelperVartime(k0, k1 *Do255sScalar, encR []byte) bool {
	// Decode encR into a point. If this fails, we can stop right away.
	var R Do255sPoint
	if R.Decode(encR) < 0 {
		return false
	}

	// Edge case: if P == N or k1 == 0, then we can use MulGen().
	if P.IsNeutral() != 0 || k1.IsZero() != 0 {
		var M Do255sPoint
		M.MulGen(k0)
		return M.Equal(&R) != 0
	}

	// We want to check that k0*G + k1*P = R. We know that k1 != 0.
	// To speed things up, we first find shorter values c0 and c1
	// such that k1 = c0/c1 mod r; the equation is then equivalent to:
	//   (k0*c1)*G + c0*P - c1*R = N
	// with c0 and c1 being half-size (at most 128 bits in absolute
	// value).
	// We negate points conditionally to ensure that c0 and c1 are
	// nonnegative.
	var winP, winR [16]Do255sPoint
	var c0, c1 [2]uint64
	negc0, negc1 := k1.ReduceBasisVartime(&c0, &c1)
	if negc0 {
		winP[0].Neg(P)
	} else {
		winP[0] = *P
	}
	if negc1 {
		winR[0] = R
	} else {
		winR[0].Neg(&R)
	}
	winP[1].Double(&winP[0])
	for i := 3; i <= 15; i += 2 {
		winP[i-1].Add(&winP[i-2], &winP[0])
		winP[i].Double(&winP[((i+1)>>1)-1])
	}
	winR[1].Double(&winR[0])
	for i := 3; i <= 15; i += 2 {
		winR[i-1].Add(&winR[i-2], &winR[0])
		winR[i].Double(&winR[((i+1)>>1)-1])
	}

	// Recode scalars. We need k0*c1 mod r. Do not forget the sign
	// of c1 (array c1[] contains |c1|, sign is in negc1).
	var sd0 [52]byte
	var kt Do255sScalar
	kt[0] = c1[0]
	kt[1] = c1[1]
	kt.Mul(&kt, k0)
	if negc1 {
		kt.Neg(&kt)
	}
	kt.recode5(&sd0)
	var sdP, sdR [26]byte
	scalar.Recode5Small(&sdP, &c0)
	scalar.Recode5Small(&sdR, &c1)

	// Initialize accumulator using the top digit of sdP (which is
	// non-negative) and add points from the other windows.
	var M Do255sPoint
	if sdP[25] != 0 {
		M = winP[int(sdP[25])-1]
	} else {
		M = do255sNeutral
	}
	M.addFromWindowVartime(&winR, sdR[25])
	M.addFromWindowAffineVartime(&do255sWin_G_xu, sd0[25])
	M.addFromWindowAffineVartime(&do255sWin_G130_xu, sd0[51])

	// Process all other digits.
	for i := 24; i >= 0; i-- {
		M.DoubleX(&M, 5)
		M.addFromWindowVartime(&winP, sdP[i])
		M.addFromWindowVartime(&winR, sdR[i])
		M.addFromWindowAffineVartime(&do255sWin_G_xu, sd0[i])
		M.addFromWindowAffineVartime(&do255sWin_G130_xu, sd0[26+i])
	}

	// Signature is valid if and only if the final output is the
	// neutral point.
	return M.IsNeutral() != 0
}

// Alternate helper function for signature verification; this returns
// 1 if k0*G + k1*P yields a point whose encoding is exactly equal to
// the first 32 bytes of encR; otherwise, it returns 0. This function
// is usually slower than VerifyHelperVartime(), but it is constant-time.
func (P *Do255sPoint) VerifyHelper(k0, k1 *Do255sScalar, encR []byte) int {
	var P1, P2 Do255sPoint
	P1.MulGen(k0)
	P2.Mul(P, k1)
	P1.Add(&P1, &P2)
	var encT [32]byte
	P1.Encode(encT[:0])
	var tt int = 0
	for i := 0; i < 32; i++ {
		tt |= int(encR[i] ^ encT[i])
	}
	return 1 - ((tt + 0xFF) >> 8)
}

// Map a sequence of bytes into a curve element. The mapping is not
// injective or surjective, and not uniform among possible outputs;
// however, any given point has only a limited number of possible
// pre-images by the map. A hash-to-curve process can be built on top
// of this map, as follows:
//  - Hash some input data in 64 bytes, with a secure hash function or
//    XOF (e.g. SHAKE).
//  - Split these 64 bytes into two halves, and map each of them to a
//    point with this map.
//  - Add the two points together.
func (P *Do255sPoint) MapBytes(bb []byte) *Do255sPoint {
	// Decode source bytes as a field element. This applies modular
	// reduction.
	var e gf.GF255s
	e.DecodeReduce(bb)

	// Formulas: we map into a point on the dual curve
	// y^2 = x*(x^2 + aa*x + bb) with aa = -2*a and bb = a^2 - 4*b.
	// We use a fixed non-square constant d.
	//
	// Then, the two candidates for x are:
	//   x1 = -aa / (1 + d*e^2)
	//   x2 = -x1 - aa
	//      = -aa*d*e^2 / (1 + d*e^2)
	// The candidates for y^2 are then yy1num/(yden^2) and
	// yy2num/(yden^2), with:
	//   yy1num = -aa*bb*d^3*e^6 + (aa^3-3*aa*bb)*d^2*e^4
	//            + (aa^3-3*aa*bb)*d*e^2 - aa*bb
	//   yy2num = -aa*bb*d^4*e^8 + (aa^3-3*aa*bb)*d^3*e^6
	//            + (aa^3-3*aa*bb)*d^2*e^4 - aa*bb*d*e^2
	//          = yy1num*d*e^2
	//   yden = (1 + d*e^2)^2
	// If yy1num is a square, then we use y = sqrt(yy1num) / yden.
	// Otherwise, we use y = -sqrt(yy2num) / yden.
	//
	// Once (x,y) is obtained, we apply the theta_{1/2} isogeny to
	// get a point in the right group.
	//
	// If 1 + d*e^2 == 0, then we map to the neutral N.
	//
	// For do255s:
	//   a = -1   b = 1/2   aa = 2   bb = -1   d = -1
	//   x1num  = -2
	//   x2num  = 2*e^2
	//   xden   = 1 - e^2
	//   yy1num = -2*e^6 + 14*e^4 - 14*e^2 + 2
	//   yy2num =  2*e^8 - 14*e^6 + 14*e^4 - 2*e^2
	//   yden   = (1 - e^2)^2
	// Since we want u = x/y, we have:
	//   u = U/T = (xnum*yden) / (xden*ynum)
	//           = (xnum*xden) / ynum
	// since yden = xden^2.

	var e2, e4, e6 gf.GF255s
	e2.Sqr(&e)
	e4.Sqr(&e2)
	e6.Mul(&e2, &e4)

	// yy1num and yy2num
	var yy1num, yy2num, tt1, tt2 gf.GF255s
	tt1.Sub(&e4, &e2)
	tt2.Lsh(&tt1, 3).Sub(&tt2, &tt1)
	tt1.Sub(&tt2, &e6).Add(&tt1, &gf.GF255s_ONE)
	yy1num.Lsh(&tt1, 1)
	yy2num.Mul(&yy1num, &e2).Neg(&yy2num)

	// Assemble the point coordinates; this is still on the dual
	// curve. If e != +/- 1, then it cannot happen that yy1num == 0.
	// However, it is possible that yy2num == 0, but only if e == 0;
	// in that case, yy1num == 2, which is not a quadratic residue,
	// and yy2num is selected.
	ls := yy1num.Legendre()
	qr1 := 1 - (ls >> 63)
	var X, Z, U, T gf.GF255s
	X.Neg(&gf.GF255s_ONE).Select(&X, &e2, qr1).Lsh(&X, 1)
	Z.Sub(&gf.GF255s_ONE, &e2)
	U.Mul(&X, &Z)
	T.Select(&yy1num, &yy2num, qr1).Sqrt(&T)
	T.CondNeg(&T, 1-qr1)

	// If e == +/- 1, then (X:Z:U:T) == (0:0:0:0)
	// If e == 0, then (X:Z:U:T) == (0:1:0:0)
	// Otherwise, Z != 0 and T != 0, and the value is correct.

	// Apply the isogeny theta_{1/2} to get a point in the proper
	// group on the right curve:
	//   x' = 4*b*u^2
	//   u' = 2*x/(u*(x^2 - bb))
	P.x.Sqr(&U).Lsh(&P.x, 1)
	P.z.Sqr(&T)
	P.u.Mul(&X, &Z).Mul(&P.u, &T).Lsh(&P.u, 1)
	Z.Sqr(&Z)
	P.t.Sqr(&X).Add(&P.t, &Z).Mul(&P.t, &U)

	// If the source value e was +/-1, then we should map to the
	// neutral N. Similarly, if e == 0, then we should obtain N.
	// In both cases, all four coordinates have been set to 0,
	// so we must simply detect that situation and set Z and T to 1.
	nz := P.z.IsZero()
	P.z.Select(&gf.GF255s_ONE, &P.z, nz)
	P.t.Select(&gf.GF255s_ONE, &P.t, nz)

	return P
}

// =====================================================================
// Precomputed windows in affine (x, u) coordinates:
//   i*G         for i = 1..16
//   i*2^65*G    for i = 1..16
//   i*2^130*G   for i = 1..16
//   i*2^195*G   for i = 1..16

var do255sWin_G_xu = [16]do255sPointAffine{
	/* 1 */
	{
		x: gf.GF255s{0x4803AC7D33B156B1, 0x3EF832265840B591,
			0x213759ECCB010B9D, 0x39BD72651783FB6D},
		u: gf.GF255s{0x0000000000000003, 0x0000000000000000,
			0x0000000000000000, 0x0000000000000000},
	},
	/* 2 */
	{
		x: gf.GF255s{0x8C6318C6318C5721, 0x18C6318C6318C631,
			0x318C6318C6318C63, 0x6318C6318C6318C6},
		u: gf.GF255s{0xB3E22F8D0D1657FC, 0xE5427944219E490E,
			0x256C2BE75887FD30, 0x6F44EC749C5869ED},
	},
	/* 3 */
	{
		x: gf.GF255s{0xC92EA6A45AF542C8, 0xACA2B1F0EB2EBB62,
			0x16AA9AB4772D49BF, 0x586E82D468993FDF},
		u: gf.GF255s{0x7204233F36F06441, 0x36F07204233F36F0,
			0x233F36F07204233F, 0x7204233F36F07204},
	},
	/* 4 */
	{
		x: gf.GF255s{0xB21E437C0F24B43F, 0xE352D46A2A191529,
			0x5105F27F8E691D57, 0x203E742185ABD1BC},
		u: gf.GF255s{0x9204A59E69223E39, 0x4E645F874B12D8E7,
			0xF2A1145C9F5D3475, 0x54E64904662F1657},
	},
	/* 5 */
	{
		x: gf.GF255s{0xD4C604388FAAFD54, 0x478144CACA27D2A0,
			0x4B4B7554A3381C38, 0x78CC3E71A8C94117},
		u: gf.GF255s{0xDF0337C00667B64D, 0x58856B292FBA673A,
			0xC17C3E9333A6D7CE, 0x52932B9A9F0CC65D},
	},
	/* 6 */
	{
		x: gf.GF255s{0x0F8CCA9A5D8AEC75, 0x2337B2C5F00CFF35,
			0x4C38E25760F3A520, 0x5B8B19C1B38EC25A},
		u: gf.GF255s{0x2378FCE7659F8304, 0xF71199A7E92DA598,
			0x05D59CECEE1D7E76, 0x7B55CE1D18E39726},
	},
	/* 7 */
	{
		x: gf.GF255s{0xE45744799916AE56, 0xAB719188444F94C0,
			0x8BCF3B53F4727AB0, 0x239A8F014119EE4B},
		u: gf.GF255s{0xBB70A3099712F248, 0x62AE8CABB5C7CED6,
			0x150E35D2C3D060D0, 0x6EB76B2647D95DA4},
	},
	/* 8 */
	{
		x: gf.GF255s{0x486D0418D4089F3F, 0x415C30647E21EC49,
			0xDB2F693C6C060489, 0x475F06F55E1CF582},
		u: gf.GF255s{0x2DE04D6F93F6EEA0, 0x13A2A2DE4BD9AB93,
			0x4BA8485FACF9CC03, 0x26DCD74FEC331ED1},
	},
	/* 9 */
	{
		x: gf.GF255s{0x6F20F664F35EBE0B, 0x72C56CBCBBE70025,
			0x1B8A888E7A1BC6B3, 0x3AE1B93B4536C750},
		u: gf.GF255s{0xF293B01428B2AEF7, 0x6544BF679F64AADF,
			0xC25DB7DB3DC0038B, 0x641BC3EAE8B16348},
	},
	/* 10 */
	{
		x: gf.GF255s{0x0B77022F05D906A1, 0x36553987A9763928,
			0x90C0593CD952E579, 0x601889E506ECA253},
		u: gf.GF255s{0xAC2450CD35E0E33E, 0xA5FBAFCA4BE098FF,
			0x15C5CFFBF6CC3634, 0x11EB002E4DD0CC8A},
	},
	/* 11 */
	{
		x: gf.GF255s{0xFC96E106BF390F8E, 0x87F4FE12D17AC934,
			0x20DD4A1C455D2425, 0x0F489EA0CF96C239},
		u: gf.GF255s{0x4845EE2F9FDD13DE, 0x343EAC1284D73CA8,
			0xD7761D7533C3A9ED, 0x44204FCFD6AF1B44},
	},
	/* 12 */
	{
		x: gf.GF255s{0x116E4C0180B4DDCF, 0xD965D1AE088C285B,
			0xCE5A14DBF0E43F69, 0x2BCC8378CC7BF398},
		u: gf.GF255s{0xEDFD8B4BAD8160D1, 0x9BBB9D365457C4A7,
			0xC2A3F9623BA7F9DB, 0x1EA722FF296F8D92},
	},
	/* 13 */
	{
		x: gf.GF255s{0xBB832D11959A60E5, 0xF3AF5A147660C5E6,
			0x1037C50FE1FAFA73, 0x251739F159A18E7B},
		u: gf.GF255s{0x4DEC233876049D71, 0x12E5410DA8EA7661,
			0xA8C80C9E1E78A79B, 0x56187141071FC1D5},
	},
	/* 14 */
	{
		x: gf.GF255s{0x238B24170CD243C2, 0xF20B2474CC7B0C22,
			0xD1AF972D83B667A6, 0x1BC20D835A45EA4E},
		u: gf.GF255s{0x0693069387390957, 0xC3E0CD34567E3C9A,
			0x2EDB82538E0CEEAE, 0x60EF0D216D3123EB},
	},
	/* 15 */
	{
		x: gf.GF255s{0x725096F4FAC9FD25, 0x4072199CA8F661BC,
			0x1E3D756A54E4B489, 0x36A82923B7C5D81C},
		u: gf.GF255s{0x3010093957EA5D5A, 0x99E3DE25AE95E20B,
			0x5594A4E0449B9EFD, 0x100615BCCFB9A5A5},
	},
	/* 16 */
	{
		x: gf.GF255s{0x6ACF66F7D508BB00, 0xFC3A8030127F8808,
			0x70481A4C6305A579, 0x30664839BD04AA0E},
		u: gf.GF255s{0xFE9991C8BEE08226, 0x70631241B6132ED0,
			0xF17032C1930CC3DF, 0x1880F6E65E61EDCA},
	},
}

var do255sWin_G65_xu = [16]do255sPointAffine{
	/* 1 */
	{
		x: gf.GF255s{0xE3AD260DF7B2D34D, 0xEF4C98E564E44AEC,
			0x41B749AF25234018, 0x0722EF995D77E307},
		u: gf.GF255s{0xC504DF70301870EC, 0xC3F5757559BEB30B,
			0x59EF9CAA0E041627, 0x3AA6C1241FDF29EA},
	},
	/* 2 */
	{
		x: gf.GF255s{0x2ADA4E59EECB8F77, 0xD577FCEA354A3FA3,
			0xE3E9AA4AC41ABDB8, 0x0FA4D0D00D3E7390},
		u: gf.GF255s{0xE395C342A49AB090, 0xBEEAEA12D994BBBA,
			0xE555C4B2316756DE, 0x0C4B1FDF2A2835EA},
	},
	/* 3 */
	{
		x: gf.GF255s{0x83437B3B44A76584, 0x85337A5B5EA37019,
			0x0AF40A3CBACBCA6E, 0x6129FB017ECE3AE0},
		u: gf.GF255s{0xDD2736E064D47F35, 0x8235DE8C4CCCB7AD,
			0x84993635937FCD8B, 0x05DADA9DFEACCC60},
	},
	/* 4 */
	{
		x: gf.GF255s{0x417B5BE5A8BA33DF, 0x1146B966C729E5BE,
			0xF19DEF75C9837F83, 0x4B20562C0C99B2FC},
		u: gf.GF255s{0x20CCC8942E0DD4AF, 0xB05A144640E6146F,
			0x39FB3501747B5584, 0x637DD68BAA0DCC5D},
	},
	/* 5 */
	{
		x: gf.GF255s{0x683D5AAFCFBCBF77, 0xB0B9F38491FC818E,
			0x8DDC79CDFA3B9EB3, 0x406A80209ECB5A2E},
		u: gf.GF255s{0x8E3257DDE0B61BA4, 0x034912BE87333A80,
			0x355389746FE3860D, 0x0430DD4DB72FDAE3},
	},
	/* 6 */
	{
		x: gf.GF255s{0xE7DCEDC9830F3BCC, 0x3190BB8BF13D83C2,
			0xB348EAB92E361E5B, 0x21A00B1A4B6B4CCB},
		u: gf.GF255s{0x41C6037317995DE4, 0xEB9F7DAEFD0E5177,
			0x3C89F25AEBEFEA23, 0x6BC51220E2A6FEA1},
	},
	/* 7 */
	{
		x: gf.GF255s{0x118E74003C58146D, 0x4E183A0D59D45A86,
			0x88CF58B85D8DD179, 0x1FCC2FD8D59414BA},
		u: gf.GF255s{0x6E0C6EC00B45CB36, 0xDE22614DA09D7B7D,
			0x68532DB386C5311B, 0x6C76A36690BF3721},
	},
	/* 8 */
	{
		x: gf.GF255s{0xAB758A963FCB67D3, 0xBDBFF8BDCDDC7179,
			0x63A4D2686CFABF15, 0x63CC5473DD3F5036},
		u: gf.GF255s{0xCF5B38B898027EEE, 0x2E0E5F286C358DFB,
			0x77F9F481B54C83A0, 0x320289A4D4AD38D8},
	},
	/* 9 */
	{
		x: gf.GF255s{0x9AFED28D850C23A5, 0x4E28F84F5436536D,
			0x493B602C85203364, 0x0A73A163C1F86147},
		u: gf.GF255s{0xE9F73ACE30D70239, 0x2D6EB28373FFD2AB,
			0xE362CB1B14264D84, 0x6A58C4C447F58EDA},
	},
	/* 10 */
	{
		x: gf.GF255s{0x6D73DF10948A18DB, 0xCE6F829C7226AFF4,
			0x27E81B57F40362F2, 0x37BFD38BEC948817},
		u: gf.GF255s{0xBD4ABD3A1D23F58F, 0xFDF2BE7E13D824B7,
			0xFDF4EBAD6DA9D7C5, 0x6F6BF0FD0A760CE6},
	},
	/* 11 */
	{
		x: gf.GF255s{0x511B8BB2C326919C, 0x36D06F274827CD94,
			0xFA874F760A089D6F, 0x796CE68B71B33782},
		u: gf.GF255s{0xAC82BC835BEF5D82, 0x1ED56E4F26768EEA,
			0x7CB1DF78AE0F6520, 0x3F27081B1A793D69},
	},
	/* 12 */
	{
		x: gf.GF255s{0x5C39D1FD89DD2FE7, 0x56D5F374BD997B36,
			0xCBC4F497125C9FBD, 0x2EF79073F850836A},
		u: gf.GF255s{0x3C2CD45E00488C1B, 0x744209C1D5AA1E98,
			0xB0A6FCE83F628FFF, 0x770D858C1138F9B1},
	},
	/* 13 */
	{
		x: gf.GF255s{0x77ACC4B81316F86D, 0x69DD38B0F257E592,
			0x0F58657B866FB0D8, 0x3C39314BB1A846EE},
		u: gf.GF255s{0x2B0FECC1159B7C7E, 0x5176B2243AD75522,
			0xF578275E6C3BD8ED, 0x3EA45CDA4A4335B5},
	},
	/* 14 */
	{
		x: gf.GF255s{0xCE74D3864A444AF9, 0x0C0F7567E10119D1,
			0x1A156CF7026B694A, 0x33D0491E8BCF686E},
		u: gf.GF255s{0xDCABCD1C2B792081, 0xC0B470A297317B12,
			0xB2C96A085D26A03D, 0x52A2BF2E7D7202E9},
	},
	/* 15 */
	{
		x: gf.GF255s{0x5D0726CD57A859A3, 0x39BB1C2D84086050,
			0x25E41A2638CFB70C, 0x455E456EBF8F2E0B},
		u: gf.GF255s{0xCB4E8E8682ED1F56, 0x95E7B3DAE6F438C6,
			0x436574C311E4EAE5, 0x5B29A9D46B6557D1},
	},
	/* 16 */
	{
		x: gf.GF255s{0x492AFE097556579D, 0x13B540A021705957,
			0xFDFA5CC01C264FF7, 0x1053932D99A6BBCD},
		u: gf.GF255s{0xA0A662EF07947CAE, 0xCFD4B9B4407ADF90,
			0x9698A42E4CDA26D4, 0x00D98818FA1C5F43},
	},
}

var do255sWin_G130_xu = [16]do255sPointAffine{
	/* 1 */
	{
		x: gf.GF255s{0xD9AB2AB4313E150A, 0xAF3680D8923D8F48,
			0x8FE12E041CE8DFAF, 0x254561A8AEB11CAB},
		u: gf.GF255s{0xC17AB82D247C18A0, 0x95542A3E6973F13C,
			0xB14CDFC79E957BD2, 0x661229C7BADE4F32},
	},
	/* 2 */
	{
		x: gf.GF255s{0x9C35A5637F6C8053, 0x5706387FD6021602,
			0x0AF3C1470C1697F6, 0x7847E5A5A7420D24},
		u: gf.GF255s{0x72C27DF187129578, 0xF443CA1C94A3CAEB,
			0xEA24368F35C3A22D, 0x17A619CA7283DADF},
	},
	/* 3 */
	{
		x: gf.GF255s{0xFCB2125561D1D07F, 0x75909948DD63D9AC,
			0x423F60D1CF853002, 0x7CA49E6CC1F4B430},
		u: gf.GF255s{0x31D729A12CB400ED, 0x7A195520CEEAA6A9,
			0x86F51D53192BE4B6, 0x22A23BF61CDCF306},
	},
	/* 4 */
	{
		x: gf.GF255s{0xF4E0408E636E5ABA, 0x9CD7416D582A44FA,
			0x942DC6F5FEA2C3C6, 0x157D068AE6F5A315},
		u: gf.GF255s{0x767FEA8208A28FA3, 0xA2DC84A037DCA7EE,
			0xEBD2AC8233E1DE6F, 0x188F708521AD21CA},
	},
	/* 5 */
	{
		x: gf.GF255s{0x5AC5803793462B8A, 0x1260408E42FAD501,
			0x04E6353CE57B6118, 0x2EC28A315514EE0D},
		u: gf.GF255s{0xB820181764875239, 0x5B4DD1BB342A7459,
			0xC66B48D80BE03B56, 0x7F08F18F636EA208},
	},
	/* 6 */
	{
		x: gf.GF255s{0xE3AB9795E746B94F, 0x838E9EFFD0E5D4BB,
			0xBF8C188E36E0B9CB, 0x22F1F02E91674AD5},
		u: gf.GF255s{0x9FFA28B92E2ECC30, 0xFB480B5678DA5C00,
			0xE1C94F9AB14804D7, 0x7D453CC948C9FA1D},
	},
	/* 7 */
	{
		x: gf.GF255s{0xC6E25DCD6722A7DB, 0x41174D150E1112E4,
			0x3E4484A32C15C12C, 0x58FBCE703CBFCC0E},
		u: gf.GF255s{0x8784811955F5DEF5, 0x48D2662133271AC2,
			0xAF3F2A3FA6C81237, 0x1A086A6FE01B2E7D},
	},
	/* 8 */
	{
		x: gf.GF255s{0x76867E220E1C7387, 0x6876C29D05E6E1C0,
			0xBFD4B391E87110E1, 0x613A84AA71E6E57F},
		u: gf.GF255s{0x1DC9A2A859DBD92B, 0xBF15C84C6245BBB3,
			0x2A8D65D0B2CAB143, 0x3DE43A8017A5525E},
	},
	/* 9 */
	{
		x: gf.GF255s{0x897F189BF41FA232, 0x31DC1F6609C0B63D,
			0x8F33869EB9BECA07, 0x6314FEF02F399003},
		u: gf.GF255s{0xDB52105D891C8CF1, 0xA1E811D5F01D372F,
			0x5CF867DE0ADED951, 0x17052800E0B4DA8B},
	},
	/* 10 */
	{
		x: gf.GF255s{0x2D9B581DC5B282B0, 0xDA006B6A93A610DE,
			0xFECD738894842EB5, 0x209F00E9867FED68},
		u: gf.GF255s{0x126C847746496B31, 0x23D37316626487D7,
			0x46505E02982B39EB, 0x234345752143582E},
	},
	/* 11 */
	{
		x: gf.GF255s{0x011BE82B2CE8661D, 0xB271BF01FFE6A835,
			0x58CA42537E335CBA, 0x55DE7CA1C28EC9D8},
		u: gf.GF255s{0xD51EDD52BECE2651, 0x792F3E2ED6FA4957,
			0x8349ED268E6F750E, 0x6BB94CDFEF6F6D88},
	},
	/* 12 */
	{
		x: gf.GF255s{0x38933FA3CD39261E, 0x3C3AF9A02BA7FB63,
			0x28F0D9086ECFFB4C, 0x7337FC83C49CEF8B},
		u: gf.GF255s{0x89340045DE2129A6, 0x10CF7923E46000DB,
			0xA17D1BEB2C69ABE2, 0x1F3592DE203C2DA4},
	},
	/* 13 */
	{
		x: gf.GF255s{0x4CD4A2285435BD02, 0x44E5A2A8E9FCB255,
			0x80DD71B6AD674422, 0x1B1EC39B8DDC6546},
		u: gf.GF255s{0x0F4767AD57355C90, 0xCAD900D819124EA8,
			0xB4B045E5E702318A, 0x7AB6CB353A7E1058},
	},
	/* 14 */
	{
		x: gf.GF255s{0xC42905109B1A192C, 0x8D33854826E55283,
			0xAA6FD4F2FA9558A9, 0x35F24A95D568423A},
		u: gf.GF255s{0xE62B26D301A95CF8, 0x84BC9EEC48398A95,
			0x08F76F6E267875A4, 0x32D3B1A49A0A50A1},
	},
	/* 15 */
	{
		x: gf.GF255s{0x63F85CD360842C43, 0xB125F90D87643CE9,
			0xC9FCE0470F0D206F, 0x1C7E6FBE705B59B8},
		u: gf.GF255s{0x9A0B63FF58819185, 0xF789BF221086D125,
			0xAD341ED8B1E1776A, 0x082617E146EBF733},
	},
	/* 16 */
	{
		x: gf.GF255s{0xB0C3DE6BEBEAEC95, 0x4D695AB6C0E9A3E4,
			0x0CA501DA9EAD2ACA, 0x7402B1938309B1CC},
		u: gf.GF255s{0xF50D8973C1409964, 0x7E7D988FECBBDB1D,
			0x64BAB1680DE07B83, 0x588C179EE277A32E},
	},
}

var do255sWin_G195_xu = [16]do255sPointAffine{
	/* 1 */
	{
		x: gf.GF255s{0xA5B332AE8B607110, 0xC6A1F5888C5C2A7B,
			0x4686BCEF5F176F78, 0x2EEA2C710771051D},
		u: gf.GF255s{0xD5386C024F6950E1, 0xB46F6A19F2C711D7,
			0xD13612F77F90F84F, 0x02B1FE40C35B1108},
	},
	/* 2 */
	{
		x: gf.GF255s{0xD08B3673DD36C5D6, 0x42875443DFF50CA5,
			0x4472B3B100AD977A, 0x71BBE2F8AAAD796C},
		u: gf.GF255s{0x4D20522801D9BE69, 0x491F6D9530005684,
			0x0394879561348B5B, 0x207807A3C7DA4C9E},
	},
	/* 3 */
	{
		x: gf.GF255s{0x91F7BAAD9ECAFD0C, 0xDF7832A52AEBB19D,
			0x72D70D6D53A736B8, 0x4053B80041C9765F},
		u: gf.GF255s{0x62329BF642399CCB, 0xD15830A754BDD8E5,
			0x052C83B745F76715, 0x2C2FEC4277A075E5},
	},
	/* 4 */
	{
		x: gf.GF255s{0x8CFBE62352C0321C, 0x35317478163A7506,
			0xB07FB5CF2F0D350C, 0x674624CC5B63BEAB},
		u: gf.GF255s{0x33F1C936823F9793, 0x45EB4CD78646E9E2,
			0xB28C3C8D49110514, 0x20B293F739F3AA65},
	},
	/* 5 */
	{
		x: gf.GF255s{0x0EFD79AE51CFE627, 0x8D6F3A49D66034F7,
			0x4D2C2C01A87347FE, 0x380D6B9FEB08A3F5},
		u: gf.GF255s{0x3F09E8FE4D4F0DAA, 0xE47EE2E2B553CA96,
			0x9FFA3C3AAC578533, 0x356B0675A08113C1},
	},
	/* 6 */
	{
		x: gf.GF255s{0xD179D52E0F60ECF2, 0x98C18C7C19CE5942,
			0xAC43B21539169363, 0x373829F1D540D16F},
		u: gf.GF255s{0xA143C257C8E7B387, 0x62D889136CC108D4,
			0xE32F52E0BB72837C, 0x076CCDEE12F67001},
	},
	/* 7 */
	{
		x: gf.GF255s{0x88F2BF91C052AA4C, 0xBD6C284630F6316B,
			0xA916918E68DB91DB, 0x29BAD78DBC302139},
		u: gf.GF255s{0x971F761E8AEDB101, 0xD10EB0D73B802862,
			0x100B543BDE5AD3F6, 0x2D2D355F27445B8E},
	},
	/* 8 */
	{
		x: gf.GF255s{0xEFF0E12764994499, 0xE0F075AC08336324,
			0xEE9F6C3330900A58, 0x12CFF18E3B68E16E},
		u: gf.GF255s{0x6AC9953979CB6EA2, 0xB754075F43904E3A,
			0x0F5C5D706C51FC54, 0x72F0811D6A945237},
	},
	/* 9 */
	{
		x: gf.GF255s{0x3B36DDD303A56217, 0x56B840E1E5485F22,
			0xE35974691ED67A31, 0x14D0AD4AA2B7708D},
		u: gf.GF255s{0xF2B5C06FCAFEB583, 0x2D0EFBE6A1A2C97B,
			0x466067CE1442F46D, 0x13349ABFC999DD67},
	},
	/* 10 */
	{
		x: gf.GF255s{0x1487FE1D0AAD152B, 0xD389875D9B57CB58,
			0x00086A583427927B, 0x10FD1E1939F33683},
		u: gf.GF255s{0x1F87C746453A8A8A, 0xDDCF1D005DA27586,
			0xC77308C7C7664BF3, 0x5774724E6436D94B},
	},
	/* 11 */
	{
		x: gf.GF255s{0xBA92EB23B8AB9BE1, 0x369244277FDE0909,
			0x603EBECE9890D5AD, 0x640E212BBA557898},
		u: gf.GF255s{0x08A8CB0C9987EBED, 0x05428BBDC4D7168A,
			0x5FDD560F0207E6AA, 0x4CD603F72497FC6D},
	},
	/* 12 */
	{
		x: gf.GF255s{0x0B9FCE4299AFD948, 0x7CA4BCD2EAA2E0EC,
			0xE3FA30C5AC388BA7, 0x606DBE6798462842},
		u: gf.GF255s{0xB4754DE27F72D537, 0xC0F01E802AC5B1C3,
			0x6492E0F04B64FC09, 0x14C0F255E1CDCAD1},
	},
	/* 13 */
	{
		x: gf.GF255s{0x8C25E45F7D80CA97, 0xCC8CB15146620FFB,
			0xA4ECCCCE54B3AA7D, 0x5725DEC68B4D2C1D},
		u: gf.GF255s{0xCF158917F580DF1C, 0x0E14A880D3E90C45,
			0x60D780AFF882B1C7, 0x420C0614ED165A4B},
	},
	/* 14 */
	{
		x: gf.GF255s{0x277645546DA86DED, 0x0306A489B4168BAF,
			0x70DC32CE22762D79, 0x49450B2E315499A1},
		u: gf.GF255s{0xE3DF626AA182F658, 0x0AEB0A55D410C765,
			0x645809F39DD03BA1, 0x40C860215AE057DB},
	},
	/* 15 */
	{
		x: gf.GF255s{0x7BBF7B434F5F8ED2, 0x41BAB11E135E5808,
			0x666BECEDBDB00CD1, 0x5E6B345CAD1A0A15},
		u: gf.GF255s{0xBA4D3B438436B829, 0xFD932027B0906610,
			0x43EFF53756D1C683, 0x0C078789FF8301C5},
	},
	/* 16 */
	{
		x: gf.GF255s{0x3AC75F22F410E610, 0xC5752B476401C98F,
			0x1B6F4E6E58028B0B, 0x794FB7F20AF770CB},
		u: gf.GF255s{0x2C94910518A991C5, 0xAAC1F71E18F00C14,
			0x2E1D217E9768A122, 0x50A69209BB0B3E0B},
	},
}
