package do255e

import (
	gf "github.com/doubleodd/go-do255/internal/field"
	"github.com/doubleodd/go-do255/internal/scalar"
)

// This file implements operations on curve points for do255e
// (specifically, on elements of the prime order group defined over
// do255e).
//
// API: a point is represented in memory by a Do255ePoint structure.
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
// Unless explicitly documented, all functions here are constant-time.

// Do255ePoint is the type for a do255e point.
//
// Default value for a point structure is not valid. The NewDo255ePoint()
// function makes sure to return only initialized structures. If allocating
// a point structure manually, make sure to properly set it to a valid point
// before using it as source.
type Do255ePoint struct {
	// Internally, we use fractional (x,u) coordinates, which have
	// complete and efficient formulas.
	x, z, u, t gf.GF255e
}

// Preallocated neutral point. Do not modify.
var do255eNeutral = Do255ePoint{
	x: gf.GF255e{0, 0, 0, 0},
	z: gf.GF255e{1, 0, 0, 0},
	u: gf.GF255e{0, 0, 0, 0},
	t: gf.GF255e{1, 0, 0, 0},
}

// Preallocated conventional generator point. Do not modify.
var do255eGenerator = Do255ePoint{
	x: gf.GF255e{2, 0, 0, 0},
	z: gf.GF255e{1, 0, 0, 0},
	u: gf.GF255e{1, 0, 0, 0},
	t: gf.GF255e{1, 0, 0, 0},
}

// Create a new point. The point is set to the group neutral element (N).
func NewDo255ePoint() *Do255ePoint {
	P := new(Do255ePoint)
	*P = do255eNeutral
	return P
}

// Set the point P to the neutral element (N).
// A pointer to this structure is returned.
func (P *Do255ePoint) Neutral() *Do255ePoint {
	*P = do255eNeutral
	return P
}

// Set the point P to the conventional generator (G).
// A pointer to this structure is returned.
func (P *Do255ePoint) Generator() *Do255ePoint {
	*P = do255eGenerator
	return P
}

// Encode a point into exactly 32 bytes. The bytes are appended to the
// provided slice; the new slice is returned. The extension is done in
// place if the provided slice has enough capacity.
func (P *Do255ePoint) Encode(dst []byte) []byte {
	// Encoded value is the representation of w = 1/u. If u == 0,
	// we encode 0.
	var w gf.GF255e
	w.Inv(&P.u).Mul(&w, &P.t)
	return w.Encode(dst)
}

// Encode a point into exactly 32 bytes.
func (P *Do255ePoint) Bytes() [32]byte {
	var d [32]byte
	P.Encode(d[:0])
	return d
}

// Encode the square of the w coordinate of a point into exactly 32 bytes.
// The bytes are appended to the provided slice; the new slice is returned.
// The extension is done in place if the provided slice has enough capacity.
// This function is meant to support ECDH.
func (P *Do255ePoint) EncodeSquaredW(dst []byte) []byte {
	// Encoded value is the representation of w^2 = 1/u^2. If u == 0,
	// we encode 0.
	var w gf.GF255e
	w.Inv(&P.u).Mul(&w, &P.t).Sqr(&w)
	return w.Encode(dst)
}

// Test whether a given chunk of 32 bytes is a valid representation of
// a do255e group element. This is faster than actually decoding it.
// Returned value is:
//    1   valid encoding of a non-neutral group element
//    0   valid encoding of the neutral point N
//   -1   invalid encoding
func Do255eCheckPoint(src []byte) int {
	var d gf.GF255e
	r := d.Decode(src)
	zz := d.IsZero()

	// d <- sqrt((w^2 - a)^2 - 4*b)   (with a = 0 and b = -2 for do255e)
	// Note: d can never be 0, because b is not a quadratic residue.
	d.SqrX(&d, 2)
	d.Add(&d, &gf.GF255e_EIGHT)

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
func (P *Do255ePoint) Decode(src []byte) int {
	var w, x, d gf.GF255e

	// Decode field element and test for zero. We want zz == 1 if
	// and only if the field element decoding worked AND the value
	// is zero.
	r := w.Decode(src)
	zz := r & w.IsZero()

	// x <- w^2 - a   (with a = 0 for do255e)
	x.Sqr(&w)

	// d <- sqrt((w^2 - a)^2 - 4*b)   (with b = -2 for do255e)
	d.Sqr(&x)
	d.Add(&d, &gf.GF255e_EIGHT)
	r &= d.Sqrt(&d)

	// x <- ((w^2 - a) + d)/2
	// We cannot get x == 0 here.
	x.Add(&x, &d)
	x.Half(&x)

	// If x is a square, then we must use the other solution, i.e.
	// ((w^2 - a) - d)/2, which we obtain by subtracting d.
	qr := x.Legendre()
	d.Select(&d, &gf.GF255e_ZERO, (qr+1)>>1)
	x.Sub(&x, &d)

	// If decoding failed, or if the value was 0, then r == 0 at
	// this point. In that case, we want to set the point to the
	// neutral (0:1:0:1). Otherwise, we set it to (x:1:1:w)
	// (because u = 1/w).
	P.x.Select(&x, &gf.GF255e_ZERO, r)
	P.z.Set(&gf.GF255e_ONE)
	P.u.Select(&gf.GF255e_ONE, &gf.GF255e_ZERO, r)
	P.t.Select(&w, &gf.GF255e_ONE, r)

	// If the point was the neutral, then r == 0 and zz == 1.
	// Otherwise, zz == 0, and r == 0 or 1, depending on point
	// validity.
	return int(int64((zz - 1) & ((r << 1) - 1)))
}

// Test whether a point is the neutral element N.
// Returned value is 1 for the neutral, 0 otherwise.
func (P *Do255ePoint) IsNeutral() int {
	return int(P.u.IsZero())
}

// Test whether this structure (P) represents the same point as the
// provided other structure (Q).
// Returned value is 1 if both points are the same, 0 otherwise.
func (P *Do255ePoint) Equal(Q *Do255ePoint) int {
	var t1, t2 gf.GF255e
	t1.Mul(&P.u, &Q.t)
	t2.Mul(&P.t, &Q.u)
	return int(t1.Eq(&t2))
}

// Copy a point structure into another.
// A pointer to this structure is returned.
func (P *Do255ePoint) Set(Q *Do255ePoint) *Do255ePoint {
	P.x.Set(&Q.x)
	P.z.Set(&Q.z)
	P.u.Set(&Q.u)
	P.t.Set(&Q.t)
	return P
}

// If ctl == 1, then copy point Q1 into P.
// If ctl == 0, then copy point Q2 into P.
// ctl MUST be 0 or 1. This is a constant-time selection primitive.
func (P *Do255ePoint) Select(P1, P2 *Do255ePoint, ctl int) {
	P.x.Select(&P1.x, &P2.x, uint64(ctl))
	P.z.Select(&P1.z, &P2.z, uint64(ctl))
	P.u.Select(&P1.u, &P2.u, uint64(ctl))
	P.t.Select(&P1.t, &P2.t, uint64(ctl))
}

// Set this point to the sum of the two provided points.
// A pointer to this structure (P) is returned.
func (P *Do255ePoint) Add(P1, P2 *Do255ePoint) *Do255ePoint {
	var t1, t2, t3, t4, t5, t6, t7, t8, t9, t10 gf.GF255e

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

	// t7 <- t1 + b*t2  (with b = -2)
	t7.Sub(&t1, &t2)
	t7.Sub(&t7, &t2)

	// t8 <- t4*t7
	t8.Mul(&t4, &t7)

	// t9 <- t3*(2*b*t5 + a*t7) = -4*t3*t5  (since a = 0 and b = -2)
	t9.Mul(&t3, &t5)
	t9.Lsh(&t9, 2)
	t9.Neg(&t9)

	// t10 <- (t4 + alpha*t3)*(t5 + t7)  (with alpha = 2)
	t5.Add(&t5, &t7)
	t3.Lsh(&t3, 1)
	t3.Add(&t3, &t4)
	t10.Mul(&t3, &t5)

	// U3 <- -t6*(t1 - b*t2) = -t6*(t1 + 2*t2)
	t6.Neg(&t6)
	t2.Lsh(&t2, 1)
	t1.Add(&t1, &t2)
	P.u.Mul(&t6, &t1)

	// Z3 <- t8 - t9
	P.z.Sub(&t8, &t9)

	// T3 <- t8 + t9
	P.t.Add(&t8, &t9)

	// X3 <- b*(t10 - t8 + beta*t9)  (with b = -2 and beta = 1/2)
	//       = 2*(t8 - t10) - t9
	t8.Sub(&t8, &t10)
	t8.Lsh(&t8, 1)
	P.x.Sub(&t8, &t9)

	return P
}

// Set this point to the difference of the two provided points (P1 - P2).
// A pointer to this structure (P) is returned.
func (P *Do255ePoint) Sub(P1, P2 *Do255ePoint) *Do255ePoint {
	var P2n Do255ePoint
	P2n.x = P2.x
	P2n.z = P2.z
	P2n.u.Neg(&P2.u)
	P2n.t = P2.t
	return P.Add(P1, &P2n)
}

// Internal type for a point in affine (x, u) coordinates. This is used
// to speed up some computations but we do not make it public so that
// the API remains simple.
type do255ePointAffine struct {
	x, u gf.GF255e
}

// Set this point to the sum of the two provided points, the second of
// which being in affine (x, u) coordinates.
// A pointer to this structure (P) is returned.
func (P *Do255ePoint) addMixed(P1 *Do255ePoint, P2 *do255ePointAffine) *Do255ePoint {
	var t1, t3, t5, t6, t7, t8, t9, t10 gf.GF255e

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

	// t7 <- t1 + b*t2  (with b = -2 and t2 = Z1)
	t7.Sub(&t1, &P1.z)
	t7.Sub(&t7, &P1.z)

	// t8 <- t4*t7 = T1*t7
	t8.Mul(&P1.t, &t7)

	// t9 <- t3*(2*b*t5 + a*t7) = -4*t3*t5  (since a = 0 and b = -2)
	t9.Mul(&t3, &t5)
	t9.Lsh(&t9, 2)
	t9.Neg(&t9)

	// t10 <- (t4 + alpha*t3)*(t5 + t7)  (with alpha = 2 and t4 = T1)
	t5.Add(&t5, &t7)
	t3.Lsh(&t3, 1)
	t3.Add(&t3, &P1.t)
	t10.Mul(&t3, &t5)

	// U3 <- -t6*(t1 - b*t2) = -t6*(t1 + 2*Z1)
	t6.Neg(&t6)
	t5.Lsh(&P1.z, 1)
	t1.Add(&t1, &t5)
	P.u.Mul(&t6, &t1)

	// Z3 <- t8 - t9
	P.z.Sub(&t8, &t9)

	// T3 <- t8 + t9
	P.t.Add(&t8, &t9)

	// X3 <- b*(t10 - t8 + beta*t9)  (with b = -2 and beta = 1/2)
	//       = 2*(t8 - t10) - t9
	t8.Sub(&t8, &t10)
	t8.Lsh(&t8, 1)
	P.x.Sub(&t8, &t9)

	return P
}

// Set this point to the difference of the two provided points, the second of
// which being in affine (x, u) coordinates.
// A pointer to this structure (P) is returned.
func (P *Do255ePoint) subMixed(P1 *Do255ePoint, P2 *do255ePointAffine) *Do255ePoint {
	var P2n do255ePointAffine
	P2n.x = P2.x
	P2n.u.Neg(&P2.u)
	return P.addMixed(P1, &P2n)
}

// Set this point (P) to the double of the provided point Q.
// A pointer to this structure (P) is returned.
func (P *Do255ePoint) Double(Q *Do255ePoint) *Do255ePoint {
	// 3M+6S formulas

	// X' = (a^2-4*b)*X*Z
	// Z' = X^2 + a*X*Z + b*Z^2
	var t1, t2, xp, zp gf.GF255e
	t1.Sqr(&Q.x)
	t2.Sqr(&Q.z)
	xp.Add(&Q.x, &Q.z).Sqr(&xp).Sub(&xp, &t1).Sub(&xp, &t2).Lsh(&xp, 2)
	zp.Lsh(&t2, 1).Sub(&t1, &zp)

	// t1 = X^2
	// t2 = Z^2

	// X'' = 4*b*X'*Z'
	// Z'' = X'^2 - 2*a*X'*Z' + (a^2-4*b)*Z'^2
	var t3, t4, t5 gf.GF255e
	t3.Sqr(&xp)
	t4.Sqr(&zp)
	t5.Add(&xp, &zp).Sqr(&t5).Sub(&t5, &t3).Sub(&t5, &t4).Lsh(&t5, 2)
	t4.Lsh(&t4, 3)
	P.x.Neg(&t5)
	P.z.Add(&t3, &t4)

	// t3 = X'^2
	// t4 = (a^2-4*b)*Z'^2

	// U'' = 2*(a^2-4*b)*(X^2 - b*Z^2)*Z'*U
	// T'' = (X'^2 - (a^2-4*b)*Z'^2)*T
	zp.Mul(&zp, &Q.u)
	t2.Lsh(&t2, 1).Add(&t1, &t2).Mul(&t2, &zp)
	P.u.Lsh(&t2, 4)
	t3.Sub(&t3, &t4)
	P.t.Mul(&Q.t, &t3)

	return P
}

// Set this point (P) to (2^n)*Q (i.e. perform n successive doublings).
// This function is constant-time with regard to the point values, but
// not to the number of doublings (n); computation time is proportional
// to n.
// A pointer to this structure (P) is returned.
func (P *Do255ePoint) DoubleX(Q *Do255ePoint, n uint) *Do255ePoint {
	var tX, tW, tZ, t1, t2 gf.GF255e

	if n == 0 {
		P.Set(Q)
		return P
	}

	// First half-doubling, combined with conversion from fractional
	// (x,u) to Jacobian (x,w); output in E[r](-2*a,a^2-4*b).
	//   X' = Z^2*T^4
	//   W' = Z*T^2 - (2*X + a*Z)*U^2
	//   Z' = Z*U*T
	// Note that a = 0 for curve do255e.
	// Cost: 4M+2S
	tW.Sqr(&Q.u)       // tW <- U^2
	t1.Mul(&Q.z, &Q.t) // t1 <- Z*T
	tW.Lsh(&tW, 1)     // tW <- 2*U^2
	t2.Mul(&t1, &Q.t)  // t2 <- Z*T^2
	tW.Mul(&tW, &Q.x)  // tW <- 2*X*U^2
	tX.Sqr(&t2)        // tX <- Z^2*T^4
	tW.Sub(&t2, &tW)   // tW <- Z*T^2 - 2*X*U^2
	tZ.Mul(&t1, &Q.u)  // tZ <- Z*U*T

	// For n-1 doublings, apply psi_1/2() then psi_1().
	for n--; n > 0; n-- {
		// Internal doublings are 1M+5S each.
		t1.Sqr(&tW)                    // t1 <- W^2
		t2.Lsh(&tX, 1).Sub(&t1, &t2)   // t2 <- t1 - 2*X
		tW.Add(&tW, &t2).Sqr(&tW)      // tW <- (W + t2)^2
		t2.Sqr(&t2)                    // t2 <- t2^2
		tW.Sub(&tW, &t1).Sub(&tW, &t2) // tW <- tW - t1 - t2
		tZ.Mul(&tZ, &tW)               // Z <- Z*tW
		t1.Sqr(&t1).Lsh(&t1, 1)        // t1 <- 2*t1^2
		tW.Sub(&t2, &t1)               // W <- t2 - t1
		tX.Sqr(&t2)                    // X <- t2^2
	}

	// Final half-doubling, combined with conversion back to
	// fractional (x,u) coordinates.
	//   X' = 4*b*Z^2
	//   Z' = W^2
	//   U' = 2*W*Z
	//   T' = 2*X - 2*a*Z^2 - W^2
	// Note that a = 0 and b = -2 for curve do255e.
	// Cost: 3S  (with 2*W*Z = (W+Z)^2 - W^2 - Z^2)
	t1.Sqr(&tZ)                       // t1 <- Z^2
	P.z.Sqr(&tW)                      // Z' <- W^2
	P.t.Lsh(&tX, 1).Sub(&P.t, &P.z)   // T' <- 2*X - W^2
	t2.Add(&tW, &tZ).Sqr(&t2)         // t2 <- (W + Z)^2
	P.u.Sub(&t2, &t1).Sub(&P.u, &P.z) // U' <- 2*W*Z
	P.x.Lsh(&t1, 3).Neg(&P.x)         // X' <- -8*Z^2

	return P
}

// Set P to the opposite of point Q.
// A pointer to this structure (P) is returned.
func (P *Do255ePoint) Neg(Q *Do255ePoint) *Do255ePoint {
	P.x.Set(&Q.x)
	P.z.Set(&Q.z)
	P.u.Neg(&Q.u)
	P.t.Set(&Q.t)
	return P
}

// do255e endomorphism:
// ====================
//
// We use one of the cases described by Gallant, Lambert and Vanstone in
// their 2001 article:
//   https://www.iacr.org/archive/crypto2001/21390189.pdf
//
// Modulus is p = 2^255 - 18651. Curve has order 2*r, with:
//   r = 2^254 - 131528281291764213006042413802501683931
//
// Let eta = sqrt(-1) mod p. There are two such roots, we use this one:
// 7656063742463026568679823572395325799027601838558345258426535816504372595438
//
// Let phi(x, w) = (-x, -eta*w)
// In (x, u) coordinates, this is expressed as:
//    phi(x, u) = (-x, eta*u)
// because 1/eta = -eta.
//
// phi() is an endomorphism over our group G or order r. For any group
// element P, phi(P) = mu*P for a constant mu which is a square root of -1
// modulo r. With our choice of eta, we have the following mu:
// 23076176648693837106500022901799924463072024427516564762134831823525232195341
//
// Let k be a scalar (integer modulo r). We decompose k into two half-size
// values k0 and k1 such that k = k0 + k1*mu mod r.
//
// Since r = 1 mod 4, it can be written as the sum of two squares. Let u
// and v such that r = u^2 + v^2; the choice is unique up to permutation and
// change of sign; these values can be found by using Lagrange's algorithm
// on the lattice ((mu, 1), (r, 0)). We choose the following:
//
//   u =  34978546233976132960203755786038370577
//   v = 166506827525740345966246169588540045182
//
// Since (u/v)^2 = -1 mod r, value mu is equal to either u/v or v/u (mod r).
// With our choices, mu = u/v mod r.
//
// It can be verified that:
//   r = u^2 + v^2
//   v + mu*u = 0 mod r
//   -u + mu*v = 0 mod r
//
// Given k, we compute integers c and d as:
//   c = round(k*v / r)
//   d = round(k*u / r)
// Note that c and d are nonnegative. Moreover, c and d fit on 127 bits each.
//
// We then have:
//   k0 = k - d*u - c*v
//   k1 = d*v - c*u
//
// It can be shown (see GLV article) that k0^2 + k1^2 <= u^2 + v^2. Since
// u^2 + v^2 = r < 2^254, this implies that |k0| < 2^127 and |k1| < 2^127.
// Thus, k0 and k1 (which are signed integers) can fit on 128 bits each,
// including their sign bit.
//
//
// Rounded division:
// =================
//
// To compute c and d, we need to do a rounded division. We use the
// fact that r = 2^254 - r0 with r0 < 2^127.
//
// Suppose that we have x = k*u or k*v, and we want to compute y = round(x/r).
// Since r is odd, we have:
//   y = round(x/r) = floor((x + (r-1)/2) / r)
// Let z = x + (r-1)/2. We can split z at the 254-bit index:
//   z = z0 + 2^254*z1
// with 0 <= z0 < 2^254, and z1 >= 0. Since k < r and given the values of u
// and v, the maximum value for z1 is about 2^126.97, i.e it fits on 127 bits.
//
// We thus have:
//   z1*r = z1*2^254 - z1*r0 < z1*2^254 < z
// and
//   (z1+2)*r = z1^254 - z1*r0 + 2*2^254 - 2*r0
//            = (z1+1)*2^254 + (2^254 - (z1+2)*r0)
// Since (z1+2)*r0 < 2^254, we have (z1+2)*r > (z1+1)*2^254 > z.
//
// It follows that the rounded division result is necessarily either z1
// or z1+1. We can thus compute that rounded division with the following
// algorithm:
//
//   Input: integer x such that 0 <= x <= (r-1)*max(u,v)
//   Output: round(x / r)
//     1. z <- x + (r-1)/2
//     2. y <- floor(z / 2^254) + 1
//     3. if y*r > z, then y <- y-1
//     4. return y
//
// The most expensive operation is the product y*r. However, we only need
// the sign of z - y*r. We can do that computation as follows:
//    z - y*r = (y-1)*2^254 + z0 - y*2^254 + y*r0
//            = z0 + y*r0 - 2^254
// We thus need to subtract 1 from y if and only if z0 + y*r0 is strictly
// lower than 2^254. y and r0 both fit on 127 bits each, and z0 is less
// than 2^254; we can thus do that computation over 255 bits.

var do255eEta = gf.GF255e{
	0xD99E0F1BAA938AEE, 0xA60D864FB30E6336,
	0xE414983FE53688E3, 0x10ED2DB33C69B85F}
var do255eMinusEta = gf.GF255e{
	0x2661F0E4556C2C37, 0x59F279B04CF19CC9,
	0x1BEB67C01AC9771C, 0x6F12D24CC39647A0}

// Multiply a point Q by a given scalar n.
// A pointer to this structure (P) is returned.
func (P *Do255ePoint) Mul(Q *Do255ePoint, n *Do255eScalar) *Do255ePoint {
	// Split input scalar into k0 and k1.
	var k0, k1 [2]uint64
	n.SplitMu(&k0, &k1)

	// Recode k0 and fill M with a copy of Q or -Q, depending on
	// whether k0 is negative or not.
	var sd0 [26]byte
	sg := scalar.Recode5SmallSigned(&sd0, &k0)
	M := *Q
	M.u.CondNeg(&M.u, sg)

	// Initialize the low window with i*M for i in 1..16
	// (point i*Q is stored at offset i-1).
	var win0 [16]Do255ePoint
	win0[0] = M
	win0[1].Double(&M)
	for i := 3; i <= 15; i += 2 {
		win0[i-1].Add(&win0[i-2], &M)
		win0[i].Double(&win0[((i+1)>>1)-1])
	}

	// Recode k1 and remember whether the sign of k1 differs from
	// that of k0.
	var sd1 [26]byte
	sg ^= scalar.Recode5SmallSigned(&sd1, &k1)

	// Apply the endomorphism on all points of the low window to
	// fill the high window. If k0 and k1 have different signs,
	// then an extra negation is performed.
	var endo gf.GF255e
	endo.CondNeg(&do255eEta, sg)
	var win1 [16]Do255ePoint
	for i := 0; i < 16; i++ {
		win1[i].x.Neg(&win0[i].x)
		win1[i].z = win0[i].z
		win1[i].u.Mul(&win0[i].u, &endo)
		win1[i].t = win0[i].t
	}

	// Lookup points corresponding to the top digits, and add them
	// to get the initial value of P. The top digits are both
	// nonnegative.
	P.do255eLookupWindow(&win0, uint(sd0[25]))
	M.do255eLookupWindow(&win1, uint(sd1[25]))
	P.Add(P, &M)

	// Process other digits from top to bottom.
	for i := 24; i >= 0; i-- {
		P.DoubleX(P, 5)
		M.do255eLookupWindow(&win0, uint(sd0[i]&0x1F))
		M.u.CondNeg(&M.u, uint64(sd0[i]>>7))
		P.Add(P, &M)
		M.do255eLookupWindow(&win1, uint(sd1[i]&0x1F))
		M.u.CondNeg(&M.u, uint64(sd1[i]>>7))
		P.Add(P, &M)
	}

	return P
}

// Multiply the conventional generator by a given scalar n. This is
// functionally equivalent (but faster) to P.Generator().Mul(&P, n).
// A pointer to this structure (P) is returned.
func (P *Do255ePoint) MulGen(n *Do255eScalar) *Do255ePoint {
	// Recode input scalar into 5-bit Booth encoding.
	var sd [52]byte
	n.recode5(&sd)

	// Lookup initial accumulator by using the top digit (which is
	// guaranteed nonnegative).
	var Ma do255ePointAffine
	Ma.do255eLookupWindowAffine(&do255eWin_G195_xu, uint(sd[51]))
	P.x = Ma.x
	P.z = gf.GF255e_ONE
	P.u = Ma.u
	P.t = gf.GF255e_ONE

	// Add points corresponding to top digits of the three other
	// quarter-scalars.
	Ma.do255eLookupWindowAffine(&do255eWin_G_xu, uint(sd[12]&0x1F))
	Ma.u.CondNeg(&Ma.u, uint64(sd[12]>>7))
	P.addMixed(P, &Ma)
	Ma.do255eLookupWindowAffine(&do255eWin_G65_xu, uint(sd[25]&0x1F))
	Ma.u.CondNeg(&Ma.u, uint64(sd[25]>>7))
	P.addMixed(P, &Ma)
	Ma.do255eLookupWindowAffine(&do255eWin_G130_xu, uint(sd[38]&0x1F))
	Ma.u.CondNeg(&Ma.u, uint64(sd[38]>>7))
	P.addMixed(P, &Ma)

	// Process all other digits from high to low. We process the
	// four quarter-scalars in parallel.
	for i := 11; i >= 0; i-- {
		P.DoubleX(P, 5)

		Ma.do255eLookupWindowAffine(&do255eWin_G_xu, uint(sd[i]&0x1F))
		Ma.u.CondNeg(&Ma.u, uint64(sd[i]>>7))
		P.addMixed(P, &Ma)
		Ma.do255eLookupWindowAffine(&do255eWin_G65_xu, uint(sd[i+13]&0x1F))
		Ma.u.CondNeg(&Ma.u, uint64(sd[i+13]>>7))
		P.addMixed(P, &Ma)
		Ma.do255eLookupWindowAffine(&do255eWin_G130_xu, uint(sd[i+26]&0x1F))
		Ma.u.CondNeg(&Ma.u, uint64(sd[i+26]>>7))
		P.addMixed(P, &Ma)
		Ma.do255eLookupWindowAffine(&do255eWin_G195_xu, uint(sd[i+39]&0x1F))
		Ma.u.CondNeg(&Ma.u, uint64(sd[i+39]>>7))
		P.addMixed(P, &Ma)
	}

	return P
}

// Constant-time lookup of a point in a window. Provided window has
// 16 elements. Input offset ('index') is in the 0..16 range. This
// function sets P to a copy of win[index - 1] if index != 0, or to
// the neutral if index == 0.
func (P *Do255ePoint) do255eLookupWindow(win *[16]Do255ePoint, index uint) {
	// Initialize P to all-zeros.
	P.x = gf.GF255e_ZERO
	P.z = gf.GF255e_ZERO
	P.u = gf.GF255e_ZERO
	P.t = gf.GF255e_ZERO

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
	P.z.CondOrFrom(&gf.GF255e_ONE, mz)
	P.t.CondOrFrom(&gf.GF255e_ONE, mz)
}

// Constant-time lookup of a point in a window. This is similar to
// do255eLookupWindow(), except that this function works on points in
// affine (x, u) coordinates.
func (P *do255ePointAffine) do255eLookupWindowAffine(win *[16]do255ePointAffine, index uint) {
	// Initialize P to all-zeros (which is the valid representation
	// of the neutral element).
	P.x = gf.GF255e_ZERO
	P.u = gf.GF255e_ZERO

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
func (P *Do255ePoint) addFromWindowVartime(win *[16]Do255ePoint, ej byte) {
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
func (P *Do255ePoint) addFromWindowAffineVartime(win *[16]do255ePointAffine, ej byte) {
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
func (P *Do255ePoint) VerifyHelperVartime(k0, k1 *Do255eScalar, encR []byte) bool {
	// Decode encR as a field element. If it is not a valid field
	// element encoding then we can reject it right away. We do not
	// need to decode it as a point.
	var Rw gf.GF255e
	if Rw.Decode(encR) != 1 {
		return false
	}

	// If P is neutral, we can ignore k1 and use MulGen().
	if P.IsNeutral() != 0 {
		var P2 Do255ePoint
		P2.MulGen(k0)
		if P2.IsNeutral() != 0 {
			return Rw.IsZero() == 1
		} else {
			Rw.Mul(&Rw, &P2.u)
			return Rw.Eq(&P2.t) == 1
		}
	}

	// Handle scalars: we want to obtain four half-width recoded
	// scalars (26 digits each) operating on G, 2^130*G, P and
	// phi(P) (with phi() being the endomorphism on the curve). For
	// G and 2^130*G, we can simply cut k0 in two halves; for P,
	// we use the endomorphism split (see above for details). We
	// also build two 5-bit windows for P and phi(P).

	// Recode scalar k0.
	var sd0 [52]byte
	k0.recode5(&sd0)

	// Split scalar k1.
	var k1_lo, k1_hi [2]uint64
	k1.SplitMu(&k1_lo, &k1_hi)

	// Recode k1_lo and fill M with a copy of P or -P, depending on
	// whether 1_lo0 is negative or not.
	var sd1_lo [26]byte
	sg := scalar.Recode5SmallSigned(&sd1_lo, &k1_lo)
	M := *P
	M.u.CondNeg(&M.u, sg)

	// Initialize the low window with i*M for i in 1..16
	// (point i*P is stored at offset i-1).
	var win_lo [16]Do255ePoint
	win_lo[0] = M
	win_lo[1].Double(&M)
	for i := 3; i <= 15; i += 2 {
		win_lo[i-1].Add(&win_lo[i-2], &M)
		win_lo[i].Double(&win_lo[((i+1)>>1)-1])
	}

	// Recode k1_hi and remember whether the sign of k1_hi differs from
	// that of k1_lo.
	var sd1_hi [26]byte
	sg ^= scalar.Recode5SmallSigned(&sd1_hi, &k1_hi)

	// Apply the endomorphism on all points of the low window to
	// fill the high window. If k1_lo and k1_hi have different signs,
	// then an extra negation is performed.
	var endo gf.GF255e
	endo.CondNeg(&do255eEta, sg)
	var win_hi [16]Do255ePoint
	for i := 0; i < 16; i++ {
		win_hi[i].x.Neg(&win_lo[i].x)
		win_hi[i].z = win_lo[i].z
		win_hi[i].u.Mul(&win_lo[i].u, &endo)
		win_hi[i].t = win_lo[i].t
	}

	// We now have the half-width scalars (low and high half of
	// sd0, sd1_lo, sd1_hi) and the windows (win_lo and win_hi).
	// We do a simple combined ppint multiplications. We use a
	// regular schedule (additions done every 5 point doublings)
	// and not NAF-w recoding, because our formulas benefit from
	// making several successive doublings.

	// Lookup initial accumulator by using the top digit of sd0
	// (which is nonnegative).
	j := int(sd0[51])
	if j != 0 {
		M.x = do255eWin_G130_xu[j-1].x
		M.u = do255eWin_G130_xu[j-1].u
	} else {
		M.x = gf.GF255e_ZERO
		M.u = gf.GF255e_ZERO
	}
	M.z = gf.GF255e_ONE
	M.t = gf.GF255e_ONE

	// Add points from the top digits of the three other half-width
	// recoded scalars.
	M.addFromWindowAffineVartime(&do255eWin_G_xu, sd0[25])
	M.addFromWindowVartime(&win_lo, sd1_lo[25])
	M.addFromWindowVartime(&win_hi, sd1_hi[25])

	// Process all other digits in top to bottom order.
	for i := 24; i >= 0; i-- {
		M.DoubleX(&M, 5)
		M.addFromWindowAffineVartime(&do255eWin_G_xu, sd0[i])
		M.addFromWindowAffineVartime(&do255eWin_G130_xu, sd0[26+i])
		M.addFromWindowVartime(&win_lo, sd1_lo[i])
		M.addFromWindowVartime(&win_hi, sd1_hi[i])
	}

	// We want to check that the w = 1/u coordinate of the result
	// (in M) matches that which was decoded from encR. We need
	// to make a special case for u = 0, since the undefined w
	// is then encoded as a zero.
	if M.IsNeutral() != 0 {
		return Rw.IsZero() == 1
	} else {
		Rw.Mul(&Rw, &M.u)
		return Rw.Eq(&M.t) == 1
	}
}

// Alternate helper function for signature verification; this returns
// 1 if k0*G + k1*P yields a point whose encoding is exactly equal to
// the first 32 bytes of encR; otherwise, it returns 0. This function
// is usually slower than VerifyHelperVartime(), but it is constant-time.
func (P *Do255ePoint) VerifyHelper(k0, k1 *Do255eScalar, encR []byte) int {
	var P1, P2 Do255ePoint
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
func (P *Do255ePoint) MapBytes(bb []byte) *Do255ePoint {
	// Decode source bytes as a field element. This applies modular
	// reduction.
	var e gf.GF255e
	e.DecodeReduce(bb)

	// Set a flag if the value is zero. We use it at the end to
	// set the result to N in that case. This allows us to assume
	// that e != 0 in the rest of the code.
	ez := e.IsZero()

	// Formulas: we map into a point on the dual curve y^2 = x^3 + bb*x,
	// with bb = -4*b:
	//   x1 = e + (1-bb)/(4*e)
	//   x2 = d*(e - (1-bb)/(4*e))
	// with d = sqrt(-1) (same square root as in the endomorphism).
	// Then, at least one of x1, x2 and x1*x2 is a proper x
	// coordinate for a curve point. We use the following:
	//   x1^3 + bb*x1 = (64*e^7 + 16*(3+bb)*e^5
	//                   + 4*(3-2*bb-bb^2)*e^3 + (1-bb)^3*e) / (64*e^4)
	//   x2^3 + bb*x2 = -d*(64*e^7 - 16*(3+bb)*e^5
	//                      + 4*(3-2*bb-bb^2)*e^3 - (1-bb)^3*e) / (64*e^4)
	//   ((x1*x2)^3 + bb*(x1*x2)) = (x1^3 + bb*x1)*(x2^3 + bb*x2)
	// For a properly deterministic result, we use the square root
	// of the numerator whose least significant bit (as an integer
	// in the 0..p-1 range) is zero; this is what Sqrt() returns.

	var e2, e3, e4, e5, e7, e8 gf.GF255e
	e2.Sqr(&e)
	e4.Sqr(&e2)
	e8.Sqr(&e4)
	e3.Mul(&e, &e2)
	e5.Mul(&e3, &e2)
	e7.Mul(&e5, &e2)

	// Compute x1 and x2, each as a fraction; the denominator is the
	// same for both (4*e). With bb = 8, we have:
	//   x1num  = 4*e^2 - 7
	//   x2num  = d*(4*e^2 + 7)
	//   x12den = 4*e
	var x1num, x2num, x12den gf.GF255e
	x1num.Lsh(&e2, 2)
	x1num.Sub(&x1num, &gf.GF255e_SEVEN)
	x2num.Lsh(&e2, 2)
	x2num.Add(&x2num, &gf.GF255e_SEVEN)
	x2num.Mul(&x2num, &do255eEta)
	x12den.Lsh(&e, 2)

	// Compute the values of y1^2 and y2^2 (candidates).
	// With bb = 8, numerator polynomials are:
	//   yy1num = 64*e^7 + 176*e^5 - 308*e^3 - 343*e
	//   yy2num = -d*(64*e^7 - 176*e^5 - 308*e^3 + 343*e)
	// for the common denominator:
	//   yyden  = 64*e^4
	// We directly use the square root of the denominator:
	//   y12den = 8*e^2
	// Note that since 1-bb = -7, and 7 and -7 are not quadratic
	// residues in the field, then neither x1 nor x2 can be zero;
	// moreover, x^2 + bb cannot be 0 for any x, since -bb is not a
	// quadratic residue. Therefore, x1^3 + bb*x1 and x2^3 + bb*x2
	// are both non-zero, and the Legendre symbol is then 1 or -1
	// for each.
	var yy1num, yy2num, y12den, tt gf.GF255e
	yy1num.Lsh(&e7, 6)
	yy2num.Set(&yy1num)
	tt.Mul(&e5, &gf.GF255e_HUNDREDSEVENTYSIX)
	yy1num.Add(&yy1num, &tt)
	yy2num.Sub(&yy2num, &tt)
	tt.Mul(&e3, &gf.GF255e_THREEHUNDREDEIGHT)
	yy1num.Sub(&yy1num, &tt)
	yy2num.Sub(&yy2num, &tt)
	tt.Mul(&e, &gf.GF255e_THREEHUNDREDFORTYTHREE)
	yy1num.Sub(&yy1num, &tt)
	yy2num.Add(&yy2num, &tt)
	yy2num.Mul(&yy2num, &do255eMinusEta)
	y12den.Lsh(&e2, 3)

	// x3 = x1*x2
	// y3^2 = (y1^2)*(y2^2)
	var x3num, x3den, yy3num, y3den gf.GF255e
	x3num.Mul(&x1num, &x2num)
	x3den.Lsh(&e2, 4)
	yy3num.Mul(&yy1num, &yy2num)
	y3den.Lsh(&e4, 6)

	// Get the Legendre symbols for y1^2 and y2^2.
	ls1 := yy1num.Legendre()
	ls2 := yy2num.Legendre()

	// If ls1 == 1, then we use x1 and yy1. Otherwise, if ls1 == -1
	// and ls2 == 1, then we use x2 and yy2. Otherwise, we use x1*x2,
	// and yy1*yy2.
	qr1 := 1 - (ls1 >> 63)
	qr2 := 1 - (ls2 >> 63)

	var xnum, xden, yynum, yden gf.GF255e
	xnum.Select(&x1num, &x2num, qr1)
	xnum.Select(&xnum, &x3num, qr1|qr2)
	xden.Select(&x12den, &x3den, qr1|qr2)
	yynum.Select(&yy1num, &yy2num, qr1)
	yynum.Select(&yynum, &yy3num, qr1|qr2)
	yden.Select(&y12den, &y3den, qr1|qr2)

	// Extract the square root of yynum; it cannot fail. The setSqrt()
	// function ensures that we get a least significant bit equal to 0.
	var ynum gf.GF255e
	ynum.Sqrt(&yynum)

	// We now have the point in fractional (x,y) coordinates. We
	// compute the coordinate u = x/y.
	var unum, uden gf.GF255e
	unum.Mul(&xnum, &yden)
	uden.Mul(&xden, &ynum)

	// Apply the isogeny theta'_{1/2} to get a point in the proper
	// group on the right curve:
	//   x' = 4*b*u^2
	//   u' = 2*x/(u*(x^2 - bb))
	P.x.Sqr(&unum).Lsh(&P.x, 3).Neg(&P.x)
	P.z.Sqr(&uden)
	P.u.Mul(&xnum, &xden).Mul(&P.u, &uden).Lsh(&P.u, 1)
	xden.Sqr(&xden).Lsh(&xden, 3)
	P.t.Sqr(&xnum).Sub(&P.t, &xden).Mul(&P.t, &unum)

	// If the source value e was zero, then all of the above is
	// invalid, and we force the point to N, i.e. x = 0 and u = 0.
	P.x.Select(&gf.GF255e_ZERO, &P.x, ez)
	P.z.Select(&gf.GF255e_ONE, &P.z, ez)
	P.u.Select(&gf.GF255e_ZERO, &P.u, ez)
	P.t.Select(&gf.GF255e_ONE, &P.t, ez)

	return P
}

// =====================================================================
// Precomputed windows in affine (x, u) coordinates:
//   i*G         for i = 1..16
//   i*2^65*G    for i = 1..16
//   i*2^130*G   for i = 1..16
//   i*2^195*G   for i = 1..16

var do255eWin_G_xu = [16]do255ePointAffine{
	/* 1 */
	do255ePointAffine{
		x: gf.GF255e{0x0000000000000002, 0x0000000000000000,
			0x0000000000000000, 0x0000000000000000},
		u: gf.GF255e{0x0000000000000001, 0x0000000000000000,
			0x0000000000000000, 0x0000000000000000},
	},
	/* 2 */
	do255ePointAffine{
		x: gf.GF255e{0xE38E38E38E38AAE3, 0x8E38E38E38E38E38,
			0x38E38E38E38E38E3, 0x638E38E38E38E38E},
		u: gf.GF255e{0xB6DB6DB6DB6D97A3, 0xDB6DB6DB6DB6DB6D,
			0x6DB6DB6DB6DB6DB6, 0x36DB6DB6DB6DB6DB},
	},
	/* 3 */
	do255ePointAffine{
		x: gf.GF255e{0x0000000000000152, 0x0000000000000000,
			0x0000000000000000, 0x0000000000000000},
		u: gf.GF255e{0xC2F21347C4043E79, 0x6B1CEBA6066D4156,
			0xAB617909A3E20224, 0x12358E75D30336A0},
	},
	/* 4 */
	do255ePointAffine{
		x: gf.GF255e{0xEA5E1BA07D1A36D1, 0xD0AF77073CCEA916,
			0xC68A748FF2D92037, 0x43554439A8DE7571},
		u: gf.GF255e{0x65A29F71130DB4AD, 0x9F71130DFA47C8BB,
			0x130DFA47C8BB65A2, 0x7A47C8BB65A29F71},
	},
	/* 5 */
	do255ePointAffine{
		x: gf.GF255e{0x9B7D88CD74D7D3CA, 0x0E31B461193896AC,
			0x93464506E97D44DB, 0x0ABC8AC61DCD9949},
		u: gf.GF255e{0x1F2B6B08DA5B43EE, 0xE40F8B8BC44A0C63,
			0x5866F1F8B35FB70C, 0x185034D250F768D7},
	},
	/* 6 */
	do255ePointAffine{
		x: gf.GF255e{0x19348B5724DF12CF, 0xAB7572F66232945A,
			0xFFAA89042A63BF4B, 0x3024ED85633EF2E2},
		u: gf.GF255e{0x0BD0C5F1F91D6B18, 0xBB4A410D263610A7,
			0xA1AB0B9D98F35F00, 0x4FA6D8B6AFDDC92B},
	},
	/* 7 */
	do255ePointAffine{
		x: gf.GF255e{0x33E381251F43C3D5, 0x49BF0E2E71C6FE8F,
			0x6AF69CC116BAEF18, 0x36199FBAD0C9585E},
		u: gf.GF255e{0x7EB52414159EF4EA, 0xB885C9D1EB4CC9E1,
			0x350914B3EE64BF7F, 0x6DD8CDFA520AED5A},
	},
	/* 8 */
	do255ePointAffine{
		x: gf.GF255e{0xDDBCA65ECDCBDFAE, 0xD19C744A47F38FBD,
			0x15099FB55653455E, 0x1FD0DC2A961EC01F},
		u: gf.GF255e{0x99DA8C93EB513A8B, 0x0706B8B95DEDFC87,
			0xC54D8F471F778CE9, 0x4766315BFA2E63E5},
	},
	/* 9 */
	do255ePointAffine{
		x: gf.GF255e{0x8DF7F2C9CE799A9C, 0x7AC7F7C3AFDD04F9,
			0x915FE4A27D833740, 0x1ED67871986F29BA},
		u: gf.GF255e{0xA84A27A9D0A08E61, 0x27E9084D132CCAC1,
			0x498C7D8B01F68C40, 0x6957FDFF940E4159},
	},
	/* 10 */
	do255ePointAffine{
		x: gf.GF255e{0xE944946059A9387C, 0x32BB0464A75287B8,
			0x122E571F46C8845D, 0x0D05AC0126E0A481},
		u: gf.GF255e{0x3AA366BBB889903E, 0x55838146CC140A37,
			0x4AA37581A9B6AD5E, 0x7B37113C916F803C},
	},
	/* 11 */
	do255ePointAffine{
		x: gf.GF255e{0x6562064C442E3709, 0xD013EB4D114A7267,
			0x166892C716D5320A, 0x2824BCCA3B493396},
		u: gf.GF255e{0xA9A8911D864E7F82, 0x65CF6B9CAB741725,
			0x8C133221E772B327, 0x158521078CD1F209},
	},
	/* 12 */
	do255ePointAffine{
		x: gf.GF255e{0x4AE1AE2876FFE733, 0x55A43A11F9D28845,
			0xBAACD8A3E4990483, 0x37B39256440F5C21},
		u: gf.GF255e{0xE4C1C725087640AA, 0xFB902D6A3EF5D5E0,
			0x53EF35932E1297EB, 0x67E65CF7E1787343},
	},
	/* 13 */
	do255ePointAffine{
		x: gf.GF255e{0x90DF642868789634, 0x267A28B9CB72C6CA,
			0x27BE4B2B937625B5, 0x62003971A89B844F},
		u: gf.GF255e{0x90F8839881061965, 0x67D0394FF2BFCB98,
			0x913200FCCD1396D8, 0x17F96D76306A3580},
	},
	/* 14 */
	do255ePointAffine{
		x: gf.GF255e{0xC76258B805A821DA, 0xDC3C29F024B5765F,
			0xB646CAAD30897EED, 0x46F594DEF8D35CB9},
		u: gf.GF255e{0x05B49673D2AC4172, 0xA016A6890D77E4E6,
			0x7C6DAA970635E1C0, 0x42C8034547A6A04A},
	},
	/* 15 */
	do255ePointAffine{
		x: gf.GF255e{0x496055CDC3DFB745, 0x0673F992F547B770,
			0x8EFAE8F99B3E5BB3, 0x33C76A13E12C07DD},
		u: gf.GF255e{0x7FFA4AF719120727, 0x705D12571BF74984,
			0x4AD1FA649FAE1F07, 0x2F4CA2B6265D7456},
	},
	/* 16 */
	do255ePointAffine{
		x: gf.GF255e{0xE1DAC53644F243F9, 0xE1154ECDDFFD59A1,
			0x1731585F8B8C6649, 0x1BFB93F1365D2CAF},
		u: gf.GF255e{0x2D808316E1227049, 0x15064C9132683177,
			0x706D8A1F41E90ED8, 0x251A19311A6DB76E},
	},
}

var do255eWin_G65_xu = [16]do255ePointAffine{
	/* 1 */
	{
		x: gf.GF255e{0xABEF504D87FDEB41, 0x3A2D867D250A6B59,
			0x5906F21AAB3E16A0, 0x58327D52081B8A67},
		u: gf.GF255e{0xAD5F8FBFC596FA71, 0x415893549DE223FF,
			0x395D2181E50A4384, 0x1B313D36A8A7626E},
	},
	/* 2 */
	{
		x: gf.GF255e{0x15356A7D39689FC6, 0xF52A5DE2A0967CE6,
			0xAD2738B8D02C8707, 0x75B28A23988AA077},
		u: gf.GF255e{0x2B7C2F71889E3F33, 0x4CA4C049A5E65CF4,
			0x5CD27D909976BFE7, 0x0BE56F359985D602},
	},
	/* 3 */
	{
		x: gf.GF255e{0xDC3DE207745F399F, 0xD7C5D084CA47E2C5,
			0xDFA47415EF678E9B, 0x5F4DD2E4049C2FE8},
		u: gf.GF255e{0x4A504A5DED61CB7F, 0xAF41508342D7801D,
			0x0519A68AAB4295EB, 0x098D3AB90B09C2B4},
	},
	/* 4 */
	{
		x: gf.GF255e{0x977017D6C9BF2E76, 0x17CE586894D3EC82,
			0xA07F1106A288282D, 0x1334610CA14E17A8},
		u: gf.GF255e{0xD98BB32E81B4D89F, 0x5A0FCD3AE1067F08,
			0x3B845EF3AFAD3191, 0x3540EB32AA6E2C23},
	},
	/* 5 */
	{
		x: gf.GF255e{0x592F8CC3E5028371, 0x45BC1E36D117EB2E,
			0x531CFBF930BEEB8F, 0x2829EA75A68A6C5E},
		u: gf.GF255e{0xDB2EB15DD80EB7AB, 0x695C5863633AC0BB,
			0x9DD36D925FB45810, 0x4B36B3BE374E74CB},
	},
	/* 6 */
	{
		x: gf.GF255e{0x94156556D9FB736E, 0x3FCA76C281E2244C,
			0x7A11BA953BA95F28, 0x2E6BCF6A3DD91BBF},
		u: gf.GF255e{0x2473D6DB1425C001, 0x763067186E58108F,
			0x0C0776C9F6ED32DE, 0x143CF8A090ACB085},
	},
	/* 7 */
	{
		x: gf.GF255e{0xE24BA37DA32EFB19, 0xB319B89EBE6E4852,
			0xB99C333F17938861, 0x053979BEC7A8BA46},
		u: gf.GF255e{0xE95E218127D26209, 0x33088969E8FB612B,
			0xBAC06821F2BF788C, 0x7EE7B8C1D61C83DD},
	},
	/* 8 */
	{
		x: gf.GF255e{0x6B70CE277F74345E, 0x5FB603F5B5D2FF08,
			0x9A28F63978F82349, 0x1E7BF7585E894EE6},
		u: gf.GF255e{0x4C033B2301944CB8, 0x3260DA08A5D337CB,
			0x34DFEFB2A6C39FFC, 0x235872299B5B142E},
	},
	/* 9 */
	{
		x: gf.GF255e{0xEDC150B81CDB24C5, 0xCB26D5421BE70FF6,
			0xF19BB5BF01AA63A4, 0x316C4EA24603AF0D},
		u: gf.GF255e{0x6742B67C6086DF92, 0x53193A718D75331F,
			0xADD356D1CA14E352, 0x07DF0CC5D7C9E88D},
	},
	/* 10 */
	{
		x: gf.GF255e{0xE36D90729B572DFF, 0xC4F6284A7B9548E7,
			0x0FDCACBCFEDA82E0, 0x40A3C68531C84199},
		u: gf.GF255e{0x7982C54051C4DF08, 0x476A1D01EBBFE2BB,
			0x3F6064CC7E287D25, 0x79E82E420C17C457},
	},
	/* 11 */
	{
		x: gf.GF255e{0xC10298FEAF5791A3, 0x36AD705BB82FDC4C,
			0xA8764B80FAC163BB, 0x1388EA98CCAF5744},
		u: gf.GF255e{0x00FF128613403E58, 0x854CBABF2D02F041,
			0xB97069DD65593890, 0x3CC8BC3168B5A376},
	},
	/* 12 */
	{
		x: gf.GF255e{0x230E07753775607D, 0xDA9A2205BF1C2AE5,
			0x975954E5C2DD10BB, 0x28652C57B5420892},
		u: gf.GF255e{0x3688782C6588D284, 0xA50B5985163B51DB,
			0x23447D42F7311FA5, 0x51CBB2B65B751D79},
	},
	/* 13 */
	{
		x: gf.GF255e{0x4CBEDE49F7C4233B, 0x4EC8BC06D529F82D,
			0x75FB83AA2630BCB9, 0x7F3F2E391D0902B6},
		u: gf.GF255e{0x6AF277E67A9FAA68, 0xFDEE555C9009703F,
			0xAA36C89F8D5B502A, 0x052051EDAFC87131},
	},
	/* 14 */
	{
		x: gf.GF255e{0x2C9FA67ECEC074F3, 0x077855CE1E217BAC,
			0xAE612A7D92BCDFF8, 0x7B5805F9F2BC9375},
		u: gf.GF255e{0x86A494FC08EB56DD, 0xF694DF6ED62219A0,
			0x214282B022098B86, 0x7182E92FA95620EB},
	},
	/* 15 */
	{
		x: gf.GF255e{0xD48702C2FB7A4ACE, 0x9009655147735647,
			0xD98938B066322FCE, 0x08B7C64434681715},
		u: gf.GF255e{0xD25E4076E71417FD, 0xAE53B9A6D617E037,
			0xEDF6131ABB08D620, 0x406E88577206C6C4},
	},
	/* 16 */
	{
		x: gf.GF255e{0xD86B084C27BD12CA, 0x772D4452961D7D94,
			0x109ACCAAED9EC998, 0x575C5B51EF95AAE4},
		u: gf.GF255e{0x7D37BBB6D0781C73, 0x1ADF4CFE82DAF6AB,
			0x6CB90D1DDFBA601C, 0x7B916A823C090D35},
	},
}

var do255eWin_G130_xu = [16]do255ePointAffine{
	/* 1 */
	{
		x: gf.GF255e{0xC796B2402F4B5577, 0xB848DE143824131E,
			0x35479390B17AF31F, 0x0EB9DD4F3DFA2F98},
		u: gf.GF255e{0x6A6F8767CAF58BCB, 0xEF9A2E5ACD520CD9,
			0x2B998E19EE40437C, 0x1E3A7692F3E02AB1},
	},
	/* 2 */
	{
		x: gf.GF255e{0x9B55C0CC7A0BD2F8, 0x73D4BA09593578E5,
			0x89A107F02B82E751, 0x23C4CFE49D8C1DB2},
		u: gf.GF255e{0x97643EA7C28259E8, 0x64C33BBEA0416456,
			0xAC5EBA85AFFFBCEB, 0x1DE0359F8936CEEA},
	},
	/* 3 */
	{
		x: gf.GF255e{0x29BC60641F6F9439, 0xA4470A0D45E0AB6F,
			0x3CF550B4CBA43D48, 0x098DBD0C76126C5C},
		u: gf.GF255e{0x4A8F6858185A63C2, 0xC29D0227105E6338,
			0x4FA122A357313D72, 0x62BAF3EF3C842009},
	},
	/* 4 */
	{
		x: gf.GF255e{0xE5AFFE37730A6B2D, 0x4AD5E971E27ED8D1,
			0xCD5F1FF510C3D246, 0x05D71BC1B236B80B},
		u: gf.GF255e{0x17625F88FFC35397, 0xEF697901DA783099,
			0xDCCCD3E1EF458EDC, 0x37EDC360CFBBEDE2},
	},
	/* 5 */
	{
		x: gf.GF255e{0x15EE9577D26B0A12, 0x796D2F22498E2395,
			0xEEBF418FA5DC2FFF, 0x6EFC71DFB61B892C},
		u: gf.GF255e{0x055D6CCDE9EC95DF, 0xEB86F24ADFCADBF5,
			0x6DDAC4AB9F7AE17C, 0x282E209409ADD692},
	},
	/* 6 */
	{
		x: gf.GF255e{0x1D2065C52CC05E06, 0xE136A646AD8CCE04,
			0x9EEE9C6F51A1AD2C, 0x29BD9CB831800380},
		u: gf.GF255e{0xBA9DA77EDCAD7528, 0xB6FA24A8892C93FF,
			0xF406AD1FFA1ACD6B, 0x7E4F46D5131E195C},
	},
	/* 7 */
	{
		x: gf.GF255e{0x3764EFB34F88E075, 0x932413347C635C0D,
			0xC9DA4905A98FA573, 0x3B23E4AB52D89A2F},
		u: gf.GF255e{0x9E343B6FEBEF1413, 0xC1293137B7F95D67,
			0x6FB5672D05CA0B5B, 0x372239482AEC172D},
	},
	/* 8 */
	{
		x: gf.GF255e{0xC2D4B896B99ECE2D, 0x12CFC73515B85FA6,
			0x68977BE4A000B16A, 0x1A7EDDB0A3FD9E15},
		u: gf.GF255e{0x21232B19EA446499, 0x0F110F3BEDDC580B,
			0x7BFEE5167F928B59, 0x4E1CF6CCC48289F4},
	},
	/* 9 */
	{
		x: gf.GF255e{0x66DD5CBDC8683D52, 0x0E57F934D5717796,
			0x06E7AFEFE940A0DE, 0x5E68DFB004FFD454},
		u: gf.GF255e{0xFFFA0532A28FF940, 0x06134380E791A9D0,
			0x24120A75934E696A, 0x6F671B022C83BC57},
	},
	/* 10 */
	{
		x: gf.GF255e{0x382654D53F6219C8, 0x5B068A05F61C4345,
			0xC106F1243C0212B6, 0x7C5492368C243A54},
		u: gf.GF255e{0x0EE1B5003B071E37, 0x3EC320448D80BCA7,
			0x108D5A1C5EC7534A, 0x64ECAAC0EBA1A511},
	},
	/* 11 */
	{
		x: gf.GF255e{0xBFB4842BF8CB93CF, 0xF994337DE38CC098,
			0xB9BB68C291517A22, 0x56B3E07FE03B3A5A},
		u: gf.GF255e{0x6CC3477AD4AB70A1, 0x0BC92225712BA4D4,
			0xD51C733D4564D8F4, 0x029CAFA5599B372A},
	},
	/* 12 */
	{
		x: gf.GF255e{0xA0A8890960C40BEC, 0xCE8ECA1A58BBBBC9,
			0x076AB31BA9BBF6DE, 0x0E1B23A02743D405},
		u: gf.GF255e{0xCB9E1D5CB7854025, 0x960067A5C064493E,
			0x41F420B5A91ED717, 0x6E81F9FA7E661D8D},
	},
	/* 13 */
	{
		x: gf.GF255e{0x56380CDCFD534861, 0xBF23AAE1C53EC995,
			0xEA53DCFF88E21C0D, 0x44019A922EE00C9D},
		u: gf.GF255e{0x1A722773E489DDB8, 0xF94983893D4AABD6,
			0x59F3D4C5BB3DFDCC, 0x653EE371D2801E6A},
	},
	/* 14 */
	{
		x: gf.GF255e{0x1AEE789343714BC2, 0x5027559B5A242FA0,
			0xD35C383123797CC4, 0x7E1BA15EC794EA34},
		u: gf.GF255e{0x460E164D09693F50, 0x96FABF7744D22EC2,
			0x216A1928595E868E, 0x50E1BEE9AC402680},
	},
	/* 15 */
	{
		x: gf.GF255e{0x0C6BCDDDE880E3FE, 0x9D52DF1A80CA1D01,
			0x6D28A55FD9814394, 0x07CD186F1749AC3B},
		u: gf.GF255e{0x442C914C64C6EE61, 0x5486463AAA3D41CD,
			0x2323BA05744FB271, 0x5CE94782B63D2983},
	},
	/* 16 */
	{
		x: gf.GF255e{0xF633589461CE1D8F, 0x43656EE4CD988663,
			0x601D9906BC58B752, 0x585EFD1D9E157C20},
		u: gf.GF255e{0x20DDE2D9560BD063, 0x68337F979386B815,
			0x9CAE33A6B5F9B94C, 0x0F2ED8418B17674E},
	},
}

var do255eWin_G195_xu = [16]do255ePointAffine{
	/* 1 */
	{
		x: gf.GF255e{0xC9D74325CBE23870, 0x3C8F799792EB5981,
			0xB4D3EBF20E4ACC9D, 0x5D505019F9CD4639},
		u: gf.GF255e{0xC6C4EA52864B8022, 0x3FACF03027F2FE05,
			0x5A78F8FDAFE0F2B2, 0x7A2058682117A352},
	},
	/* 2 */
	{
		x: gf.GF255e{0xB5FFEA9757F8DCE6, 0x992EAADD0950F49E,
			0xC30EE566764DB296, 0x77BEE9FA736A26DC},
		u: gf.GF255e{0xE7A9C455FDCCE69F, 0xB043C24E23C52866,
			0xCBC1DD8A3179B032, 0x597FE7EC4E366C38},
	},
	/* 3 */
	{
		x: gf.GF255e{0xF4762F3FF4FB8115, 0x60D22515C8F29371,
			0xE64A746FA6B9C81A, 0x107BB7D6EED2E10E},
		u: gf.GF255e{0xEBF1C192782AD7E7, 0xDAC867CE0228990B,
			0xB0C0AEB839C2A9BB, 0x5D529C2B2E3222F2},
	},
	/* 4 */
	{
		x: gf.GF255e{0xA0C19926FC389D68, 0xD93A4C18F4C2CCCD,
			0x09EB7E080DF8E02C, 0x5CD950BF71F691C5},
		u: gf.GF255e{0xE5AD1FDA050CE7CF, 0xD5179DCB398FE9E2,
			0x880F0F9CA2B23DE9, 0x73E9DA1D7C583AB6},
	},
	/* 5 */
	{
		x: gf.GF255e{0x97A1BBD1E5AA8841, 0x09B24BDC0DC8FFB6,
			0xA57657C8DCE4DE79, 0x228ECFC1B5307822},
		u: gf.GF255e{0x81533FD0CED00FE0, 0x2B41B323457375B0,
			0x3428954D0B0B6412, 0x3FB05C6B656FCDE7},
	},
	/* 6 */
	{
		x: gf.GF255e{0xF03AEB994C6B5021, 0xB0156B2F0414CE7E,
			0x64B75C8B8346FCE3, 0x4F54E4B6B9C3FD25},
		u: gf.GF255e{0xABECE8DEAA4DEFF3, 0xA6B25F5370FF8BED,
			0xC70C1F018B95875D, 0x1EEE7F4380019FB8},
	},
	/* 7 */
	{
		x: gf.GF255e{0x9854F37D39978A1E, 0xA8401C2863D1E85D,
			0x021EDD635FDF6914, 0x317884D08D246053},
		u: gf.GF255e{0xC2513D53CE5A6CD2, 0x8AF5B5BD4C9ADB58,
			0xDF748C1856292D78, 0x1C54D437C147EB47},
	},
	/* 8 */
	{
		x: gf.GF255e{0x09F1B77F26BC9F8E, 0x219F33E838C90E61,
			0x320CDCEF213ECF00, 0x67131909DEA4A881},
		u: gf.GF255e{0x961441BFA4853698, 0x76396E28425429D3,
			0x23187D9C49399AB4, 0x47EC89341C754A72},
	},
	/* 9 */
	{
		x: gf.GF255e{0x1BE2A51DB2720321, 0x645B8A7DD6376B36,
			0x686F85695A8133B5, 0x07E6F607B6D91397},
		u: gf.GF255e{0x790DDCB656DA5EBB, 0x899CA30F8A6C1157,
			0xB055E943E160FF52, 0x0C4E4B67E97A3F02},
	},
	/* 10 */
	{
		x: gf.GF255e{0x5940A4043221C587, 0x5FA59E201799740F,
			0x73197D2BCF61AC84, 0x303417596A2CF352},
		u: gf.GF255e{0x7CBAD6A204D1F6B9, 0xEB725D50999A4399,
			0x7C80807D104B0670, 0x44B8942C5DF07889},
	},
	/* 11 */
	{
		x: gf.GF255e{0x97C6A80A64F51D16, 0xEAC3C45ABBB9912B,
			0x72FB0626EAAD16C4, 0x5DA3C2E5773227DE},
		u: gf.GF255e{0x30FE96716E7DF796, 0xB6A214969A98317F,
			0xEB5423DB7543C3F6, 0x7DD81F0AD475BB65},
	},
	/* 12 */
	{
		x: gf.GF255e{0x0FE33479735C7A13, 0xA6F8DF8C4AC20F15,
			0xD19D08B9A74DDF08, 0x7D923FBD990B2A82},
		u: gf.GF255e{0x041881DC4BB15593, 0x0620628690F070A8,
			0xF6647FFFA6239BFD, 0x60406418E0A8E484},
	},
	/* 13 */
	{
		x: gf.GF255e{0x9D7AE8C91CE6917B, 0xB5DB6A7F7DEA58CC,
			0x9B1C372C39EEA275, 0x58742A848CC4EB79},
		u: gf.GF255e{0x02A5C58CE250DF40, 0x80A9A62960313C1E,
			0xEDC81102A9A286D9, 0x03C8B0610E5DE932},
	},
	/* 14 */
	{
		x: gf.GF255e{0xDD44633AC9F7D040, 0x2632460F3A278E77,
			0x092C7F41C75E959E, 0x00D4C6E23ED1D27A},
		u: gf.GF255e{0x2359C265188EE74D, 0x3498D65BDB7611FC,
			0x97D9FFD6286A6BD1, 0x63F775224E36165F},
	},
	/* 15 */
	{
		x: gf.GF255e{0x515E61CB24DAF0E2, 0x606D332C7B076CF7,
			0x4792A1E865D47BC0, 0x5A839AAC1191AB83},
		u: gf.GF255e{0x14BD04EFD34FC573, 0x58F89E267B42BEA3,
			0xDD4D47FA083CC9BF, 0x5C69FCC38DA29629},
	},
	/* 16 */
	{
		x: gf.GF255e{0xF04AA6935B6381F5, 0xA05C72D30E2E31F2,
			0x5414DCBFD3642F59, 0x75A8FE2604285EC3},
		u: gf.GF255e{0x805E028481B3D10D, 0x3E6A069F39FCEFCB,
			0xDA636B907FED771B, 0x162581D9B675A4E1},
	},
}
