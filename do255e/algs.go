package do255e

import (
	"crypto"
	cryptorand "crypto/rand"
	"encoding/binary"
	"errors"
	"golang.org/x/crypto/sha3"
	"io"
)

// This file implements high-level operations over do255e:
//
//   - Key pair generation
//   - Key exchange (ECDH)
//   - Signature generation and verification
//   - Hash-to-curve

// A private key structure contains a private key, i.e. a non-zero
// scalar for do255e. For efficiency reasons, it internally caches a
// copy of the public key as well.
type Do255ePrivateKey struct {
	d    Do255eScalar
	pub  Do255ePoint
	epub [32]byte
}

// A public key structure contains a non-neutral group element.
type Do255ePublicKey struct {
	pub  Do255ePoint
	epub [32]byte
}

// Test whether a public key is equal to another.
func (pk Do255ePublicKey) Equal(other crypto.PublicKey) bool {
	pk2, ok := other.(Do255ePublicKey)
	if !ok {
		return false
	}
	var t byte = 0
	for i := 0; i < 32; i++ {
		t |= pk.epub[i] ^ pk2.epub[i]
	}
	return t == 0
}

// Decode a private key from bytes. This function expects exactly
// 32 bytes. If the provided slice does not have length exactly 32,
// or if what it contains is not the canonical encoding of a valid
// non-zero scalar for do255e, then this function returns nil and an
// error.
func Do255eDecodePrivateKey(src []byte) (*Do255ePrivateKey, error) {
	if len(src) != 32 {
		return nil, errors.New("Invalid private key")
	}
	sk := new(Do255ePrivateKey)
	if sk.d.Decode(src) != 1 {
		return nil, errors.New("Invalid private key")
	}
	sk.pub.MulGen(&sk.d)
	sk.pub.Encode(sk.epub[:0])
	return sk, nil
}

// Encode a private key into bytes. The private key (exactly 32 bytes)
// is appended to the provided slice. If 'dst' has enough capacity, then
// it is returned; otherwise, a new slice is allocated, and receives
// the concatenation of the current contents of 'dst' and the encoded
// private key.
func (sk *Do255ePrivateKey) Encode(dst []byte) []byte {
	return sk.d.Encode(dst)
}

// Get the public key corresponding to a given private key.
func (sk *Do255ePrivateKey) Public() *Do255ePublicKey {
	pk := new(Do255ePublicKey)
	pk.pub.Set(&sk.pub)
	copy(pk.epub[:], sk.epub[:])
	return pk
}

// Decode a public key from bytes. This function expects exactly
// 32 bytes. If the provided slice does not have length exactly 32,
// or if what it contains is not the canonical encoding of a valid
// non-neutral do255e element, then this function returns nil and an
// error.
func Do255eDecodePublicKey(src []byte) (*Do255ePublicKey, error) {
	if len(src) != 32 {
		return nil, errors.New("Invalid public key")
	}
	pk := new(Do255ePublicKey)
	if pk.pub.Decode(src) != 1 {
		return nil, errors.New("Invalid public key")
	}
	copy(pk.epub[:], src)
	return pk, nil
}

// Encode a public key into bytes. The public key (exactly 32 bytes)
// is appended to the provided slice. If 'dst' has enough capacity, then
// it is returned; otherwise, a new slice is allocated, and receives
// the concatenation of the current contents of 'dst' and the encoded
// public key.
func (pk *Do255ePublicKey) Encode(dst []byte) []byte {
	n := len(dst)
	n2 := n + 32
	var dst2 []byte
	if cap(dst) >= n2 {
		dst2 = dst[:n2]
	} else {
		dst2 = make([]byte, n2)
		copy(dst2, dst)
	}
	copy(dst2[n:], pk.epub[:])
	return dst2
}

// Key pair generation with do255e: from a random source 'rand', a
// private key (a scalar value) and the corresponding public key (group
// element) are generated. The random source MUST be cryptographically
// secure. If 'rand' is nil, then crypto/rand.Reader is used (this is
// the recommended way).
func Do255eGenerateKeyPair(rand io.Reader) (*Do255ePrivateKey, error) {
	// We obtain 32 bytes from the random source and decode them
	// as a scalar with modular reduction. The group order is very
	// close to 2^254 (difference is less than 2^127) so the bias
	// is negligible.
	if rand == nil {
		rand = cryptorand.Reader
	}
	var bb [32]byte
	if _, err := io.ReadFull(rand, bb[:]); err != nil {
		return nil, err
	}
	sk := new(Do255ePrivateKey)
	sk.d.DecodeReduce(bb[:])

	// In the utterly improbable case that we got a zero here, we
	// replace it with 1.
	var bb2 [32]byte
	bb2[0] = byte(sk.d.IsZero())
	var d2 Do255eScalar
	d2.Decode(bb2[:])
	sk.d.Add(&sk.d, &d2)

	// Compute public key.
	sk.pub.MulGen(&sk.d)
	sk.pub.Encode(sk.epub[:0])

	return sk, nil
}

// Key exchange with do255e: given our private key, and the public key
// from the peer, a shared secret of length 'len' bytes is produced. The
// peer's public key is provided encoded; it should have length exactly
// 32 bytes. If the provided sequence of bytes has not length exactly 32
// bytes, or if it is not otherwise a valid do255e point encoding, then
// the key exchange fails. On failure, a byte sequence of the requested
// length is still generated; that byte sequence is not predictable by
// outsiders, and cannot be distinguished from the output of a
// successful ECDH exchange by outsiders. This is meant to support rare
// protocols where exchanged keys are not public, and the exchange
// should not have any validation semantics. The 'ok' returned value has
// value 1 on success, 0 on error (an 'int' is used to promote
// constant-time processing).
func Do255eKeyExchange(sk *Do255ePrivateKey, peer_pk []byte, secretLen int) (secret []byte, ok int) {
	// Decode peer point; remember if it worked.
	var P Do255ePoint
	var venc int
	var eppk [32]byte
	if len(peer_pk) == 32 {
		venc = P.Decode(peer_pk)
		copy(eppk[:], peer_pk)
	} else {
		P.Neutral()
		venc = -1
	}
	// Input point is valid if and only if it decoded properly AND it
	// was not the neutral (i.e. venc == 1, not 0 or -1).
	ok = int(-uint32(venc) >> 31)

	// ECDH
	var P2 Do255ePoint
	P2.Mul(&P, &sk.d)

	// Generate secret as SHAKE256 over the concatenation of:
	//  - a domain separation string
	//  - a byte with value 0x00 (success) or 0xFF (failure)
	//  - the encoded shared point (on success) or our private key
	//    (on failure)
	//  - our public key and the peer public key, both encoded;
	//    first one is the numerically lowest, when values are
	//    interpreted as integers in unsigned little-endian convention
	sh := sha3.NewShake256()
	sh.Write([]byte("do255e-ecdh:"))

	// Get bb[0] == 0x00 on success, 0xFF otherwise.
	var bb [33]byte
	var esk [32]byte
	bb[0] = byte(ok) - 1
	// bb[1..32] is either the shared point, or our own private key.
	P2.EncodeSquaredW(bb[1:1])
	sk.Encode(esk[:0])
	for i := 0; i < 32; i++ {
		bb[1+i] ^= bb[0] & (bb[1+i] ^ esk[i])
	}
	sh.Write(bb[:])

	// Get the two public keys in numerical order. Our public key is
	// in sk.epub; the peer public key is in eppk.
	var cc uint = 0
	for i := 0; i < 32; i++ {
		x := uint(sk.epub[i]) - uint(eppk[i]) - cc
		cc = (x >> 8) & 1
	}
	// If cc == 1, then sk.epub < eppk; otherwise, cc == 0 and
	// sk.epub >= eppk.
	var bb2 [64]byte
	m := byte(-cc)
	for i := 0; i < 32; i++ {
		bb2[i] = (sk.epub[i] & m) | (eppk[i] & ^m)
		bb2[32+i] = (sk.epub[i] & ^m) | (eppk[i] & m)
	}
	sh.Write(bb2[:])

	// Generate and return secret.
	secret = make([]byte, secretLen)
	io.ReadFull(sh, secret)
	return
}

// Get hash OID string for the provided hash identifier.
func getHashOID(opts crypto.SignerOpts) (oid []byte, err error) {
	switch opts.HashFunc() {
	case crypto.Hash(0):
		oid = []byte("")
	case crypto.SHA224:
		oid = []byte("2.16.840.1.101.3.4.2.4")
	case crypto.SHA256:
		oid = []byte("2.16.840.1.101.3.4.2.1")
	case crypto.SHA384:
		oid = []byte("2.16.840.1.101.3.4.2.2")
	case crypto.SHA512:
		oid = []byte("2.16.840.1.101.3.4.2.3")
	case crypto.SHA512_224:
		oid = []byte("2.16.840.1.101.3.4.2.5")
	case crypto.SHA512_256:
		oid = []byte("2.16.840.1.101.3.4.2.6")
	case crypto.SHA3_224:
		oid = []byte("2.16.840.1.101.3.4.2.7")
	case crypto.SHA3_256:
		oid = []byte("2.16.840.1.101.3.4.2.8")
	case crypto.SHA3_384:
		oid = []byte("2.16.840.1.101.3.4.2.9")
	case crypto.SHA3_512:
		oid = []byte("2.16.840.1.101.3.4.2.10")
	default:
		return nil, errors.New("Unknown hash identifier")
	}
	return oid, nil
}

// Compute the "challenge" in a signature process.
func mkChallenge(dk *Do255eScalar, Renc []byte, epub []byte, hashOID []byte, data []byte) {
	sh := sha3.NewShake256()
	sh.Write([]byte("do255e-sign-e:"))
	sh.Write(Renc)
	sh.Write(epub)
	sh.Write(hashOID)
	sh.Write([]byte(":"))
	sh.Write(data)
	var challenge [32]byte
	io.ReadFull(sh, challenge[:])
	dk.DecodeReduce(challenge[:])
}

// Schnorr signature with do255e. The data to sign ('data') may be either
// raw data, or a hash value. The 'opts' parameter specifies the hash
// function that was used to pre-hash the data (use crypto.Hash(0) for
// raw data).
//
// The signature process is deterministic, for a given 'seed' value. The
// seed may be nil, which is equivalent to a seed of length 0. How the
// seed contents were chosen does not impact the algorithmic security of
// the signature (using a non-repeating seed may increase implementation
// robustness against some classes of physical attacks).
//
// The signature is returned as a newly allocated slice. Its length is
// exactly 64 bytes. An error is reported if the hash function
// identified by 'opts' is not known.
func (sk *Do255ePrivateKey) signWithSeed(seed []byte, data []byte, opts crypto.SignerOpts) (signature []byte, err error) {
	// Get hash OID.
	var hashOID []byte
	hashOID, err = getHashOID(opts)
	if err != nil {
		return nil, err
	}

	// Generate the per-signature k value.
	sh := sha3.NewShake256()
	sh.Write([]byte("do255e-sign-k:"))
	var bb [32]byte
	sk.d.Encode(bb[:0])
	sh.Write(bb[:32])
	var tt [8]byte
	if seed == nil {
		sh.Write(tt[:8])
	} else {
		binary.LittleEndian.PutUint64(tt[:], uint64(len(seed)))
		sh.Write(tt[:8])
		sh.Write(seed)
	}
	sh.Write(hashOID)
	sh.Write([]byte(":"))
	sh.Write(data)
	io.ReadFull(sh, bb[:])
	var k Do255eScalar
	k.DecodeReduce(bb[:])

	// R = k*G
	var R Do255ePoint
	R.MulGen(&k)
	var R_enc [32]byte
	R.Encode(R_enc[:0])

	// Compute challenge e.
	var e Do255eScalar
	mkChallenge(&e, R_enc[:], sk.epub[:], hashOID, data)

	// s = k + e*d
	var s Do255eScalar
	s.Mul(&e, &sk.d).Add(&s, &k)

	// Signature is (R,s)
	signature = make([]byte, 64)
	copy(signature[:32], R_enc[:])
	s.Encode(signature[:32])
	return
}

// Schnorr signature with do255e. The data to sign ('data') may be either
// raw data, or a hash value. The 'opts' parameter specifies the hash
// function that was used to pre-hash the data (use crypto.Hash(0) for
// raw data).
//
// If 'rand' is nil, then the signature is deterministic (this is safe).
// If 'rand' is not nil, then 32 bytes are read from it and used to
// complement the internal per-signature nonce generation process,
// making the signature non-deterministic, in case a specific protocol
// requires this property. Non-deterministic signatures might also
// improve implementation robustness against some kinds of physical
// attacks (in particular fault attacks). It is not necessary that the
// extra randomness returned by 'rand' has high quality; the security of
// the signature will be maintained in all case, even if that data is
// fully predictable.
//
// The signature is returned as a newly allocated slice. Its length is
// exactly 64 bytes. An error is reported if 'rand' is not nil but a
// read attempt returns an error. An error is also reported if the hash
// function identified by 'opts' is not known.
func (sk *Do255ePrivateKey) Sign(rand io.Reader, data []byte, opts crypto.SignerOpts) (signature []byte, err error) {
	if rand == nil {
		return sk.signWithSeed(nil, data, opts)
	} else {
		var seed [32]byte
		_, err = io.ReadFull(rand, seed[:])
		if err != nil {
			return nil, err
		}
		return sk.signWithSeed(seed[:], data, opts)
	}
}

// Verify a signature on a message, relatively to a public key.
//
// The message data is provided in 'data'. This is interpreted as raw data
// if opts is crypto.Hash(0)); otherwise, it will be considered to be
// pre-hashed data, processed with the hash function identified by opts.
//
// Returned value is true if the hash function is recognized, and the
// signature is valid relatively to the provided public key. In all other
// cases, false is returned.
//
// This function is not constant-time, under the assumption that public
// keys and signatures are public data.
func (pk *Do255ePublicKey) VerifyVartime(data []byte, opts crypto.SignerOpts, sig []byte) bool {
	// Valid signatures have length 64 bytes exactly.
	if len(sig) != 64 {
		return false
	}

	// Sanity check: public keys cannot be the neutral point.
	if pk.pub.IsNeutral() != 0 {
		return false
	}

	// Get the hash function OID string. If unrecognized, then we
	// have to reject the signature.
	hashOID, _ := getHashOID(opts)
	if hashOID == nil {
		return false
	}

	// The signature splits into an encoded point (Renc) and a
	// scalar. Decode the scalar; reject the signature if that value
	// is not canonical. Note that we accept zero.
	Renc := sig[:32]
	var s Do255eScalar
	if s.Decode(sig[32:]) < 0 {
		return false
	}

	// Compute the challenge.
	var e Do255eScalar
	mkChallenge(&e, Renc, pk.epub[:], hashOID, data)

	// Signature is valid if and only if s*G - e*pk = R.
	e.Neg(&e)
	return pk.pub.VerifyHelperVartime(&s, &e, Renc)
}

// Hash some input data into a curve point. The input data ('data') is
// either raw or pre-hashed, as identified by the opts parameter. If
// the hash function identifier is not recognized, then an error is
// returned. Otherwise, the output point is returned.
func Do255eHashToCurve(data []byte, opts crypto.SignerOpts) (*Do255ePoint, error) {
	// Get internal identifier for the hash function.
	oid, err := getHashOID(opts)
	if err != nil {
		return nil, err
	}

	// Hash the input into 64 bytes with SHAKE.
	sh := sha3.NewShake256()
	sh.Write([]byte("do255e-hash-to-curve:"))
	sh.Write(oid)
	sh.Write([]byte(":"))
	sh.Write(data)
	var bb [64]byte
	io.ReadFull(sh, bb[:])

	// Map each half into a curve point, and add them together.
	var P1, P2 Do255ePoint
	P1.MapBytes(bb[:32])
	P2.MapBytes(bb[32:])
	return NewDo255ePoint().Add(&P1, &P2), nil
}
