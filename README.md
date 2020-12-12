# Go Implementation of do255e and do255s

This is a plain Go implementation of do255e and do255s. It is
considered secure; all relevant functions should be constant-time
on modern architectures.

The `do255e` and `do255s` packages provide APIs for:

  - the prime order group itself;
  - scalars (i.e. integers modulo the primer order of the group);
  - high-level functionalities: key pair generation, key exchange
    (Diffie-Hellman), signature generation and verification, and
    hash-to-curve.

The implementation requires Go-1.12 or later. It is portable; there
is no assembly or other architecture-dependent feature.
