# Understanding ECC cryptography

This a pure Python implementation of a generic multiplication and addition
operations used key generation, document signing, and
[Diffie–Hellman key exchange](https://en.wikipedia.org/wiki/Diffie%E2%80%93Hellman_key_exchange).

## Acknowledgment

This code would be imposable without
[Understanding Cryptography](https://www.springer.com/gp/book/9783642041006)
by Christof Paar and Jan Pelzl.
The validity of the generic implementation would be hard to prove without
NIST test vectors.

## Usage

```python

x3, y3 = point_add(x1, y1, x2, y2, prime, a)
x2, y2 = point_multiply(x1, y1, multiplicant, prime, a)
```

Above the X and Y coordinates are assumed to be on the curve and the prime and
a to belong to the curve domain parameters.

NIST domain parameters may be used ass follows:

```python
from functional_ecc import domain_params

p_p21_dp = domain_params["P-521"]


x3, y3 = point_add(x1, y1, x2, y2, p_p21_dp.prime, p_p21_dp.a)
```

New domain parameters and curves may be defined as follows.

```
from functional_ecc import domain_param


paar_17 = domain_param(prime=17, a=2, b=2, g_x=5, g_y=1)
```

## Supported curves

This implementation should work for any prime field elliptic curve, however the
following NIST curve domain parameters can be found in the domain_params
dictionary:

* P-192
* P-224
* P-256
* P-384
* P-521
