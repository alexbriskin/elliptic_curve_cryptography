# -*- coding: utf-8 -*-
"""alex_ecc_experiment.ecc

    My attempt att ecc
"""
from collections import namedtuple
from functools import partial
from itertools import chain
from operator import itemgetter
from parse import parse
from sympy import mod_inverse
import pytest


# @pytest.fixture
def ecc_test_vector(vector_file):
    def parse_vectors_for_curve(lines):
        def parse_gen_info(lines):
            for line in lines:
                if not line.strip() or line.startswith("#  "):
                    continue
                result = parse("[{}]", line)
                assert result
                return result[0]

        def parse_num_vectors(lines):
            result = parse("N = {}", next(lines))
            assert result
            return result[0]

        def parse_test_vector(lines, variables):
            assert not next(lines).strip()
            vector = {}
            for expected_name, line in zip(variables, lines):
                name, value = parse("{} = {}", line)
                assert name == expected_name
                vector[name] = int(value, 16)
            return vector

        while True:
            curve_name = parse_gen_info(lines)
            curve_descr = parse_gen_info(lines)
            num_vectors = parse_num_vectors(lines)
            for _ in range(int(num_vectors)):
                test_vector = parse_test_vector(lines, "d Qx Qy".split())
                assert test_vector
                yield curve_name, curve_descr, test_vector

    with open(vector_file) as raw_vectors:
        yield from parse_vectors_for_curve(raw_vectors)


def calculate_slope(x1, y1, x2, y2, prime, a):
    assert x1 == x1 % prime
    assert y1 == y1 % prime
    assert x2 == x2 % prime
    assert y2 == y2 % prime
    return (
        mod_inverse(2 * y1, prime) * (3 * x1 ** 2 + a) % prime
        if (x1, y1) == (x2, y2)
        else mod_inverse(x2 - x1, prime) * (y2 - y1) % prime
    )


def point_add(x1, y1, x2, y2, prime, a):
    """ Semantics of point at infinity is y coordinate is None"""
    if (x1, y1) == (x2, -1 * y2 % prime):
        return x1, None
    if y1 is None:
        return x2, y2
    if y2 is None:
        return x1, y1
    slope = calculate_slope(x1, y1, x2, y2, prime, a)
    x3 = (slope ** 2 - x1 - x2) % prime
    y3 = (slope * (x1 - x3) - y1) % prime
    return x3, y3


def point_multiply(x1, y1, multiplicant, prime, a):
    assert multiplicant >= 1
    x2, y2 = x1, y1
    # Skip first '1' digit and start double and add sequence
    for digit in bin(multiplicant)[3:]:
        x2, y2 = point_add(x2, y2, x2, y2, prime, a)
        if digit == "1":
            x2, y2 = point_add(x1, y1, x2, y2, prime, a)

    return x2, y2


domain_param = namedtuple("domain_param", "prime a b G_x G_y n")
domain_params = {
    "P-192": domain_param(
        prime=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF,
        b=0x64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1,
        a=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFC,
        G_x=0x188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012,
        G_y=0x07192B95FFC8DA78631011ED6B24CDD573F977A11E794811,
        n=0xFFFFFFFFFFFFFFFFFFFFFFFF99DEF836146BC9B1B4D22831,
    ),
    "P-224": domain_param(
        prime=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF000000000000000000000001,
        a=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFE,
        b=0xB4050A850C04B3ABF54132565044B0B7D7BFD8BA270B39432355FFB4,
        G_x=0xB70E0CBD6BB4BF7F321390B94A03C1D356C21122343280D6115C1D21,
        G_y=0xBD376388B5F723FB4C22DFE6CD4375A05A07476444D5819985007E34,
        n=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFF16A2E0B8F03E13DD29455C5C2A3D,
    ),
    "P-256": domain_param(
        prime=0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF,
        a=0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC,
        b=0x5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B,
        G_x=0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296,
        G_y=0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5,
        n=0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551,
    ),
    "P-384": domain_param(
        prime=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFF,
        a=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFF0000000000000000FFFFFFFC,
        b=0xB3312FA7E23EE7E4988E056BE3F82D19181D9C6EFE8141120314088F5013875AC656398D8A2ED19D2A85C8EDD3EC2AEF,
        G_x=0xAA87CA22BE8B05378EB1C71EF320AD746E1D3B628BA79B9859F741E082542A385502F25DBF55296C3A545E3872760AB7,
        G_y=0x3617DE4A96262C6F5D9E98BF9292DC29F8F41DBD289A147CE9DA3113B5F0B8C00A60B1CE1D7E819D7A431D7C90EA0E5F,
        n=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFC7634D81F4372DDF581A0DB248B0A77AECEC196ACCC52973,
    ),
    "P-521": domain_param(
        prime=2 ** 521 - 1,
        a=2 ** 521 - 4,
        b=0x0051953EB9618E1C9A1F929A21A0B68540EEA2DA725B99B315F3B8B489918EF109E156193951EC7E937B1652C0BD3BB1BF073573DF883D2C34F1EF451FD46B503F00,
        G_x=0xC6858E06B70404E9CD9E3ECB662395B4429C648139053FB521F828AF606B4D3DBAA14B5E77EFE75928FE1DC127A2FFA8DE3348B3C1856A429BF97E7E31C2E5BD66,
        G_y=0x11839296A789A3BC0045C8A5FB42C7D1BD998F54449579B446817AFBD17273E662C97EE72995EF42640C550B9013FAD0761353C7086A272C24088BE94769FD16650,
        n=0x01FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA51868783BF2F966B7FCC0148F709A5D03BB5C9B8899C47AEBB6FB71E91386409,
    ),
}


def test_key_gen_vectors():
    for curve_name, curve_description, vector in ecc_test_vector(
        "/home/alex/Documents/Crypto/TestVectors/ECDSA/KeyPair.rsp"
    ):
        if curve_name not in domain_params:
            continue
        prime, a, _, G_x, G_y, n = domain_params[curve_name]
        _, y_res = point_multiply(G_x, G_y, n, prime, a)
        assert y_res is None
        d, Qx, Qy = itemgetter("d", "Qx", "Qy")(vector)
        calculated_point = point_multiply(G_x, G_y, d, prime, a)
        assert calculated_point == (
            Qx,
            Qy,
        ), f"{curve_name} - {calculated_point} != {(Qx, Qy)}"


def test_sign_gen_vectors():
    Msg = 0x699325D6FC8FBBB4981A6DED3C3A54AD2E4E3DB8A5669201912064C64E700C139248CDC19495DF081C3FC60245B9F25FC9E301B845B3D703A694986E4641AE3C7E5A19E6D6EDBF1D61E535F49A8FAD5F4AC26397CFEC682F161A5FCD32C5E780668B0181A91955157635536A22367308036E2070F544AD4FFF3D5122C76FAD5D
    d = 0x16797B5C0C7ED5461E2FF1B88E6EAFA03C0F46BF072000DFC830D615
    Qx = 0x605495756E6E88F1D07AE5F98787AF9B4DA8A641D1A9492A12174EAB
    Qy = 0xF5CC733B17DECC806EF1DF861A42505D0AF9EF7C3DF3959B8DFC6669
    k = 0xD9A5A7328117F48B4B8DD8C17DAE722E756B3FF64BD29A527137EEC0
    R = 0x2FC2CFF8CDD4866B1D74E45B07D333AF46B7AF0888049D0FDBC7B0D6
    S = 0x8D9CC4C8EA93E0FD9D6431B9A1FD99B88F281793396321B11DAC41EB
    prime, a, _, G_x, G_y, n = domain_params["P-224"]
    x_res, y_res = point_multiply(G_x, G_y, d, prime, a)
    assert (x_res, y_res) == (Qx, Qy)
    assert 1 <= k < n
    r, _ = point_multiply(G_x, G_y, k, prime, a)
    r = r % n
    assert r
    assert r == R
    from hashlib import sha224

    msg = Msg.to_bytes((Msg.bit_length() + 7) // 8, byteorder="big")
    z = int.from_bytes(sha224(msg).digest(), byteorder="big")
    s = (z + r * d) * mod_inverse(k, n) % n
    assert s
    assert s == S, f"{s:x} != {S:x}"
    assert (r, s) == (R, S)
    ################################################################################
    w = mod_inverse(s, n)
    u = z * w % n
    v = r * w % n
    x, y = point_add(
        *point_multiply(G_x, G_y, u, prime, a),
        *point_multiply(Qx, Qy, v, prime, a),
        prime,
        a,
    )
    assert r == x % n, f"{r:x} != {x % n:x}"


if __name__ == "__main__":
    pass
