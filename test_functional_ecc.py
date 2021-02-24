# -*- coding: utf-8 -*-
"""test_functional_ecc

    Test the functional implementation of ECC multiplication and addition.
"""
from operator import itemgetter
from parse import parse
from sympy import mod_inverse

from functional_ecc import point_add, point_multiply, at_infinity, domain_params


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


def test_key_gen_vectors():
    for curve_name, curve_description, vector in ecc_test_vector("KeyPair.rsp"):
        prime, a, _, G_x, G_y, n = domain_params[curve_name]
        _, y_res = point_multiply(G_x, G_y, n, prime, a)
        assert y_res is at_infinity
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
