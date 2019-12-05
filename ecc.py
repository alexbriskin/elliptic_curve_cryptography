#!/usr/bin/env python3
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


def calculate_sum(x1, y1, x2, y2, prime, a):
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
    _calculate_sum = partial(calculate_sum, prime=prime, a=a)
    x2, y2 = x1, y1
    # Skip first '1' digit and start double and add sequence
    for digit in bin(multiplicant)[3:]:
        x2, y2 = _calculate_sum(x2, y2, x2, y2)
        if digit == "1":
            x2, y2 = _calculate_sum(x1, y1, x2, y2)

    return x2, y2


domain_param = namedtuple("domain_param", "prime a b G_x G_y n")
domain_params = {
    "p-192": domain_param(
        prime=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF,
        b=0x64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1,
        a=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFC,
        G_x=0x188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012,
        G_y=0x07192B95FFC8DA78631011ED6B24CDD573F977A11E794811,
        n=0xFFFFFFFFFFFFFFFFFFFFFFFF99DEF836146BC9B1B4D22831,
    )
}


if __name__ == "__main__":
    # value = ecc_test_vector('/home/alex/Documents/Crypto/TestVectors/ECDSA/KeyPair.rsp')
    for curve_name, curve_description, vector in ecc_test_vector(
        "/home/alex/Documents/Crypto/TestVectors/ECDSA/KeyPair.rsp"
    ):
        _, y_res = point_multiply(G_x, G_y, n, prime_192, a_192)
        assert y_res is None
        d, Qx, Qy = itemgetter("d", "Qx", "Qy")(vector)
        calculated_point = point_multiply(G_x, G_y, d, prime_192, a_192)
        assert calculated_point == (
            Qx,
            Qy,
        ), f"{curve_name} - {calculated_point} != {(Qx, Qy)}"
    # d = 0xE5CE89A34ADDDF25FF3BF1FFE6803F57D0220DE3118798EA
    # Qx = 0x8ABF7B3CEB2B02438AF19543D3E5B1D573FA9AC60085840F
    # Qy = 0xA87F80182DCD56A6A061F81F7DA393E7CFFD5E0738C6B245
    # expected = point_multiply(G_x, G_y, d, prime_192, a_192)
    # assert expected == (Qx, Qy)

#     for multiplicant in range(1, 23):
#         print(f"{multiplicant:-2d}: {point_multiply(5, 1, multiplicant, 17, 2)}")
#         Q = Q + P
#         print(Q)
