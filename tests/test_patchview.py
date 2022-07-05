# -*- coding: utf-8 -*-
"""Tests for `patchview` package."""

import random

import pytest

from patchview import patchview


@pytest.fixture
def generate_numbers():
    """Sample pytest fixture. Generates list of random integers.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """

    return random.sample(range(100), 10)


def test_sum_numbers(generate_numbers):
    """Sample test function for sum_numbers, using pytest fixture."""

    our_result = patchview.sum_numbers(generate_numbers)
    assert our_result == sum(generate_numbers)


def test_max_number(generate_numbers):
    """Sample test function for max_number, using pytest fixture."""

    our_result = patchview.max_number(generate_numbers)
    assert our_result == max(generate_numbers)


# def test_max_number_bad(generate_numbers):
#     """Sample test function that fails. Uncomment to see."""
#
#     our_result = patchview.max_number(generate_numbers)
#     assert our_result == max(generate_numbers) + 1
