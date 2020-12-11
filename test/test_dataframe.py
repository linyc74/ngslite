import pandas as pd
from ngslite.dataframe import left_join, outer_join, inner_join
from .setup import TestCase


class TestLeftJoin(TestCase):

    def test_success(self):
        left = pd.DataFrame({
            'C1': [1, 2, 3],
            'C2': [4, 5, 6],
        })

        right = pd.DataFrame({
            'C1': [1, 2, 3],
            'C3': [7, 8, 9],
        })

        merged = left_join(
            left=left,
            right=right,
            on=None
        )

        expected = pd.DataFrame({
            'C1': [1, 2, 3],
            'C2': [4, 5, 6],
            'C3': [7, 8, 9],
        })

        self.assertTrue(merged.equals(expected))

    def test_catch_length_error(self):
        left = pd.DataFrame({
            'C1': [1, 2, 3],
            'C2': [4, 5, 6],
        })

        right = pd.DataFrame({
            'C1': [1, 2, 3, 3],
            'C3': [7, 8, 9, 0],
        })

        with self.assertRaises(AssertionError) as context:
            left_join(
                left=left,
                right=right,
                on='C1'
            )
            msg = 'merged length (4) != left length 3'
            self.assertEqual(msg, str(context.exception))


class TestOuterJoin(TestCase):

    def test_success(self):
        left = pd.DataFrame({
            'C1': [1, 2],
            'C2': [3, 4],
        })

        right = pd.DataFrame({
            'C1': [2, 3],
            'C3': [5, 6],
        })

        merged = outer_join(
            left=left,
            right=right,
            on=None
        )

        expected = pd.DataFrame({
            'C1': [1, 2, 3],
            'C2': [3, 4, None],
            'C3': [None, 5, 6],
        })

        self.assertTrue(merged.equals(expected))


class TestInnerJoin(TestCase):

    def test_success(self):
        left = pd.DataFrame({
            'C1': [1, 2],
            'C2': [3, 4],
        })

        right = pd.DataFrame({
            'C1': [2, 3],
            'C3': [5, 6],
        })

        merged = inner_join(
            left=left,
            right=right,
            on=None
        )

        expected = pd.DataFrame({
            'C1': [2],
            'C2': [4],
            'C3': [5],
        })

        self.assertTrue(merged.equals(expected))
