import unittest


from ngslite.arrfunc import *


class TestArrfunc(unittest.TestCase):
    def test_replace_zero_with(self):
        arr_in = np.array([0, 1, 0, 2, 0, 3, 0])
        arr_in = replace_zero_with(arr_in, 100)
        arr_out = np.array([100, 1, 100, 2, 100, 3, 100])
        for i, j in zip(arr_in, arr_out):
            self.assertEqual(i, j)

