import unittest
import numpy as np
from ngslite.arrfunc import replace_zero_with


class TestArrfunc(unittest.TestCase):

    def test_replace_zero_with(self):
        arr1 = replace_zero_with(
            np.array([0, 1, 0, 2, 0, 3, 0]), 100
        )
        arr2 = np.array([100, 1, 100, 2, 100, 3, 100])
        diff = arr1 - arr2
        self.assertEqual(np.sum(diff), 0)


if __name__ == '__main__':
    unittest.main()
