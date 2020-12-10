import random
import numpy as np
from typing import Union, List


def replace_zero_with(
        np_array: np.ndarray,
        value: Union[int, float]) -> np.ndarray:
    """
    Replaces the 0s in an numpy array with a specified value

    Args:
        np_array:
            1D numpy.ndarray, dtype = np.int or np.float

        value

    Returns: 1D numpy.ndarray, dtype = np.float
    """
    bool_arr = np.equal(np_array, 0)
    val_arr = np.ones((len(np_array), ), dtype=np.float) * value
    return np_array + bool_arr * val_arr


def jitter(
        array: Union[List[Union[int, float]],
                     np.ndarray],
        value: float) -> np.ndarray:
    """
    Add the input <array> with some random numbers
    The random numbers are uniformly distributed between -<value> to +<value>

    Args:
        array:
            An array of numbers, could be list, np.ndarray with int or float data types

        value:
            The max amount of value to be jittered

    Returns: np.ndarray, dtype=np.float32
    """
    # Create a new instance of the input array with np.float32 data type
    array = np.array(array, dtype=np.float32)
    noise = [random.uniform(-value, value) for _ in range(len(array))]
    noise = np.array(noise, dtype=np.float32)
    return array + noise


def densCols(
        x: np.ndarray,
        y: np.ndarray,
        bins: int) -> np.ndarray:
    """
    This is an implementation of the densCols() function in R
    Python just doesn't have this feature!

    Args:
        x:
            x coordinates

        y:
            y coordinates

        bins:
            number of bins

    Returns: 1D numpy array, length = len(x), dtype = np.float64
        Local densities (i.e. counts) of each input 2D point
    """
    # 2D histogram matrix, dtype = np.float64
    hist2D = np.histogram2d(x, y, bins)[0]

    # The bin size is a scalar
    bin_size_x = (max(x) - min(x)) / bins
    bin_size_y = (max(y) - min(y)) / bins

    # For each x and y coordinate, get the position with regard to which bin
    x_pos = (x - min(x)) / bin_size_x
    y_pos = (y - min(y)) / bin_size_y

    # Round it down to integer, so that it could be indexed
    x_pos = x_pos.astype(np.uint, copy=False)
    y_pos = y_pos.astype(np.uint, copy=False)

    # The max of x_pos and y_pos == bins, which could not be indexed (out of bound)
    # Decrease by 1 so that it can be indexed as the last bin
    x_pos[np.argmax(x_pos)] = bins - 1
    y_pos[np.argmax(y_pos)] = bins - 1

    # For each point, find the local density, which is encoded in the hist2D matrix
    densities = hist2D[x_pos, y_pos]

    return densities
