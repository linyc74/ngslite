import pandas as pd
from typing import List, Union, Optional


def left_join(
        left: pd.DataFrame,
        right: pd.DataFrame,
        on: Optional[Union[str, List[str]]] = None) -> pd.DataFrame:

    if on is None:
        on = list(set(left.columns).intersection(set(right.columns)))

    n = len(left)

    df = left.merge(
        right=right,
        on=on,
        how='left'
    )

    assert len(df) == n, f'merged length ({len(df)}) != left length {n}'

    return df


def outer_join(
        left: pd.DataFrame,
        right: pd.DataFrame,
        on: Optional[Union[str, List[str]]] = None) -> pd.DataFrame:

    if on is None:
        on = list(set(left.columns).intersection(set(right.columns)))

    df = left.merge(
        right=right,
        on=on,
        how='outer'
    )

    return df


def inner_join(
        left: pd.DataFrame,
        right: pd.DataFrame,
        on: Optional[Union[str, List[str]]] = None) -> pd.DataFrame:

    if on is None:
        on = list(set(left.columns).intersection(set(right.columns)))

    df = left.merge(
        right=right,
        on=on,
        how='inner'
    )

    return df
