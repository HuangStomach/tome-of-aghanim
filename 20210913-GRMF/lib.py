import numpy as np
from typing import TypeVar

P = TypeVar('P', bound='KNKN')
class KNKN:
    def __init__(self, k) -> None:
        self._k = k

    def fit(self, _X) -> P:
        self.X = _X.copy()
        return self
    
    def neighbors(self, i) -> np.ndarray:
        drugs_sim = self.X.iloc[i].to_numpy()
        values = np.sort(drugs_sim)[::-1][1:self._k+1]
        indexes = np.argsort(drugs_sim)[::-1][1:self._k+1]
        return (indexes, values)
