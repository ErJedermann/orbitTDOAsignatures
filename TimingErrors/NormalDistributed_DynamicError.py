from __future__ import annotations
from TimingErrors.TimingErrorInterface import TimingError

import numpy as np


class NormalDistributed_DynamicError(TimingError):
    # A gaussian distributed error with a given mean and a given standard deviation (both in sec).
    # The error is dynamic (changes with every call of get_error()).

    def __init__(self, mean: float = 0.0, std: float = 10e-9):
        self.standard_deviation = std
        self.mean = mean

    def __str__(self):
        return f"NormalErr(mean:{self.mean},std:{self.standard_deviation})"

    def get_error(self, amount):
        return np.random.normal(self.mean, self.standard_deviation, amount)

    def get_mean_std(self):
        return self.mean, self.standard_deviation
