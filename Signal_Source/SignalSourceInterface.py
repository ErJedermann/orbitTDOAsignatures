import abc
from astropy.time import Time
import astropy.coordinates as coord


class SignalSource(abc.ABC):
    @abc.abstractmethod
    def __init__(self, auth: bool):
        pass

    @abc.abstractmethod
    def get_positions_TEME(self, receivers: [coord.ITRS], measurement_times: [Time]):
        pass

    @abc.abstractmethod
    def is_authentic(self):
        pass

    @abc.abstractmethod
    def set_deviation(self, range: float):
        pass

    @abc.abstractmethod
    def get_names(self):
        pass
