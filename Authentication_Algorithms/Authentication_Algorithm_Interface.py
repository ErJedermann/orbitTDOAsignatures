from TLEcalculator import TLEcalculator
from astropy import coordinates as coord
import abc

class Authentication_Algorithm(abc.ABC):
    @abc.abstractmethod
    def __init__(self, valid_sat_names: [str]):
        pass

    @abc.abstractmethod
    def minimal_receivers(self):
        pass

    @abc.abstractmethod
    # returns (a,b,c)
    # a: True/False if an satellite is authenticated.
    # b: The index (i) of the satellite, that is authenticated.
    # c: The estimated position (in TEME) of the signal source, at the last time. (Interesting for attacks.)
    def authenticate_satellite(self, receivers: [coord.ITRS], receiving_jDay: [(float, float)], TDOA_measurements: [[float]], tle_calculator: TLEcalculator) -> (bool, int, (float, float, float)):
        pass