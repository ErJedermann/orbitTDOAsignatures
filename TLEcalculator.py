from sgp4.api import Satrec, SatrecArray, SGP4_ERRORS
from sgp4.api import days2mdhms, jday

import numpy as np

# The SGP4 propagator returns raw x,y,z Cartesian coordinates in a “True Equator Mean Equinox” (TEME)
# reference frame that’s centered on the Earth but does not rotate with it — an “Earth centered inertial” (ECI)
# reference frame.
# The satellite deviate from the ideal orbits described in TLE files about 1–3 km/day.

# The purpose of this class is to load a set of satellites from TLEs and calculate their position at a given point in
# time. The returned positions are in TEME frame.
class TLEcalculator:
    def __init__(self, tleFile: str, warnings: bool=True, verbose: bool=True):
        self.tleFile = tleFile
        self.warnings = warnings
        self.verbose = verbose
        self.valid_days = 7
        self.__parseFile()

    def __parseFile(self):
        file = open(self.tleFile, "r")
        lines = file.readlines()
        elements = int (len(lines) / 3)
        sat_list = []
        sat_names = []
        for i in range(elements):
            line_name = lines[i*3]
            line_one = lines[i*3 + 1]
            line_two = lines[i*3 + 2]
            if line_name[-1] is "\n":
                line_name = line_name[:-1]
                line_name = line_name.lstrip()
                line_name = line_name.rstrip()
            if line_one[-1] is "\n":
                line_one = line_one[:-1]
            if line_two[-1] is "\n":
                line_two = line_two[:-1]
            tempSat = Satrec.twoline2rv(line_one, line_two)
            sat_list.append(tempSat)
            sat_names.append(line_name)
        self.sat_names = sat_names
        self.satList = sat_list
        self.satrec = SatrecArray(sat_list)
        if self.verbose:
            print(f"INFO:TLEcalculator: {len(self.satList)} sats parsed")

    def get_min_max_epoch(self) -> ((int, float), (int, float)):
        minYr = 999
        minDay = 999
        maxYr = 000
        maxDay = 0
        for sat in self.satList:
            tempYr = sat.epochyr
            tempDay = sat.epochdays
            if tempYr < minYr:
                minYr = tempYr
                minDay = tempDay
            elif tempYr == minYr and tempDay < minDay:
                minDay = tempDay
            if tempYr > maxYr:
                maxYr = tempYr
                maxDay = tempDay
            elif tempYr == maxYr and tempDay > maxDay:
                maxDay = tempDay
        return (minYr, minDay), (maxYr, maxDay)

    def calculate_position_single(self, satellite: Satrec, jDay: float, jDayF: float):
        err, pos, vel = satellite.sgp4(jDay, jDayF)
        if err is not 0 and self.verbose:
            print(f"INFO:TLEcalculator.calc_pos_single: error {err}: '{SGP4_ERRORS.get(err)}'. pos:{pos}, "
                  f"sat_epoch:{satellite.jdsatepoch+satellite.jdsatepochF}, target_epoch:{jDay+jDayF}")
        return err, pos, vel

    def calculate_positions_all(self, jTimeDay: float, jTimeFr: float=0.0):
        if self.warnings and abs(jTimeDay - self.satList[0].jdsatepoch) > self.valid_days:
            print(f"WARNING: TLEcalculator.calculate_position: Large difference between TLE-epoch and given time.")
        jd = np.array([jTimeDay])
        fr = np.array([jTimeFr])
        err, pos, vel = self.satrec.sgp4(jd, fr)
        return pos, vel

    def calculate_multiple_positions_all(self, jTimeDay: np.ndarray, jTimeFr: np.ndarray):
        err, pos, vel = self.satrec.sgp4(np.array(jTimeDay), np.array(jTimeFr))
        return pos, vel

    def abs_time_to_jDay(self, year: int, month: int, day: int, hour:int, minute: int, second: int) -> (float, float):
        return jday(year, month, day, hour, minute, second)

    def relative_time_to_jDay(self, seconds_to_add: float) -> (float, float):
        # The seconds are added to the first epoch-entry of the sat-list
        epochJday = self.satList[0].jdsatepoch
        epochjDayF = self.satList[0].jdsatepochF
        day2add = seconds_to_add / (24*3600)
        day2add += epochjDayF
        epochJday += int(day2add)
        epochjDayF = day2add - int(day2add)
        return epochJday, epochjDayF

    def epoch_2_jDay(self, epochYr: int, epochDay: float) -> (float, float):
        # epoch[month, day, hour, minute, second]
        epoch = days2mdhms(epochYr, epochDay)
        # add constant 2000 to the year. This is not the best solution (epoches before 2000 will be corrupted)
        jDate = jday(epochYr+2000, epoch[0], epoch[1], epoch[2], epoch[3], epoch[4])
        return jDate





