from ITRSconverter import LatLonHeight_2_ITRS, Cartesian_2_ITRS
from astropy import coordinates as coord
from astropy import units as u
import random
import numpy as np

class BenchmarkSetupReceiver:
    def __init__(self):
        oxford = LatLonHeight_2_ITRS(51.750411, -1.282607, 72)
        thun = LatLonHeight_2_ITRS(46.759, 7.63, 560)
        kaiserslautern = LatLonHeight_2_ITRS(49.4268, 7.608899, 251)
        madrid = LatLonHeight_2_ITRS(40.4378698, -3.8196193, 667)
        rome = LatLonHeight_2_ITRS(41.909986, 12.3959154, 37)
        warsaw = LatLonHeight_2_ITRS(52.232855, 20.9211129, 113)
        paris = LatLonHeight_2_ITRS(48.8588377, 2.2770206, 130)
        athens = LatLonHeight_2_ITRS(37.990832, 23.70332, 277)
        stockholm = LatLonHeight_2_ITRS(59.3260668, 17.8419723, 0)
        budapest = LatLonHeight_2_ITRS(47.4811281, 18.9902212, 102)
        vienna = LatLonHeight_2_ITRS(48.2205998, 16.2399782, 157)
        marseilles = LatLonHeight_2_ITRS(43.2803051, 5.2404141, 12)
        edinburgh = LatLonHeight_2_ITRS(55.9411885, -3.2753779, 5)
        lisbon = LatLonHeight_2_ITRS(38.7436214, -9.1952226, 100)
        palermo = LatLonHeight_2_ITRS(38.1405228, 13.2872489, 14)
        sarajevo = LatLonHeight_2_ITRS(43.8937019, 18.312952, 511)
        prague = LatLonHeight_2_ITRS(50.0595854, 14.3255424, 192)
        oslo = LatLonHeight_2_ITRS(59.8937806, 10.6450364, 23)
        amsterdam = LatLonHeight_2_ITRS(52.3545983, 4.8339212, 0)
        barcelona = LatLonHeight_2_ITRS(41.3947688, 2.0787283, 12)
        bucharest = LatLonHeight_2_ITRS(44.43792, 26.0245, 70)
        self.europe_receivers = [oxford, thun, kaiserslautern, madrid, rome, warsaw, paris, athens, stockholm, budapest,
                                 vienna, marseilles, edinburgh, lisbon, palermo, sarajevo, prague, oslo, amsterdam,
                                 barcelona, bucharest]

    def get_random_receivers(self, amount: int, min_diameter: float = None, max_diameter: float = None) -> [coord.ITRS]:
        center_location = self.__select_random_europe_receivers(1)
        return self.__create_new_receiver_locations(amount, min_diameter, max_diameter, center_location[0], 10)

    def __select_random_europe_receivers(self, amount: int) -> [coord.ITRS]:
        output = []
        while len(output) < amount:
            index = int(random.random() * len(self.europe_receivers))
            candidate = self.europe_receivers[index]
            if candidate not in output:
                output.append(candidate)
        return output

    def __get_pos_in_m(self, input: coord.ITRS):
        if input.cartesian.xyz.unit is u.km:
            return input.cartesian.xyz.value * 1000
        else:
            return input.cartesian.xyz.value

    # distributes receivers on a circle line
    def __create_new_receiver_locations(self, amount: int, diameter_min: int, diameter_max: int,
                                        center_location: coord.ITRS, height_diff: float=0) -> [coord.ITRS]:
        radius_min = diameter_min/2
        radius_max = diameter_max/2
        radius_range = radius_max - radius_min
        angle_step = 360/amount
        angle_ramdomness = angle_step * 0.2 # the angle is randomly influenced by ±20% of the step-size
        c_lon = center_location.earth_location.geodetic.lon.value
        c_lat = center_location.earth_location.geodetic.lat.value
        c_height = center_location.earth_location.geodetic.height.value
        c_pos = center_location.cartesian.xyz.value
        if center_location.cartesian.xyz.unit is u.km:
            c_pos = c_pos * 1000 # km in m
        output = []
        while len(output) < amount:
            # select a random point in 2D based on  min < r [m] < max  and  between 0° to 360° (0° is north)
            temp_radius = random.random() * radius_range + radius_min
            temp_angle = len(output) * angle_step + ((random.random() * angle_ramdomness * 2) - angle_ramdomness)
            # calculate the shift in north and east in meter
            shift_north = np.cos(np.radians(temp_angle)) * temp_radius
            shift_east = np.sin(np.radians(temp_angle)) * temp_radius
            if height_diff is 0:
                up_shift = 0
            else:
                up_shift = random.random() * height_diff - height_diff/2
            # convert this east-north-up shift to ITRS coordinates:
            x_shift, y_shift, z_shift = self.__east_north_up_shift_2_xyz_shift(center_location, shift_north, shift_east, up_shift)
            # use this to create a new location
            temp_pos = [c_pos[0]+x_shift, c_pos[1]+y_shift, c_pos[2]+z_shift]
            temp_pos = np.array(temp_pos)
            temp_pos = temp_pos / 1000 # m in km
            temp_location = Cartesian_2_ITRS(temp_pos, (0,0,0))
            output.append(temp_location)
        return output

    def __east_north_up_shift_2_xyz_shift(self, start_pos: coord.ITRS, north_shift: float, east_shift: float, up_shift):
        # north_shift and east_shift in [m]
        # from https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
        phi = start_pos.earth_location.geodetic.lat.value
        lam = start_pos.earth_location.geodetic.lon.value
        delta_x = - np.sin(np.radians(lam)) * east_shift \
                  - np.cos(np.radians(lam)) * np.sin(np.radians(phi)) * north_shift \
                  + np.cos(np.radians(lam)) * np.cos(np.radians(phi)) * up_shift
        delta_y = + np.cos(np.radians(lam)) * east_shift \
                  - np.sin(np.radians(lam)) * np.sin(np.radians(phi)) * north_shift \
                  + np.sin(np.radians(lam)) * np.cos(np.radians(phi)) * up_shift
        delta_z = 0 + np.cos(np.radians(phi)) * north_shift + np.sin(np.radians(phi)) * up_shift
        return delta_x, delta_y, delta_z # in [m]

    def get_centroid_pos(self, locations: [coord.ITRS]) -> (float, float, float):
        size1 = len(locations)
        positions = np.zeros((size1, 3))
        for i in range(len(locations)):
            pos = locations[i].cartesian.xyz.value
            positions[i] = pos
        sum = np.sum(positions, 0)
        middle_point = sum / size1
        # [x,y,z] values in km
        return middle_point

