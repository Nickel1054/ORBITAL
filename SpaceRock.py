from math import pi, sqrt, sin, cos, sinh, cosh


class SpaceRock:
    def __init__(self, mass=2e30):
        self.__mu = 6.678e-11 * 86400 * 86400 / 1e18 * mass   # [a.u.^3 / JD^2]
        self.__name = ''
        self.__a = None  # semimajor axis
        self.__b = self.__a
        self.__e = None  # Eccentricity
        self.__i = None  # Inclination
        self.__w = None  # Argument of periapsis
        self.__lo = None  # Longitude of the ascending node
        self.__t0 = None  # epoch
        self.__toICS = None # matrix orbital -> inertial

    def set(self, name, a, e, i, w, lo, t0):
        self.__name = name
        self.__a = a                  # semimajor axis
        self.__e = e                  # Eccentricity
        if self.__e < 1:
            self.__b = a*sqrt(1 - self.__e*self.__e)
        elif self.__e > 1:
            self.__b = a * sqrt(self.__e * self.__e - 1)
        self.__i = i * pi / 180       # Inclination
        self.__w = w * pi / 180       # Argument of periapsis
        self.__lo = lo * pi / 180     # Longitude of the ascending node
        self.__t0 = t0                # epoch
        self.__toICS = [
            [cos(self.__lo)*cos(self.__w) - sin(self.__lo)*cos(self.__i)*sin(self.__w),
             -cos(self.__lo)*sin(self.__w) - sin(self.__lo)*cos(self.__i)*cos(self.__w),
             sin(self.__lo) * sin(self.__i)],
            [sin(self.__lo) * cos(self.__w) + cos(self.__lo) * cos(self.__i) * sin(self.__w),
             -sin(self.__lo) * sin(self.__w) + cos(self.__lo) * cos(self.__i) * cos(self.__w),
             -cos(self.__lo) * sin(self.__i)],
            [sin(self.__i)*sin(self.__w),
             sin(self.__i)*cos(self.__w),
             cos(self.__i)]
                      ]

    def __mat_vec(self, mat, vec):
        r = []
        for i in range(3):
            z = 0
            for j in range(3):
                z += mat[i][j] * vec[j]
            r.append(z)
        return r

    def __str__(self):
        return (self.__name + ' is an object with orbital parameters:\n'
                            'a: ' + str(self.__a) + '\n'
                            'e: ' + str(self.__e) + '\n'
                            'i: ' + str(self.__i * 180 / pi) + '\n'
                            'w: ' + str(self.__w * 180 / pi) + '\n'
                            'W: ' + str(self.__lo * 180 / pi) + '\n'
                            't0: ' + str(self.__t0) + '\n')

    def julian_date(self, year, month, day, hour=12, minute=0, sec=0):
        a = (14 - month) // 12
        y = year + 4800 - a
        m = month + 12 * a - 3
        jdn = day + (153 * m + 2) // 5 + 365 * y + y // 4 - y // 100 + y // 400 - 32045
        jd = jdn + (hour - 12) / 24 + minute / 1440 + sec / 86400
        return jd

    def __get_E(self, t, eps=1e-7):
        n = sqrt(self.__mu/(self.__a*self.__a*self.__a))
        Es = n*(t-self.__t0)
        En = n*(t-self.__t0) + self.__e*sin(Es)
        while abs(Es - En) > eps:
            Es = En
            En = n*(t-self.__t0) + self.__e*sin(Es)
        return (En + Es)/2

    def __G(self, t, n, x):     # function, which is being put in the bisection method
        return self.__e * sinh(x) - x - n * (t - self.__t0)

    def sign(self, x):
        if x > 0:
            return 1
        elif x < 0:
            return -1
        else:
            return 0

    def get_G2(self, t, x1, x2, eps=1e-7):      # bisection method
        n = sqrt(self.__mu / (self.__a * self.__a * self.__a))
        xi = 0
        while abs(x1 - x2) > eps:
            xi = x1 + (x2 - x1)/2
            if self.sign(self.__G(t, n, x1)) != self.sign(self.__G(t, n, xi)):
                x2 = xi
            else:
                x1 = xi
        return xi

    def get_r(self, t):
        if self.__e < 1:
            E = self.__get_E(t)
            r0 = [self.__a * (cos(E) - self.__e),
                  self.__b * sin(E),
                  0]
            return self.__mat_vec(self.__toICS, r0)

        elif self.__e > 1:                              # something wrong
            G = self.get_G2(t, -100., 100.)
            r0 = [self.__a * (cosh(G) - self.__e),
                  self.__b * sinh(G),
                  0]
            return self.__mat_vec(self.__toICS, r0)


class SpaceRockGEO(SpaceRock):
    def __init__(self):
        super().__init__(6e24)


class SpaceRockJOV(SpaceRock):
    def __init__(self):
        super().__init__(1.8982e27)

