# ORBITAL
This project contains a class for calculating an orbit of an object around the Sun or any other massive object. 
File SpaceRock.py uses 7 Keplerian elements of orbit to calculate and plot the orbit around the central object, which mass is written as a parameter (kg).

mat_vec(self, mat, vec) is used to multiply matrix and vector (for further transformation of frame of reference to inertial - toICS)

julian_date(self, year, month, day, hour=12, minute=0, sec=0) - is used to transfer Gregorian date to Julian date. 

get_E(self, t, eps) is used to calculate the eccentric anomaly for a moment in time using the iteration method.

get_r(self, t) calculates the r vector for a moment in time in the inertial frame of reference.

G and get_G2 are methods to calculate hiperbolic, not eccliptic orbit. Are still in progress, NOT WORKING.


main.py contains Keplerian parameters written for 4 inner planets of Solar System and code to plot and show it. 

Jovian.py makes the same for Galilean moons of Jupiter.
