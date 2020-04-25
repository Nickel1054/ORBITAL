from SpaceRock import SpaceRockJOV
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# (name, a, e, i, w, lo, t0)


Io = SpaceRockJOV()
Io.set('Io', 421.8, 0.0041, 0.036, 84.129, 43.977, 0.)

Europa = SpaceRockJOV()
Europa.set('Europa', 671.1, 0.094, 0.466, 88.970, 219.106, 0.)

Gan = SpaceRockJOV()
Gan.set('Ganymede', 1070.4, 0.0013, 0.177, 192.417, 63.552, 0.)

Call = SpaceRockJOV()
Call.set('Callisto', 1882.7, 0.0074, 0.192, 52.643, 298.848, 0.)

xi = []
yi = []
zi = []

xe = []
ye = []
ze = []

xg = []
yg = []
zg = []

xc = []
yc = []
zc = []

for i in range(24*17):
    ri = Io.get_r(i)
    re = Europa.get_r(i)
    rg = Gan.get_r(i)
    rc = Call.get_r(i)

    xi.append(ri[0])
    yi.append(ri[1])
    zi.append(ri[2])
    xe.append(re[0])
    ye.append(re[1])
    ze.append(re[2])
    xg.append(rg[0])
    yg.append(rg[1])
    zg.append(rg[2])
    xc.append(rc[0])
    yc.append(rc[1])
    zc.append(rc[2])


fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')
ax1.scatter(xi, yi, zi, s=1, c='r', marker="o", label='Io')
ax1.scatter(xe, ye, ze, s=1, c='b', marker="o", label='Europa')
ax1.scatter(xg, yg, zg, s=1, c='y', marker="o", label='Ganymede')
ax1.scatter(xc, yc, zc, s=1, c='#999999', marker="o", label='Callisto')
ax1.scatter(0, 0, 0, s=10, c='#000000', marker="o", label='Jupiter')
ax1.set(xlim=(-2000, 2000), ylim=(-2000, 2000), zlim=(-2000, 2000))
plt.legend(loc='upper left')
plt.show()
