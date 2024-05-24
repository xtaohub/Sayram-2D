import numpy as np
import matplotlib.pyplot as plt
import meshio
from scipy.interpolate import griddata

E_MIN = 0.2  # MeV
E_MAX = 5  # MeV
E0 = 0.511875
c = 1.0
dlogE = (np.log(E_MAX) - np.log(E_MIN)) / (80 - 1)


def calP(Ec):
    return np.sqrt(Ec ** 2 + 2 * Ec * E0) / c


def locateE(Ec):
    return (np.log(Ec) - np.log(E_MIN)) / dlogE


p1 = calP(0.5)
p2 = calP(2.0)
P_MIN = calP(E_MIN)
P_MAX = calP(E_MAX)

M = 50
N = 50
x1 = np.linspace(5, 90, M)
x2 = np.linspace(5 + 85 / (2 * N), 90 - 85 / (2 * N), N)

# TODO Tao's data
with open("p80x80/p80x802") as dT01:
    dataT01_ = dT01.readlines()
del (dataT01_[0])

with open("p80x80/p80x8020") as dT10:
    dataT10_ = dT10.readlines()
del (dataT10_[0])

alphav = np.linspace(5, 90, 80)
y_05_0 = np.exp(-(0.5 - 0.2) / 0.1) * (np.sin(alphav * np.pi / 180) - np.sin(5.0 * np.pi / 180))
y_20_0 = np.exp(-(2.0 - 0.2) / 0.1) * (np.sin(alphav * np.pi / 180) - np.sin(5.0 * np.pi / 180))


with open("output/SMPPFV_10s/smppfv0") as raw0:
    d0 = raw0.readlines()

data0 = np.zeros((50, 50))

for i in range(50):
    for j in range(50):
        data0[i, j] = float(d0[50*i + j])

pos_p1 = ((p1-P_MIN - (P_MAX - P_MIN) / 100) / (P_MAX - P_MIN) * 50)
loc_p1 = int(pos_p1)
w1 = pos_p1 - loc_p1
w2 = 1 - w1
y_0s05 = (data0[loc_p1] * 1/w1 / (1/w1 + 1/w2) + data0[loc_p1 + 1] * 1/w2 / (1/w1 + 1/w2)) * p1**2

pos_p2 = ((p2-P_MIN - (P_MAX - P_MIN) / 100) / (P_MAX - P_MIN) * 50)
loc_p2 = int(pos_p2)
w1 = pos_p2 - loc_p2
w2 = 1 - w1
y_0s20 = (data0[loc_p2] * 1/w1 / (1/w1 + 1/w2) + data0[loc_p2 + 1] * 1/w2 / (1/w1 + 1/w2)) * p2**2


e1 = 0.5
pos1 = int((calP(e1) - P_MIN) / (P_MAX - P_MIN) * 400)
e2 = 2.0
pos2 = int((calP(e2) - P_MIN) / (P_MAX - P_MIN) * 400)

xx_range = np.linspace(5, 90, 400)


dataT01 = np.zeros((80, 80), dtype=float)
dataT10 = np.zeros((80, 80), dtype=float)
for i in range(80):
    temp01 = dataT01_[i].split()
    temp10 = dataT10_[i].split()
    for j in range(80):
        dataT01[i, j] = float(temp01[j])
        dataT10[i, j] = float(temp10[j])

loc05 = int(np.floor(locateE(0.5)))
w05 = locateE(0.5) - loc05
loc20 = int(np.floor(locateE(2.0)))
w20 = locateE(2.0) - loc20


y_05_t01 = dataT01[:, loc05] * (1/w05 / (1/w05 + 1/(1-w05))) + dataT01[:, loc05 + 1] * (1 - 1/w05 / (1/w05 + 1/(1-w05)))
y_05_t10 = dataT10[:, loc05] * (1/w05 / (1/w05 + 1/(1-w05))) + dataT10[:, loc05 + 1] * (1 - 1/w05 / (1/w05 + 1/(1-w05)))
y_20_t01 = dataT01[:, loc20] * (1/w20 / (1/w20 + 1/(1-w20))) + dataT01[:, loc20 + 1] * (1 - 1/w20 / (1/w20 + 1/(1-w20)))
y_20_t10 = dataT10[:, loc20] * (1/w20 / (1/w20 + 1/(1-w20))) + dataT10[:, loc20 + 1] * (1 - 1/w20 / (1/w20 + 1/(1-w20)))


fig = plt.figure(figsize=(14, 5))

ax1 = fig.add_subplot(1, 2, 1)

l1, = plt.semilogy(alphav, y_05_0, color="black", label='T = 0.0day')
l2, = plt.semilogy(alphav, y_05_t01, color="blue", label='T = 0.1day')
l3, = plt.semilogy(alphav, y_05_t10, color="red", label='T = 1.0day')
l4, = plt.semilogy(x2, y_0s05, "b--", label='0.0d')


plt.title("0.5MeV", fontsize=16)
plt.ylabel('Structured Mesh \n\n flux', fontsize=16)
plt.legend(handles=[l1, l2, l3],
           labels=['T = 0.0 day', 'T = 0.1 day', 'T = 1.0 day'], loc='best', fontsize=16)
plt.tick_params(pad=10, labelsize=16)
plt.xlim(0, 90)
plt.ylim(1e-3, 1)

ax2 = fig.add_subplot(1, 2, 2)


l1, = plt.semilogy(alphav, y_20_0, color="black", label='T = 0.0 day')
l2, = plt.semilogy(alphav, y_20_t01, color="blue", label='layer method')
l3, = plt.semilogy(alphav, y_20_t10, color="red", label='T = 1.0day')
l4, = plt.semilogy(x2, y_0s20, "b--", label='PPFV,SM')


plt.title("2MeV", fontsize=16)
plt.legend(handles=[l2, l4],
           labels=['layer method', 'PPFV'], loc='best', fontsize=16)
plt.tick_params(pad=10, labelsize=16)
plt.xlim(0, 90)
plt.ylim(1e-10, 1e-2)


plt.tight_layout()
plt.show()
