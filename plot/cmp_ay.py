import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.interpolate import interpn


gE0 = 0.511875


def fname_base(run_id):
    return '../output/' + run_id + '/' + run_id


def read_ay():
# Tao's data
    ay01 = np.loadtxt("../output/p80x80/p80x802", skiprows=1)
    ay10 = np.loadtxt("../output/p80x80/p80x8020", skiprows=1)

    return ay01, ay10


def ay_coord():
    alpha_lc = 5
    alpha_max = 90

    logemin = np.log(0.2)
    logemax = np.log(5)

    alphav = np.linspace(alpha_lc, alpha_max, 80)
    logev = np.linspace(logemin, logemax, 80)

    return alphav, logev 


def ay_init():
    alphav = np.linspace(5, 90, 80)
    y05_0 = np.exp(-(0.5 - 0.2) / 0.1) * (np.sin(np.deg2rad(alphav)) - np.sin(np.deg2rad(5)))
    y20_0 = np.exp(-(2.0 - 0.2) / 0.1) * (np.sin(np.deg2rad(alphav)) - np.sin(np.deg2rad(5)))

    return y05_0, y20_0


def p2e(p, E0=gE0, cv=1):
    return np.sqrt(p**2 * cv**2 + E0**2) - E0


def e2p(E, E0=gE0, cv=1):
    return np.sqrt(E * (E + 2 * E0)) / cv


def read_xy(data):
    alphav = data['alpha0'][:]
    logEN = data['logEN'][:]

    return alphav, logEN


def f1d(fmat, alphav, ev, e):
    points = (alphav, ev)
    xi = np.zeros((alphav.size, 2))
    xi[:, 0] = alphav[:]
    xi[:, 1] = e

    return interpn(points, fmat, xi)


if __name__ == '__main__':

    run_id = 'AlbertYoung'

    data = h5py.File(fname_base(run_id) + '_data.h5', 'r')
    alphav, logEN = read_xy(data)

    f01_2d = data['f/1'][:]
    f10_2d = data['f/10'][:]

    f0501 = f1d(f01_2d, alphav, logEN + np.log(gE0), np.log(0.5)) * e2p(0.5)**2
    f0510 = f1d(f10_2d, alphav, logEN + np.log(gE0), np.log(0.5)) * e2p(0.5)**2

    f2001 = f1d(f01_2d, alphav, logEN + np.log(gE0), np.log(2.0)) * e2p(2.0)**2
    f2010 = f1d(f10_2d, alphav, logEN + np.log(gE0), np.log(2.0)) * e2p(2.0)**2

    ay_alphav, ay_logev = ay_coord()
    ay_01, ay_10 = read_ay()

    f0500_ay, f2000_ay = ay_init()
    f0501_ay = f1d(ay_01, ay_alphav, ay_logev, np.log(0.5))
    f0510_ay = f1d(ay_10, ay_alphav, ay_logev, np.log(0.5))

    f2001_ay = f1d(ay_01, ay_alphav, ay_logev, np.log(2.0))
    f2010_ay = f1d(ay_10, ay_alphav, ay_logev, np.log(2.0))

    fig = plt.figure(figsize=(6, 3))
    ax1 = fig.add_subplot(1, 2, 1)

    ax1.plot(ay_alphav, f0500_ay, color="black", label='T = 0.0day')
    ax1.plot(ay_alphav, f0501_ay, "o", mfc='none', color="C0")
    ax1.plot(ay_alphav, f0510_ay, "o", mfc='none', color="C1")
    ax1.plot(alphav, f0501, color="C0", label='T = 0.1day')
    ax1.plot(alphav, f0510, color="C1", label='T = 1.0day')
    ax1.set(yscale='log')

    ax1.set_title("0.5 MeV")
    ax1.set_xlabel(r'$\alpha_0$ $(^\mathrm{o})$')
    ax1.set_ylabel(r'flux (arbitrary units)')
    ax1.legend(loc='best')
    ax1.set_xlim(0, 90)
    ax1.set_ylim(1e-3, 1)

    ax2 = fig.add_subplot(1, 2, 2)

    ax2.plot(ay_alphav, f2000_ay, color="black")
    ax2.plot(ay_alphav, f2001_ay, "o", mfc='none', color="C0", label='Albert Young 2005')
    ax2.plot(ay_alphav, f2010_ay, "o", mfc='none', color="C1")
    ax2.plot(alphav, f2001, color="C0", label='Sayram')
    ax2.plot(alphav, f2010, color="C1")
    ax2.set(yscale='log')

    ax2.set_title("2 MeV")
    ax2.set_xlabel(r'$\alpha_0$ $(^\mathrm{o})$')
    ax2.set_ylabel(r'flux (arbitrary units)')

    ax2.legend(loc='best')
    ax2.set_xlim(0, 90)
    ax2.set_ylim(1e-10, 1e-2)

    plt.tight_layout()
    plt.show()
