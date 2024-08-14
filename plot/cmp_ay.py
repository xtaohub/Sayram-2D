import numpy as np
import matplotlib.pyplot as plt
import configparser
import sys
from scipy.interpolate import interpn

def fname_base(run_id):
    return '../output/' + run_id + '/' + run_id

def parse_ini(run_id):
    config = configparser.ConfigParser()
    ini_fname= fname_base(run_id) + '.ini'

    config.read(ini_fname)

    paras = {}

    try:
        paras['E_MIN'] = float(config.get('basic', 'Emin'))  # MeV
        paras['E_MAX'] = float(config.get('basic', 'Emax'))  # MeV
        paras['nx'] = int(config.get('basic', 'nalpha'))
        paras['ny'] = int(config.get('basic', 'nE'))
        paras['ALPHA_LC'] = float(config.get('basic', 'alpha0_lc'))
        paras['ALPHA_MAX']=90
    except Exception as e:
        print("section_name or option_name wrong, check the input file.")
        sys.exit(1)

    return paras

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

def p2e(p, E0=0.511875, cv=1):
    return np.sqrt(p**2* cv**2 + E0**2) - E0

def e2p(E, E0=0.511875, cv=1):
    return np.sqrt(E * (E + 2 * E0)) / cv 

def read_xy(run_id):
    xfname = fname_base(run_id) + '_x.dat'
    yfname = fname_base(run_id) + '_y.dat'
    
    alphav = np.loadtxt(xfname)
    p = np.loadtxt(yfname)

    alphav = np.rad2deg(alphav)
    ev = p2e(p)

    return alphav, ev

def f1d(fmat, alphav, ev, e):
    points = (alphav, ev)
    xi = np.zeros((alphav.size, 2))
    xi[:,0] = alphav[:]
    xi[:,1] = e 

    return interpn(points, fmat, xi)

if __name__ == '__main__':

    run_id = 'albert_young'
    alphav, ev = read_xy(run_id)

    fnamed01 = fname_base(run_id) + '1'   
    fnamed10 = fname_base(run_id) + '10'

    f01 = np.loadtxt(fnamed01)
    f10 = np.loadtxt(fnamed10)

    f0501 = f1d(f01, alphav, ev, 0.5) * e2p(0.5)**2
    f0510 = f1d(f10, alphav, ev, 0.5) * e2p(0.5)**2

    f2001 = f1d(f01, alphav, ev, 2.0) * e2p(2.0)**2
    f2010 = f1d(f10, alphav, ev, 2.0) * e2p(2.0)**2

    ay_alphav, ay_logev = ay_coord()
    ay_01, ay_10 = read_ay()

    f0500_ay, f2000_ay = ay_init()
    f0501_ay = f1d(ay_01, ay_alphav, ay_logev, np.log(0.5))
    f0510_ay = f1d(ay_10, ay_alphav, ay_logev, np.log(0.5))

    f2001_ay = f1d(ay_01, ay_alphav, ay_logev, np.log(2.0))
    f2010_ay = f1d(ay_10, ay_alphav, ay_logev, np.log(2.0))

    plt.close()
    fig,axes = plt.subplots(1, 2, figsize = (9, 4))

    axes[0].plot(ay_alphav, f0500_ay, color='k')
    axes[0].plot(ay_alphav, f0501_ay, color='C0')
    axes[0].plot(ay_alphav, f0510_ay, color='C1')
    axes[0].plot(alphav, f0501, ls = '--', color='C0')
    axes[0].plot(alphav, f0510, ls = '--', color='C1')

    axes[0].set(yscale='log', ylim=(1e-5, 1e1))
    # axes[0].set(yscale='log')

    axes[1].plot(ay_alphav, f2000_ay, color='k')
    axes[1].plot(ay_alphav, f2001_ay)
    axes[1].plot(ay_alphav, f2010_ay)

    axes[1].plot(alphav, f2001, ls = '--', color='C0')
    axes[1].plot(alphav, f2010, ls = '--', color='C1')

    axes[1].set(yscale='log', ylim=(1e-10, 1e-2))
    # axes[1].set(yscale='log')
   
    plt.show()

