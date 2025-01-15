import numpy as np
import matplotlib.pyplot as plt
import configparser
import sys
import h5py
from scipy.interpolate import interpn
from matplotlib.lines import Line2D


gE0 = 0.511875


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
        paras['nalpha0'] = int(config.get('basic', 'nalpha0'))
        paras['nE'] = int(config.get('basic', 'nE'))
        # paras['ALPHA_MAX']=90
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

def read_xy(data):
    alphav = data['alpha0'][:]
    logEN = data['logEN'][:]

    return alphav, logEN

def f1d(fmat, alphav, ev, e):
    points = (alphav, ev)
    xi = np.zeros((alphav.size, 2))
    xi[:,0] = alphav[:]
    xi[:,1] = e 

    return interpn(points, fmat, xi)

if __name__ == '__main__':

    run_id = 'AlbertYoung'
    paras = parse_ini(run_id)
    
    data = h5py.File(fname_base(run_id) + '_data.h5', 'r')
    alphav, logEN = read_xy(data)

    
    f01_2d = data['f/1'][:]
    f10_2d = data['f/10'][:]
    
    f0501 = f1d(f01_2d, alphav, logEN + np.log(gE0), np.log(0.5)) * e2p(0.5)**2
    f0510 = f1d(f10_2d, alphav, logEN + np.log(gE0), np.log(0.5)) * e2p(0.5)**2
    
    zero_ = np.array([0])
    
    f0501 = np.concatenate((zero_, f0501))
    f0510 = np.concatenate((zero_, f0510))

    f2001 = f1d(f01_2d, alphav, logEN + np.log(gE0), np.log(2.0)) * e2p(2.0)**2
    f2010 = f1d(f10_2d, alphav, logEN + np.log(gE0), np.log(2.0)) * e2p(2.0)**2
    
    f2001 = np.concatenate((zero_, f2001))
    f2010 = np.concatenate((zero_, f2010))

    ay_alphav, ay_logev = ay_coord()
    ay_01, ay_10 = read_ay()

    f0500_ay, f2000_ay = ay_init()
    f0501_ay = f1d(ay_01, ay_alphav, ay_logev, np.log(0.5))
    f0510_ay = f1d(ay_10, ay_alphav, ay_logev, np.log(0.5))

    f2001_ay = f1d(ay_01, ay_alphav, ay_logev, np.log(2.0))
    f2010_ay = f1d(ay_10, ay_alphav, ay_logev, np.log(2.0))
    
    alphalc_ = np.array([5])
    alphav = np.concatenate((alphalc_, alphav))

    fig = plt.figure(figsize=(6, 3))
    ax1 = fig.add_subplot(1, 2, 1)

    l1, = plt.semilogy(ay_alphav, f0500_ay, color="black", label='T = 0.0day')
    l2, = plt.semilogy(ay_alphav, f0501_ay, "x", color="C0", label='T = 0.1day', markevery=2, markeredgewidth=2)
    l3, = plt.semilogy(ay_alphav, f0510_ay, "x", color="C1", label='T = 1.0day', markevery=2, markeredgewidth=2)

    l4, = plt.semilogy(alphav, f0501, color="C0", label='T = 0.1day')
    l5, = plt.semilogy(alphav, f0510, color="C1", label='T = 1.0day')

    plt.title("0.5 MeV")
    plt.xlabel(r'$\alpha_0$ $(^\mathrm{o})$')
    plt.ylabel(r'flux (arbitrary units)')
    plt.legend(handles=[l1, l4, l5],
               labels=['T = 0.0 day', 'T = 0.1 day', 'T = 1.0 day'], loc='best')
    plt.xlim(0, 90)
    plt.ylim(1e-3, 1)

    ax2 = fig.add_subplot(1, 2, 2)

    l1, = plt.semilogy(ay_alphav, f2000_ay, color="black", label='T = 0.0 day')
    l2, = plt.semilogy(ay_alphav, f2001_ay, "x", color="C0", label='layer method', markevery=2, markeredgewidth=2)
    l3, = plt.semilogy(ay_alphav, f2010_ay, "x", color="C1", label='T = 1.0day', markevery=2, markeredgewidth=2)
    l4, = plt.semilogy(alphav, f2001, color="C0", label='3DPPFV')
    l5, = plt.semilogy(alphav, f2010, color="C1", label='1d')

    plt.title("2 MeV")
    plt.xlabel(r'$\alpha_0$ $(^\mathrm{o})$')
    plt.ylabel(r'flux (arbitrary units)')

    proxy_lines = [Line2D([0], [0], linestyle='none', marker='x', color='black', markevery=2, markeredgewidth=2),
                   Line2D([0], [0], color='black')]

    plt.legend(proxy_lines, ['Albert and Young', 'Sayram'], loc='best')
    plt.xlim(0, 90)
    plt.ylim(1e-10, 1e-2)

    plt.tight_layout()
    plt.show()

