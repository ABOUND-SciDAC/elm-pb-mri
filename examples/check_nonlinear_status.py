#!/usr/bin/env python

"""Check status of the nonlinear simulation:
"""

from __future__ import (absolute_import, division,
                        print_function) #, unicode_literals)
import matplotlib as mpl
import matplotlib.pyplot as plt
from cycler import cycler
from subprocess import check_output
import os
import numpy as np

from boututils import file_import, parser_config, elmsize, deriv, save2nc
import boututils.functions as boutf
from boutdata import boutgrid, collect
import visualization as boutv

plt.ion()

# control parameters
loaddata = 1
savefig = 1
savedata = 1
calculation = loaddata
filter_n0 = 1   # True

var = 'P'
var0 = "{}0".format(var)    # equilibrium profiles
case = 'crash'  # case name
path = case     # case dmp file path
parts = (['rms_cs', 'rms_vs_t', 'data0_dc', 'fn', 'elmsize']
         if var == 'P' else ['rms_cs', 'rms_vs_t', 'fn'])

## elmsize calculation
norm = 1 # ['3f', '6f', float], normalization factor for `var`, used in elmsize

## plot setting
nlevels = 40    # nlevels for contourf plot
tevery = 100    # plot fn & p0_dc every t step

db_path = 'db/{}'.format(case.split('/')[-1])
if (savedata or savefig) and (not os.path.exists(db_path)):
    os.makedirs(db_path)

gridpath = "./"

## read BOUT.inp
boutinp = parser_config(path + '/BOUT.inp')
timestep = float(boutinp['timestep'])
zperiod = int(boutinp['zperiod'])
lowpass = int(boutinp['highbeta']['low_pass_z'])

grid = os.path.join(gridpath, boutinp['grid'].split('/')[-1])
if not os.path.exists(grid):
    raise ValueError('gridfile NOT exists:\n  {}'.format(grid))

grid = boutgrid(grid)

timeout = 0
strcase = '_' + case if len(case) else ''

## confirm controls
print("path: ", path)
print("var: ", var)

if loaddata:
    print("collecting data ...")
    t_array = collect('t_array', path=path)
    print("**P0**")
    #p0 = collect('P0', path=path)
    p0 = collect('P0', tind=[0, len(t_array)-1], path=path)
    if var0 != 'P0':
        if var0 == 'Ni0':
            var0 = 'N0'
        print("**{}**".format(var0))
        data0 = collect(var0, path=path)
    else:
        data0 = p0
    print("**{}**".format(var))
    data = collect(var, tind=[0, len(t_array)-1], path=path)
    nx, ny, nz, nt = data.shape

# get some time points, including last step
t_arr = boutf.nrange(1, nt-1, tevery)
t_array = t_array[:nt]

psi, iyind, oyind = grid.get_psin(yind='omp', index=True)
xpeak = deriv(psi, p0[:, oyind]).argmin()
print("peak position of Grad(p0): ind={}, psin={:.2f}".format(
    xpeak, psi[xpeak]))

if calculation:
    print("calculating rms & dc parts ...")
    rms = data.std(axis=-2).squeeze()
    dc = data.mean(axis=-2).squeeze()    
    elm_size = elmsize(dc, data0, grid, norm=norm)
    np.save('{}/p0.npy'.format(db_path),data0)
    np.save('{}/p.npy'.format(db_path),data)
    np.save('{}/dcP.npy'.format(db_path),dc)
    np.save('{}/rmsP.npy'.format(db_path),rms)
    np.save('{}/elmsize.npy'.format(db_path),elm_size['s3']*100.)

    if 'fn' in parts:
        print("calculating mode amplitude ...")
        if filter_n0:
            print("filter_n0: True")
            fp = np.fft.fft(data-dc.reshape(nx, ny, 1, nt), axis=-2)
            fpa = np.abs(fp)
        else:
            fp = np.fft.fft(data, axis=-2)
            fpa = np.abs(fp)

if 'rms_cs' in parts:
    print("========rms in contour/surface")
    tind = (nt - 1)*timestep
    title = "rms{}{}_t{}".format(var, strcase, tind)
    print("t = {} Ta".format(tind))
    print("title = ", title)
    plt.figure('rms_c')
    plt.clf()
    plt.contourf(psi, np.arange(ny), rms[:, :, -1].T, nlevels)
    plt.xlabel('$\psi_n$')
    plt.ylabel('poloidal')
    plt.title('{}$\\tau_A$'.format(title))
    # psi = 1.0
    plt.axvline(1.0, lw=1.5, color='w', ls='--')
    # outer mid-plane
    plt.axhline(oyind, lw=1.5, color='w', ls='--')
    plt.tight_layout()
    plt.show()
    if boutf.get_yesno('update fig', timeout=timeout):
        boutf.run_command()
    if savefig:
        boutv.savefig('{}/{}_c'.format(db_path, title))

    plt.figure('rms_s')
    plt.clf()
    boutv.surface(psi, np.arange(ny), rms[:, :, -1])
    plt.ylabel('$\psi_n$')
    plt.xlabel('poloidal')
    plt.title('{}$\\tau_A$'.format(title))
    plt.tight_layout()
    plt.show()
    if boutf.get_yesno('update fig', timeout=timeout):
        boutf.run_command()
    if savefig:
        boutv.savefig('{}/{}_s'.format(db_path, title))

if 'rms_vs_t' in parts:
    print("========plot rms_vs_t")
    plt.figure("rms_vs_t")
    plt.clf()
    plt.plot(np.arange(nt)*timestep, rms[xpeak, oyind, :])
    plt.title("rms{}, x{}y{}".format(var, xpeak, oyind))
    plt.xlabel(r'time/$\tau_A$')
    plt.yscale('log')
    plt.tight_layout()
    plt.show()
    if savedata:
        print("saving rms data ...")
        save2nc("{}/rms{}_x{}y{}.nc".format(db_path, var, xpeak, oyind),
                "w", rms=rms[xpeak, oyind, :].squeeze())
        save2nc("{}/dc{}.nc".format(db_path, var),
                "w", dc=dc[:, :, :].squeeze())
    if boutf.get_yesno('update fig', timeout=timeout):
        boutf.run_command()
    if savefig:
        boutv.savefig('{}/rms{}{}_vs_t'.format(db_path, var, strcase))

if 'data0_dc' in parts:
    print("========plot data0_dc")
    plt.figure("data0_dc", figsize=(30, 10))
    plt.clf()
    mpl.rcParams['axes.prop_cycle'] = cycler(
        color=boutv.color_list(len(t_arr)))
    for i in t_arr: plt.plot(psi, data0[:, oyind]+dc[:, oyind, i].T)
    plt.legend(t_arr, loc='center left', bbox_to_anchor=(1.02, 0.5), title=r'time/{}$\tau_A$'.format(
        timestep if timestep != 1 else ''), ncol=2, fontsize=16)
    plt.xlabel('$\psi_n$')
    plt.title('$\hat P$, y{}'.format(oyind))
    mpl.rcParams['axes.prop_cycle'] = cycler(color=boutv.colors_default)
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.show()
    if boutf.get_yesno('update fig', timeout=timeout):
        boutf.run_command()
    if savefig:
        boutv.savefig('{}/{}{}_y{}'.format(
            db_path, var, strcase, oyind))

if 'fn' in parts:
    print('========plot Amp. vs n')
    plt.figure('fn', figsize=(30, 10))
    plt.clf()
    mpl.rcParams['axes.prop_cycle'] = cycler(
        color=boutv.color_list(len(t_arr)))
    plt.plot(np.arange(lowpass+1)*zperiod,
             fpa[xpeak, oyind, :lowpass+1, t_arr].T)
    plt.xlabel('toroidal mode number')
    plt.title('${}_n$'.format(var))
    plt.legend(t_arr, loc='center left', bbox_to_anchor=(1.02, 0.5), title=r'time/{}$\tau_A$'.format(
        timestep if timestep != 1 else ''), ncol=2, fontsize=16)
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    if boutf.get_yesno('update fig', timeout=timeout):
        boutf.run_command()
    if savefig:
        boutv.savefig('{}/{}n{}'.format(db_path, var, strcase))

    print('========plot Amp. vs time')
    plt.figure('fn vs t')
    plt.clf()
    mpl.rcParams['axes.prop_cycle'] = cycler(
        color=boutv.color_list(lowpass+1))
    plt.plot(t_array, fpa[xpeak, oyind, :lowpass+1, :].T)
    plt.xlabel(r'time/$\tau_A$')
    plt.yscale('log')
    plt.ylabel('Amp. of ${}_n$'.format(var))
    plt.title('{}x{}y{}'.format(var, xpeak, oyind))
    plt.legend(np.arange(lowpass+1)*zperiod, ncol=3, fontsize=18)
    mpl.rcParams['axes.prop_cycle'] = cycler(color=boutv.colors_default)
    plt.tight_layout()
    if boutf.get_yesno('update fig', timeout=timeout):
        boutf.run_command()
    if savefig:
        boutv.savefig('{}/{}n_vs_time_x{}y{}'.format(
            db_path, var, xpeak, oyind))

if 'elmsize' in parts:
    print('========elmsize')
    plt.figure('elmsize')
    plt.clf()
    elm_size = elmsize(dc, data0, grid, norm=norm)
    np.save('{}/elmsize_1ne.npy'.format(db_path),elm_size['s3']*100.)
    plt.plot(t_array, elm_size['s3']*100.)
    plt.ylabel('elmsize$_{{{}}}$ (%)'.format(var))
    plt.xlabel(r'time/$\tau_A$')
    plt.tight_layout()
    if savedata:
        print("saving elmsize data ...")
        save2nc("{}/elmsize_{}".format(db_path, var), 'w',
                timestep=timestep, zperiod=zperiod, norm=norm, **elm_size)
    if boutf.get_yesno('update fig', timeout=timeout):
        boutf.run_command()
    if savefig:
        boutv.savefig('{}/{}_elmsize'.format(db_path, var))

plt.show()

