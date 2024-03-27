import numpy as np
from atlas_plot_funcs import *
import glob


print('Starting')
formatPlots(usetex=False, show=False, figSize=(3,3))

IDs = [
        "SJD-500x600-Uni4E16",
        "SJD-500x600-Ferro-Uni4E16",
        "SJD-500x600-Ferro-Uni4E16-pshift",
        "SJD-500x600-Uni4E16-eps100",
]
for ID in IDs:
    # ID = "SJD-500x600-Ferro-Uni2E16"
    dataDir = "Z:/Atlas/Vertical_SJ_FET/" + ID + "/TonyPlot/"
    extractDir = "Z:/Atlas/Vertical_SJ_FET/" + ID + "/extract/"
    saveDir = "Z:/Atlas/Vertical_SJ_FET/" + ID + '/fullPlots/'


    xrange = [0,(0.5+0.6)/2]
    yrange = [0,10.25]
    copies = 3

    for file in glob.glob(dataDir + '*.png'):
        # if 'forward' not in file:
        #     continue
        try:
            if 'EFieldLin' in file:
                print(file)
                plot_2d(file, xrange, yrange, clabel=r'Electric Field [MV/cm]', crange=[0, 4], saveDir=saveDir, copies=copies)
            if 'EFieldLin_zoom' in file:
                print(file)
                plot_2d(file, xrange, [0,1.25], clabel=r'Electric Field [MV/cm]', crange=[0, 4], saveDir=saveDir, copies=1, fig_size=(3,3))
            if 'Permittivity' in file:
                print(file)
                plot_2d(file, xrange, yrange, clabel=r'Relative Permittivity', crange=[0, 150], saveDir=saveDir, copies=copies)
            if 'EFieldLog' in file:
                print(file)
                plot_2d(file, xrange, yrange, clabel=r'Electric Field [V/cm]', crange=[1E0, 1E7], c_scale='log', saveDir=saveDir, copies=copies)
            if 'eConc' in file:
                print(file)
                plot_2d(file, xrange, yrange, clabel=r'Electron Conc. [$\mathrm{cm^{-3}}$]', crange=[1E4, 1E20], c_scale='log', saveDir=saveDir, copies=copies)
            if 'curDensity_lin' in file:
                print(file)
                plot_2d(file, xrange, yrange, clabel=r'Current Density [$\mathrm{A/cm^{2}}$]', crange=[0, 2E4], c_scale='linear', saveDir=saveDir, copies=copies)
            if 'curDensity_log' in file:
                print(file)
                plot_2d(file, xrange, yrange, clabel=r'Current Density [$\mathrm{A/cm^{2}}$]', crange=[1E0, 3.16E4], c_scale='log', saveDir=saveDir, copies=copies)
            if 'curDen_lin_z' in file:
                print(file)
                plot_2d(file, xrange, [8.6,10.0], clabel=r'Current Density [$\mathrm{A/cm^{2}}$]', crange=[1E0, 2E4], c_scale='linear', saveDir=saveDir, copies=1)
            if 'carrierConc' in file:
                print(file)
                plot_2d(file, xrange, yrange, clabel=r'Carrier Conc. [$\mathrm{cm^{-3}}$]', crange=[1E-1, 1E19], c_scale='log', saveDir=saveDir, copies=copies)
            if 'Potential' in file:
                print(file)
                nameparts = file.split('_')
                VD = float(nameparts[-2])
                plot_2d(file, xrange, yrange, clabel=r'Potential [$\mathrm{V}$]', crange=[-1, 1000], c_scale='lin', saveDir=saveDir, copies=copies)
        except Exception as e:
            print(e)
            pass
        plt.close('all')

    for file in glob.glob(extractDir + '*.dat'):
        print(file)
        try:
            if 'EFieldCenter' in file:
                print(file)
                name = file.split('/')[-1].split('\\')[-1].strip('.dat')
                print(name)
                plot_cutline(file, name, x_label='y-depth [$\mathrm{\mu m}$]', y_label='Electric Field [$\mathrm{MV/cm}$]',  saveDir=saveDir, xlim=None, yfactor=1E-6)
            if 'EFieldEdge' in file:
                print(file)
                name = file.split('/')[-1].split('\\')[-1].strip('.dat')
                print(name)
                plot_cutline(file, name, x_label='y-depth [$\mathrm{\mu m}$]', y_label='Electric Field [$\mathrm{MV/cm}$]',  saveDir=saveDir, yfactor=1E-6)
            if 'PotentialCenter' in file:
                print(file)
                name = file.split('/')[-1].split('\\')[-1].strip('.dat')
                print(name)
                plot_cutline(file, name, x_label='y-depth [$\mathrm{\mu m}$]', y_label='Potential [$\mathrm{V}$]',  saveDir=saveDir)
            if 'forward_IV' in file:
                print(file)
                name = file.split('/')[-1].split('\\')[-1].strip('.dat')
                print(name)
                plot_IV(file, name, xrange, x_label='Voltage [V]', y_label='Current Density [$\mathrm{kA/cm^{-2}}$]',  saveDir=saveDir)
            if 'forward_IV' in file:
                print(file)
                name = file.split('/')[-1].split('\\')[-1].strip('.dat') + '_log'
                print(name)
                plot_IV(file, name, xrange, x_label='Voltage [V]', y_label='Current Density [$\mathrm{kA/cm^{-2}}$]',  saveDir=saveDir, y_scale='log', xlim=[0,3])
            if 'reverse_IV' in file:
                print(file)
                name = file.split('/')[-1].split('\\')[-1].strip('.dat') + '_log'
                print(name)
                plot_IV(file, name, xrange, x_label='Voltage [V]', y_label='Current Density [$\mathrm{kA/cm^{-2}}$]',  saveDir=saveDir, y_scale='log', xlim=[0,3])
            if 'doping.dat' in file:
                print(file)
                name = file.split('/')[-1].split('\\')[-1].strip('.dat')
                print(name)
                plot_cutline(file, name, x_label='y-depth [$\mathrm{\mu m}$]', y_label='Doping Concentration [$\mathrm{cm^{-3}}$]', y_scale='log', ylim=[1E15, 1E20], saveDir=saveDir)
        except Exception as e:
                pass
        plt.close('all')






























# %%
