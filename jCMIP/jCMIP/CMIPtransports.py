# Computes ocean based transports and streamfunctions:

import numpy as np
#import glob
#import sys
import os
from netCDF4 import Dataset

from . import CMIPobject
from . import CMIPread

# Compute barotropic streamfunction:
def computePSI(Model,EXP,ENS,outfile):
    print('Computing Barotropic Streamfunction for ' + Model.name + ', ' + EXP + ', ' + ENS)
    
    Ufiles = Model.getFiles('uo',EXP,ENS,vtype='Omon',gtype='gn')
    nf = len(Ufiles)
    Tfiles = Model.getFiles('thkcello',EXP,ENS,vtype='Omon',gtype='gn')
    Vfiles = Model.getFiles('volcello',EXP,ENS,vtype='Omon',gtype='gn')
    
    # Read in mesh information:
    meshmask = Model.Omeshmask
    ncid  = Dataset(meshmask,'r')
    dxu   = ncid.variables['dxu'][:,:]
    dyu   = ncid.variables['dyu'][:,:]
    dzu   = ncid.variables['dzu'][:,:,:]
    umask = np.squeeze(ncid.variables['umask'][:,:,:,:])
    if ((len(Tfiles) == 0) & (len(Vfiles) != 0)):
        dxt   = ncid.variables['dxt'][:,:]
        dyt   = ncid.variables['dyt'][:,:]     
    ncid.close()
    
    # Get file information:
    dims  = CMIPread.getDims(Ufiles[0],'uo')
    dtime = dims[0].name
    dlev  = dims[1].name
    dlat  = dims[2].name
    dlon  = dims[3].name
    
    nk  = dims[1].size
    nj  = dims[2].size
    ni  = dims[3].size - np.sum(Model.OextraWE)
    
    lon,lat = CMIPread.Olatlon(Model,Ufiles[0],'uo')
    
    # Get calendar information:
    ncid   = Dataset(Ufiles[0],'r')
    cal    = ncid.variables[dtime].calendar
    units  = ncid.variables[dtime].units
    bounds = ncid.variables[dtime].bounds
    ncid.close()
    
    # Initialize output file:
    if os.path.isfile(outfile):
        # Determine how much has already been computed
        ncid = Dataset(outfile,'r')
        nn = ncid.variables[dtime].size
        ncid.close()
    else:
        # Create file
        nn = 0 # Number of months computed is zero

        ncid = Dataset(outfile, 'w', format='NETCDF4')
        # coordinates:
        ncid.createDimension(dlon,ni)
        ncid.createDimension(dlat,nj)
        ncid.createDimension('bnds',2)
        ncid.createDimension(dtime,None)
        # variables:
        ncid.createVariable(Model.Olon,'f8',(dlat,dlon,))
        ncid.createVariable(Model.Olat,'f8',(dlat,dlon,))
        ncid.createVariable(dtime,'f8',(dtime,))
        ncid.createVariable(bounds,'f8',(dtime,'bnds',))
        ncid.createVariable('psi','f8',(dtime,dlat,dlon,))

        # fill variables:
        ncid.variables[Model.Olon][:,:] = lon
        ncid.variables[Model.Olat][:,:] = lat

        # Add Calendar info:
        ncid.variables[dtime].calendar = cal
        ncid.variables[dtime].units    = units
        ncid.variables[dtime].bounds   = bounds

        # close:
        ncid.close()
    
    
    # Loop through each uncomputed month:
    nm = nn
    # Loop through each file zonal velocity file:
    for ff in range(0,nf):
        print(('computing file ' + str(ff+1) + ' of ' + str(nf)))
        # Determine how many months of data are in the current file:
        dims = CMIPread.getDims(Ufiles[ff],'uo')
        nt = dims[0].size
        ncid    = Dataset(Ufiles[ff],'r')
        units2  = ncid.variables[dtime].units
        ttime   = ncid.variables[dtime][:]
        tbounds = ncid.variables[bounds][:,:]
        ncid.close()
    
        # Determine which months have been computed:
        if nm >= nt:
            # All months from this file have been computed, move on to next file:
            nm = nm - nt
        else:
            # Loop through and compute all the missing months:
            for mm in range(nm,nt):  # mm is month in the file
                print(('month ' + str(mm+1) + ' of ' + str(nt)))
                
                # Velocity:
                uo = CMIPread.Oread3Ddata(Model,Ufiles[ff],'uo',time=mm)*umask
                
                # Take varying depth into account:
                if len(Tfiles) != 0:
                    # Find matching month:
                    fff,nnn = CMIPread.timeFile(nn,Tfiles)
                    dzt     = CMIPread.Oread3Ddata(Model,Tfiles[fff],'thkcello',time=nnn)
                    # Move data to u point:
                    dzu     = CMIPread.moveData(Model,'T','U',dzt,'min')*umask
                elif len(Vfiles) != 0:
                    # Find matching month:
                    fff,nnn = CMIPread.timeFile(nn,Vfiles)
                    dzt     = CMIPread.Oread3Ddata(Model,Vfiles[fff],'volcello',time=nnn)/np.tile(dxt*dyt,(nk,1,1))
                    # Move data to u point:
                    dzu     = CMIPread.moveData(Model,'T','U',dzt,'min')*umask               
                
                # Time:
                tt = ttime[mm]
                tb = tbounds[mm,:]
                # Fix time:
                tt,tb = CMIPread.fixTime(units,units2,cal,tb,tb)
            
                # Main computation:
                PSI    = -np.cumsum(np.sum(uo*dzu,axis=0)*dyu,axis=0)
                PSI    = (PSI - np.tile(PSI[-1,:],(nj,1)))*umask[0,:,:]

                ncido = Dataset(outfile, 'a', format='NETCDF4')
                ncido.variables[dtime][nn]    = tt
                ncido.variables[bounds][nn,:] = tb
                ncido.variables['psi'][nn,:,:] = PSI/1e6
                ncido.close()
                nn = nn + 1
            nm = 0
                
    print('-= DONE =-')