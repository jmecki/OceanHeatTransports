##############################################################################
#
# Computes means over specified time periods:
#
# #############################################################################

import os
import numpy as np
from netCDF4 import Dataset
import cftime
import copy
import sys

from . import CMIPobject
from . import CMIPread

# Mean of seasonal cycle over several years:
def meanSC(Model,EXP,ENS,var,vtype,styr,fnyr,outfile,gtype='gn'):
    print('computing ' + Model.name + ' ' + EXP + ' ' + ENS)
    modeldir = os.path.dirname(outfile)
    
    # Make directory if it doesn't exist:
    if not os.path.exists(modeldir):
        os.makedirs(modeldir)
    
    # Find all files needed for the computations:
    Files = Model.getFiles(var,EXP=EXP,ENS=ENS,vtype=vtype,gtype=gtype)
    nf = len(Files)
    print(nf)
    print(gtype)
    
    # Get model information:
    dims  = CMIPobject.getDims(Files[0],var)
    nd    = len(dims)   
    nyr   = fnyr - styr + 1
    dtime = dims[0].name
    dlon  = dims[nd-1].name
    nx    = dims[nd-1].size
    dlat  = dims[nd-2].name
    ny    = dims[nd-2].size
    if nd == 4:
        dlev = dims[1].name
        nz   = dims[1].size
    if Model.OextraT:  # Add more as required:
        if ((var == 'vo') | (var == 'uo') | (var == 'tauuo') | (var == 'tauvo')):
            ny = ny + 1
    # Decrease dimension size if data is W-E periodic:
    nx = nx - np.sum(Model.OextraWE)
            
                    
    # Get latitude and longitude:
    extra = False
    if vtype == 'Omon':
        nlon = Model.Olon
        nlat = Model.Olat
        reg  = Model.Oreg
        lon,lat = CMIPread.Olatlon(Model,Files[0],var)
    elif vtype == 'Amon':
        nlon = Model.Alon
        nlat = Model.Alat
        reg  = Model.Areg
        lon,lat = CMIPread.Alatlon(Model,Files[0],var)
    else:
        print((vtype + ' is an invalid/uncoded type'))
        
    # Loop through all files and add together the relavent years:
    if nd == 3:
        data = np.zeros((12,ny,nx),'float')
    elif nd == 4:
        data = np.zeros((12,nz,ny,nx),'float')  
    days = np.zeros((12),'float')
    bnds = np.zeros((12,2),'float')
    nn = 0

    # Loop through all files:
    for ff in range(0,nf):
        # Get time information from file:
        print(Files[ff])
        ncid = Dataset(Files[ff],'r')
        time   = ncid.variables[dtime][:]
        bounds = ncid.variables[dtime].bounds
        time_bnds = ncid.variables[bounds][:,:]
        cal   = ncid.variables[dtime].calendar
        units = ncid.variables[dtime].units
        ncid.close()
        if Model.name == 'FGOALS-g2':
            units = (units + '-01')  
            
        nt = np.size(time)
            
        for tt in range(0,nt):
            yr = cftime.num2date(time[tt],units,cal).year
            mm = cftime.num2date(time[tt],units,cal).month
            
            if ((yr >= styr) & (yr <= fnyr)):
                nn = nn + 1
                if vtype == 'Omon':
                    if nd == 3:
                        tmp = CMIPread.Oread2Ddata(Model,Files[ff],var,time=tt)
                    elif nd == 4:
                        tmp = CMIPread.Oread3Ddata(Model,Files[ff],var,time=tt) 
                elif vtype == 'Amon':
                    if nd == 3:
                        tmp = CMIPread.Aread2Ddata(Model,Files[ff],var,time=tt)
                    elif nd == 4:
                        print('Need to code 3D atmosphere reading in of data') 
                else:
                    print((vtype + ' is not coded yet!'))
                
                data[mm-1,:,:] = data[mm-1,:,:] + tmp
                days[mm-1]     = days[mm-1]     + time[tt]
                bnds[mm-1,:]   = bnds[mm-1,:]   + time_bnds[tt,:]

    data = data/nyr
    days = days/nyr
    bnds = bnds/nyr
    
    if np.sum(days) != 0:
        # Save data to file:
        ncid = Dataset(outfile, 'w', format='NETCDF4')
        # coordinates:
        ncid.createDimension(dlon,nx)
        ncid.createDimension(dlat,ny)
        if nd == 4:
            ncid.createDimension(dlev,nz)
        ncid.createDimension(dtime,None)
        ncid.createDimension('bnds',2)
        # variables:
        if reg:
            ncid.createVariable(nlon,'f8',(dlon,))
            ncid.createVariable(nlat,'f8',(dlat,))
        else:
            ncid.createVariable(nlon,'f8',(dlat,dlon,))
            ncid.createVariable(nlat,'f8',(dlat,dlon,))
        ncid.createVariable(dtime,'f8',(dtime,))
        ncid.createVariable(bounds,'f8',(dtime,'bnds',))
        if nd == 4:
            ncid.createVariable(dlev,'f8',(dlev,))
            ncid.createVariable(var,'f8',(dtime,dlev,dlat,dlon,))
        else:
            ncid.createVariable(var,'f8',(dtime,dlat,dlon,))


        ncid.variables[dtime].calendar = cal
        ncid.variables[dtime].units    = units
        ncid.variables[dtime].bounds   = bounds

        # fill variables:
        if reg:
            ncid.variables[nlon][:] = lon
            ncid.variables[nlat][:] = lat
        else:
            ncid.variables[nlon][:,:] = lon
            ncid.variables[nlat][:,:] = lat
        ncid.variables[dtime][0:12]   = days
        ncid.variables[bounds][:,:]   = bnds
        if nd == 3:
            ncid.variables[var][0:12,:,:] = data
        elif nd == 4:
            ncid.variables[var][0:12,:,:,:] = data

        # close:
        ncid.close()
    else:
        print('No data for the specified years')
        sys.exit()
    
# Compute seasonal mean:
def seasonal_means(Model,EXP,ENS,var,vtype,mons,outfile,gtype='gn'):
    print('computing ' + Model.name + ' ' + EXP + ' ' + ENS)
    # Number of months for mean:
    nm = len(mons)

    # List all files:
    files = Model.getFiles(var,EXP=EXP,ENS=ENS,vtype=vtype,gtype=gtype)
    ns = len(files)
    files.sort()
    print(files)

    # Read information from first file:
    ncid  = Dataset(files[0],'r')
    cal   = ncid.variables['time'].calendar
    units = ncid.variables['time'].units
    ncid.close()
    # Check if model has a regular grid:
    if vtype == 'Omon':
        lon,lat = CMIPread.Olatlon(Model,files[0],var)
        if Model.Oreg:
            lon,lat = np.meshgrid(lon,lat)
    elif vtype == 'Amon':
        lon,lat = CMIPread.Alatlon(Model,files[0],var)
        if Model.Areg:
            lon,lat = np.meshgrid(lon,lat)
        
    else:
        print('need to code')
        sys.exit()

    ni = np.size(lon,axis=1)
    nj = np.size(lon,axis=0)

    # Specify and initiate output file:
    modeldir = os.path.dirname(outfile)
    # Make directory if it doesn't exist:
    if not os.path.exists(modeldir):
        os.makedirs(modeldir)
        
    if not os.path.isfile(outfile):
        ncid = Dataset(outfile,'w')
        # Dimensions:
        ncid.createDimension('year',None)
        ncid.createDimension('x',ni)
        ncid.createDimension('y',nj)

        # Variables:
        ncid.createVariable('year','f8',('year',))
        ncid.createVariable('lon', 'f8',('y','x',))
        ncid.createVariable('lat', 'f8',('y','x',))
        ncid.createVariable(var, 'f8',('year','y','x',))

        # Data:
        ncid.variables['lon'][:,:] = lon
        ncid.variables['lat'][:,:] = lat

        ncid.close()
        nyr = 0
    else:
        # Find out how much has been computed:
        ncid = Dataset(outfile,'r')
        nyr  = np.size(ncid.variables['year'][:],axis=0)
        ncid.close()

    # Loop through each file:
    yy = 0
    for ss in range(0,ns):
        infile = files[ss]
        print(infile)
        ncid      = Dataset(infile,'r')
        time      = ncid.variables['time'][:]
        time_bnds = ncid.variables[ncid.variables['time'].bounds][:,:]
        units2    = ncid.variables['time'].units
        ncid.close()

        # Fix dates if they are different:
        if units2 != units:
            print('need to fix dates! - probably not but check!')

        nt = np.size(time,axis=0)
        # Initialize if first file:
        if ss == 0:
            mm   = 0
            days = 0
            tmp  = np.zeros((nj,ni),'float')

        for tt in range(0,nt):
            # Check if one of the months going into the average:
            if (cftime.num2date(time[tt],units2,cal).month == int(mons[mm])):
                if yy >= nyr:
                    days = days + time_bnds[tt,1] - time_bnds[tt,0]
                    if vtype == 'Omon':
                        tmp = tmp + CMIPread.Oread2Ddata(Model,infile,var,tt)*(time_bnds[tt,1] - time_bnds[tt,0])
                    elif vtype == 'Amon':
                        tmp = tmp + CMIPread.Aread2Ddata(Model,infile,var,tt)*(time_bnds[tt,1] - time_bnds[tt,0])
                    else:
                        print('need to code')
                mm = mm + 1
                # Save if it is last month and if so save:
                if mm == nm:
                    year  = cftime.num2date(time[tt],units2,cal).year
                    print(year)
                    if yy >= nyr:
                        ncido = Dataset(outfile, 'a', format='NETCDF4')
                        ncido.variables['year'][yy]  = year
                        ncido.variables[var][yy,:,:] = tmp/days
                        ncido.close()

                    mm   = 0
                    days = 0
                    tmp  = np.zeros((nj,ni),'float')
                    yy = yy + 1 # Count years

    print('-= DONE =-')
    
# Compute means over a specified lat-lon box (currently only works for data on T-points):
def box_means(Model,EXP,ENS,var,vtype,imin,imax,jmin,jmax,outfile,gmask='',gtype='gn'):
    print('computing ' + Model.name + ' ' + EXP + ' ' + ENS)
    # List all files:
    Tfiles = Model.getFiles(EXP=EXP,ENS=ENS,var=var,vtype=vtype)
    nf = len(Tfiles)
    
    # Read information from first file:
    ncid  = Dataset(Tfiles[0],'r')
    cal   = ncid.variables['time'].calendar
    units = ncid.variables['time'].units
    ncid.close()
    # Check if model has a regular grid:
    if vtype == 'Omon':
        lon,lat = CMIPread.Olatlon(Model,Tfiles[0],var)
        if Model.Oreg:
            lon,lat = np.meshgrid(lon,lat)
    elif vtype == 'Amon':
        lon,lat = CMIPread.Alatlon(Model,Tfiles[0],var)
        if Model.Areg:
            lon,lat = np.meshgrid(lon,lat)
    else:
        print('need to code')
        sys.exit()    
    
    dims  = CMIPobject.getDims(Tfiles[0],var)
    nd    = len(dims)   
    dtime = dims[0].name
    dlon  = dims[nd-1].name
    ni    = dims[nd-1].size
    dlat  = dims[nd-2].name
    nj    = dims[nd-2].size
    if nd == 4:
        dlev = dims[1].name
        nk   = dims[1].size
    elif nd != 3:
        print('Unclear resoulation')
        
    # Set-up masks:
    if vtype == 'Omon':
        # Setup grid:
        meshmask = Model.Omeshmask

        # Set-up mask and known weights for computations:
        ncid   = Dataset(meshmask,'r')
        dxt    = np.squeeze(ncid.variables['dxt'][:,:])
        dyt    = np.squeeze(ncid.variables['dyt'][:,:])
        if nd == 3:
            tmask  = np.squeeze(ncid.variables['tmask'][:,0,:,:])
            amask = copy.deepcopy(tmask)
        elif nd == 4:
            tmask  = np.squeeze(ncid.variables['tmask'][:,:,:,:])
            amask  = copy.deepcopy(tmask[0,:,:])
            dzt    = np.squeeze(ncid.variables['dzt'][:,:,:])
            lev    = np.squeeze(ncid.variables[dlev][:])
        lon    = ncid.variables['tlon'][:,:]
        lat    = ncid.variables['tlat'][:,:]
        ncid.close()
    else:
        print('need to code for atmosphere')
        
    # Mask out the region of interest:
    # Mask data:
    amask[np.where(lat > jmax)] = 0
    amask[np.where(lat < jmin)] = 0
    # Rearrange to satisfy range of input data and fix input data if needed:

    # Make longitude between -180 and 180 (needs to be updated for areas that cross 180):
    lon[np.where(lon >  180)] = lon[np.where(lon >  180)] - 360
    lon[np.where(lon < -180)] = lon[np.where(lon < -180)] + 360
    amask[np.where(lon > imax)] = 0
    amask[np.where(lon < imin)] = 0
    amask = amask*dxt*dyt
        
    # Combine masks:
    if gmask != []:
        amask = amask*gmask

    wmask = copy.deepcopy(tmask)
    if nd == 3:
        wmask = wmask*amask
    elif nd == 4:
        vol   = copy.deepcopy(dzt)
        for kk in range(0,nk):
            wmask[kk,:,:] = wmask[kk,:,:]*amask
            vol[kk,:,:]   = vol[kk,:,:]*amask
        
    
    areaxy = np.sum(np.sum(wmask,axis=-1),axis=-1)
    
    # Initialize output file:
    modeldir = os.path.dirname(outfile)
    # Make directory if it doesn't exist:
    if not os.path.exists(modeldir):
        os.makedirs(modeldir)
        
    
    if os.path.isfile(outfile):    
        # Determine how much has already been computed
        ncid = Dataset(outfile,'r')
        nn = ncid.variables[dtime].size
        ncid.close()
    else:
        # Create file (need to edit for 3D means)
        nn = 0 # Number of months computed is zero

        ncid = Dataset(outfile, 'w', format='NETCDF4')
        # coordinates:
        ncid.createDimension(dtime,None)
        ncid.createDimension('bnds',2)
        if nd == 4:
            ncid.createDimension(dlev,nk)
        # variables:
        ncid.createVariable(dtime,'f8',(dtime,))
        ncid.createVariable('time_bnds','f8',(dtime,'bnds',))
        if nd == 4:
            ncid.createVariable(dlev,'f8',(dlev,))
            ncid.createVariable('vol','f8',(dlev,))
            ncid.createVariable(var,'f8',(dtime,dlev,))
        else:
            ncid.createVariable(var,'f8',(dtime,))

        # Add Calendar info:
        ncid.variables[dtime].calendar = cal
        ncid.variables[dtime].units    = units
        ncid.variables[dtime].bounds   = 'time_bnds'

        # close:
        ncid.close()
        
    # Loop through each uncomputed month:
    nm = nn
    # Loop through each file temperature file:
    for ff in range(0,nf):
        print(('computing file ' + str(ff+1) + ' of ' + str(nf)))
        # Determine how many months of data are in the current temperature file:
        dims  = CMIPobject.getDims(Tfiles[ff],var)
        ncid.close
        nt = dims[0].size

        # Determine which months have been computed:
        if nm >= nt:
            # All months from this file have been computed, move on to next file:
            nm = nm - nt
        else:
            # Loop through and compute all the missing months:
            for mm in range(nm,nt):  # mm is month in the file
                print(('month ' + str(mm+1) + ' of ' + str(nt)))

                if nd == 3:
                    tos = CMIPread.Oread2Ddata(Model,Tfiles[ff],var,time=mm)*wmask
                elif nd == 4:
                    tos = CMIPread.Oread3Ddata(Model,Tfiles[ff],var,time=mm)*wmask
                ncid = Dataset(Tfiles[ff],'r')
                tt   = ncid.variables[dtime][mm]
                # Get Calendar information:
                units2 = ncid.variables[dtime].units
                bounds = ncid.variables[dtime].bounds
                tb   = ncid.variables[bounds][mm,:]
                ncid.close()

                # Fix time:
                time,bounds = CMIPobject.fixTime(units,units2,cal,tt,tb)

                # Main computation:
                sst = np.sum(np.sum(tos,axis=-1),axis=-1)/areaxy

                ncido = Dataset(outfile, 'a', format='NETCDF4')
                ncido.variables[dtime][nn]         = time
                ncido.variables['time_bnds'][nn,:] = bounds
                if nd == 3:
                    ncido.variables[var][nn]       = sst
                elif nd == 4:
                    ncido.variables[var][nn,:]     = sst
                ncido.close()

                # Tidy:
                nn = nn + 1
            nm = 0

    print('-= DONE =-')
