# Reads in data from CMIP files, making adjustments.

from netCDF4 import Dataset
import numpy as np
import copy
import calendar
import cftime

# Determine dimensions using a given file and variable:
def getDims(infile,var):
    ncid = Dataset(infile,'r')
    dims = ncid.variables[var].get_dims()
    ncid.close
    
    return dims

# Checks and fixes time if calendars don't match:
def fixTime(units,units2,cal,time,time_bnds):
    if units2 != units:
        yr_init = int(units.split(' ')[2].split('-')[0])
        #yr_init = int(units.split(' ')[2][0:4])
        yr_new  = int(units2.split(' ')[2].split('-')[0])
        #yr_new  = int(units2.split(' ')[2][0:4])
        nleap   = 0 
        if ((cal == 'standard') | (cal == 'gregorian')):
            days = 365
            for yy in range(yr_init,yr_new):
                if calendar.isleap(yy):
                    nleap = nleap + 1
        elif cal == 'noleap':
            days = 365
        else:
            days = int(cal[0:3])
        offset = days*(yr_new-yr_init) + nleap
        time      = time      + offset
        time_bnds = time_bnds + offset
        
    return time, time_bnds

# Determine starting file and position:
def timeFile(Model,cal,units2,tt,Files):
    ff    = 0
    nn    = np.nan
    found = False
    dd = cftime.num2date(tt,units2,cal)
    
    nf = len(Files)
    while(((ff < nf) & (found == False))):
        ncid = Dataset(Files[ff],'r')
        ttf = ncid.variables['time'][:]
        units = ncid.variables['time'].units
        ncid.close()
        if ((Model.name == 'FGOALS-g2')):
            units = (units + '-01')
        
        dd2 = cftime.date2num(dd,units,cal)
        
        if(((dd2 >= ttf[0]) & (dd2 <= ttf[-1]))):
            nn = 0
            while(((nn < len(ttf)) & (found == False))):
                if ((cftime.num2date(ttf[nn],units,cal).year  == cftime.num2date(tt,units2,cal).year) &
                    (cftime.num2date(ttf[nn],units,cal).month == cftime.num2date(tt,units2,cal).month)):
                    found = True
                else:
                    nn = nn + 1
        else:
            ff = ff + 1
            
    return ff, nn

# Reads in lat and lon data from ocean:
def Olatlon(Model,infile,var):
    ncid = Dataset(infile,'r')
    if (Model.name == 'MPI-ESM1-2-XR'):
        if (var == 'uo'):
            if Model.Oreg:
                lon = ncid.variables['lon_2'][:]
                lat = ncid.variables['lat_2'][:]
            else:
                lon = ncid.variables['lon_2'][:,:]
                lat = ncid.variables['lat_2'][:,:]
        elif  (var == 'vo'):
            if Model.Oreg:
                lon = ncid.variables['lon_3'][:]
                lat = ncid.variables['lat_3'][:]
            else:
                lon = ncid.variables['lon_3'][:,:]
                lat = ncid.variables['lat_3'][:,:]
        else:
            if Model.Oreg:
                lon = ncid.variables[Model.Olon][:]
                lat = ncid.variables[Model.Olat][:]
            else:
                lon = ncid.variables[Model.Olon][:,:]
                lat = ncid.variables[Model.Olat][:,:]           
    else:
        if Model.Oreg:
            lon = ncid.variables[Model.Olon][:]
            lat = ncid.variables[Model.Olat][:]
        else:
            lon = ncid.variables[Model.Olon][:,:]
            lat = ncid.variables[Model.Olat][:,:]
    ncid.close()
    
    # Flip North-South:
    if Model.OflipNS:
        lat  = np.flip(lat,axis=0)
        if not Model.Oreg:
            lon  = np.flip(lon,axis=0)
            
        
    # Extra row in T fields (coded only for regular grid):
    if Model.OextraT:
        if ((var == 'vo') | (var == 'uo') | (var == 'tauuo') | (var == 'tauvo') | (var == 'hfy')):
            if Model.Oreg:
                lat = np.concatenate((lat,[-90,]),0)
            else:
                lat = np.concatenate((lat,-90*np.ones((1,np.size(lat,axis=1)),'float')),0)
                lon = np.concatenate((lon,lon[-1:,:]),0)
                
                
    # Remove extra W-E columns:
    if Model.Oreg:
        ni  = np.size(lon,axis=0)
        lon = lon[Model.OextraWE[0]:(ni-Model.OextraWE[1])]
    else:
        ni  = np.size(lon,axis=1)
        lon = lon[:,Model.OextraWE[0]:(ni-Model.OextraWE[1])]
        lat = lat[:,Model.OextraWE[0]:(ni-Model.OextraWE[1])]
        
    return lon,lat

# Reads in lat and lon data from atmosphere:
def Alatlon(Model,infile,var):
    ncid = Dataset(infile,'r')
    if Model.Areg:
        lon = ncid.variables[Model.Alon][:]
        lat = ncid.variables[Model.Alat][:]
    else:
        lon = ncid.variables[Model.Alon][:,:]
        lat = ncid.variables[Model.Alat][:,:]
    ncid.close()
    
    # Flip North-South:
    if Model.AflipNS:
        lat  = np.flip(lat,axis=0)
        if not Model.Areg:
            lon  = np.flip(lon,axis=0)
        
    # Extra row in u and v fields (coded only for regular grid):
    if Model.AextraT:
        print('Need to code for AextraUV')
                
    # Remove extra W-E columns:
    if Model.Areg:
        ni  = np.size(lon,axis=0)
        lon = lon[Model.AextraWE[0]:(ni-Model.AextraWE[1])]
    else:
        ni  = np.size(lon,axis=1)
        lon = lon[:,Model.AextraWE[0]:(ni-Model.AextraWE[1])]
        lat = lat[:,Model.AextraWE[0]:(ni-Model.AextraWE[1])]
        
    return lon,lat

# Reads in data from a 2D ocean field:
def Oread2Ddata(Model,infile,var,time=None,lev=None,mask=False):
    ncid = Dataset(infile,'r')
    # Flip Up-Down:
    if ((lev != None) & ((Model.name == 'CFSv2-2011') | (Model.name == 'FGOALS-gl') | (Model.name == 'HadGEM2-AO'))):
        nk  = len(ncid.variables['lev'][:])
        lev = nk - 1 - lev
    if mask:
        if time == None:
            if lev == None:
                data = 1-np.squeeze(ncid.variables[var][:,:]).mask
            else:
                data = 1-np.squeeze(ncid.variables[var][lev,:,:]).mask
        else:
            if lev == None:
                data = 1-np.squeeze(ncid.variables[var][time,:,:]).mask
            else:
                data = 1-np.squeeze(ncid.variables[var][time,lev,:,:]).mask
    else:
        if time == None:
            if lev == None:
                data = np.squeeze(ncid.variables[var][:,:]).data
            else:
                data = np.squeeze(ncid.variables[var][lev,:,:]).data
        else:
            if lev == None:
                data = np.squeeze(ncid.variables[var][time,:,:]).data
            else:
                data = np.squeeze(ncid.variables[var][time,lev,:,:]).data
    ncid.close()
    
    # Flip North-South:
    if Model.OflipNS:
        data = np.flip(data,axis=0)
        
    # Extra row in u and v fields:
    if Model.OextraT:
        if ((var == 'vo') | (var == 'uo') | (var == 'tauvo') | (var == 'tauuo') | (var == 'hfy')):
            data = np.concatenate((data,np.expand_dims(data[-1,:],0)),0)
                
    # Remove extra W-E columns:
    ni   = np.size(data,axis=1)
    data = data[:,Model.OextraWE[0]:(ni-Model.OextraWE[1])]
        
    return data

# Reads in data from a 3D ocean field:
def Oread3Ddata(Model,infile,var,time=None,mask=False):
    ncid = Dataset(infile,'r')
    if mask:
        if time == None:
            data = 1-np.squeeze(ncid.variables[var][:,:,:]).mask
        else:
            data = 1-np.squeeze(ncid.variables[var][time,:,:,:]).mask
    else:
        if time == None:
            data = np.squeeze(ncid.variables[var][:,:,:]).data
        else:
            data = np.squeeze(ncid.variables[var][time,:,:,:]).data
    ncid.close()
    
    # Flip North-South:
    if Model.OflipNS:
        data = np.flip(data,axis=1)
        
    # Extra row in u and v fields:
    if Model.OextraT:
        if ((var == 'vo') | (var == 'uo') | (var == 'tauvo') | (var == 'tauuo') | (var == 'hfy')):
            data = np.concatenate((data,np.expand_dims(data[:,-1,:],1)),1)
            
    # Remove extra W-E columns:
    ni   = np.size(data,axis=2)
    data = data[:,:,Model.OextraWE[0]:(ni-Model.OextraWE[1])]
    
    # Flip Up-Down:
    if ((Model.name == 'CFSv2-2011') | (Model.name == 'FGOALS-gl') | (Model.name == 'HadGEM2-AO')):
        data = data[::-1,:,:]
        
    return data

# Reads in data from a 2D atmosphere field:
def Aread2Ddata(Model,infile,var,time=None,lev=None,mask=False):
    ncid = Dataset(infile,'r')
    if mask:
        if time == None:
            if lev == None:
                data = 1-np.squeeze(ncid.variables[var][:,:]).mask
            else:
                data = 1-np.squeeze(ncid.variables[var][lev,:,:]).mask
        else:
            if lev == None:
                data = np.squeeze(ncid.variables[var][time,:,:]).mask
            else:
                data = np.squeeze(ncid.variables[var][time,lev,:,:]).mask
    else:
        if time == None:
            if lev == None:
                data = np.squeeze(ncid.variables[var][:,:]).data
            else:
                data = np.squeeze(ncid.variables[var][lev,:,:]).data
        else:
            if lev == None:
                data = np.squeeze(ncid.variables[var][time,:,:]).data
            else:
                data = np.squeeze(ncid.variables[var][time,lev,:,:]).data
    ncid.close()
    
    # Flip North-South:
    if Model.AflipNS:
        data = np.flip(data,axis=0)
                
    # Remove extra W-E columns:
    ni   = np.size(data,axis=1)
    data = data[:,Model.AextraWE[0]:(ni-Model.AextraWE[1])]
        
    return data

# Move data onto different grid points:
def moveData(Model,grid1,grid2,data,computation='mean',dyt=[]):
    if ((grid1 == 'T') & (grid2 == 'U')):
        if Model.Ogrid[0] == 'B':
            if np.size(np.shape(data)) == 2:
                tmp = np.tile(data,(4,1,1))
                ncid  = Dataset(Model.Omeshmask,'r')
                umask = ncid.variables['umask'][0,:,:]
                ncid.close()
                if Model.Ogrid[1] == 'b':
                    tmp[2,:-1,:] = data[1:,:]
                    tmp[3,:-1,:] = data[1:,:]
                elif Model.Ogrid[1] == 't':
                    tmp[2,1:,:]  = data[:-1,:]
                    tmp[3,1:,:]  = data[:-1,:]
                if Model.Ogrid[2] == 'l':
                    tmp[1,:,:] = np.roll(tmp[0,:,:],1,axis=1)
                    tmp[3,:,:] = np.roll(tmp[2,:,:],1,axis=1)
                elif Model.Ogrid[2] == 'r':
                    tmp[1,:,:] = np.roll(tmp[0,:,:],-1,axis=1)
                    tmp[3,:,:] = np.roll(tmp[2,:,:],-1,axis=1)                
            elif np.size(np.shape(data)) == 3:
                tmp = np.tile(data,(4,1,1,1))
                ncid  = Dataset(Model.Omeshmask,'r')
                umask = ncid.variables['umask'][:,:,:]
                ncid.close()
                if Model.Ogrid[1] == 'b':
                    tmp[2,:,:-1,:] = data[:,1:,:]
                    tmp[3,:,:-1,:] = data[:,1:,:]
                elif Model.Ogrid[1] == 't':
                    tmp[2,:,1:,:]  = data[:,:-1,:]
                    tmp[3,:,1:,:]  = data[:,:-1,:]
                if Model.Ogrid[2] == 'l':
                    tmp[1,:,:,:] = np.roll(tmp[0,:,:,:],1,axis=2)
                    tmp[3,:,:,:] = np.roll(tmp[2,:,:,:],1,axis=2)
                elif Model.Ogrid[2] == 'r':
                    tmp[1,:,:,:] = np.roll(tmp[0,:,:,:],-1,axis=2)
                    tmp[3,:,:,:] = np.roll(tmp[2,:,:,:],-1,axis=2)
        elif Model.Ogrid[0] == 'C':
            if np.size(np.shape(data)) == 2:
                tmp = np.tile(data,(2,1,1))
                ncid  = Dataset(Model.Omeshmask,'r')
                umask = ncid.variables['umask'][0,:,:]
                ncid.close()
                if Model.Ogrid[2] == 'l':
                    tmp[1,:,:] = np.roll(tmp[0,:,:],1,axis=1)
                elif Model.Ogrid[2] == 'r':
                    tmp[1,:,:] = np.roll(tmp[0,:,:],-1,axis=1)                
            elif np.size(np.shape(data)) == 3:
                tmp = np.tile(data,(2,1,1,1))
                ncid  = Dataset(Model.Omeshmask,'r')
                umask = ncid.variables['umask'][:,:,:]
                ncid.close()
                if Model.Ogrid[2] == 'l':
                    tmp[1,:,:,:] = np.roll(tmp[0,:,:,:],1,axis=2)
                elif Model.Ogrid[2] == 'r':
                    tmp[1,:,:,:] = np.roll(tmp[0,:,:,:],-1,axis=2)
        elif Model.Ogrid[0] == 'A':
            if np.size(np.shape(data)) == 2:
                tmp = np.tile(data,(1,1,1)) 
                ncid  = Dataset(Model.Omeshmask,'r')
                umask = ncid.variables['umask'][0,:,:]
                ncid.close()              
            elif np.size(np.shape(data)) == 3:
                tmp = np.tile(data,(1,1,1,1))
                ncid  = Dataset(Model.Omeshmask,'r')
                umask = ncid.variables['umask'][:,:,:]
                ncid.close()
                
        if computation == 'mean':
            datanew = np.squeeze(np.mean(tmp,axis=0)*umask)
        elif computation == 'min':
            datanew = np.squeeze(np.min(tmp,axis=0)*umask)
        elif computation == 'max':
            datanew = np.squeeze(np.max(tmp,axis=0)*umask)
            
    elif ((grid1 == 'T') & (grid2 == 'VT')):
        # Data won't be masked:
        ncid  = Dataset(Model.Omeshmask,'r')
        dyt   = ncid.variables['dyt'][:,:]
        ncid.close()
        if ((Model.Ogrid[0] == 'A') | (Model.Ogrid[1] == 'b')):
            if np.size(np.shape(data)) == 2:
                tmp = np.tile(data,(2,1,1))
                dyt = np.tile(dyt,(2,1,1))
                tmp[1,1:,:] = tmp[0,:-1,:]
                dyt[1,1:,:] = dyt[0,:-1,:]
            elif np.size(np.shape(data)) == 3:
                tmp = np.tile(data,(2,1,1,1))
                dyt = np.tile(dyt,(2,np.size(tmp,1),1,1))
                tmp[1,:,1:,:] = tmp[0,:,:-1,:]
                dyt[1,:,1:,:] = dyt[0,:,:-1,:] 
        else:
            if np.size(np.shape(data)) == 2:
                tmp = np.tile(data,(2,1,1))
                dyt = np.tile(dyt,(2,1,1))
                tmp[1,:-1,:] = tmp[0,1:,:]
                dyt[1,:-1,:] = dyt[0,1:,:]
            elif np.size(np.shape(data)) == 3:
                tmp = np.tile(data,(2,1,1,1))
                dyt = np.tile(dyt,(2,np.size(tmp,1),1,1))
                tmp[1,:,:-1,:] = tmp[0,:,1:,:]
                dyt[1,:,:-1,:] = dyt[0,:,1:,:] 
                
            
        if computation == 'mean':
            datanew = np.squeeze(np.sum(tmp*dyt,axis=0)/np.sum(dyt,axis=0))
        elif computation == 'min':
            datanew = np.squeeze(np.min(tmp,axis=0))
        elif computation == 'max':
            datanew = np.squeeze(np.max(tmp,axis=0))
        
    elif ((grid1 == 'V') & (grid2 == 'VT')):
        # Data won't be masked:
        if Model.Ogrid[0] == 'A':
            ncid  = Dataset(Model.Omeshmask,'r')
            dyv   = ncid.variables['dyv'][:,:]
            ncid.close()
            if np.size(np.shape(data)) == 2:
                tmp = np.tile(data,(2,1,1))
                dyv = np.tile(dyv,(2,1,1))
                tmp[1,1:,:] = tmp[0,:-1,:]
                dyv[1,1:,:] = dyv[0,:-1,:]
            elif np.size(np.shape(data)) == 3:
                tmp = np.tile(data,(2,1,1,1))
                dyv = np.tile(dyv,(2,np.size(tmp,1),1,1))
                tmp[1,:,1:,:] = tmp[0,:,:-1,:]
                dyv[1,:,1:,:] = dyv[0,:,:-1,:] 
            
            if computation == 'mean':
                datanew = np.squeeze(np.sum(tmp*dyv,axis=0)/np.sum(dyv,axis=0))
            elif computation == 'min':
                datanew = np.squeeze(np.min(tmp,axis=0))
            elif computation == 'max':
                datanew = np.squeeze(np.max(tmp,axis=0))
                
        elif Model.Ogrid[0] == 'B':
            if Model.Ogrid[2] == 'r':
                if np.size(np.shape(data)) == 2:
                    tmp = np.tile(data,(2,1,1))
                    tmp[1,:,:] = np.roll(tmp[0,:,:],1,axis=1)
                elif np.size(np.shape(data)) == 3:
                    tmp = np.tile(data,(2,1,1,1))
                    tmp[1,:,:,:] = np.roll(tmp[0,:,:,:],1,axis=2)
            if Model.Ogrid[2] == 'l':
                if np.size(np.shape(data)) == 2:
                    tmp = np.tile(data,(2,1,1))
                    tmp[1,:,:] = np.roll(tmp[0,:,:],-1,axis=1)
                elif np.size(np.shape(data)) == 3:
                    tmp = np.tile(data,(2,1,1,1))
                    tmp[1,:,:,:] = np.roll(tmp[0,:,:,:],-1,axis=2)
            
            if computation == 'mean':
                datanew = np.squeeze(np.mean(tmp,axis=0))
            elif computation == 'min':
                datanew = np.squeeze(np.min(tmp,axis=0))
            elif computation == 'max':
                datanew = np.squeeze(np.max(tmp,axis=0))
        elif Model.Ogrid[0] == 'C':
            datanew = tmp
            
    elif ((grid1 == 'T') & (grid2 == 'V')):
        # Data won't be masked:
        if (Model.Ogrid[0] == 'A'):
            datanew = data
        elif (Model.Ogrid[0] == 'B'):
            if np.size(np.shape(data)) == 2:
                tmp = np.tile(data,(4,1,1))
            elif np.size(np.shape(data)) == 3:
                tmp = np.tile(data,(4,1,1,1))
                
            if (Model.Ogrid[1] == 't'):
                if np.size(np.shape(data)) == 2:
                    tmp[2,:-1,:] = tmp[0,1:,:]
                    tmp[3,:-1,:] = tmp[0,1:,:]
                elif np.size(np.shape(data)) == 3:
                    tmp[2,:,:-1,:] = tmp[0,:,1:,:]
                    tmp[3,:,:-1,:] = tmp[0,:,1:,:]
            elif (Model.Ogrid[1] == 'b'):
                if np.size(np.shape(data)) == 2:
                    tmp[2,1:,:] = tmp[0,:-1,:]
                    tmp[3,1:,:] = tmp[0,:-1,:]
                elif np.size(np.shape(data)) == 3:
                    tmp[2,:,1:,:] = tmp[0,:,:-1,:]
                    tmp[3,:,1:,:] = tmp[0,:,:-1,:]
            
            if (Model.Ogrid[2] == 'r'):
                if np.size(np.shape(data)) == 2:
                    tmp[1,:,:] = np.roll(tmp[0,:,:],-1,axis=1)
                    tmp[2,:,:] = np.roll(tmp[3,:,:],-1,axis=1)
                elif np.size(np.shape(data)) == 3:
                    tmp[1,:,:,:] = np.roll(tmp[0,:,:,:],-1,axis=2)
                    tmp[2,:,:,:] = np.roll(tmp[3,:,:,:],-1,axis=2)
            elif (Model.Ogrid[2] == 'l'):
                if np.size(np.shape(data)) == 2:
                    tmp[1,:,:] = np.roll(tmp[0,:,:],1,axis=1)
                    tmp[2,:,:] = np.roll(tmp[3,:,:],1,axis=1)
                elif np.size(np.shape(data)) == 3:
                    tmp[1,:,:,:] = np.roll(tmp[0,:,:,:],1,axis=2)
                    tmp[2,:,:,:] = np.roll(tmp[3,:,:,:],1,axis=2)
            
            if computation == 'mean':
                datanew = np.squeeze(np.mean(tmp,axis=0))
            elif computation == 'min':
                datanew = np.squeeze(np.min(tmp,axis=0))
            elif computation == 'max':
                datanew = np.squeeze(np.max(tmp,axis=0))
        else:
            if (np.size(dyt) == 0):
                ncid  = Dataset(Model.Omeshmask,'r')
                dyt   = ncid.variables['dyt'][:,:]
                ncid.close()
            if ((Model.Ogrid[1] == 'b')):
                if np.size(np.shape(data)) == 2:
                    tmp = np.tile(data,(2,1,1))
                    dyt = np.tile(dyt,(2,1,1))
                    tmp[1,1:,:] = tmp[0,:-1,:]
                    dyt[1,1:,:] = dyt[0,:-1,:]
                elif np.size(np.shape(data)) == 3:
                    tmp = np.tile(data,(2,1,1,1))
                    dyt = np.tile(dyt,(2,np.size(tmp,1),1,1))
                    tmp[1,:,1:,:] = tmp[0,:,:-1,:]
                    dyt[1,:,1:,:] = dyt[0,:,:-1,:] 
            else:
                if np.size(np.shape(data)) == 2:
                    tmp = np.tile(data,(2,1,1))
                    dyt = np.tile(dyt,(2,1,1))
                    tmp[1,:-1,:] = tmp[0,1:,:]
                    dyt[1,:-1,:] = dyt[0,1:,:]
                elif np.size(np.shape(data)) == 3:
                    tmp = np.tile(data,(2,1,1,1))
                    dyt = np.tile(dyt,(2,np.size(tmp,1),1,1))
                    tmp[1,:,:-1,:] = tmp[0,:,1:,:]
                    dyt[1,:,:-1,:] = dyt[0,:,1:,:] 
                
            
            if computation == 'mean':
                datanew = np.squeeze(np.sum(tmp*dyt,axis=0)/np.sum(dyt,axis=0))
            elif computation == 'min':
                datanew = np.squeeze(np.min(tmp,axis=0))
            elif computation == 'max':
                datanew = np.squeeze(np.max(tmp,axis=0))
            
    else:
        print('Need to code')
    
    return datanew
