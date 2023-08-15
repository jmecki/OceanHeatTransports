import numpy as np
import jCMIP as jc
from netCDF4 import Dataset
import sys
import os
import shutil
import copy
import cftime

# Get information:
inputs = sys.argv

Clist = jc.readList('../../CMIP' + inputs[1] + 'list')
Model   = Clist[inputs[2]]
EXP     = inputs[3]
ENS     = inputs[4]
grid    = inputs[5]
gtype   = inputs[6]
ctype   = inputs[7]
outfile = inputs[8]
runfile = inputs[9]

# Parameters:
cp  = 4000 # Specific heat capacity J/kg/K
rho = 1026 # Reference density  kg/m^3
datadir = '/gws/pw/j07/wishbone/jmecki/'

if not os.path.isfile(outfile):
    # Find files for computation:
    if ctype == 'None':
        Tfiles = Model.getFiles('thetao'  ,EXP=EXP,ENS=ENS,vtype='Omon',gtype=gtype)
        Sfiles = Model.getFiles('so'      ,EXP=EXP,ENS=ENS,vtype='Omon',gtype=gtype)
        Vfiles = Model.getFiles('vo'      ,EXP=EXP,ENS=ENS,vtype='Omon',gtype=gtype)
        Dfiles = Model.getFiles('thkcello',EXP=EXP,ENS=ENS,vtype='Omon',gtype=gtype)
        nf = len(Tfiles)
        
        if ((Model.name == 'CSIRO-Mk3-6-0') & (ENS == 'r1i1p1')):
            # There is a seasonal cycle thkcello but only ensemble member to have it
            Dfiles = []

        if len(Dfiles) == 0:
            VOLfiles = Model.getFiles('volcello',EXP=EXP,ENS=ENS,vtype='Omon',gtype=gtype)
            
            if ((Model.name == 'ACCESS-ESM1-5') & (EXP == 'piControl') & (ENS == 'r1i1p1f1')):
                # The years don't line up...
                VOLfiles = []
    elif ctype == 'vbar':
        Tfiles = Model.getFiles('thetao'  ,EXP=EXP,ENS=ENS,vtype='Omon',gtype=gtype)
        Sfiles = Model.getFiles('so'      ,EXP=EXP,ENS=ENS,vtype='Omon',gtype=gtype)
        Dfiles = Model.getFiles('thkcello',EXP=EXP,ENS=ENS,vtype='Omon',gtype=gtype)
        nf = len(Tfiles)
        
        if ((Model.name == 'CSIRO-Mk3-6-0') & (ENS == 'r1i1p1')):
            # There is a seasonal cycle thkcello but only ensemble member to have it
            Dfiles = []

        if len(Dfiles) == 0:
            VOLfiles = Model.getFiles('volcello',EXP=EXP,ENS=ENS,vtype='Omon',gtype=gtype)
            
            if ((Model.name == 'ACCESS-ESM1-5') & (EXP == 'piControl') & (ENS == 'r1i1p1f1')):
                # The years don't line up...
                VOLfiles = []
        if (EXP == 'piControl'):
            Vfile  = (datadir + 'CMIP' + Model.cmip + '/' + Model.name + '/piControl/' + ENS +
                      '/vo_' + Model.name + '_piControl_' + ENS + '_' + gtype + '_first500_SC.nc')
        else:
            Vfile  = (datadir + 'CMIP' + Model.cmip + '/' + Model.name + '/historical/' + ENS +
                      '/vo_' + Model.name + '_historical_' + ENS + '_' + gtype + '_197001-199912_SC.nc')
        ncid   = Dataset(Vfile,'r')
        vo_mon = np.squeeze(ncid.variables['vo'][:,:,:,:])
        ncid.close()
    elif ctype == 'Tbar':
        Tfiles = Model.getFiles('thetao'  ,EXP=EXP,ENS=ENS,vtype='Omon',gtype=gtype)
        Vfiles = Model.getFiles('vo'      ,EXP=EXP,ENS=ENS,vtype='Omon',gtype=gtype)
        Dfiles = Model.getFiles('thkcello',EXP=EXP,ENS=ENS,vtype='Omon',gtype=gtype)
        nf = len(Vfiles)
        
        if ((Model.name == 'CSIRO-Mk3-6-0') & (ENS == 'r1i1p1')):
            # There is a seasonal cycle thkcello but only ensemble member to have it
            Dfiles = []

        if len(Dfiles) == 0:
            VOLfiles = Model.getFiles('volcello',EXP=EXP,ENS=ENS,vtype='Omon',gtype=gtype)
            
            if ((Model.name == 'ACCESS-ESM1-5') & (EXP == 'piControl') & (ENS == 'r1i1p1f1')):
                # The years don't line up...
                VOLfiles = []
        if (EXP == 'piControl'):
            Tfile  = (datadir + 'CMIP' + Model.cmip + '/' + Model.name + '/piControl/' + ENS +
                      '/thetao_' + Model.name + '_piControl_' + ENS + '_' + gtype + '_first500_SC.nc')
            Sfile  = (datadir + 'CMIP' + Model.cmip + '/' + Model.name + '/piControl/' + ENS +
                      '/so_' + Model.name + '_piControl_' + ENS + '_' + gtype + '_first500_SC.nc')
        else:
            Tfile  = (datadir + 'CMIP' + Model.cmip + '/' + Model.name + '/historical/' + ENS +
                      '/thetao_' + Model.name + '_historical_' + ENS + '_' + gtype + '_197001-199912_SC.nc')
            Sfile  = (datadir + 'CMIP' + Model.cmip + '/' + Model.name + '/historical/' + ENS +
                      '/so_' + Model.name + '_historical_' + ENS + '_' + gtype + '_197001-199912_SC.nc')
        ncid   = Dataset(Tfile,'r')
        thetao_mon = np.squeeze(ncid.variables['thetao'][:,:,:,:])
        ncid.close()
        ncid   = Dataset(Sfile,'r')
        so_mon = np.squeeze(ncid.variables['so'][:,:,:,:])
        ncid.close()
    
    # Get Model information:
    print(Model.name + ' ' + ENS + ' ' + EXP)
    print(Model.notes)

    # Information about files:
    dims = jc.getDims(Tfiles[0],'thetao')

    # Dimention sizes:
    nx = dims[3].size - np.sum(Model.OextraWE)
    ny = dims[2].size
    nz = dims[1].size

    # Dimension names:
    dlon = dims[3].name
    dlat = dims[2].name
    dlev = dims[1].name
    dtime = dims[0].name

    # Get calendar information:
    ncid = Dataset(Tfiles[0],'r')
    cal    = ncid.variables[dtime].calendar
    units  = ncid.variables[dtime].units
    bounds = ncid.variables[dtime].bounds
    ncid.close()
    if ((Model.name == 'FGOALS-g2')):
        units = (units + '-01')
    
    # Setup grid for transports:

    # Read in mask data:
    if grid == 'Atlantic':
        basinfile = (datadir + 'CMIP' + Model.cmip + '/basins/' + Model.name + '_basins.nc')
        ncid  = Dataset(basinfile,'r')
        bmask = np.squeeze(ncid.variables['Atlantic'][:,:])
        ncid.close()
    elif grid == 'IndoPacific':
        basinfile = (datadir + 'CMIP' + Model.cmip + '/basins/' + Model.name + '_basins.nc')
        ncid  = Dataset(basinfile,'r')
        bmask = np.squeeze(ncid.variables['IndoPacific'][:,:])
        ncid.close()
    elif grid == 'Global':
        basinfile = (datadir + 'CMIP' + Model.cmip + '/basins/' + Model.name + '_basins.nc')
        ncid  = Dataset(basinfile,'r')
        bmask = np.squeeze(ncid.variables['Atlantic'][:,:]) + np.squeeze(ncid.variables['IndoPacific'][:,:])
        ncid.close()
    elif grid == 'BayOfBengal':
        basinfile = (datadir + 'CMIP' + Model.cmip + '/basins/' + Model.name + '_basins.nc')
        ncid  = Dataset(basinfile,'r')
        bmask = np.squeeze(ncid.variables['BayOfBengal'][:,:])
        ncid.close()
    elif grid == 'ArabianSea':
        basinfile = (datadir + 'CMIP' + Model.cmip + '/basins/' + Model.name + '_basins.nc')
        ncid  = Dataset(basinfile,'r')
        bmask = np.squeeze(ncid.variables['ArabianSea'][:,:])
        ncid.close()
    else:
        print((region + ' is not yet specified'))
        sys.exit()

    # Read in non-varying grid data:
    meshmask = Model.Omeshmask
    ncid  = Dataset(meshmask,'r')
    tmask = np.squeeze(ncid.variables['tmask'][:,:,:,:])
    vmask = np.squeeze(ncid.variables['vmask'][:,:,:,:])
    dxv   = np.squeeze(ncid.variables['dxv'][:,:])
    dyv   = np.squeeze(ncid.variables['dyv'][:,:])
    dxt   = np.squeeze(ncid.variables['dxt'][:,:])
    dyt   = np.squeeze(ncid.variables['dyt'][:,:])
    dzt   = np.squeeze(ncid.variables['dzt'][:,:,:])
    tlat  = np.squeeze(ncid.variables['tlat'][:,:])
    vlat  = np.squeeze(ncid.variables['vlat'][:,:])
    depth = np.squeeze(ncid.variables[dlev][:])
    ncid.close()

    print(np.shape(tmask))
    print(np.shape(bmask))
    
    print(nz)
    tmask = copy.deepcopy(tmask)*np.tile(bmask,(nz,1,1))
    vmask = copy.deepcopy(vmask)*np.tile(bmask,(nz,1,1))

    tlenx   = np.sum(dxt*tmask[0,:,:],axis=1)
    vlenx   = np.sum(dxv*vmask[0,:,:],axis=1)

    tareaxy   = np.sum(np.tile(dyt*dxt,(nz,1,1))*tmask,axis=2)

    # Compute the latitude used for the output data:
    Tlat = np.sum(tlat*dxt*tmask[0,:,:],axis=1)/tlenx
    Vlat = np.sum(vlat*dxv*vmask[0,:,:],axis=1)/vlenx
    
    # Initialize output file:
    if os.path.isfile(runfile):
        # Determine how much has already been computed
        ncid = Dataset(runfile,'r')
        nn = ncid.variables[dtime].size
        ncid.close()
    else:
        # Create file
        nn = 0 # Number of months computed is zero
        
        # Make directory if it doesn't exist:
        rundir = os.path.dirname(runfile)
        if not os.path.exists(rundir):
            os.makedirs(rundir)
        outdir = os.path.dirname(outfile)
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        ncid = Dataset(runfile, 'w', format='NETCDF4')
        # coordinates:
        ncid.createDimension(dlat,ny)
        ncid.createDimension(dlev,nz)
        ncid.createDimension('bnds',2)
        ncid.createDimension(dtime,None)
        # variables:
        ncid.createVariable('tlat'        ,'f8',(dlat,))
        ncid.createVariable('vlat'        ,'f8',(dlat,))
        ncid.createVariable('lev'         ,'f8',(dlev,))
        ncid.createVariable(dtime         ,'f8',(dtime,))
        ncid.createVariable(bounds        ,'f8',(dtime,'bnds',))
        ncid.createVariable('tarea'       ,'f8',(dlev,dlat,))
        ncid.createVariable('tvol'        ,'f8',(dtime,dlev,dlat,))
        ncid.createVariable('thetao_zonal','f8',(dtime,dlev,dlat,))
        ncid.createVariable('so_zonal'    ,'f8',(dtime,dlev,dlat,))
        ncid.createVariable('vo_zonal'    ,'f8',(dtime,dlev,dlat,))
        ncid.createVariable('thetao_ref'  ,'f8',(dtime,dlat,))
        ncid.createVariable('so_ref'      ,'f8',(dtime,dlat,))
        ncid.createVariable('vo_ref'      ,'f8',(dtime,dlat,))
        ncid.createVariable('moc_topdown' ,'f8',(dtime,dlev,dlat,))
        ncid.createVariable('moc_bottomup','f8',(dtime,dlev,dlat,))
        ncid.createVariable('moc_section' ,'f8',(dtime,dlev,dlat,))
        ncid.createVariable('moc_depth'   ,'f8',(dtime,dlev,dlat,))
        ncid.createVariable('Ttrans'      ,'f8',(dtime,dlat,))
        ncid.createVariable('Strans'      ,'f8',(dtime,dlat,))
        ncid.createVariable('Vtrans'      ,'f8',(dtime,dlat,))
        ncid.createVariable('Hov'         ,'f8',(dtime,dlat,))
        ncid.createVariable('Mov'         ,'f8',(dtime,dlat,))
        ncid.createVariable('Haz'         ,'f8',(dtime,dlat,))
        ncid.createVariable('Maz'         ,'f8',(dtime,dlat,))
        ncid.createVariable('H'           ,'f8',(dtime,dlat,))
        ncid.createVariable('M'           ,'f8',(dtime,dlat,))

        # fill variables:
        ncid.variables['tlat'][:]  = Tlat
        ncid.variables['vlat'][:]  = Vlat
        ncid.variables['tarea'][:] = tareaxy
        ncid.variables['lev'][:]   = depth

        # Add Calendar info:
        ncid.variables[dtime].calendar = cal
        ncid.variables[dtime].units    = units
        ncid.variables[dtime].bounds   = bounds

        # close:
        ncid.close()
        
    # Loop through each uncomputed month:
    nm = nn
    # Loop through each file temperature file:
    for ff in range(0,nf):
        print(('computing file ' + str(ff+1) + ' of ' + str(nf)))
        # Determine how many months of data are in the current temperature file:
        if (ctype == 'Tbar'):
            dims = jc.getDims(Vfiles[ff],'vo')
            nt = dims[0].size
        else:
            dims = jc.getDims(Tfiles[ff],'thetao')
            nt = dims[0].size

        # Determine which months have been computed:
        if nm >= nt:
            # All months from this file have been computed, move on to next file:
            nm = nm - nt
        else:
            # Loop through and compute all the missing months:
            for mm in range(nm,nt):  # mm is month in the file
                print(('month ' + str(mm+1) + ' of ' + str(nt)))
                
                if (ctype == 'Tbar'):                    
                    # Velocity:
                    vo = jc.Oread3Ddata(Model,Vfiles[ff],'vo',mm)*vmask
                    ncid = Dataset(Vfiles[ff],'r')
                    tt   = ncid.variables[dtime][mm]
                    # Get Calendar information:
                    units2 = ncid.variables[dtime].units
                    tb   = ncid.variables[bounds][mm,:]
                    ncid.close()
                    if ((Model.name == 'FGOALS-g2')):
                        units2 = (units2 + '-01')

                    # Temperature & Salinity:   
                    mon = cftime.num2date(tt,units2,cal).month
                    thetao = thetao_mon[mon-1,:,:,:]*tmask
                    so = so_mon[mon-1,:,:,:]*tmask
                    
                else:
                    # Temperature:
                    thetao = jc.Oread3Ddata(Model,Tfiles[ff],'thetao',mm)*tmask
                    ncid = Dataset(Tfiles[ff],'r')
                    tt   = ncid.variables[dtime][mm]
                    # Get Calendar information:
                    units2 = ncid.variables[dtime].units
                    tb   = ncid.variables[bounds][mm,:]
                    ncid.close()
                    if ((Model.name == 'FGOALS-g2')):
                        units2 = (units2 + '-01')

                    # Salinity:
                    fs,ns = jc.timeFile(Model,cal,units2,tt,Sfiles)
                    if ((np.isnan(ns))):
                        so = np.nan*tmask
                    else:
                        so    = jc.Oread3Ddata(Model,Sfiles[fs],'so',ns)*tmask

                    # Velocity:
                    if (ctype == 'vbar'):
                        mon = cftime.num2date(tt,units2,cal).month
                        vo = vo_mon[mon-1,:,:,:]*vmask
                    else:
                        fv,nv = jc.timeFile(Model,cal,units2,tt,Vfiles)
                        if ((np.isnan(nv))):
                            vo = np.nan*tmask
                        else:
                            vo    = jc.Oread3Ddata(Model,Vfiles[fv],'vo',nv)*vmask

                # Grid box depth information:
                if ((len(Dfiles) != 0) & (Model.name != 'CSIRO-Mk3-6-0')):
                    if ((Model.name == 'GFDL-CM4') | ((Model.name == 'ACCESS-CM2') & (EXP == 'historical'))):
                        # volume not thickness is saved
                        fd,nd  = jc.timeFile(cal,units2,tt,Dfiles)
                        dzt    = jc.Oread3Ddata(Model,Dfiles[fd],'volcello',nd)*tmask
                        dzt    = dzt/np.tile(dxt*dyt,(np.size(dzt,axis=0),1,1))
                    else:
                        fd,nd = jc.timeFile(Model,cal,units2,tt,Dfiles)
                        if ((np.isnan(nd))):
                            dzt    = dzt*tmask
                        else:
                            dzt   = jc.Oread3Ddata(Model,Dfiles[fd],'thkcello',nd)*tmask
                elif (len(VOLfiles) != 0):
                    fd,nd  = jc.timeFile(Model,cal,units2,tt,VOLfiles)
                    dzt    = jc.Oread3Ddata(Model,VOLfiles[fd],'volcello',nd)*tmask
                    dzt    = dzt/np.tile(dxt*dyt,(np.size(dzt,axis=0),1,1))
                else:
                    dzt    = dzt*tmask
                    
                    
                if (Model.name == 'CSIRO-Mk3-6-0'):
                    thetao[np.where(thetao > 1e10)] = 0
                
                # Fix time:
                tt,tb = jc.fixTime(units,units2,cal,tt,tb)
                
                # Fix units:
                if (np.max(thetao) > 200):
                    thetao = thetao - 273.15
                if (np.max(so) < 30):
                    so = so*1000

                # Move data to transport grid:
                dzv     = jc.moveData(Model,'T','V',dzt,computation='min')*vmask
                thetaov = jc.moveData(Model,'T','V',thetao    ,computation='mean')
                sov     = jc.moveData(Model,'T','V',so        ,computation='mean')

                tarea   = dxt*dzt*tmask
                tareax  = np.sum(tarea,axis=2)
                tareaxz = np.sum(tareax,axis=0)

                varea   = dxv*dzv*vmask
                vareax  = np.sum(varea,axis=2)
                vareaz  = np.sum(varea,axis=0)
                vareaxz = np.sum(vareax,axis=0)

                tvol    = dxt*dyt*dzt*tmask
                tvolx   = np.sum(tvol,axis=2)

                # Compute zonal means:
                Tzon  = np.sum(thetao *tarea,axis=2)/tareax
                Szon  = np.sum(so     *tarea,axis=2)/tareax
                Tzonv = np.sum(thetaov*varea,axis=2)/vareax
                Szonv = np.sum(sov    *varea,axis=2)/vareax
                Vzon  = np.sum(vo     *varea,axis=2)/vareax

                # Compute section means:
                Tref = np.sum(np.sum(thetaov*varea,axis=2),axis=0)/vareaxz
                Sref = np.sum(np.sum(sov    *varea,axis=2),axis=0)/vareaxz
                Vref = np.sum(np.sum(vo     *varea,axis=2),axis=0)/vareaxz

                # Compute total transports:
                Vtrans =    np.sum(np.sum(        vo*varea,axis=2),axis=0)/1e6
                Ttrans =    np.sum(np.sum(thetaov*vo*varea,axis=2),axis=0)*rho*cp/1e15
                Strans = -1*np.sum(np.sum(sov    *vo*varea,axis=2),axis=0)/Sref/1e6

                # Compute MOC-topdown:
                moc_topdown =    np.cumsum(np.sum(vo*varea,axis=2),axis=0)/1e6/(vareax > 0) 
                # Compute MOC-bottomup:
                moc_bottomup = -(np.cumsum(np.sum(vo*varea,axis=2)[::-1,:],axis=0))[::-1,:]/1e6/(vareax > 0)
                # Compute MOC-section:
                vmean       = np.sum(np.sum(vo*varea,axis=2),axis=0)/vareaxz
                vstar       = vo - np.swapaxes(np.tile(vmean,(nz,1,1)),2,1)
                moc_section =    np.cumsum(np.sum(vstar*varea,axis=2),axis=0)/1e6/(vareax > 0)     
                # Compute MOC-depth:
                vmean      = np.sum(vo*varea,axis=0)/vareaz
                vbar       = vo - np.tile(vmean,(nz,1,1))
                moc_depth  = np.cumsum(np.sum(vbar*varea,axis=2),axis=0)/1e6/(vareax > 0)

                # Compute overturing and gyre transports:
                Vszon    = np.sum(vstar*varea,axis=2)/vareax
                vp       = vstar   - np.moveaxis(np.tile(Vszon,(nx,1,1)),0,2)*vmask
                sp       = sov     - np.moveaxis(np.tile(Szonv,(nx,1,1)),0,2)*vmask
                tp       = thetaov - np.moveaxis(np.tile(Tzonv,(nx,1,1)),0,2)*vmask
                Mov      = -1*np.sum(np.sum(np.moveaxis(np.tile(Szon*Vszon,(nx,1,1)),0,2)*varea,axis=2),axis=0)/Sref/1e6
                Hov      =    np.sum(np.sum(np.moveaxis(np.tile(Tzon*Vszon,(nx,1,1)),0,2)*varea,axis=2),axis=0)*rho*cp/1e15
                Maz      = -1*np.sum(np.sum(vp*sp*varea,axis=2),axis=0)/Sref/1e6
                Haz      =    np.sum(np.sum(vp*tp*varea,axis=2),axis=0)*rho*cp/1e15
                M        = Mov+Maz
                H        = Hov+Haz

                # Save data to file:
                ncido = Dataset(runfile, 'a', format='NETCDF4')
                ncido.variables[dtime][nn]              = tt
                ncido.variables[bounds][nn,:]           = tb
                ncido.variables['tvol'][nn,:,:]         = tvolx
                ncido.variables['thetao_zonal'][nn,:,:] = Tzon
                ncido.variables['so_zonal'][nn,:,:]     = Szon
                ncido.variables['vo_zonal'][nn,:,:]     = Vzon
                ncido.variables['thetao_ref'][nn,:]     = Tref
                ncido.variables['so_ref'][nn,:]         = Sref
                ncido.variables['vo_ref'][nn,:]         = Vref
                ncido.variables['moc_topdown'][nn,:,:]  = moc_topdown
                ncido.variables['moc_bottomup'][nn,:,:] = moc_bottomup
                ncido.variables['moc_section'][nn,:,:]  = moc_section
                ncido.variables['moc_depth'][nn,:,:]    = moc_depth
                ncido.variables['Ttrans'][nn,:]         = Ttrans
                ncido.variables['Strans'][nn,:]         = Strans
                ncido.variables['Vtrans'][nn,:]         = Vtrans
                ncido.variables['Hov'][nn,:]            = Hov
                ncido.variables['Mov'][nn,:]            = Mov
                ncido.variables['Haz'][nn,:]            = Haz
                ncido.variables['Maz'][nn,:]            = Maz
                ncido.variables['H'][nn,:]              = H
                ncido.variables['M'][nn,:]              = M
                ncido.close()

                # Tidy:
                nn = nn + 1
            nm = 0

    shutil.move(runfile, outfile)