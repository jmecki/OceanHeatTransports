##############################################################################
#
# Generates list of CMIP model objects:
#
# #############################################################################

import os
import _pickle as pickle

from . import CMIPobject

# Read in list of models:
def readList(filename):
    infile = open(filename,'rb')
    Clist  = pickle.load(infile)
    infile.close()
    
    return Clist

# Save list of models:
def saveList(filename,Clist):
    outfile = open(filename,'wb')
    pickle.dump(Clist,outfile)
    outfile.close()

# If the list doesn't exist yet, generates a list of models for the specified
# cmip otherwise adds missing models to list.
def generateList(cmip,filename):
    Models = CMIPobject.getModels(cmip=cmip)

    if os.path.isfile(filename):
        Clist = readList(filename)
    else:
        Clist = dict()
        
    print('starting models on list: ' + str(len(Clist)))

    for mm in Models:
        if not mm in Clist:
            Clist[mm] = CMIPobject.CMIPmodel(mm,cmip=cmip)

            # Save to file:
            saveList(filename,Clist)
            print('added ' + mm)
            
    print('total models on list: ' + str(len(Clist)))
    
# Updates information in the CMIP models list:
def updateList(filename,Model,cmip='',datadir='',savedir='',\
               Ogrid='',Omeshmask='',Oreg='',Olon='',Olat='',OflipNS='',OextraT='',OextraWE='',\
               Agrid='',Ameshmask='',Areg='',Alon='',Alat='',AflipNS='',AextraT='',AextraWE='',\
               Igrid='',Imeshmask='',Ireg='',Ilon='',Ilat='',IflipNS='',IextraT='',IextraWE='',\
               Lgrid='',Lmeshmask='',Lreg='',Llon='',Llat='',LflipNS='',LextraT='',LextraWE='',\
               notes=''):
    
    Clist = readList(filename)
    
    # Generic Model:
    if (cmip != ''):
        Clist[Model].cmip      = cmip
    if (datadir != ''):
        Clist[Model].datadir   = datadir
    if (savedir != ''):
        Clist[Model].savedir   = savedir
    
    # Ocean:
    if (Ogrid != ''):
        Clist[Model].Ogrid     = Ogrid
    if (Omeshmask != ''):
        Clist[Model].Omeshmask = Omeshmask
    if (Oreg != ''):
        Clist[Model].Oreg      = Oreg
    if (Olon != ''):
        Clist[Model].Olon      = Olon
    if (Olat != ''):
        Clist[Model].Olat      = Olat
    if (OflipNS != ''):
        Clist[Model].OflipNS   = OflipNS
    if (OextraT != ''):
        Clist[Model].OextraT   = OextraT
    if (OextraWE != ''):
        Clist[Model].OextraWE  = OextraWE
    
    # Atmosphere:
    if (Agrid != ''):
        Clist[Model].Agrid     = Agrid
    if (Ameshmask != ''):
        Clist[Model].Ameshmask = Ameshmask
    if (Areg != ''):
        Clist[Model].Areg      = Areg
    if (Alon != ''):
        Clist[Model].Alon      = Alon
    if (Alat != ''):
        Clist[Model].Alat      = Alat
    if (AflipNS != ''):
        Clist[Model].AflipNS   = AflipNS
    if (AextraT != ''):
        Clist[Model].AextraT  = AextraT
    if (AextraWE != ''):
        Clist[Model].AextraWE  = AextraWE
    
    # Sea Ice:
    if (Igrid != ''):
        Clist[Model].Igrid     = Igrid
    if (Imeshmask != ''):
        Clist[Model].Imeshmask = Imeshmask
    if (Ireg != ''):
        Clist[Model].Ireg      = Ireg
    if (Ilon != ''):
        Clist[Model].Ilon      = Ilon
    if (Ilat != ''):
        Clist[Model].Ilat      = Ilat
    if (IflipNS != ''):
        Clist[Model].IflipNS   = IflipNS
    if (IextraT != ''):
        Clist[Model].IextraT  = IextraT
    if (IextraWE != ''):
        Clist[Model].IextraWE  = IextraWE
    
    # Land:
    if (Lgrid != ''):
        Clist[Model].Lgrid     = Lgrid
    if (Lmeshmask != ''):
        Clist[Model].Lmeshmask = Lmeshmask
    if (Lreg != ''):
        Clist[Model].Lreg      = Lreg
    if (Llon != ''):
        Clist[Model].Llon      = Llon
    if (Llat != ''):
        Clist[Model].Llat      = Llat
    if (LflipNS != ''):
        Clist[Model].LflipNS   = LflipNS
    if (LextraT != ''):
        Clist[Model].LextraT  = LextraT
    if (LextraWE != ''):
        Clist[Model].LextraWE  = LextraWE
        
    # General useful information:
    if (notes != ''):
        Clist[Model].notes  = notes
    
    # Save to file:
    saveList(filename,Clist)
    
# Updates information for all members in list:
def updateAllList(filename,cmip='',datadir='',savedir='',\
               Ogrid='',Omeshmask='',Oreg='',Olon='',Olat='',OflipNS='',OextraT='',OextraWE='',\
               Agrid='',Ameshmask='',Areg='',Alon='',Alat='',AflipNS='',AextraT='',AextraWE='',\
               Igrid='',Imeshmask='',Ireg='',Ilon='',Ilat='',IflipNS='',IextraT='',IextraWE='',\
               Lgrid='',Lmeshmask='',Lreg='',Llon='',Llat='',LflipNS='',LextraT='',LextraWE='',\
               notes=''):
    
    Clist = readList(filename)
    Models = list(Clist.keys())
    for Model in Models:    
        # Generic Model:
        if (cmip != ''):
            Clist[Model].cmip      = cmip
        if (datadir != ''):
            Clist[Model].datadir   = datadir
        if (savedir != ''):
            Clist[Model].savedir   = savedir

        # Ocean:
        if (Ogrid != ''):
            Clist[Model].Ogrid     = Ogrid
        if (Omeshmask != ''):
            Clist[Model].Omeshmask = Omeshmask
        if (Oreg != ''):
            Clist[Model].Oreg      = Oreg
        if (Olon != ''):
            Clist[Model].Olon      = Olon
        if (Olat != ''):
            Clist[Model].Olat      = Olat
        if (OflipNS != ''):
            Clist[Model].OflipNS   = OflipNS
        if (OextraT != ''):
            Clist[Model].OextraT  = OextraT
        if (OextraWE != ''):
            Clist[Model].OextraWE  = OextraWE

        # Atmosphere:
        if (Agrid != ''):
            Clist[Model].Agrid     = Agrid
        if (Ameshmask != ''):
            Clist[Model].Ameshmask = Ameshmask
        if (Areg != ''):
            Clist[Model].Areg      = Areg
        if (Alon != ''):
            Clist[Model].Alon      = Alon
        if (Alat != ''):
            Clist[Model].Alat      = Alat
        if (AflipNS != ''):
            Clist[Model].AflipNS   = AflipNS
        if (AextraT != ''):
            Clist[Model].AextraT  = AextraT
        if (AextraWE != ''):
            Clist[Model].AextraWE  = AextraWE

        # Sea Ice:
        if (Igrid != ''):
            Clist[Model].Igrid     = Igrid
        if (Imeshmask != ''):
            Clist[Model].Imeshmask = Imeshmask
        if (Ireg != ''):
            Clist[Model].Ireg      = Ireg
        if (Ilon != ''):
            Clist[Model].Ilon      = Ilon
        if (Ilat != ''):
            Clist[Model].Ilat      = Ilat
        if (IflipNS != ''):
            Clist[Model].IflipNS   = IflipNS
        if (IextraT != ''):
            Clist[Model].IextraT  = IextraT
        if (IextraWE != ''):
            Clist[Model].IextraWE  = IextraWE

        # Land:
        if (Lgrid != ''):
            Clist[Model].Lgrid     = Lgrid
        if (Lmeshmask != ''):
            Clist[Model].Lmeshmask = Lmeshmask
        if (Lreg != ''):
            Clist[Model].Lreg      = Lreg
        if (Llon != ''):
            Clist[Model].Llon      = Llon
        if (Llat != ''):
            Clist[Model].Llat      = Llat
        if (LflipNS != ''):
            Clist[Model].LflipNS   = LflipNS
        if (LextraT != ''):
            Clist[Model].LextraT  = LextraT
        if (LextraWE != ''):
            Clist[Model].LextraWE  = LextraWE
        
        # General useful information:
        if (notes != ''):
            Clist[Model].notes  = notes

    
    # Save to file:
    saveList(filename,Clist)