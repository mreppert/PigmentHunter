import numpy as np
import parmed as pmd

import os,sys
sys.path.append('./main/')

import pigment

# Pigments are required to contain only four atoms: NA, NB, NC, and ND
def calculate_coupling(PigList, ChainList):

    CoupTraj = []
    DipTraj = []
    RotTraj = []
    
    # First create list of chain-selected pigments
    SelPigs = []
    for pig in PigList:
        if ChainList.count(pig.residue.chain)>0:
            SelPigs.append(pig)
    
    Npigs = len(SelPigs)
    h = 6.62607015e-34 #J*s
    c = 2.998e10 #cm/s
    eo = 4.80320451e-10 # esu
    
    Erg2J = 1.0e-7
    ang2cm = 1.0e-8
    
    # Monitor errors
    error = False
    
    # All that's needed for a coupling calculation are the coordinates of the 
    # NA, NB, NC, and ND atoms. The oscillator strength is determined from
    # the pig.species.diplength value (transition dipole magnitude, in statC*cm). 

    # Check how many frames there are in the structure.
    # We use the first selected pigment as a proxy since all
    # should have the same number of frames. 
    # First dimension of atcoords is number of frames.
    Nframes = np.shape(SelPigs[0].atcoords)[0]

    # Now calculate centers and dipoles
    CentMat = np.zeros((Nframes,Npigs,3))
    Dips = np.zeros((Nframes,Npigs,3))
    CentAtoms = ['NA', 'NB', 'NC', 'ND']
    for p in range(0, Npigs):
        pig = SelPigs[p]
        for name in CentAtoms:
            if pig.atnames.count(name)==0:
                print('Error: Could not locate center atom ' + name + " in pigment " + pig.residue.name + " " + pig.residue.chain + " " + str(pig.residue.number))
                error = True
                print('Aborting PDC calculation.')
                break
            else:
                ndx = pig.atnames.index(name)
                
                # CentMat is average of all N atom positions
                CentMat[:,p,:] += pig.atcoords[:,ndx,:]/float(len(CentAtoms))
                
                # Dips gets positive contribution from ND and negative from NB
                if name=='NB':
                    Dips[:,p,:] -= pig.atcoords[:,ndx,:]
                elif name=='ND':
                    Dips[:,p,:] += pig.atcoords[:,ndx,:]
        
    # Note that Dips dipoles are *not* normalized!
    
    # If no errors: 
    if error==False:
        
        # And then calculate the couplings
        for fr in range(0, Nframes):

            CoupMat = np.zeros((Npigs,Npigs))
            RotMat = np.zeros((Npigs,Npigs))
            DipMat = np.zeros((Npigs,3))
            
            # First, add normalized dipoles to DipMat
            for p in range(0, Npigs):
                DipMat[p,:] = Dips[fr,p,:]/np.linalg.norm(Dips[fr,p,:])
            
            # Now calculate interactions
            for p1 in range(0, Npigs):
                pig1 = SelPigs[p1]
                
                # Note we only calculate half of the elements: the other half are symmetric
                for p2 in range(0, p1):
                    pig2 = SelPigs[p2]
                    
                    # NB: here we convert from Ang to cm:
                    Rmn = ang2cm*(CentMat[fr,p1,:] - CentMat[fr,p2,:])

                    # rmn is the distance
                    rmn = np.linalg.norm(Rmn)

                    # The oscillator strength is the product of p1 and p2 dipole lengths:
                    osc = pig1.species.diplength*pig2.species.diplength

                    # Check units: 
                    # -- Osc ==> (statC*cm)^2 = dyne*cm^4
                    # -- DipMat ==> unitless
                    # rmn and Rmn ==> cm
                    # Both terms scale as osc/(rmn^3) ==> dyne*cm^4 / cm^3 = dyne*cm =  erg.
                    # We multiply by Erg2J to get to Joules and then divide by (h*c) to get to 1/cm. 
                    prefac = Erg2J/(h*c)

                    # And PDC value: 
                    CoupMat[p1,p2] = prefac*osc*(np.dot(DipMat[p1,:],DipMat[p2,:]) / (rmn**3) - 3.0*np.dot(DipMat[p1,:], Rmn)*np.dot(DipMat[p2,:], Rmn) / (rmn**5))
                    
            CoupMat = CoupMat + np.transpose(CoupMat)
            
            # Calculate rotation matrix
            # Note that RotMat is *not* scaled by dipole length. 
            # This is done when the spectrum is calculated. 
            for m in range(0, Npigs):
                for n in range(0, Npigs):
                    Rmn = CentMat[fr,n,:] - CentMat[fr,m,:]
                    
                    # MER changed sign on RotMat on 12/16/2021
                    RotMat[m,n] = - np.dot(Rmn, np.cross(DipMat[m,:], DipMat[n,:]))

            CoupTraj.append(CoupMat)
            DipTraj.append(DipMat)
            RotTraj.append(RotMat)
        
    # NB: oscillator strengths are all finally scaled to unity!
    # This is true in both DipTraj and RotTraj.
    # The oscillator strengths must be introduced during spectrum calculation
    # based on the values stored in diplengths.txt.
    return CoupTraj, DipTraj, RotTraj
