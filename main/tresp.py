import numpy as np
import parmed as pmd

import os,sys
sys.path.append('./main/')

import pigment

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
    
    # Now check whether TrESP parameters are available for *all*
    # selected pigments and (if so) store them for reference.
    ListAtoms = []
    ListQ00 = []
    ListQ10 = []
    ListQ11 = []
    error = False
    
    # Loop over pigments
    for p in range(0, Npigs):
        pig = SelPigs[p]
        
        # TrESP parameters should be in file 'misc/TrESP/' + pig.species.stdname + '.txt'
        fname = 'misc/TrESP/' + pig.species.stdname + '.txt'
        if os.path.isfile(fname)==False:
            print('Error: Could not locate TrESP parameter file for species ' + pig.species.stdname)
            error = True
            break
        else:
            atoms=  []
            q00 = []
            q10 = []
            q11 = []
            with open(fname) as fd:
                for line in fd:
                    if len(line)>0 and line[0]!="#":
                        items = line.split()
                        format_error = False
                        try:
                            atoms.append(items[0])
                            
                            q00.append(float(items[1]))
                            q10.append(float(items[2]))
                            q11.append(float(items[3]))
                            
                        except:
                            print('Error reading TrESP input file ' + fname + '.')
                            print('Aborting TrESP calculation.')
                            error = True
                            break
                            
            ListAtoms.append(atoms)
            ListQ00.append(q00)
            ListQ10.append(q10)
            ListQ11.append(q11)
            
        if(error):
            break
                            
    # If no errors, we've located all relevant TrESP files and imported the data. 
    # Now check whether all necessary atoms are available in the structure.
    if error==False:
        
        for p in range(0, Npigs):
            
            # Added 10/5/2021 by MER
            pig = SelPigs[p]
        
            # Loop over TrESP atoms
            for n in range(0, len(ListAtoms[p])):

                # name for the nth TrESP atom in pth pigment
                name = ListAtoms[p][n]

                # If we can't find it in the structure, throw an error. 
                # pig.atnames stores the names of all atoms associated with pigment pig
                if pig.atnames.count(name)==0:
                    print('Error: Could not locate TrESP atom ' + name + " in pigment " + pig.residue.name + " " + pig.residue.chain + " " + str(pig.residue.number))
                    print('Aborting TrESP calculation.')
                    error = True
                    break

            # If there's already been an error, don't continue. 
            if error:
                break
                
    # Check how many frames there are in the structure.
    # We use the first selected pigment as a proxy since all
    # should have the same number of frames. 
    # First dimension of atcoords is number of frames.
    Nframes = np.shape(SelPigs[0].atcoords)[0]

    # If no errors, calculate centers.
    CentMat = np.zeros((Nframes,Npigs,3))
    CentAtoms = ['NA', 'NB', 'NC', 'ND']
    for p in range(0, Npigs):
        pig = SelPigs[p]
        for name in CentAtoms:
            if pig.atnames.count(name)==0:
                print('Error: Could not locate center atom ' + name + " in pigment " + pig.residue.name + " " + pig.residue.chain + " " + str(pig.residue.number))
                error = True
                print('Aborting TrESP calculation.')
                break
            else:
                ndx = pig.atnames.index(name)
                CentMat[:,p,:] += pig.atcoords[:,ndx,:]/float(len(CentAtoms))

    # If no errors, all pigments have TrESP parameters and necessary atoms, 
    # and pigment centers have been calculated (for all frames) and stored in CentMat.
    if error==False:
        # Build list of TrESP coordinates for each pigment
        # NB: Units are converted to cm! (cgs units)
        CoordList = []
        for p in range(0, Npigs):
            pig = SelPigs[p]
            coords = np.zeros((Nframes, len(ListAtoms[p]), 3))

            # Find index associated with each TrESP atom in the structure.
            # Coordinate list will be sorted in atom order specified in TrESP 
            # parameter file *not* order of crystal structure. 
            for n in range(0, len(ListAtoms[p])):
                name = ListAtoms[p][n]
                ndx = pig.atnames.index(name)
                coords[:,n,:] = pig.atcoords[:,ndx,:].copy()
            CoordList.append(coords*ang2cm)

    ###############################################
    # Here is where we begin looping over frames
    ###############################################


        for fr in range(0, Nframes):

            CoupMat = np.zeros((Npigs,Npigs))
            RotMat = np.zeros((Npigs,Npigs))
            Dips = []

            # Calculate pigment dipole lengths (as per default TrESP charges)
            # We'll divide charges for each pigment by the ratio of the *calculated*
            # dipole length to the *standard* dipole length. 
            # Charges are in units of cm^3/2 g^1/2 s^−1.
            DipLengths = []
            for p in range(0, Npigs):
                pig = SelPigs[p]
                dip = np.zeros((3,))
                for n in range(0, len(ListAtoms[p])):
                    dip += eo*ListQ10[p][n]*CoordList[p][fr,n,:]
                    
                DipLengths.append(np.linalg.norm(dip))
                Dips.append(dip/np.linalg.norm(dip))
            Dips = np.array(Dips)
            
            # Calculate interactions
            for p1 in range(0, Npigs):
                pig1 = SelPigs[p1]
                df1 = pig1.species.diplength/DipLengths[p1]
                for p2 in range(0, p1):
                    pig2 = SelPigs[p2]
                    df2 = pig2.species.diplength/DipLengths[p2]
                    for atm in range(0, len(ListAtoms[p1])):
                        for atn in range(0, len(ListAtoms[p2])):
                            Rmn = CoordList[p1][fr,atm,:] - CoordList[p2][fr,atn,:]
                            rmn = np.linalg.norm(Rmn)
                            CoupMat[p1,p2] += eo*eo*df1*df2*(ListQ10[p1][atm]*ListQ10[p2][atn]/rmn)*(Erg2J/h)/c
                            
                
            CoupMat = CoupMat + np.transpose(CoupMat)
            
            # Calculate rotation matrix
            for m in range(0, Npigs):
                for n in range(0, Npigs):
                    Rmn = CentMat[fr,n,:] - CentMat[fr,m,:]
                    #osc = SelPigs[m].species.diplength*SelPigs[n].species.diplength
                    osc = 1.0
                    RotMat[m,n] = np.dot(Rmn, np.cross(Dips[m,:], Dips[n,:]))*osc

            CoupTraj.append(CoupMat)
            DipTraj.append(Dips)
            RotTraj.append(RotMat)
        
    # NB: oscillator strengths are all finally scaled to unity!
    return CoupTraj, DipTraj, RotTraj


def calculate_shift(PigList, ChainList, instruc):
    
    # First create list of chain-selected pigments
    SelPigs = []
    for pig in PigList:
        if ChainList.count(pig.residue.chain)>0:
            SelPigs.append(pig)
    
    Npigs = len(SelPigs)
    FreqTraj = []
        
    h = 6.62607015e-34 #J*s
    c = 2.998e10 #cm/s
    eo = 4.80320451e-10 # esu
    eps_eff = 2.5
    
    Erg2J = 1.0e-7
    ang2cm = 1.0e-8
    
    # Now check whether TrESP parameters are available for *all*
    # selected pigments and (if so) store them for reference.
    ListAtoms = []
    ListQ00 = []
    ListQ10 = []
    ListQ11 = []
    error = False
    
    # Loop over pigments
    for p in range(0, Npigs):
        pig = SelPigs[p]
        
        # TrESP parameters should be in file 'misc/TrESP/' + pig.species.stdname + '.txt'
        fname = 'misc/TrESP/' + pig.species.stdname + '.txt'
        if os.path.isfile(fname)==False:
            print('Error: Could not locate TrESP parameter file for species ' + pig.species.stdname)
            error = True
            break
        else:
            atoms=  []
            q00 = []
            q10 = []
            q11 = []
            with open(fname) as fd:
                for line in fd:
                    if len(line)>0 and line[0]!="#":
                        items = line.split()
                        format_error = False
                        try:
                            atoms.append(items[0])
                            q00.append(float(items[1]))
                            q10.append(float(items[2]))
                            q11.append(float(items[3]))
                        except:
                            print('Error reading TrESP input file ' + fname + '.')
                            print('Aborting TrESP calculation.')
                            error = True
                            break
                            
            ListAtoms.append(atoms)
            ListQ00.append(q00)
            ListQ10.append(q10)
            ListQ11.append(q11)
            
        if(error):
            break
                            
    # If no errors, we've located all relevant TrESP files and imported the data. 
    # Now check whether all necessary atoms are available in the structure.
    if error==False:
        for p in range(0, Npigs):
            
            # Added 10/5/2021 by MER
            pig = SelPigs[p]
            
            # Loop over TrESP atoms
            for n in range(0, len(ListAtoms[p])):

                # name for the nth TrESP atom in pth pigment
                name = ListAtoms[p][n]

                # If we can't find it in the structure, throw an error. 
                # pig.atnames stores the names of all atoms associated with pigment pig
                if pig.atnames.count(name)==0:
                    print('Error: Could not locate TrESP atom ' + name + " in pigment " + pig.residue.name + " " + pig.residue.chain + " " + str(pig.residue.number))
                    print('Aborting TrESP calculation.')
                    error = True
                    break

            # If there's already been an error, don't continue. 
            if error:
                break
                
    # If no errors, all pigments have TrESP parameters and necessary atoms. 
    if error==False:

        # Number of frames in trajectory
        Nframes = np.shape(SelPigs[0].atcoords)[0]

        # Build list of TrESP coordinates for each pigment
        # NB: Units are converted to cm! (cgs units)
        CoordList = []
        for p in range(0, Npigs):
            pig = SelPigs[p]
            coords = np.zeros((Nframes, len(ListAtoms[p]), 3))
            # Find index associated with each TrESP atom in the structure.
            # Coordinate list will be sorted in atom order specified in TrESP
            # parameter file *not* order of crystal structure. 
            for n in range(0, len(ListAtoms[p])):
                name = ListAtoms[p][n]
                ndx = pig.atnames.index(name)
                coords[:,n,:] = pig.atcoords[:,ndx,:].copy()
            CoordList.append(coords*ang2cm)


        for fr in range(0, Nframes):

            tFreqs = np.zeros((Npigs,))

            # Coordinates for local frame
            frcoords = ang2cm*instruc.get_coordinates(fr)
        
            for p in range(0, Npigs):
                
                pig = SelPigs[p]
                shift = 0.0
                
                data_frame = instruc.to_dataframe()
                
                # We exclude the pigment itself from calculation
                ndcs = (data_frame.resid != pig.residue.idx)
                                
                for atm in range(0, len(ListAtoms[p])):
                    dQ = ListQ11[p][atm] - ListQ00[p][atm]
                    Rvec = CoordList[p][fr,atm,:] - frcoords[ndcs,:]
                    Dvec = np.sqrt(np.sum(np.power(Rvec,2),1))
                    Qvec = (data_frame[ndcs].charge).to_numpy()
                    shift += eo*eo*dQ*np.sum(np.divide(Qvec,Dvec))*Erg2J/(h*c*eps_eff)
                    
                tFreqs[p] = shift

            FreqTraj.append(tFreqs)
        
    return FreqTraj, error
