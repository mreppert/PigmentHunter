import numpy as np
import parmed as pmd

import os,sys
sys.path.append('./main/')

import pigment

def calculate_coupling(PigList, ChainList):
    
    # First create list of chain-selected pigments
    SelPigs = []
    for pig in PigList:
        if ChainList.count(pig.residue.chain)>0:
            SelPigs.append(pig)
    
    Npigs = len(SelPigs)
    CoupMat = np.zeros((Npigs,Npigs))
    RotMat = np.zeros((Npigs,Npigs))
    Dips = []
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
                            
                            # TrESP file charges are in units of 1000*eo.
                            # We convert here to eo. 
                            q00.append(float(items[1])*1.0e-3)
                            q10.append(float(items[2])*1.0e-3)
                            q11.append(float(items[3])*1.0e-3)
                            
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
                
    # If no errors, calculate centers
    CentMat = np.zeros((Npigs,3))
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
                CentMat[p,:] += pig.atcoords[ndx]/float(len(CentAtoms))

#     for p in range(0, Npigs):
#         pig = SelPigs[p]
#         MGndx = pig.atnames.index('MG')
#         print('Calculated: ' + str(CentMat[p]))
#         print('MG: ' + str(pig.atcoords[MGndx]))
#         print('')
    
    # If no errors, all pigments have TrESP parameters and necessary atoms, 
    # and pigment centers have been calculated and stored in CentMat.
    if error==False:
        # Build list of TrESP coordinates for each pigment
        # NB: Units are converted to cm! (cgs units)
        CoordList = []
        for p in range(0, Npigs):
            pig = SelPigs[p]
            coords = []
            # Find index associated with each TrESP atom in the structure.
            # Coordinate list will be sorted in atom order specified in TrESP 
            # parameter file *not* order of crystal structure. 
            for name in ListAtoms[p]:
                ndx = pig.atnames.index(name)
                coords.append(pig.atcoords[ndx])
            CoordList.append(np.array(coords)*ang2cm)
            
        # Calculate pigment dipole lengths (as per default TrESP charges)
        # We'll divide charges for each pigment by the ratio of the *calculated*
        # dipole length to the *standard* dipole length. 
        # Charges are in units of cm^3/2 g^1/2 s^−1.
        DipLengths = []
        for p in range(0, Npigs):
            pig = SelPigs[p]
            dip = np.zeros((3,))
            for n in range(0, len(ListAtoms[p])):
                dip += eo*ListQ10[p][n]*CoordList[p][n]
                
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
                        Rmn = CoordList[p1][atm,:] - CoordList[p2][atn,:]
                        rmn = np.linalg.norm(Rmn)
                        CoupMat[p1,p2] += eo*eo*df1*df2*(ListQ10[p1][atm]*ListQ10[p2][atn]/rmn)*(Erg2J/h)/c
                        
            
        CoupMat = CoupMat + np.transpose(CoupMat)
        
        # Calculate rotation matrix
        for m in range(0, Npigs):
            for n in range(0, Npigs):
                Rmn = CentMat[n,:] - CentMat[m,:]
                #osc = SelPigs[m].species.diplength*SelPigs[n].species.diplength
                osc = 1.0
                RotMat[m,n] = np.dot(Rmn, np.cross(Dips[m,:], Dips[n,:]))*osc
        
    # NB: oscillator strengths are all finally scaled to unity!
    return CoupMat, Dips, RotMat

# Picks the residue to which the atom closest to MG belongs
def locate_ligand(pig, instruc):
    # Only do this if an MG atom is listed in the structure
    if pig.atnames.count('MG')!=0:
        MGndx = pig.atnames.index('MG')
        data_frame = instruc.to_dataframe()
        bgndcs = (data_frame.resid != pig.residue.idx)
        dvec = np.sqrt(np.sum(np.power(pig.atcoords[MGndx] - instruc.coordinates[bgndcs],2),1))
        ndx0 = np.argmin(dvec)
        at0 = instruc[bgndcs].atoms[ndx0]
        return at0.residue
    else:
        return []


def calculate_shift(PigList, ChainList, instruc):
    
    # First create list of chain-selected pigments
    SelPigs = []
    for pig in PigList:
        if ChainList.count(pig.residue.chain)>0:
            SelPigs.append(pig)
    
    Npigs = len(SelPigs)
    tFreqs = np.zeros((Npigs,))
        
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
        # Build list of TrESP coordinates for each pigment
        # NB: Units are converted to cm! (cgs units)
        CoordList = []
        for p in range(0, Npigs):
            pig = SelPigs[p]
            coords = []
            # Find index associated with each TrESP atom in the structure.
            # Coordinate list will be sorted in atom order specified in TrESP
            # parameter file *not* order of crystal structure. 
            for name in ListAtoms[p]:
                ndx = pig.atnames.index(name)
                coords.append(pig.atcoords[ndx])
            CoordList.append(np.array(coords)*ang2cm)
            
        
        for p in range(0, Npigs):
            
            pig = SelPigs[p]
            shift = 0.0
            
            #lig = locate_ligand(pig, instruc)
            data_frame = instruc.to_dataframe()
            
            # We exclude the pigment itself from calculation
            #ndcs = np.logical_and(data_frame.resid != pig.residue.idx, data_frame.resid != lig.idx)
            ndcs = data_frame.resid != pig.residue.idx
            
            for atm in range(0, len(ListAtoms[p])):
                dQ = ListQ11[p][atm] - ListQ00[p][atm]
                Rvec = CoordList[p][atm,:] - instruc.coordinates[ndcs,:]*ang2cm
                Dvec = np.sqrt(np.sum(np.power(Rvec,2),1))
                Qvec = data_frame[ndcs].charge
                shift += eo*eo*dQ*np.sum(np.divide(Qvec,Dvec))*Erg2J/(h*c*eps_eff)
            tFreqs[p] = shift
        
    return tFreqs, error
