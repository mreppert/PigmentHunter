import numpy as np
import parmed as pmd

from scipy import spatial

import os,sys
sys.path.append('./main/')

import pigment


def nsd_transformation(refxyz, coordinates):
        
    # non weighted center of mass
    rcm = np.sum(coordinates,0)/np.shape(refxyz)[0]
    
    # shift
    xyz = coordinates - rcm
    
    # Each elements of U should contain a sum over atom indices. 
    # a and b iterate over coordinates (x,y,z).
    # The sum over atoms is accomplished using np.dot(). 
    U = np.zeros((3,3))
    for a in range(0,3):
        for b in range(0,3):
            U[a,b] = np.dot(xyz[:,a], xyz[:,b])
    
    # Passing two output arguments means that the first will be eigenvalues and the second eigenvectors
    evals, evecs = np.linalg.eig(U)
    
    # This orders all eigenvalues from smallest to largest. 
    # This isn't really needed here, but it's often useful. 
    # Note that the nth *column* (not row) of evecs is the 
    # eigenvector corresponding to the nth eigenvalue. 
    idx = evals.argsort()
    sevals = evals[idx]
    sevecs = evecs[:,idx]
    m = sevecs[:,0]
    if np.dot(np.cross(xyz[1,:], xyz[0,:]), m)<0:
        m = -m

        
    ##############
    # Rotation 1 #
    ##############
    
    # Rotation axis unit vector
    u = np.cross(m, np.array([0,0,1]))/np.linalg.norm(np.cross(m, np.array([0,0,1])))
    
    # angle
    b = np.arccos(m[2])
    
    # Transformation object
    R = spatial.transform.Rotation.from_rotvec(u*(b))
    newxyz = R.apply(xyz)
    
    
    ##############
    # Rotation 2 #
    ##############
    
    # Calculate theta
    top = (np.dot(refxyz[:,0], newxyz[:,1])-np.dot(refxyz[:,1], newxyz[:,0]))
    bot = (np.dot(refxyz[:,0], newxyz[:,0])+np.dot(refxyz[:,1], newxyz[:,1]))
    theta0 = np.arctan(top/bot)
    
    # Keep between 0 and 2*pi
    if theta0 < 0:
        theta1 = theta0 + np.pi
        theta2 = theta0 + 2.0*np.pi
    else:
        theta1 = theta0
        theta2 = theta0 + np.pi
    
    # Now distinguish maximum vs. minimum. (Shifted by pi from each other.)
    # Note typo in paper!
    discrim1 = (np.dot(refxyz[:,0], newxyz[:,1])-np.dot(refxyz[:,1], newxyz[:,0]))*np.sin(theta1) + (np.dot(refxyz[:,0], newxyz[:,0])+np.dot(refxyz[:,1], newxyz[:,1]))*np.cos(theta1)
    discrim2 = (np.dot(refxyz[:,0], newxyz[:,1])-np.dot(refxyz[:,1], newxyz[:,0]))*np.sin(theta2) + (np.dot(refxyz[:,0], newxyz[:,0])+np.dot(refxyz[:,1], newxyz[:,1]))*np.cos(theta2)
    
    # b1 is the rotation angle
    if discrim1>0:
        b1 = theta1
    elif discrim2>0:
        b1 = theta2
    
    # check if values make sense
    if (discrim1>0) and (discrim2)>0:
        print('Error! Both rotation angles were positive!')
    if (discrim1<0) and (discrim2)<0:
        print('Error! Neither rotation angle was positive!')
        
    # Rotation axis unit vector
    u1 = np.array([0,0,1])
    
    # Transformation object
    R1 = spatial.transform.Rotation.from_rotvec(u1*(-b1))
    newnewxyz = R1.apply(newxyz)
    
    return newnewxyz-refxyz




# Dobs is the difference between refxyz coordinates and pigment coordinates
# after rotation into the NSD reference frame. 
def nsd_calc(Dobs, Dgamma):
    
    Nrows = np.shape(Dgamma)[0]
    Ncols = np.shape(Dgamma)[1]
    
    # Deformation 
    D0 = []
    OOP = []
    
    # G loops through symmetry groups. 
    # The first few symmetry groups have three independent vectors
    for G in range(0, Ncols-2, 3):
        Dvecs = Dgamma[:,G:G+3]
        
        # Minimal basis calculation
        P = np.zeros((3,))
        for j in range(0,3):
            P[j] = np.dot(Dvecs[:,j], Dobs[:,2])
        D0.append(P[0])
        
        # Full basis calculation
        B = np.zeros((3,3))
        for j in range(0,3):
            for m in range(0,3):
                B[j,m] = np.dot(Dvecs[:,j], Dvecs[:,m])
        
        Binv = np.linalg.inv(B)
        dvec = Binv@P
        dtot = np.sqrt(dvec.transpose()@B@dvec)
        OOP.append(dtot)
        
    # The last symmetry group has only two vectors
    # NB: Note hard-coded column numbers. :(
    Dvecs = Dgamma[:,15:17]
    
    # Mimimal basis
    P = np.zeros((2,))
    for j in range(0,2):
        P[j] = np.dot(Dvecs[:,j], Dobs[:,2])
    D0.append(P[0])
    
    # Full basis
    B = np.zeros((2,2))
    for j in range(0,2):
        for m in range(0,2):
            B[j,m] = np.dot(Dvecs[:,j], Dvecs[:,m])
            
    Binv = np.linalg.inv(B)
    dvec = Binv@P
    dtot = np.sqrt(dvec.transpose()@B@dvec)
    OOP.append(dtot)
    
    return np.array(D0)


def siteenergy(d0, pars):
    
    # other constants
    h = 6.626e-37 #kJ*s
    c = 2.998e10 #cm/sec

    #out of plane displacements from table
    ED = 0

    #force constant (1/((cm*A^2))
    k = pars[0:6]
    
    # H --> L gap
    a = pars[6]
    
    # H-1 --> L+1 gap
    b = pars[7]
    
    # Coupling
    C = pars[8]
    
    for n in range(0,6):
        ED += k[n] * d0[n]**2
        
    X1 = a + k[3]*(d0[3]**2) - k[5]*(d0[5]**2) # frequency of H to L transition (cm-1)
    X2 = b + k[4]*(d0[4]**2) - k[2]*(d0[2]**2) # frequency of H-1 to L+1 transition (cm-1)
    E = 0.5*( (X1+X2) - np.sqrt( (X2-X1)**2 + (4*(C**2)) )) # site energy of deformation energy in kJ/mol 
    
    return E
        
def calculate_shift(PigList, ChainList, instruc):
    
    # First create list of chain-selected pigments
    SelPigs = []
    for pig in PigList:
        if ChainList.count(pig.residue.chain)>0:
            SelPigs.append(pig)
    
    if len(SelPigs)==0:
        return [], True, 'No pigments in selection.'

    
    Npigs = len(SelPigs)
    FreqTraj = []
    D0Traj = []
    DOOPTraj = []
    
    h = 6.62607015e-34 #J*s
    c = 2.998e10 #cm/s
    eo = 4.80320451e-10 # esu
    eps_eff = 2.5
    
    Erg2J = 1.0e-7
    
    # Now check whether NSD parameters are available for *all*
    # selected pigments and (if so) store them for reference. 
    
    # Loop over pigments
    error = False
    msg = ''
    ListNames = []
    ListPar = []
    ListRefXYZ = []
    ListDGamma = []
    for p in range(0, Npigs):
        pig = SelPigs[p]
        
        # NSD parameters should be in two separate files:
        #    'misc/NSD/' + pig.species.stdname + '/names.txt' stores atom names
        #         -- the first column is the PDB standard, second is NSD convention
        #    'misc/NSD/' + pig.species.stdname + '/dgamma.txt' stores the DGamma matrix
        #    'misc/NSD/' + pig.species.stdname + '/xyz.txt' stores xyz reference structure
        #    'misc/NSD/' + pig.species.stdname + '/par.txt' stores mixing parameters
        #         -- Columns 0 - 5 are force constants. Columns 6 - 8 are H --> L gap, H-1 --> L+1 gap, and coupling term. 
        fpar = 'misc/NSD/' + pig.species.stdname + '/par.txt'
        fxyz = 'misc/NSD/' + pig.species.stdname + '/xyz.txt'
        fdgamma = 'misc/NSD/' + pig.species.stdname + '/dgamma.txt'
        fnsdnms = 'misc/NSD/' + pig.species.stdname + '/names.txt'
        
        flist = [fnsdnms, fdgamma, fxyz, fpar]
        for fnm in flist:
            if os.path.isfile(fnm)==False:
                msg = 'Error: Could not locate NSD file ' + fnm + '. Aborting NSD calculation.'
                error = True
                break
                
        # If we found all the files, load the parameters and store in the appropriate List
        if not error:
            
            ListPar.append(np.loadtxt(fpar))
            ListRefXYZ.append(np.loadtxt(fxyz))
            ListDGamma.append(np.loadtxt(fdgamma)*(1.0e-4))
            
            atnms = []
            with open(fnsdnms) as fd:
                for line in fd:
                    if len(line)>0 and line[0]!="#":
                        lst = line.split()
                        if len(lst)<2:
                            msg = 'Error reading NSD input file ' + fnsdnms + '. Aborting NSD calculation.'
                            error = True
                            break
                        else:
                            atnms.append([lst[0], lst[1]])
                            
            # Check that new XYZ and DGamma additions are of the same length as new name list.
            # Number of rows in DGamma and RefXYZ should be the same as len(atnms). 
            natoms = len(atnms)
            if natoms==np.shape(ListRefXYZ[p])[0] and natoms==np.shape(ListDGamma[p])[0]:
                
                # Note that this is a list of tuple-lists, not a NumPy array
                ListNames.append(atnms)
            else:
                error = True
                msg = 'Error: NSD name file ' + fnsdnms + ' appears to reference a different number (' + str(natoms) + ') than the corresponding DGamma (' + str(np.shape(ListDGamma[p])[0]) + ') and/or refxyz file (' + str(np.shape(ListRefXYZ[p])[0]) + '). Aborting NSD calculation.'
                
        # If there's an error, we break out of the loop over pigments. 
        # If not, proceed to the next pigment
        if(error):
            break
                            
    # If no errors, we've located all relevant parameter files and imported the data. 
    # Now check whether all necessary atoms are available in the structure.
    if error==False:
        for p in range(0, Npigs):
            
            pig = SelPigs[p]
            
            # Loop over atom names. PDB std names are in first column, NSD names in second. 
            for n in range(0, len(ListNames[p])):

                # PDB name for the nth NSD atom in pth pigment
                name = ListNames[p][n][0]

                # If we can't find it in the structure, throw an error. 
                # pig.atnames stores the names of all atoms associated with pigment pig
                if pig.atnames.count(name)==0:
                    msg = 'Error: Could not locate NSD atom ' + name + " in pigment " + pig.residue.name + " " + pig.residue.chain + " " + str(pig.residue.number) + '. Aborting NSD calculation.'
                    error = True
                    break

            # If there's already been an error, don't continue. 
            if error:
                break
                
    # If no errors, all pigments have NSD parameters and necessary atoms. 
    if error==False:

        # Number of frames in trajectory
        Nframes = np.shape(SelPigs[0].atcoords)[0]

        # Build list of NSD coordinates for each pigment
        # NB: Units are converted to cm! (cgs units)
        CoordList = []
        for p in range(0, Npigs):
            pig = SelPigs[p]
            coords = np.zeros((Nframes, len(ListNames[p]), 3))
            
            # Find index associated with each NSD atom in the structure.
            # Coordinate list will be sorted in atom order specified in NSD
            # name file *not* order of crystal structure. 
            for n in range(0, len(ListNames[p])):
                name = ListNames[p][n][0]
                ndx = pig.atnames.index(name)
                coords[:,n,:] = pig.atcoords[:,ndx,:].copy()
            CoordList.append(coords)
            
        # Coordinates for the nth NSD atom in fth frame are in CoordList[f,n,:]
        for fr in range(0, Nframes):

            tFreqs = np.zeros((Npigs,))
            tD0 = np.zeros((Npigs,6))
            tDOOP = np.zeros((Npigs,6))

            for p in range(0, Npigs):
                
                pig = SelPigs[p]
                shift = 0.0
                
                # These are coordinates for the NSD atoms in the current frame
                XYZ = CoordList[p][fr,:,:]
                
                # Transform to reference basis
                Dobs = nsd_transformation(ListRefXYZ[p], XYZ)
                
                # Run NSD calculation
                d0, doop = nsd_calc(Dobs, ListDGamma[p])
                tD0[p,:] = d0
                tDOOP[p,:] = doop
                
                # Calculate site energy 
                tFreqs[p] = siteenergy(d0, ListPar[p])
                
            FreqTraj.append(tFreqs)
            D0Traj.append(tD0)
            DOOPTraj.append(tDOOP)
            
    return FreqTraj, D0Traj, DOOPTraj, error, msg
