# class category is used to classify input molecules into broad groups. 
# Each instance of a category includes two traits:
#  * category.name: String used to identify the class
#
#  * category.BaseNames: List of atom names (strings) that *must*
#       be identified in the input molecule if it is to be considered a 
#       member of this category
#
# atom names that *must* be present for a molecule to be considered 
class category:
    def __init__(self, name, BaseNames):
        self.name = name
        self.BaseNames = BaseNames

# class species define a specific chemical identities for different 
# pigments, e.g., chlorophyll b or pheophytin a. The class traits include:
# 
# * pclass: The pigment.category to which the species belongs
# 
# * name: The (not necessarily unique) residue name commonly used 
#      to refer to the pigment, e.g., "Chl b" for Chlorophyll a  or "Pheo a" for
#      pheophytin a. 
# 
# * stdname: A unique, standardized residue name abbreviation used 
#      for internal referencing within PigmentHunter. E.g., CLB for Chlorophyll b
#      or PHA for pheophytin a. 
# 
# * xnames: A list of atom names used to define the pigment species. This list
#      should not include atoms on the BaseNames list used to define the
#      pigment.category, although different species may hold xnames in common. 
#
class species:
    def __init__(self, ResName, StdName, pClass, XAtNames, diplength):
        self.pclass = pClass,
        self.name = ResName
        self.stdname = StdName
        self.xnames = XAtNames
        self.diplength = diplength
        
        
# Instances of class pigment store data for a particular pigment molecule 
# described by imported structural data (e.g., from a PDB file). Its
# traits include:
# * idx: The residue index for this molecule in the parmed.struc
#     from which it was extracted. Thus struc.residues[idx] returns
#     the corresponding parmed residue object. 
# 
# * species: the instance of the pigment.species class associated with this pigment
#     (e.g., representing CLA, CLB, PHA, etc.)
# 
# * alist: A list of "allowed" species that *could* be associated
#       with this pigment, based solely on the structural date present 
#       in the parmed.struc from which it was extracted. Many pigment types
#       cannot be identified uniquely from heavy atoms only (particularly since
#       the structure could be incomplete); alist enumerates all recognized 
#       species consistent with the available structural data. 
#
# * residue: the parmed residue object corresponding to this pigment
# 
# * atnames: a list of atom name strings for each atom in the 
#       corresponding parmed residue
# 
# * atcoords: a coordinates matrix for the atoms named in atnames
#
class pigment:
    def __init__(self, idx, spec, alist, residue, atnames, atcoords):
        self.idx = idx
        self.species = spec
        self.alist = alist
        self.atnames = atnames
        self.atcoords = atcoords
        self.widget = []
        self.residue = residue

# Returns a list of all residue indices (in struc)
# that contain all required atoms
def find_pigments(NameList, instruc):

    SelStr = "(:0-"+str(len(instruc.residues)) + ")"
    for ndx in range(1, len(NameList)):
        tstruc = instruc.view[SelStr + ' & (@'+NameList[ndx] + ")"]
        SelNdx = []
        for at in tstruc:
            SelNdx.append(at.residue.idx)
        SelStr = "(:" + ",".join([str(i+1) for i in SelNdx]) + ")"
        if len(SelNdx)==0:
            break

    return SelNdx
    
# Now check which pigment types can be *excluded* based on the atoms
# that are present in each pigment. NB: We do *not* require all atoms
# to be present. 
# The return argument AList is a list of pigment types that *are* allowed. 
def eliminate_types(PigNdcs, PClass, TypeList, UNK, instruc):
    
    AList = []
    for ndx in PigNdcs:
        tstruc = instruc[':'+str(ndx+1)]
        
        # For each atom in the pigment, check whether it's allowed
        # under each pigment type. 
        TypeAllowed = [1 for i in range(0,len(TypeList))]
        for at in tstruc:
            # If not a Hydrogen...we always ignore hydrogens
            if at.name[0]!='H':
                # First check if it's a base atom. If not, we check 
                # the characteristic atoms for each pigment. 
                if (PClass.BaseNames.count(at.name)==0):
                    for tndx in range(0, len(TypeList)):
                        if TypeList[tndx].xnames.count(at.name)==0:
                            TypeAllowed[tndx] = 0
                            #print('Could not match atom ' + at.name + ' to a known ' + TypeList[tndx].stdname + ' atom')
                            
        alist = []
        for tndx in range(0, len(TypeList)):
            if TypeAllowed[tndx]:
                alist.append(TypeList[tndx])
                
        if len(alist)==0:
            alist = [UNK]
        
        AList.append(alist)
        
    return AList


# MList[p] is a list of matched types for pigment p. Here a "match"
# means that *all* the species.xnames atoms are present
# in the pigments parmed.residue.
def match_types(PigNdcs, AList, instruc):
    
    MList = []
    # Loop through pigments
    for p in range(0, len(PigNdcs)):
        ndx = PigNdcs[p]
        
        tstruc = instruc[':'+str(ndx+1)]
        
        # alist is a list of pigment types allowed for pigment p
        alist = AList[p]
        
        # mlist will store a list of types that match the pigment exactly
        mlist = []
        
        # Loop through possible types
        for t in range(0, len(alist)):
            
            # Start assuming the pigment matches the type
            match = True
            typ = alist[t]
            for name in typ.xnames:
                
                # If any atom in type.xnames is NOT found in the corresponding
                # instruc.residue set match to False and break
                if len(tstruc['@'+name].atoms)==0:
                    match = False
                    break
                
            # If we have a match, append it to the pigment's match list
            if match:
                mlist.append(typ)
                
        MList.append(mlist)
    return MList

# 1. Assign all pigments to the top entry in MList[p] if non-empty
# 2. For each pigment still unassigned, check whether any pigments with the 
#    same residue.name have already been assigned. If so, assign to the same type. 
# 3. For each pigment still unassigned, assign to the top entry in AList[p] if non-empty
def assign_pigments(PigNdcs, MList, AList, instruc):
    TList = []
    
    # 1. Assign by MList
    for p in range(0, len(PigNdcs)):
        if len(MList[p])>0:
            TList.append(MList[p][0])
        else:
            TList.append([])
            
    # 2. Assign by name
    for p in range(0, len(PigNdcs)):
        if TList[p]==[]:
            rname = instruc.residues[PigNdcs[p]].name
            
            # Check if any assigned pigment has a matching name 
            for pref in range(0, len(PigNdcs)):
                if (TList[pref]!=[]) and (rname==instruc.residues[PigNdcs[pref]].name):
                    TList[p] = TList[pref]
                    break
                    
    # 3. If pigments are still unmatched, assign to top alist entry
    for p in range(0, len(PigNdcs)):
        if TList[p]==[]:
            TList[p] = AList[p][0]
                
    return TList