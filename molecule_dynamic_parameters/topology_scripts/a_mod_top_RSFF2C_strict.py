"""
2021/01/27  fix the bug of N-term amino acids and C-term amino acids by CJN
2021/02/02  ignore residues around unnatural residues
2021/03/09  fix the bug of ignoring the residues near ACE and NME
"""
from sys import argv
from subprocess import Popen, PIPE

topBaseName = argv[1]
CMAPfileName = 'RSFF2C_CMAP_final.dat' # which is 'CMAP_decrease_a_L.dat'

dih_E_chi2 = '-0.457   0.453  -0.105   0.269'
dih_E_chi3 = ' 0.028   0.201   0.052   0.025'
dih_Q_chi2 = ' 0.603   0.672   0.341   0.311'
dih_Q_chi3 = '-0.187   0.012  -0.440   0.032'
dih_N_chi2 = '-0.343   0.484  -0.154  -0.108'

AAs_0 = ('GLY', 'ALA')
AAs_1 = ('PRO', 'GLU', 'GLN', 'LYS', 'ARG', 'MET',
         'LEU', 'PHE', 'TYR', 'TRP', 'CYS', 'SER',
         'ASP', 'ASN', 'VAL', 'ILE', 'THR', 'HID', 'HIE', 'HIP' )

Atoms = ('N', 'CA', 'C', 'CB')

top_Title = """
%FLAG FORCE_FIELD_TYPE
%FORMAT(i2,a78)
 1 Residue-Specific Force Field modified from ff14SB using CHARMM-like CMAP  
%COMMENT RSFF2C by JIANG Fan and KANG Wei
"""

top_Add = """
%FLAG CHARMM_UREY_BRADLEY_COUNT
%COMMENT  Number of Urey Bradley terms and types
%FORMAT(2i8)               
     3      1
%FLAG CHARMM_UREY_BRADLEY
%COMMENT  in each UB term: i,k,index
%FORMAT(10i8)              
       1       5       1       2       5       1       3       5       1       
%FLAG CHARMM_UREY_BRADLEY_FORCE_CONSTANT
%FORMAT(5e16.8)            
  0.00000000E+02
%FLAG CHARMM_UREY_BRADLEY_EQUIL_VALUE
%FORMAT(5e16.8)            
  1.00000000E+00
%FLAG CHARMM_NUM_IMPROPERS
%FORMAT(10i8)              
      1
%FLAG CHARMM_IMPROPERS
%COMMENT  i,j,k,l,index  i,j,k,l,index
%FORMAT(10i8)              
      1       2       3       5       1      
%FLAG CHARMM_NUM_IMPR_TYPES
%FORMAT(i8)                
      1
%FLAG CHARMM_IMPROPER_FORCE_CONSTANT
%FORMAT(5e16.8)            
  0.00000000E+03  
%FLAG CHARMM_IMPROPER_PHASE
%FORMAT(5e16.8)            
  0.00000000E+00
"""

def split( line, width ) :
    r = [ ]
    i = 0
    while i < len(line) :
        r.append( line[i:i+width].strip() )
        i += width
    if r[-1] == '' :
        return r[:-1]
    else :
        return r

def searchResidue( n_atom ) :
    for residue in residues :
        if residue.initAtomNum <= n_atom <= residue.lastAtomNum :
            return residue

def addTorsionCMD(atoms_involved, torsion_param):
    coeff = list(map(float, torsion_param.split()))
    coeffStr = list(map(str, list(map(abs, coeff))))    # ks
    phaseStr = ['180' if x > 0 else '0' for x in coeff]
    print((coeffStr, phaseStr))
    cmd = ''
    for i, k in enumerate(coeffStr):
        cmd += 'addDihedral ' + atoms_involved + k + ' ' \
             + str(i+1) + ' ' + phaseStr[i] + ' 1.0 2.0 type normal\n'
    return cmd

class Residue :
    def __init__(self, name, num) :
        self.name = name
        self.num = num      # residue number starting from 1
        self.N3 = 0         # The number of N3 in the residue
        self.Nterm = False
        self.O2 = 0
        self.Cterm = False
    def output(self) :
        s = self.name + ' ' \
             + str(self.num) + ' ' \
             + str(self.initAtomNum) + '-' \
             + str(self.lastAtomNum)
        for attr in ('N', 'CA', 'C', 'CB', 'G') :
            if hasattr(self, attr) :
                s += ' ' + attr + ':' + str( getattr(self, attr) )
        return s
    def isGamma(self, atom) :
        if self.name == 'VAL' :
            if atom == 'CG2' :
                return True
        elif atom in ( 'CG', 'OG', 'SG', 'CG1', 'OG1' ) :
            return True

        return False

add_torsions_cmd  = addTorsionCMD(':GLU@CA :GLU@CB :GLU@CG :GLU@CD  ', dih_E_chi2 ) \
                  + addTorsionCMD(':GLU@CB :GLU@CG :GLU@CD :GLU@OE1 ', dih_E_chi3 ) \
                  + addTorsionCMD(':GLU@CB :GLU@CG :GLU@CD :GLU@OE2 ', dih_E_chi3 ) \
                  + addTorsionCMD(':GLN@CA :GLN@CB :GLN@CG :GLN@CD  ', dih_Q_chi2 ) \
                  + addTorsionCMD(':GLN@CB :GLN@CG :GLN@CD :GLN@OE1 ', dih_Q_chi3 ) \
                  + addTorsionCMD(':ASN@CA :ASN@CB :ASN@CG :ASN@OD1 ', dih_N_chi2 ) \
                  + 'parmout ' + topBaseName + '.sd.top\ngo\n'
print(add_torsions_cmd)

p = Popen('parmed -p ' + topBaseName + '.hmr.top -s -O --logfile parmed.log\n',
          stdin = PIPE, stdout= PIPE, shell = True)
p.communicate(add_torsions_cmd.encode())

with open( topBaseName + '.sd.top', 'r' ) as ifile :
    top_lines = ifile.readlines()

read = False
for line in top_lines :
    if '%FLAG POINTERS' in line :
        read = True
    if read and '%' not in line :
        numOfAtoms = int( line.split()[0] )
        break
    
print(('Total number of atoms:', numOfAtoms))

residues = [ ]
read = False
resNum = 0
for line in top_lines :
    if read and '%FLAG' in line :
        break
    if '%FLAG RESIDUE_LABEL' in line :
        read = True
    if read and '%' not in line :
        for each in split( line, 4 ) :
            resNum += 1
            residues.append( Residue( each, resNum ) )

residues = residues[:min(len([each for each in residues if each.name not in ['WAT', 'Na+', 'Cl-']])+1, len(residues))]
print(('Number of residues:', len(residues)))

read = False
i_res = 0
for line in top_lines :
    if read and '%FLAG' in line :
        break
    if '%FLAG RESIDUE_POINTER' in line :
        read = True
    if read and '%' not in line :
        for each in split( line, 8 ) :
            if i_res > len(residues) :
                break
            if i_res < len(residues) :
                residues[i_res].initAtomNum = int(each)
            if i_res > 0 :
                residues[i_res-1].lastAtomNum = int(each) - 1
            i_res += 1

residues[-1].lastAtomNum = numOfAtoms

read = False
n_atom = 1
for line in top_lines :
    if read and '%FLAG' in line :
        break
    if '%FLAG ATOM_NAME' in line :
        read = True
    if read and '%' not in line :
        for each in split( line, 4 ) :
            residue = searchResidue( n_atom )
            if each in Atoms :
                #print n_atom, residue.name, each
                setattr( residue, each, n_atom )
            elif residue.isGamma( each ) :
                setattr( residue, 'G', n_atom )
            n_atom += 1

for residue in residues :
    print((residue.output()))

# count N3 and O2
read = False
n_atom = 1
for line in top_lines :
    if read and '%FLAG' in line :
        break
    if '%FLAG AMBER_ATOM_TYPE' in line :
        read = True
    if read and '%' not in line :
        for each in split( line, 4 ) :
            residue = searchResidue( n_atom )
            if residue.name in AAs_0 or residue.name in AAs_1: 
                if each == "N3" :
                    residue.N3 += 1
                    if residue.name != "LYS" and residue.N3 == 1:
                        residue.Nterm = True
                    elif residue.name == "LYS" and residue.N3 == 2:
                        residue.Nterm = True
                    elif residue.name == "LYS" and residue.N3 == 1:
                        pass
                    else:
                        print((residue.name, residue.N3))
                        raise
                if each == "O2":
                    residue.O2 += 1
                    if residue.name not in ["ASP", "GLU"] and residue.O2 == 2:
                        residue.Cterm = True
                    elif residue.name in ["ASP", "GLU"] and residue.O2 == 4:
                        residue.Cterm = True
                    elif residue.name not in ["ASP", "GLU"] and residue.O2 < 2:
                        pass
                    elif residue.name in ["ASP", "GLU"] and residue.O2 < 4:
                        pass
                    else:
                        print((residue.num, residue.name, residue.O2, n_atom))
                        raise
            n_atom += 1

read = False
for line in top_lines :
    if read and '%FLAG' in line :
        read = False
    
    if '%FLAG LENNARD_JONES_ACOEF' in line :
        read = True
        top_Add += '%FLAG LENNARD_JONES_14_ACOEF\n'
    
    if '%FLAG LENNARD_JONES_BCOEF' in line :
        read = True
        top_Add += '%FLAG LENNARD_JONES_14_BCOEF\n'
    
    if read and '%FLAG' not in line :
        top_Add += line


CMAP_terms = [ ]
for i_res in range(1, len(residues)-1) :
    name = residues[i_res].name
    if residues[i_res-1].name not in AAs_0+AAs_1+("CYX","ACE","NME") or residues[i_res].name not in AAs_0+AAs_1 or residues[i_res+1].name not in AAs_0+AAs_1+("CYX","ACE","NME"):
        print(("ignore residue:", i_res+1, name))
        continue
    if residues[i_res].Nterm or residues[i_res].Cterm:
        print(("terminal residue:", i_res+1, name))
        continue
    print(("current residue:", i_res+1, name))
    if name in AAs_0 :
        i = residues[i_res-1].C
        j = residues[i_res].N
        k = residues[i_res].CA
        l = residues[i_res].C
        m = residues[i_res+1].N
        n = 1 + AAs_0.index( name )
        CMAP_terms.append( (i,j,k,l,m, n) )
    
    elif name in AAs_1 :
        i = residues[i_res-1].C
        j = residues[i_res].N
        k = residues[i_res].CA
        l = residues[i_res].C
        m = residues[i_res+1].N
        CMAP_terms.append( (i,j,k,l,m, 2) )

        n = 3 + AAs_1.index( name ) * 2
        
        i = residues[i_res].G
        j = residues[i_res].CB
        k = residues[i_res].CA
        l = residues[i_res].N
        m = residues[i_res-1].C
        CMAP_terms.append( (i,j,k,l,m, n) )

        l = residues[i_res].C
        m = residues[i_res+1].N
        CMAP_terms.append( (i,j,k,l,m, n+1) )
    else:
        print(name)
        
for each in CMAP_terms :
    print(each)


top_Add_CMAP = """
%FLAG CHARMM_CMAP_COUNT
%COMMENT  Total number of CMAP terms, number of unique CMAP types
%FORMAT(2I8)
""" + str(len(CMAP_terms)).rjust(8) + """      42 
%FLAG CHARMM_CMAP_RESOLUTION
%COMMENT  Resolution for each CMAP type
%FORMAT(20I4)
  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24
  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24  24
  24  24
"""

# load the CMAP parameters
with open( CMAPfileName, 'r' ) as ifile :
    for line in ifile :
        top_Add_CMAP += line

top_Add_CMAP += """%FLAG CHARMM_CMAP_INDEX
%COMMENT  Atom index i,j,k,l,m for CHARMM_CMAP_PARAMETER_n
%FORMAT(6I8)
"""

for each in CMAP_terms :
    top_Add_CMAP += "%8d%8d%8d%8d%8d%8d\n" %each

with open( topBaseName + '.run.top', 'w' ) as ofile :
 ofile.write( top_Title )
 ofile.writelines( top_lines )
 ofile.write( top_Add + top_Add_CMAP )

