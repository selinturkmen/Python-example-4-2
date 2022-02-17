pdbfile =open("1ubq.pdb","r")
pdbfile.read()

def protein_parser(pdbfile):
    pdbfile =open("1ubq.pdb","r")
    atoms = []
    position = []
    
    for line in pdbfile:
        if line.startswith('ATOM'):
            line = line.split()
            element = line[11]
            dim = (line[6:9])
            dimne = []
            for x in dim:
                dimne.append(float(x))
            elements=[]
            elements.append(element)
            position = dimne + elements
            atoms.append(position)
  

          
    return atoms

                

def center_of_mass(atoms):

    mass={'C':12.01, 'O':16.00, 'H':1.008, 'N': 14.01, 'S':32.06}
    mass_cm = 0
    mass_z = 0
    mass_x = 0
    mass_y = 0
    for x in atoms:
        ri_x = x[0]
        ri_y = x[1]
        ri_z = x[2]
        mz = x[3]
        mi = mass[mz]
        mass_cm = mass_cm + mi
        mass_x = mass_x + ri_x * mi
        mass_y = mass_y + ri_y * mi
        mass_z = mass_z+ ri_z * mi
        cm_x = mass_x / mass_cm
        cm_y = mass_y / mass_cm
        cm_z = mass_z / mass_cm
        cm = [cm_x, cm_y, cm_z]
    return cm

def shift(atoms,vec):

    atomss = []
    for a in atoms:
        xdimension = a[0] + vec[0]
        ydimension = a[1] + vec[1]
        zdimension = a[2] + vec[2]
        isa = a[3]
        new = [xdimension,ydimension,zdimension,isa]
        atomss.append(new)
    return atomss


atoms=protein_parser('1ubq.pdb')

com=center_of_mass(atoms[:])
print(f"Center of mass of the protein is {com[0]:6.3f}, {com[1]:6.3f}, {com[2]:6.3f}")

for i in range(3):
    com[i] *= -1
atomsnew=shift(atoms,com)

com=center_of_mass(atomsnew[:])
print(f"New center of mass of the protein is {com[0]:6.3f}, {com[1]:6.3f}, {com[2]:6.3f}")

