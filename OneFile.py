import autode as ade
import numpy as np
orca = ade.methods.ORCA()
xtb = ade.methods.XTB()
import os



MoI1 = ade.Molecule('aBiCyDNH.xyz', solvent_name='acetonitrile', charge=0)
MoI2 = ade.Molecule('aBiCyDNH2.xyz', solvent_name='acetonitrile', charge=1)
R1 = ade.Molecule('p-fluorophenol',smiles='C1=CC(=CC=C1O)F')
P1 = ade.Molecule('p-fluorophenolate',smiles='C1=CC(=CC=C1[O-])F')

print(MoI1)
print(MoI2)

MoI1.optimise(method=xtb)
MoI2.optimise(method=xtb)
R1.optimise(method=xtb)
P1.optimise(method=xtb)

MoI1.optimise(method=orca)
MoI2.optimise(method=orca)
R1.optimise(method=orca)
P1.optimise(method=orca)

CoI1 = ade.Calculation(name=MoI1.name, molecule=MoI1, method=orca, keywords=orca.keywords.hess, n_cores=128)
CoI2 = ade.Calculation(name=MoI2.name, molecule=MoI2, method=orca, keywords=orca.keywords.hess, n_cores=128)
CoIR1 = ade.Calculation(name=R1.name, molecule=R1, method=orca, keywords=orca.keywords.hess, n_cores=128)
CoIP1 = ade.Calculation(name=P1.name, molecule=P1, method=orca, keywords=orca.keywords.hess, n_cores=128))

CoI1.output.filename = MoI1.name+'.out'
CoI2.output.filename = MoI2.name+'.out'
CoIR1.output.filename = R1.name+'.out'
CoIP1.output.filename = P1.name+'.out'

MoI1.calc_thermo(calc=CoI1, n_cores=128)
MoI2.calc_thermo(calc=CoI2, n_cores=128)
R1.calc_thermo(calc=CoIR1, n_cores=128)
P1.calc_thermo(calc=CoIP1, n_cores=128)

GibbsE1 = MoI1.free_energy
GibbsE2 = MoI2.free_energy
GibbsR1 = R1.free_energy
GibbsP1 = P1.free_energy

DelGibbs = (GibbsR1 + GibbsE1) - (GibbsE2 + GibbsP1)

print(f'delG = {DelGibbs} Ha for'+MoI1.name)

NameofFile = MoI1.name+".txt"

with open(NameofFile, 'w') as f
    f.write(f'delG = {DelGibbs} Ha for'+MoI1.name)
    f.close()