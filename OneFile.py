import autode as ade
import numpy as np
orca = ade.methods.ORCA()
xtb = ade.methods.XTB()
import os



MoI1 = ade.Molecule('aBiCyDNH.xyz', solvent_name='acetonitrile', charge=0)
MoI2 = ade.Molecule('aBiCyDNH2.xyz', solvent_name='acetonitrile', charge=1)

print(MoI1)
print(MoI2)

MoI1.optimise(method=xtb)
MoI2.optimise(method=xtb)

MoI1.optimise(method=orca)
MoI2.optimise(method=orca)

CoI1 = ade.Calculation(name=MoI1.name, molecule=MoI1, method=orca, keywords=orca.keywords.hess, n_cores=128)
CoI2 = ade.Calculation(name=MoI2.name, molecule=MoI2, method=orca, keywords=orca.keywords.hess, n_cores=128)

CoI1.output.filename = MoI1.name+'.out'
CoI2.output.filename = MoI2.name+'.out'

MoI1.calc_thermo(calc=CoI1, n_cores=128)
MoI2.calc_thermo(calc=CoI2, n_cores=128)

GibbsE1 = MoI1.free_energy
GibbsE2 = MoI2.free_energy

DelGibbs = GibbsE1 - GibbsE2

print(f'delG = {DelGibbs} Ha')