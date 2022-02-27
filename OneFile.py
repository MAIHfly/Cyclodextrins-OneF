import autode as ade
import numpy as np
orca = ade.methods.ORCA()
xtb = ade.methods.XTB()
import os
ade.Config.n_cores = 48
xtb.keywords.opt = ade.OptKeywords(['HF-3c'])
orca.keywords.opt = ade.OptKeywords(['HF-3c'])
orca.keywords.hess = ade.HessianKeywords(['HF-3c'])



MoI1 = ade.Molecule('butane.xyz', solvent_name=None, charge=0)
MoI2 = ade.Molecule('pentane.xyz', solvent_name=None, charge=0)
R1 = ade.Molecule('p-fluorophenol',smiles='C1=CC(=CC=C1O)F', solvent_name='acetonitrile')
P1 = ade.Molecule('p-fluorophenolate',smiles='C1=CC(=CC=C1[O-])F', solvent_name='acetonitrile', charge=-1)

print(MoI1)
print(MoI2)
print(R1)
print(P1)

MoI1.optimise(method=xtb, n_cores=48)
MoI2.optimise(method=xtb, n_cores=48)
R1.optimise(method=xtb, n_cores=48)
P1.optimise(method=xtb, n_cores=48)

print('molecules have been optimised with XTB')

MoI1.optimise(method=orca, n_cores=48)
MoI2.optimise(method=orca, n_cores=48)
R1.optimise(method=orca, n_cores=48)
P1.optimise(method=orca, n_cores=48)

print('molecules have been optimised with ORCA')

CoI1 = ade.Calculation(name=MoI1.name, molecule=MoI1, method=orca, keywords=orca.keywords.hess, n_cores=48)
CoI2 = ade.Calculation(name=MoI2.name, molecule=MoI2, method=orca, keywords=orca.keywords.hess, n_cores=48)
CoIR1 = ade.Calculation(name=R1.name, molecule=R1, method=orca, keywords=orca.keywords.hess, n_cores=48)
CoIP1 = ade.Calculation(name=P1.name, molecule=P1, method=orca, keywords=orca.keywords.hess, n_cores=48)

print('Calculations have been carried out')

CoI1.output.filename = MoI1.name+'_hess_orca.out'
CoI2.output.filename = MoI2.name+'_hess_orca.out'
CoIR1.output.filename = R1.name+'_hess_orca.out'
CoIP1.output.filename = P1.name+'_hess_orca.out'

print('calculations have been output to files')

MoI1.calc_thermo(calc=CoI1, n_cores=48)
MoI2.calc_thermo(calc=CoI2, n_cores=48)
R1.calc_thermo(calc=CoIR1, n_cores=48)
P1.calc_thermo(calc=CoIP1, n_cores=48)

print('thermodynamic calculations have been carried out')

GibbsE1 = MoI1.free_energy
GibbsE2 = MoI2.free_energy
GibbsR1 = R1.free_energy
GibbsP1 = P1.free_energy

print('gibbs free energies acquired')

print(f'gibbsE {MoI1.name} = {GibbsE1} Ha')
print(f'gibbsE {MoI2.name} = {GibbsE2} Ha')
print(f'gibbsE {R1.name} = {GibbsR1} Ha')
print(f'gibbsE {P1.name} = {GibbsP1} Ha')

print('gibbs free energies printed')

DelGibbs = (GibbsR1 + GibbsE1) - (GibbsE2 + GibbsP1)

print('Gibbs Delta has been Acquired')

print(f'delG = {DelGibbs} Ha for {MoI1.name} reaction')

NameofFile = f'{MoI1.name}.txt'
String = f'delG = {DelGibbs} Ha for {MoI1.name}'

with open(NameofFile, 'w') as f:
    f.write(String)
    f.close()