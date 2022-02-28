import autode as ade
import numpy as np
import os

Cores = 48
ade.Config.n_cores = Cores

orca = ade.methods.ORCA()
xtb = ade.methods.XTB()

MoI1 = ade.Molecule('butane.xyz', solvent_name=None, charge=0)
MoI2 = ade.Molecule('pentane.xyz', solvent_name=None, charge=0)
R1 = ade.Molecule('p-fluorophenol',smiles='C1=CC(=CC=C1O)F', solvent_name='acetonitrile')
P1 = ade.Molecule('p-fluorophenolate',smiles='C1=CC(=CC=C1[O-])F', solvent_name='acetonitrile', charge=-1)

print(MoI1, flush=True)
print(MoI2, flush=True)
print(R1, flush=True)
print(P1, flush=True)

MoI1.optimise(method=xtb, n_cores=Cores)
MoI2.optimise(method=xtb, n_cores=Cores)
R1.optimise(method=xtb, n_cores=Cores)
P1.optimise(method=xtb, n_cores=Cores)

print('molecules have been optimised with XTB', flush=True)

MoI1.optimise(method=orca, keywords=ade.OptKeywords(['HF-3c']), n_cores=Cores)
MoI2.optimise(method=orca, keywords=ade.OptKeywords(['HF-3c']), n_cores=Cores)
R1.optimise(method=orca, keywords=ade.OptKeywords(['HF-3c']), n_cores=Cores)
P1.optimise(method=orca, keywords=ade.OptKeywords(['HF-3c']),n_cores=Cores)

print('molecules have been optimised with ORCA', flush=True)

CoI1 = ade.Calculation(name=MoI1.name, molecule=MoI1, method=orca, keywords=ade.HessianKeywords(['PBEPBE', 'def2SVP', 'Freq']), n_cores=Cores)
CoI2 = ade.Calculation(name=MoI2.name, molecule=MoI2, method=orca, keywords=ade.HessianKeywords(['PBEPBE', 'def2SVP', 'Freq']), n_cores=Cores)
CoIR1 = ade.Calculation(name=R1.name, molecule=R1, method=orca, keywords=ade.HessianKeywords(['PBEPBE', 'def2SVP', 'Freq']), n_cores=Cores)
CoIP1 = ade.Calculation(name=P1.name, molecule=P1, method=orca, keywords=ade.HessianKeywords(['PBEPBE', 'def2SVP', 'Freq']), n_cores=Cores)

print('Calculations have been carried out', flush=True)

CoI1.output.filename = MoI1.name+'_hess_orca.hess'
CoI2.output.filename = MoI2.name+'_hess_orca.hess'
CoIR1.output.filename = R1.name+'_hess_orca.hess'
CoIP1.output.filename = P1.name+'_hess_orca.hess'

print('calculations have been output to files', flush=True)

MoI1.calc_thermo(calc=CoI1, keywords=['HF-3c'], temp=298.15, ss='1atm', sn=1, n_cores=Cores)
MoI2.calc_thermo(calc=CoI2, keywords=['HF-3c'], temp=298.15, ss='1atm', sn=1, n_cores=Cores)
R1.calc_thermo(calc=CoIR1, keywords=['HF-3c'], temp=298.15, ss='1atm', sn=1, n_cores=Cores)
P1.calc_thermo(calc=CoIP1, keywords=['HF-3c'], temp=298.15, ss='1atm', sn=1, n_cores=Cores)

print('thermodynamic calculations have been carried out', flush=True)

GibbsE1 = MoI1.g_cont
GibbsE2 = MoI2.g_cont
GibbsR1 = R1.g_cont
GibbsP1 = P1.g_cont

print('gibbs free energies acquired', flush=True)

print(f'gibbsE {MoI1.name} = {GibbsE1} Ha', flush=True)
print(f'gibbsE {MoI2.name} = {GibbsE2} Ha', flush=True)
print(f'gibbsE {R1.name} = {GibbsR1} Ha', flush=True)
print(f'gibbsE {P1.name} = {GibbsP1} Ha', flush=True)

print('gibbs free energies printed', flush=True)

DelGibbs = (GibbsR1 + GibbsE1) - (GibbsE2 + GibbsP1)

print('Gibbs Delta has been Acquired', flush=True)

print(f'delG = {DelGibbs} Ha for {MoI1.name} reaction', flush=True)

NameofFile = f'{MoI1.name}.txt'
String = f'delG = {DelGibbs} Ha for {MoI1.name}'

with open(NameofFile, 'w') as f:
    f.write(String)
    f.close()