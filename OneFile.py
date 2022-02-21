import autode as ade
import numpy as np
orca = ade.methods.ORCA()
xtb = ade.methods.XTB()
import os

# get pKa (there is an error here so result calculation will be done manually)
def GetPka(DeltGibbsE):
    DeltGibbsEkJpM = 2600*DeltGibbsE
    R = 0.008314
    T = input("What's the Temperature(if you give nothing it will be 273.15K): ")
    if T == '':
        T = 273.15
    LnK = DeltGibbsEkJpM/(-1*R*T)
    K = np.exp(LnK)
    pKa = -1*np.log10(K)
    print(f'pKa = {pKa}')
    return(pKa)

def PtFVal(DeltGibbsE, pKa):
    GibbsTxt = f'DeltG = {DeltGibbsE} Ha '
    pKaTxt = f'pKa = {pKa} '
    with open('ValsofInt.txt','w') as Valsfile:
        Valsfile.write(GibbsTxt)
        Valsfile.write(pKaTxt)
    Valsfile.close()
    return(Valsfile)

Input1 = open("butane.xyz", 'rw')
Input2 = open("propane.xyz", 'rw')

MoI1 = ade.Molecule(Input1, solvent_name=None, charge=0)
MoI2 = ade.Molecule(Input2, solvent_name=None, charge=0)

MoI1.optimise(method=xtb)
MoI2.optimise(method=xtb)

CoI1 = ade.Calculation(name=MoI1.name, molecule=MoI1, method=orca, keywords=orca.keywords.hess, n_cores=48)
CoI2 = ade.Calculation(name=MoI2.name, molecule=MoI2, method=orca, keywords=orca.keywords.hess, n_cores=48)

CoI1.output.filename = MoI1.name+'.out'
CoI2.output.filename = MoI2.name+'.out'

MoI1.calc_thermo(calc=CoI1, n_cores=48)
MoI2.calc_thermo(calc=CoI2, n_cores=48)

GibbsE1 = MoI1.free_energy
GibbsE2 = MoI2.free_energy

DelGibbs = GibbsE1 - GibbsE2

print(f'G = {DelGibbs} Ha')