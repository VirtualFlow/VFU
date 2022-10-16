from chimera import openModels, Molecule
from AddCharge import estimateFormalCharge
mol = openModels.list(modelTypes=[Molecule])[0]

fc = estimateFormalCharge(mol.atoms)
atomicSum = 0
for a in mol.atoms:
	atomicSum += a.element.number
if (atomicSum + fc) % 2 == 0:
	print  fc		#charge estimade
else:
	print  fc		#Bad charge estimade
