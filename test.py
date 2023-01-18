from OGEN.OPLS import create_cmd, create_params, create_zmat, get_bonded, get_dispersion
from rdkit import Chem
import subprocess

mol_file = '../OGEN_C6H6/C6H6.mol'
mol = Chem.MolFromMolFile(mol_file, removeHs=False, sanitize=False)
create_zmat(mol)
create_params('init')
create_cmd('init')
subprocess.run('bash bosscmd', check=True, shell=True)
create_params('internal')
create_cmd('internal')
subprocess.run('bash bosscmd', check=True, shell=True)
create_params('end')
create_cmd('end')
subprocess.run('bash bosscmd', check=True, shell=True)

# moleculeRDkit, newIndexToOriginalIndex, atomsNameOriginal, residueNameOriginal = utilities.generateRDkitMolecule(self.ifile, self.smile, self.workdir, molnameA, self.debug)
        # Зачем-то меняется порядок атомов

# moleculeA = Molecule.fromRDkitMolecule(moleculeRDkit)
        # dummyAtomsShift = 4 # dummy atoms (+3) and index to serial (+1)
        # bondsVariable, bondsAdditional, anglesVariable, anglesAdditional, torsionsVariable, torsionsAdditional = \
        #     geometry.generateBondAnglesAndDihedrals(moleculeRDkit, dummyAtomsShift)
        # cls._checkIfDiscontinuitiesInBondedLst(bondsVariable)
        # cls._checkIfDiscontinuitiesInBondedLst(anglesVariable)
        # cls._checkIfDiscontinuitiesInBondedLst(torsionsVariable)
        # shiftcoordinates = utilities.translateToceroZcoord(moleculeRDkit)  # BOSS imposes Z coordinate of the first atom in cartesians to be cero (internal coordinates)
        # atoms = cls.generateAtomsFromRDkit(moleculeRDkit, bondsVariable, anglesVariable, torsionsVariable, dummyAtomsShift)
        # return cls(atoms, bondsVariable, bondsAdditional, anglesVariable, anglesAdditional, torsionsVariable, 
        #     torsionsAdditional, None, None, 0, shiftcoordinates[0], shiftcoordinates[1], shiftcoordinates[2])
# logger.info('--- SUMMARY MOLECULE STRUCTURE ---')
# moleculeA.report()

# zmatName = zmat.write(moleculeA, molnameA, self.workdir,  writeAtomParameters = True)

# outFile, pdbFile = boss.run(zmatName, self.chargeAlgorithm, self.numberOfOptimizations, self.charge, molnameA, 
#     self.workdir, self.debug)

# moleculeA = Molecule.fromBOSS(zmatName, outFile, pdbFile, moleculeA.shiftX, moleculeA.shiftY, moleculeA.shiftZ)

# utilities.fixNonIntegerCharge(moleculeA)

# self.updateOriginalAtomIndexesAndSerials(moleculeA, newIndexToOriginalIndex, atomsNameOriginal)
# moleculeA.residueName = self.resname

# moleculeA_dual = moleculeA  # Fake possible dual topology for NAMD and others
# molnameAB_dual = molnameA