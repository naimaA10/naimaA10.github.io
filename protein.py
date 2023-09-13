from Bio.PDB.DSSP import DSSP
from Bio import PDB 
from Bio.PDB import PDBParser
from sphere import *
from vector import *
import sys

class Protein:

	def __init__(self, name):
		""" 
		Initialisation de la protéine

			Nom du fichier pdb
			Dictionnaire des résidus contient:
				Les coordonnées atomiques 
				La surface accessible au solvant de chaque 
				Le type hydrophobe ou hydrophile
				Le nom du résidu
			Les coordonnées du centre de masse (null initialement)
		"""
		self.name = name 
		self.dict_residus = dict() 
		self.center_of_mass = [0, 0, 0] 
		
	def __str__(self):
		"""	Redéfinition de la fonction print

		Return : 
			Retourne les information sur la protéine dont le nombre de résidus
			accessible au solvant et les coordonnées du centre de
			masse.
		"""
		return (f"\n\nProtéine :\nNom :{self.name[:-4]}\nNombre de résidus"
				f" accessible au solvant : {len(self.dict_residus)}\nCoordonnées"
				f" du centre de masse : {self.center_of_mass}\n\n")

	def extract_CA_access_solvant (self):
		"""
		Extraction des coordonnées atomiques des carbones alpha (CA)
		Calcul de la surface accessible au solvant pour chaque residu
		Calcul du centre de masse 
		"""
		# Parser le fichier pdb
		parser = PDBParser()
		structure = parser.get_structure(self.name[0:-4], self.name)
		
		# Objet DSSP
		dssp = DSSP(structure[0], self.name, dssp = 'mkdssp')

		residus_hydrophobe = ["PHE", "GLY", "ILE", "LEU", "MET", "VAL", "TRP", "TYR"]

		#Extraction des coordonnées atomiques des carbones alpha (CA)
		for model in structure:
			for chain in model:
				for residue in chain:
					
					# Identifiant et nom du residu
					id_residu = residue.get_id()[1]
					name_residu = residue.get_resname()
					
					if residue.get_id()[0] == ' ' and residue.get_id()[2] == ' ':
						for atom in residue:
							if atom.get_name() == 'CA':

								# Coordonnées du résidu et accessibilité au solvant à partir de dssp
								coord_residu = atom.get_coord()
								#print(f"chain id : {type(chain.id)}\n id_residu {type(id_residu)}")
								accessibility = dssp[(chain.id, id_residu)][3]
								if name_residu in residus_hydrophobe:
									type_hydrophobe = True
									
								else :
									type_hydrophobe = False

								if (accessibility > 0.3):
									self.dict_residus[id_residu] = {'name' : name_residu,
																	'coord' : coord_residu,
																	'hydrophobe' : type_hydrophobe,
																	'accessibility' : accessibility}
		print(f"Extraction des coordonnées atomiques des CA terminés.\nCalcul"
			f" de l'accessibilité au solvant terminés.")

	def calcul_center_of_mass(self):
		"""Calcul du centre de masse"""
		total_mass = 0

		for id_residu in self.dict_residus:
			total_mass += 12
			x = self.dict_residus[id_residu]['coord'][0]
			y = self.dict_residus[id_residu]['coord'][1]
			z = self.dict_residus[id_residu]['coord'][2]
			
			self.center_of_mass[0] += x * 12
			self.center_of_mass[1] += y * 12
			self.center_of_mass[2] += z * 12

		self.center_of_mass[0] /= total_mass
		self.center_of_mass[1] /= total_mass
		self.center_of_mass[2] /= total_mass
		print(f"Calcul du centre de masse terminé.")


if __name__ == "__main__":
	""" Programme principal, crée une protéine, extraire les coordonnées CA,
	calculer l'accessibilité au solvant et calculer le centre de masse.
	"""
	try :
		name_pdb = sys.argv[1]
		if name_pdb[-4:] == ".pdb" :
			print(f"Analyse du fichier {name_pdb}")
	except:
		print(f"Veuillez entrez le nom d'un fichier pdb.")

	protein = Protein(name_pdb)	
	protein.extract_CA_access_solvant()
	protein.calcul_center_of_mass()
	print(protein) # Information sur la protéine

	""" Crée une sphère de rayon 20 avec pour centre le centre de masse et 
	placer aléatoirement 10 points"""

	sphere = Sphere(protein.center_of_mass, 20, 50)
	sphere.build_sphere()
	print(sphere) # Information sur la sphere

	""" Crée un vecteur	entre le centre de masse et chaque points placer 
	aléatoirement et calculer l'équation du plan orthogonal à ce vecteur
	Ces information sont ensuite stockés dans le dictionnaire dic_vect_eq"""

	list_vectors = list()

	for coord_point in sphere.coord_points:
		vector = Vector(protein.center_of_mass, coord_point)
		list_vectors.append(vector)
	
	""" Crée une membrane à partir des équations cartésiennes représentant les
	deux plan. Calculer l'hydrophobicité stocker cette valeur et déplacer la 
	membrane de 1 angstrom et calculer de nouveau l'hydrophobicité si la valeur
	est supérieur à cette présedente enregistrer les valeurs de la nouvelle 
	membrane sinon continuer à déplacer la membrane. Faire de même pour chaque 
	membrane orthogonal au vecteur normal dont les valeurs sont dans 
	dict_vect_eq """
	#AFFICHAGE
	print(f"liste vecteurs :")
	for v in list_vectors:
		print(v)
	membrane = Membrane(list_vectors[0].equations_cartesiennes[0], list_vectors[0].equations_cartesiennes[1])
	hydrophobicity, equation1, equation2 = membrane.maximum_hydrophobicity(protein.dict_residus)


	for vector in list_vectors: # dict_vect_eq[cood_point] = {'coord_vector': ,'equation1': , 'equation2':}
		equation1 = vector.equations_cartesiennes[0]
		equation2 = vector.equations_cartesiennes[1]
		new_membrane = Membrane(equation1,equation2)
		new_hydrophobicity, new_equation1, new_equation2= membrane.maximum_hydrophobicity(protein.dict_residus)
		#print(f"CALCUL {hydrophobicity}")

		if new_hydrophobicity > hydrophobicity :
			membrane = Membrane(new_equation1, new_equation2)
			hydrophobicity, equation1, equation2 = membrane.maximum_hydrophobicity(protein.dict_residus)
			#print(f"NEW {hydrophobicity}")
	print(f"Calcul de la meilleur valeur d'hydrophobicité terminés.\nLa"
		f" membrane dont l'hydrophobicité est maximal  : {membrane}")

	# Elargir la membrane pour trouver la meilleur hydrophobicité
	#membrane.maximum_elargir(1, protein.dict_residus)
	print([protein.dict_residus[i]["coord"] for i in protein.dict_residus])
	membrane.add_pdb(name_pdb, [protein.dict_residus[i]["coord"] for i in protein.dict_residus]) # coords_points les coordonnées des carbones alphas