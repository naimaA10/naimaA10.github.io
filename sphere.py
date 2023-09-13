from protein import *
from vector import *
from membrane import *
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Sphere:

	def __init__(self, center_mass, radius, num_points):
		""" Initialisation de la sphère

		Arg : 
			center_mass : Les coordonnées du centre de masse de la protéine
			radius :Le rayon de la sphère
			num_points : Le nombre de points à placer sur la sphère
			coord_points : Liste pour stocker les coordonnées des points
		"""
		self.center_mass = center_mass
		self.radius = radius
		self.num_points = num_points
		self.coord_points = []

	def __str__(self):
		"""	Redéfinition de la fonction print

		Return :
			Chaîne de caractère décrivant la sphère avec le centre de la sphere
			son le rayon et le nombre de points à placer.
		"""
		return (f"\nSphère : \nCoordonnées du centre de la sphère : "
				f"{self.center_mass}\nRayon : {self.radius} angstroms\n"
				f"Nombre de points à placer : {self.num_points}\n")


	def build_sphere(self):
		"""	Construire la sphère avec pour centre le centre de masse et placer 
		les points sur une demi-sphère
		"""
		# Saaf-Kuijlaars
		print(f"sssssssss {self.num_points} et aussi {self.num_points + 1}")
		for k in range(1, self.num_points + 1):  # n-points pour éviter division par 0
			h = -1 + (2 * (k - 1) / (self.num_points - 1))
			print(h)
			if h != 1 :
				theta = math.acos(h)
				print(theta)
				if k == 1 or k is self.num_points:
					phi = 0
				else:
					phi = (phi + (3.6 / math.sqrt(self.num_points) * (1 / math.sqrt(1 - h * h)))) % (2 * math.pi)
				# Centering on the mass center
				x = self.center_mass[0] + math.sin(phi) * math.sin(theta)
				y = self.center_mass[1] + math.cos(theta)
				z = self.center_mass[2] + math.cos(phi) * math.sin(theta)
				print(f" z: {z} center : {self.center_mass[2]}")
				# Garder les points sur une demi-sphère
				if z > self.center_mass[2] :
					self.coord_points.append((x, y, z))
			print(f"Création de la sphère terminé.\nPlacer {self.num_points} points"
			f" alatoirement sur la sphère terminé.")
		return(self.coord_points)

if __name__ == "__main__":
	""" Programme principal, crée une protéine, extraire les coordonnées CA,
	calculer l'accessibilité au solvant et calculer le centre de masse.
	"""
	protein = Protein("1uaz_A.pdb")	
	protein.extract_CA_access_solvant()
	protein.calcul_center_of_mass()
	print(protein) # Information sur la protéine

	""" Crée une sphère de rayon 20 avec pour centre le centre de masse et 
	placer aléatoirement 10 points"""

	sphere = Sphere(protein.center_of_mass, 20,20)
	sphere.build_sphere()
	print(sphere) # Information sur la sphere

	""" Crée un vecteur	entre le centre de masse et chaque points placer 
	aléatoirement et calculer l'équation du plan orthogonal à ce vecteur
	Ces information sont ensuite stockés dans le dictionnaire dic_vect_eq"""

	dict_vect_eq = dict()

	for coord_point in sphere.coord_points: # Pour chaque point placer aléatoirement
		vector = Vector(protein.center_of_mass, coord_point) # crée vecteur et calcul équation cartésienne orthogonal au vecteur 
		dict_vect_eq[coord_point] = {'coord_vector' : vector,
							'equation1' : vector.equations_cartesiennes[0],
							'equation2' : vector.equations_cartesiennes[1]}
	
	""" Crée une membrane à partir des équations cartésiennes représentant les
	deux plan. Calculer l'hydrophobicité stocker cette valeur et déplacer la 
	membrane de 1 angstrom et calculer de nouveau l'hydrophobicité si la valeur
	est supérieur à cette présedente enregistrer les valeurs de la nouvelle 
	membrane sinon continuer à déplacer la membrane. Faire de même pour chaque 
	membrane orthogonal au vecteur normal dont les valeurs sont dans 
	dict_vect_eq """

	first_key = next(iter(dict_vect_eq)) # première clé du dictionnaire
	membrane = Membrane(dict_vect_eq[first_key]['equation1'], dict_vect_eq[first_key]['equation2'])
	hydrophobicity, equation1, equation2 = membrane.maximum_hydrophobicity(protein.dict_residus)
	print(f"INITIAL {hydrophobicity}")

	for coord_point in dict_vect_eq: # dict_vect_eq[cood_point] = {'coord_vector': ,'': , 'equation2':}
		equation1 = dict_vect_eq[coord_point]['equation1']
		equation2 = dict_vect_eq[coord_point]['equation2']
		new_membrane = Membrane(equation1,equation2)
		new_hydrophobicity, new_equation1, new_equation2= membrane.maximum_hydrophobicity(protein.dict_residus)
		print(f"CALCUL {hydrophobicity}")

		if new_hydrophobicity > hydrophobicity :
			membrane = Membrane(new_equation1, new_equation2)
			hydrophobicity, equation1, equation2 = membrane.maximum_hydrophobicityequation1(protein.dict_residus)
			print(f"NEW {hydrophobicity}")
	print(f"FINAL {hydrophobicity}")

##############################################################""
	# Créez une grille de points 3D
	x = np.linspace(0, 250, 100)
	y = np.linspace(0, 80, 100)
	x, y = np.meshgrid(x, y)

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	sphere_points = np.array(sphere.coord_points)

	# Placer les points aléatoire (en blue)
	ax.scatter(sphere_points[:, 0], sphere_points[:, 1], sphere_points[:, 2],
				c='b', marker='o')

	# Placer les CA aléatoire (en rouge)
	ca_coordinates = list()
	for id_residu in protein.dict_residus:
		x = protein.dict_residus[id_residu]['coord'][0]
		y = protein.dict_residus[id_residu]['coord'][1]
		z = protein.dict_residus[id_residu]['coord'][2]
		ca_coordinates.append([x,y,z])
	
	sc = np.array(ca_coordinates)
	scp = np.array(sc)
	ax.scatter(scp[:, 0], scp[:, 1], scp[:, 2],
				c='r', marker='o')

	# Placer le centre de masse (en vert)
	cdm = np.array(protein.center_of_mass)
	cdmp = np.array(cdm)
	ax.scatter(cdmp[0], cdmp[1], cdmp[2],
				c='g', marker='o')

	# Placer les vecteurs
	print(dict_vect_eq)
	for key in dict_vect_eq :
		vector = dict_vect_eq[key]['coord_vector']
		ax.quiver(protein.center_of_mass[0], protein.center_of_mass[1], protein.center_of_mass[2], vector.coord_vector[0], vector.coord_vector[1], vector.coord_vector[2])


	# Définissez les équations cartésiennes du plan 1 (par exemple, plan1 = [a, b, c, d])
	a, b, c, d = dict_vect_eq[first_key]['equation1']

	# Créez une grille 2D de points (x, y)
	x_range = np.linspace(0, 250, 100)
	y_range = np.linspace(0, 80, 100)
	x_mesh, y_mesh = np.meshgrid(x_range, y_range)

	# Calculez les coordonnées z pour le plan 1 en utilisant ses équations
	if c != 0:
		z_mesh1 = (-a * x_mesh - b * y_mesh - d) / c
	else:
		z_mesh1 = np.zeros_like(x_mesh)

	# Tracez le plan 1 (par exemple, en bleu)
	ax.plot_surface(x_mesh, y_mesh, z_mesh1, color='b', alpha=0.5, label='Plan 1')
	ax.plot_surface(x_mesh, y_mesh, z_mesh1+14, color='b', alpha=0.5, label='Plan 1')


	# Placer la membrane
	plt.title("Protéine (rouge), centre de masse (vert), les points (blue) et plan en 3D")
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	plt.show() # Affichez sur le terminal
