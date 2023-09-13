from protein import *
from sphere import *
import math

class Vector:
	""" Crée et manipuler les vecteurs pour déterminer la position de la membrane 
	"""

	def __init__(self, coord_origin, coord_point):
		"""Crée un Vecteur

		Parameters
		----------
		coord_origin : list 1D
		coord_point  : numpy 1D-array
	    """
		self.coord_origin = coord_origin # origin (x,y,z)
		self.coord_point = coord_point
		self.coord_vector = self.creat_vector()
		self.equations_cartesiennes = self.calcul_equations_cartesiennes() # selon point
		
	def creat_vector(self):
		""" Crée un vecteur à partir de 2 points et retourner les coordonées
		du vecteur.

		Note
		----
		La création du vecteur nécessite les coordonnées du centre de masse et
		d'un point placer aléatoirement sur la sphère

		Parameters
		----------
		self.coord_origin : 
		self.coord_point : 

		Returns
		-------
		list 1D
			Les coordonnées du vecteur
		"""
		x, y, z = self.coord_origin
		x1 = self.coord_point[0]
		y1 = self.coord_point[1]
		z1 = self.coord_point[2]
		vector = [x1 - x, y1 - y, z1 - z]

		#print("Création du vecteur terminé.")
		return(vector)

	def calcul_equations_cartesiennes(self):
		""" Calcul de l'équation cartésienne du plan : ax + by + cz + d = 0
		avec pour vecteur normale self.coord_vector et pour point appartenant
		au plan self.coord_point
		"""
		x = self.coord_point[0]
		y = self.coord_point[1]
		z = self.coord_point[2]
		#x, y, z = self.coord_origin
		a, b, c = self.coord_vector
		d = -a * x - b * y - c * z

		equation_plan1 = [a,b,c,d]
		equation_plan2 = [a,b,c,d+14] # Plan2 distant à 14 angstroms du plan1
		
		#print("Le calcul de l'équation cartésienne des plans est terminé.")
		return(equation_plan1,equation_plan2)

	def __str__(self):
		""" Redéfinition de print 
		Return : 
			Les coordonées du vecteur, les coordonnées du points d'origine et
			celle du point placer aléatoirement sur la sphère. Les équation 
			cartésienne.
		"""
		a, b, c, d = self.equations_cartesiennes[0] # Plan1
		a1, b1, c1, d1 = self.equations_cartesiennes[1] # Plan2

		return(f"\n\nVecteur :\nLes coordonnées du vecteur : {self.coord_vector}"
				f"\nLe point d'origine {self.coord_origin}\nLe point aléatoire"
				f" sur la sphere : {self.coord_point}\nEquation cartésienne "
				f"plan1 : {a}x {b:+.2f}y {c:+.2f}z {d:+.2f} = 0\nEquation"
				f" cartésienne plan2 : {a1}x {b1:+.2f}y {c1:+.2f}z "
				f"{d1:+.2f} = 0\n\n")