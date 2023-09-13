from protein import *
from vector import *
from membrane import *
import math
import numpy as np

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