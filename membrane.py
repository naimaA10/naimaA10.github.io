import pymol
import math
from vector import *
from protein import *
from sphere import *
import numpy as np

class Membrane:
	def __init__(self, eq_cartesienne1, eq_cartesienne2):
		self.membrane1 = eq_cartesienne1 # [a,b,c,d] pour équation cartésienne ax+by+cz+d=0
		self.membrane2 = eq_cartesienne2 
		self.hydrophobicity = 0 # hydrophobicité
		self.protein_in_membrane = True
		#print("Création de l'objet de type Membrane terminé.")

	def __str__(self):
		""" Redéfinition de la fonction print pour un objet de type Membrane
		Arg :
			Cette méthode ne prend pas d'argument

		Return : 
			Chaîne de caractère indiquant les équations cartésiennes des deux
			plans de l'objet Membrane et d'hydrophobicité associé
		"""
		a, b, c, d = self.membrane1	
		a1, b1, c1, d1 = self.membrane2
		return(f"\n\nMembrane :\nEquation cartésienne du plan1 : {a}x {b:+.2f}y "
				f"{c:+.2f}z {d:+.2f} = 0\nEquation cartésienne plan2 : {a1}x "
				f"{b1:+.2f}y {c1:+.2f}z {d1:+.2f} = 0\nHydrophobicité : {self.hydrophobicity}\n\n")


	def membrane_movement(self, one_movement_step):
		"""Déplacer la membrane de 1 Angstrom selon la direction de la droite
		"""
		self.membrane1[3] += one_movement_step
		self.membrane2[3] += one_movement_step
		#print(f"La membrane est déplacée de {one_movement_step} angstrom(s).")
	


	def hydrophobicity_calculation(self, dict_residus):
		"""Calcul de l'hydrophobicité de la membrane """
		#print("DANS hydrophobicity_calculation(dict_residus)")
		# Initialiser les variables
		nb_hydrophobe = 0
		nb_hydrophile = 0
		total_hydrophobe = 0
		total_hydrophile =  0

		# Nombre total de résidus hydrophobe et hydrophile
		for id_residu in dict_residus:
			if (dict_residus[id_residu]['hydrophobe']):
				total_hydrophobe+=1
			else:
				total_hydrophile+=1
		#print(f"total_hydrophobe {total_hydrophobe}\ntotal_hydrophile {total_hydrophile}")
		
		# Donnée de équation cartésienne
		a, b, c, d = self.membrane1 
		a1, b1, c1, d1 = self.membrane2
		#print(f"membrane 1 : {a}x {b:+.2f}y {c:+.2f}z {d:+.2f} = 0\n"
		#	f"membrane 2 : {a1}x {b1:+.2f}y {c1:+.2f}z {d1:+.2f} = 0")

		# Pour chaque CA exposé de la protéine
		for id_residu in dict_residus:

			# Coordonnées du CA
			x = dict_residus[id_residu]['coord'][0]
			y = dict_residus[id_residu]['coord'][1]
			z = dict_residus[id_residu]['coord'][2]

			# La distance entre le CA et le plan :
			# |(ax + by + cz + d)| / √(a^2 + b^2 + c^2)
			distance1 = abs(a * x + b * y + c * z + d) / math.sqrt(a**2 + b**2 + c**2)
			distance2 = abs(a1 * x + b1 * y + c1 * z + d1) / math.sqrt(a1**2 + b1**2 + c1**2)

			# Vérifier que le CA est entre les deux plans donc dans la membrane
			if (distance1 < 14 and distance2 < 14) : 
				# Indiquer si le résidus est hydrophobe
				if (dict_residus[id_residu]['hydrophobe']):  
					nb_hydrophobe += 1
				else :
					nb_hydrophile += 1
			# Si aucun résidus hydrophobe ou hydrophile est trouver alors la protéine est hors de la membrane			
			if (nb_hydrophobe == 0 and nb_hydrophile == 0):
				self.protein_in_membrane= False
				#print(f"nb_hydrophobe {nb_hydrophobe}, 	nb_hydrophile {nb_hydrophile}")
				#print("La protéine n'est pas dans la membrane")
			else:
				## Calcule de l'hydrophobicité : 
				# (nb résidus exposé polaire dans mb / nb résidus exposé polaire)  + (nb résidus exposé hydrophobe dans mb/ nb rédisus exposé hydrophobe)
				if (total_hydrophobe != 0):
					#print(f"nb_hydrophobe {nb_hydrophobe}, 	nb_hydrophile {nb_hydrophile}")
					self.hydrophobicity = (nb_hydrophile / total_hydrophile) + (nb_hydrophobe / total_hydrophobe)
				else:
					self.hydrophobicity = 0	
		#print(f"nb_hydrophobe {nb_hydrophobe}, 	nb_hydrophile {nb_hydrophile}, hydrophobicity {self.hydrophobicity}")
		#print("Calcul de l'hydrophobicité de la membrane terminé.\n")
		return self.hydrophobicity
	
	def maximum_hydrophobicity(self, dict_residus):
		""" Calcul l'hydrophobicité de la membrane à cette position stocke 
		cette valeur et déplace la membrane de 1 angstroms jusqu'a obtenir 
		la position de la membrane pour laquelle l'hydrophobicité est 
		maximum"""

		hydrophobicity = self.hydrophobicity_calculation(dict_residus)
		equation_plan1 = self.membrane1
		equation_plan2 = self.membrane2
		
		self.protein_in_membrane = True
		#print(f"initiale hydrophobicité : {hydrophobicity}\ndebut dans mb ? {self.protein_in_membrane}")
		while (self.protein_in_membrane): # La protéin est dans la membrane
			self.membrane_movement(1)
			new_hydrophobicity = self.hydrophobicity_calculation(dict_residus)
			if new_hydrophobicity > hydrophobicity:
				hydrophobicity = new_hydrophobicity
				equation_plan1 = self.membrane1
				equation_plan2 = self.membrane2
				
		self.protein_in_membrane = True
		while (self.protein_in_membrane): # La protéin est dans la membrane
			self.membrane_movement(-1)
			new_hydrophobicity = self.hydrophobicity_calculation(dict_residus)
			if new_hydrophobicity > hydrophobicity:
				hydrophobicity = new_hydrophobicity
				equation_plan1 = self.membrane1
				equation_plan2 = self.membrane2
		return(hydrophobicity,equation_plan1,equation_plan2)

	def movement_elargie(self,pas):
		"""Elarfie la membrane de pas*2 Angstroms selon la direction de la droite
		"""
		self.membrane1[3] += pas
		self.membrane2[3] -= pas
		#print(f"La membrane est déplacée de {one_movement_step} angstrom(s).")

	def maximum_elargir(self, pas, dict_residus):
		"""Elargir la membrane de pas*2
		
		Parametre :
			pas : int
		"""
		hydrophobicity = self.hydrophobicity

		self.protein_in_membrane = True
		while (self.protein_in_membrane):
			self.movement_elargie(1) 
			new_hydrophobicity = self.hydrophobicity_calculation(dict_residus)
			if new_hydrophobicity > hydrophobicity:
				hydrophobicity = new_hydrophobicity
				equation_plan1 = self.membrane1
				equation_plan2 = self.membrane2

		return(hydrophobicity)

	def add_pdb(self, name_pdb, coords_CA):
		""" 
		Ajouter la membrane à la pdb
		"""
		print(f"coords points {coords_CA}\n\n\n")
		print( coords_CA[0][0])
		xmin = xmax = coords_CA[0][0]
		ymin = ymax = coords_CA[0][1]

		for point in coords_CA:
			x = point[0]
			y = point[1]
			if x > xmax :
				xmax = x
			if x < xmin :
				xmin = x
			if y > ymax : 
				ymax = y
			if y < ymin :
				ymin = y

		list_x = list(np.arange(xmin, xmax, 1))
		list_y = list(np.arange(ymin, ymax, 1))
		plan1 = list()
		plan2 = list()

		# membrane 1
		a, b, c, d = self.membrane1
		
		# membrane 2
		a1, b1, c1, d1 = self.membrane2	
		print("ok")
		for x in list_x:
			for y in list_y:
				z= (-a * x - b * y - d) / c# membrane 1
				z1= (-a1 * x - b1 * y - d1) / c1 # membrane 2
				plan1.append([x,y,z])
				plan2.append([x,y,z1])
		
		# Ouvrir pdb sur pymol et ajouter les points et vecteurs
		pymol.finish_launching(['pymol', '-qc'])
		pymol.cmd.load(name_pdb, name_pdb[:-4])
		
		for i, point in enumerate(plan1):
			x, y , z= point
			pymol.cmd.pseudoatom(f'plan1_{i}', resi='CM', pos=[x,y,z], color=(1, 0, 0), label=f'plan1 {i}')
		
		for i, point in enumerate(plan2):
			x, y , z= point
			pymol.cmd.pseudoatom(f'plan1_{i}', resi='CM', pos=[x,y,z], color=(0, 0, 1), label=f'plan1 {i}')

		# save pdb
		new_pdb = f'{name_pdb[:-4]}_membrane.pdb'
		pymol.cmd.save(new_pdb)