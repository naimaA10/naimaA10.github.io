# Project assignment and detection of the transmembrane parts of a protein

The aim of this project is to identify the position of the membranes of a transmembrane protein from its PDB structural file. The approach of the method is described in the following article:

> Tusnády GE, Dosztányi Z, Simon I. Transmembrane proteins in the Protein Data Bank: 
identification and classification. Bioinformatics. 2004 Nov 22;20(17):2964-72. Epub 2004 Jun 4


## Install

to use this project, follow these installation steps :

Clone this repository :

	https://github.com/naimaA10/naimaA10.github.io.git

PyMol

	$ sudo apt-get install pymol

Create conda environment and install dependendies:

	$  conda env create -f projet_transmembranaire.yml

Load conda environment:

	$ conda activate projet_transmembranaire

## Run

To run the program
	
	$ python3 protein.py 1uaz.pdb

with 1uaz.pdb corresponding to the name of the pdb file

- This program will create a protein from the name_pdb file and the centre of mass of the protein will be calculated.

- The sphere is created from the co-ordinates of the centre of mass, a value corresponding to the radius of the sphere and the number of random points to be placed on the sphere.

- A list of vectors is then created, taking into account the centre of mass and the randomly placed points. For each vector, a corresponding member is produced.

- The membrane is moved in the direction of the vector with a step of 1 angstrom and its hydrophobicity is calculated to obtain the position of the membrane for which the hydrophobicity is maximum.

This program produces a display on the terminal indicating the progress of the process and a pdb file containing the position of the membrane. It can then be viewed on pymol.

Open pymol :

	$ pymol

Drag the pdb file into the graphical interface.