#!/usr/bin/env python

#---# Stale uzywane w programie #---#

n = 3		#liczba komorek oligofenylu polaczonych ze soba
T = 0.417	#stala translacja komorki oligofenylu wzdluz osi X (obliczone z odleglosci miedzyatomowych)
file1 = open("hamiltoniany/positions.txt","w")	#Plik positions.txt
pos = [[0.0,0.0,0.0], [0.0695,0.1204,0.0], [0.0695,-0.1204,0.0], [0.2085,0.1204,0.0], [0.2085,-0.1204,0.0], [0.278,0.0,0.0]]	#Pozycje atomow w pierwszym oligofenylu

#-----------------------------------#

def zapiszPozycje(i, ifNextLine):
	"""Generuje pozycje kolejnych atomów w łańcuchu oligofenylowym o zadanej długości.

	:param i: (int) Numer obecnie rozpatrywanego oczka w łańcuchu benzenowym.
	:param ifNextLine: (bool) Informacja o potrzebie przejścia do nowej linii.
	"""
	for j in range(6):
		file1.write(str(round(pos[j][0] + (i * T), 4)))
		file1.write(",")
		file1.write(str(round(pos[j][1], 4)))
		file1.write(",")
		file1.write(str(round(pos[j][2], 4)))
		if(ifNextLine == True or (j != 5)):
			file1.write("\n")

for i in range(n):
	if(i != (n - 1)):
		zapiszPozycje(i, True)
	else:
		zapiszPozycje(i, False)

file1.close()
