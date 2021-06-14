#---------------------------------------------------------------#
#								  #
#    generator_struktura->hamiltoniany.py 		          #
#								  #
#---------------------------------------------------------------#

import numpy as np
import matplotlib.pyplot as plt
import sys
from math import sqrt, pi, cos

#---# Wstepne operacje #--------------------------------------------#

pos = []
hd = []

with open('hamiltoniany/positions.txt') as plik:
	pos = [list(map(float, wiersz.split(','))) for wiersz in plik]
	pos = np.array(pos)

dist = 0.

if(sys.argv[1] == "benzen"):
	Vpp = -2.7
	acc = 0.139
else:
	Vpp = -2.8	#grafen
	acc = 0.142

dist_tab = np.zeros([len(pos)])
dist_tab.fill(1000.) #Bo wtedy nie policzy "pustych" miejsc w tablicy jako najmniejszych odleglosci
hd = np.zeros([len(pos), len(pos)])
hd_app = np.zeros([len(pos), len(pos)])

#--------------------------------------------------------------------#

#Definicje funkcji---------------------------------------------------#

def distSites(n, m):
	"""Oblicza odległość między danymi dwoma atomami i wyraża ją w nanometrach.

	:param n: (int) Numer pierwszego rozważanego atomu.
	:param m: (int) Numer drugiego rozważanego atomu.
	:return: (double) Odległość w nanometrach między danymi dwoma atomami.
	"""
	return sqrt((pos[n][0] - pos[m][0])**2 + (pos[n][1] - pos[m][1])**2)
	
def upperTriangle():
	"""Wypełnia górny trójkąt hamiltonianu ribbonu.

	"""
	for i in range(len(pos)):
		for j in range(len(pos)):
			temp = distSites(i, j)
			temp = round(temp, 2)
			if(temp != 0.):
				dist_tab[j] = temp
	
		min_value = dist_tab.min()
	

		for j in range(len(pos)):
			if(min_value != 1000.):
				if(dist_tab[j] == min_value):
					hd[i][j] = Vpp * cos((abs(pos[i][0] - pos[j][0])) / acc) #w nanometrach!

		dist_tab.fill(1000.)

def diagAndLowerTriangle():
	"""Uzupełnia dolny trójkąt i przekątną hamiltonianu ribbonu.

	"""
	for i in range(len(hd)):
		if(sys.argv[1] == "benzen"):
			hd[i][i] = -8.57
		else:
			hd[i][i] = 0	#grafen
		for j in range(len(hd)):
			hd[j][i] = hd[i][j]
			
def printAppHD():
	"""Wyświetla w konsoli hamiltonian ribbonu z przybliżonymi wartościami.

	"""
	for i in range(len(hd)):
		for j in range(len(hd)):
			hd_app[i][j] = round(hd[i][j], 2)

	print(hd_app)
	
def saveHDtoFile():
	"""Zapisuje hamiltonian ribbonu do pliku.

	"""
	file1 = open("hamiltoniany/pozostale/z_generatora/hamiltonian_hd.txt","w")
	for i in range(len(hd)):
		for j in range(len(hd)):
			file1.write(str(hd[i][j]))
			if(j != len(hd) - 1):
				file1.write(",")
		file1.write("\n")
	
	file1.close()
	
def plotSaveStructure():
	"""Zapisuje wizualizację struktury do pliku w postaci wykresu XY.

	"""
	x = []
	y = []

	for i in range(len(pos)):
		x.append(pos[i][0])
		y.append(pos[i][1])

	plt.figure()
	plt.scatter(x, y)
	plt.savefig("wykresy/_struktura.png")

#--------------------------------------------------------------------#

upperTriangle()
diagAndLowerTriangle()
printAppHD()
saveHDtoFile()
plotSaveStructure()
