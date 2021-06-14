#---------------------------------------------------------------#
#								  #
#    transmission_zmieniona_metoda_sigmy.py  			  #
#								  #
#---------------------------------------------------------------#

import numpy as np
import matplotlib.pyplot as plt
import sys
from math import sqrt, pi, cos

#Definicje funkcji---------------------------------------------------#

def argNowa4(new_h0):
	"""Generuje macierze h0, h1l (związaną z doczepieniem lewej elektrody) i h1r (związaną z doczepieniem prawej elektrody) dla struktury oligoacenu o szerokości 2,48 anstrema.
	
	:param new_h0: (array) Macierz h0 pobrana z pliku. W tej funkcji nadpisuje się ją i zwraca.
	:return: (array, array, array) Macierze h0, h1l i h1r dla rozważanej struktury.
	"""
	
	Vpp = -2.8
	new_h0 = [
		  [0,Vpp,Vpp,0], 
		  [Vpp,0,0,Vpp], 
		  [Vpp,0,0,0], 
		  [0,Vpp,0,0]]
	
	h1r = [
		  [0,0,Vpp,0], 
		  [0,0,0,Vpp], 
		  [0,0,0,0], 
		  [0,0,0,0]]
		  
	h1l = [
		  [0,0,0,0], 
		  [0,0,0,0], 
		  [Vpp,0,0,0], 
		  [0,Vpp,0,0]]
		  
		  
	h1l = np.asarray(h1l)
	h1r = np.asarray(h1r)
	new_h0 = np.asarray(new_h0)
	
	
	return new_h0, h1r, h1l
	
def argNowa6(new_h0, sysarg3):
	"""Generuje macierze h0, h1r (związaną z doczepieniem prawej elektrody) oraz ewentualnie h1l (związaną z doczepieniem lewej elektrody) dla struktury oligoacenu o szerokości 4,92 anstrema.
	
	:param new_h0: (array) Macierz h0 pobrana z pliku. W tej funkcji nadpisuje się ją i zwraca.
	:param sysarg3: (int) Wartość logiczna w formie liczby, określająca, czy elektrody prawa i lewa się różnią od siebie.
	:return: (array, array, array/int) Macierze h0, h1l i h1r dla rozważanej struktury.
	"""

	Vpp = -2.8

	new_h0 = [[0,Vpp,0,Vpp,0,0], 
		  [Vpp,0,0,0,Vpp,0], 
		  [0,0,0,0,0,Vpp], 
		  [Vpp,0,0,0,0,0], 
		  [0,Vpp,0,0,0,Vpp], 
		  [0,0,Vpp,0,Vpp,0]]
	
	h1r = [[0,0,0,Vpp,0,0], 
		[0,0,0,0,Vpp,0], 
		[0,0,0,0,0,Vpp], 
		[0,0,0,0,0,0], 
		[0,0,0,0,0,0], 
		[0,0,0,0,0,0]]
		  
	new_h0 = np.asarray(new_h0)
	h1r = np.asarray(h1r)
	
	if(sysarg3 == 1):
		h1l = [[0,0,0,0,0,0], 
		     [0,0,0,0,0,0], 
		     [0,0,0,0,0,0], 
		     [Vpp,0,0,0,0,0], 
		     [0,Vpp,0,0,0,0], 
		     [0,0,Vpp,0,0,0]]
		     
		h1l = np.asarray(h1l)
		    
		return new_h0, h1r, h1l	#Tak, maja byc odwrotnie, bo chyba sie pomylilem

	#Koniec elementow macierzowych
	
	return new_h0, h1r, 0

#Koniec definicji funkcji ----------------------------------------------------#

#-----------------------------------------------------------------------------#


with open('hamiltoniany/hamiltonian_h0.txt') as plik:
	h0 = [list(map(float, wiersz.split(','))) for wiersz in plik]
	h0 = np.array(h0)

with open('hamiltoniany/hamiltonian_h1.txt') as plik:
	h1 = [list(map(float, wiersz.split(','))) for wiersz in plik]
	h1 = np.array(h1)

with open('hamiltoniany/hamiltonian_hd.txt') as plik:
	hd = [list(map(float, wiersz.split(','))) for wiersz in plik]
	hd = np.array(hd)
	
with open('hamiltoniany/hamiltonian_h1L.txt') as plik:
	h1l = [list(map(float, wiersz.split(','))) for wiersz in plik]
	h1l = np.array(h1l)

#---------# Sprawdzenie argumentow startowych programu #---------------------------#

hamilt = hd
new_h0 = h0
new_hamilt = []

sysarg3 = 0

if(sys.argv[2] == "lr"):

	sysarg3 = 1

if(sys.argv[1] == "oligoacen_4"):

	new_h0, h1, h1l = argNowa4(new_h0)
	
elif(sys.argv[1] == "oligoacen_6"):

	new_h0, h1, h1l = argNowa6(new_h0, sysarg3)

#---------------------------------------
f_start = -15.0	#Musi byc ujemne, bo jeszcze w 3 linijki nizej nie dalem mozliwosci ujemnego konca albo dodatniego poczatku
f_koniec = 0.0
df = 0.1
T_args = np.linspace(f_start, f_koniec, round((abs(f_start) / df) + (f_koniec / df) + 1), endpoint=True) #200 bez zera, 201 z zerem (numerycznie nie liczy z zerem)
T_vals_re = np.zeros([round((abs(f_start) / df) + (f_koniec / df) + 1)])			#, dtype = 'complex_')
T_vals_im = np.zeros([round((abs(f_start) / df) + (f_koniec / df) + 1)])			#, dtype = 'complex_')
T_iter = 0
S_0 = np.zeros([len(new_h0), len(new_h0)])
S_1 = np.zeros([len(new_h0), len(new_h0)])
np.fill_diagonal(S_0, 1)
Sd = np.zeros([len(hamilt), len(hamilt)])
np.fill_diagonal(Sd, 1)
SIGMA_R_minus = []
SIGMA_R_plus = []
for i in range(round((abs(f_start) / df) + (f_koniec / df) + 1)):
	SIGMA_R_plus.append(np.zeros([len(hamilt), len(hamilt)], dtype = 'complex_'))
for i in range(round((abs(f_start) / df) + (f_koniec / df) + 1)):
	SIGMA_R_minus.append(np.zeros([len(hamilt), len(hamilt)], dtype = 'complex_'))
delta = 0.01
#----------------------------------------

f = f_start

while(f <= f_koniec):

	E = f

	#petla dla z+

	sigma_plus = np.zeros([len(new_h0), len(new_h0)])
	z = E + 1j * delta
	iteration = 0
	alfa = 0.5
	
	sigma_plus_minus_minus = np.zeros([len(new_h0), len(new_h0)])
	sigma_plus_minus = np.zeros([len(new_h0), len(new_h0)])
	
	while_number = 1
	while(while_number == 1):
		suma = 0
		if(iteration >= 2):
			sigma_plus = alfa * sigma_plus_minus_minus + (1 - alfa) * sigma_plus_minus
		sigma_plus = (h1 - z * S_1) @ (np.linalg.inv(z * S_0 - new_h0 - sigma_plus)) @ (h1.conj().T - z * S_1.conj().T)
		if(iteration >= 2):

			for i in range(len(new_h0)):
				for j in range(len(new_h0)):
					suma += (sigma_plus[i, j] - sigma_plus_minus[i, j])**2
			if(abs(suma) <= 0.00001):
				while_number = 0

		sigma_plus_minus_minus = sigma_plus_minus
		sigma_plus_minus = sigma_plus
		iteration = iteration + 1
	
	#koniec petli dla z+
	#petla dla z-

	sigma_mniej = np.zeros([len(new_h0), len(new_h0)])
	z_minus = E - 1j * delta
	iteration = 0
	
	sigma_mniej_minus_minus = np.zeros([len(new_h0), len(new_h0)])
	sigma_mniej_minus = np.zeros([len(new_h0), len(new_h0)])
	
	while_number = 1
	while(while_number == 1):
		suma = 0
		if(iteration >= 2):
			sigma_mniej = alfa * sigma_mniej_minus_minus + (1 - alfa) * sigma_mniej_minus
		sigma_mniej = (h1 - z_minus * S_1) @ (np.linalg.inv(z_minus * S_0 - new_h0 - sigma_mniej)) @ (h1.conj().T - z_minus * S_1.conj().T)
		if(iteration >= 2):

			for i in range(len(new_h0)):
				for j in range(len(new_h0)):
					suma += (sigma_mniej[i, j] - sigma_mniej_minus[i, j])**2
			if(abs(suma) <= 0.00001):
				while_number = 0

		sigma_mniej_minus_minus = sigma_mniej_minus
		sigma_mniej_minus = sigma_mniej
		iteration = iteration + 1

	#koniec petli dla z-

	sigma_r_plus = np.zeros([len(hamilt), len(hamilt)], dtype = 'complex_')
	sigma_l_plus = np.zeros([len(hamilt), len(hamilt)], dtype = 'complex_')
	sigma_r_minus = np.zeros([len(hamilt), len(hamilt)], dtype = 'complex_')
	sigma_l_minus = np.zeros([len(hamilt), len(hamilt)], dtype = 'complex_')

	for i in range(len(new_h0)):
		for j in range(len(new_h0)):
			sigma_r_plus[i, j] = sigma_plus[i, j]
	if(sysarg3 != 1):
		for i in range(len(new_h0)):
			for j in range(len(new_h0)):
				sigma_l_plus[len(hamilt) - len(new_h0) + i, len(hamilt) - len(new_h0) + j] = sigma_plus[i, j]
	for i in range(len(new_h0)):
		for j in range(len(new_h0)):
			sigma_r_minus[i, j] = sigma_mniej[i, j]
	if(sysarg3 != 1):
		for i in range(len(new_h0)):
			for j in range(len(new_h0)):
				sigma_l_minus[len(hamilt) - len(new_h0) + i, len(hamilt) - len(new_h0) + j] = sigma_mniej[i, j]
			
	#Jezeli nie ma sysarg3 == 1 to wtedy od razy liczy transmisje. Jezeli sysarg3 == 1 to wtedy nie liczy, tylko zapisuje sigmy do tablicy SIGMA_R i breakuje.
	
	if(sysarg3 == 1):
		SIGMA_R_plus[T_iter] = sigma_r_plus
		SIGMA_R_minus[T_iter] = sigma_r_minus
		print(E)
		
	if(sysarg3 != 1):

		gamma_r = 2 * np.imag(sigma_r_plus)
		gamma_l = 2 * np.imag(sigma_l_plus)

		G_plus = np.linalg.inv(z * Sd - hamilt - sigma_l_plus - sigma_r_plus)
		G_minus = np.linalg.inv(z_minus * Sd - hamilt - sigma_l_minus - sigma_r_minus)
		T = np.matrix.trace(gamma_l @ G_minus @ gamma_r @ G_plus)
		
		print(E)
		print(T)
		
		T_vals_re[T_iter] = np.real(T)
		T_vals_im[T_iter] = np.imag(T)

	f = f + df
	T_iter = T_iter + 1
	
if(sysarg3 == 1):

	f = f_start
	T_iter = 0

	while(f <= f_koniec):

		E = f

		#petla dla z+

		sigma_plus = np.zeros([len(new_h0), len(new_h0)])
		z = E + 1j * delta
		iteration = 0
		alfa = 0.5
		
		sigma_plus_minus_minus = np.zeros([len(new_h0), len(new_h0)])
		sigma_plus_minus = np.zeros([len(new_h0), len(new_h0)])
		
		while_number = 1
		while(while_number == 1):
			suma = 0
			if(iteration >= 2):
				sigma_plus = alfa * sigma_plus_minus_minus + (1 - alfa) * sigma_plus_minus
			sigma_plus = (h1l - z * S_1) @ (np.linalg.inv(z * S_0 - new_h0 - sigma_plus)) @ (h1l.conj().T - z * S_1.conj().T)
			if(iteration >= 2):
				for i in range(len(new_h0)):
					for j in range(len(new_h0)):
						suma += (sigma_plus[i, j] - sigma_plus_minus[i, j])**2
				if(abs(suma) <= 0.00001):
					while_number = 0

			sigma_plus_minus_minus = sigma_plus_minus
			sigma_plus_minus = sigma_plus
			iteration = iteration + 1
		
		#koniec petli dla z+
		#petla dla z-
	
		sigma_mniej = np.zeros([len(new_h0), len(new_h0)])
		z_minus = E - 1j * delta
		iteration = 0
		
		sigma_mniej_minus_minus = np.zeros([len(new_h0), len(new_h0)])
		sigma_mniej_minus = np.zeros([len(new_h0), len(new_h0)])
		
		while_number = 1
		while(while_number == 1):
			suma = 0
			if(iteration >= 2):
				sigma_mniej = alfa * sigma_mniej_minus_minus + (1 - alfa) * sigma_mniej_minus
			sigma_mniej = (h1l - z_minus * S_1) @ (np.linalg.inv(z_minus * S_0 - new_h0 - sigma_mniej)) @ (h1l.conj().T - z_minus * S_1.conj().T)
			if(iteration >= 2):
				for i in range(len(new_h0)):
					for j in range(len(new_h0)):
						suma += (sigma_mniej[i, j] - sigma_mniej_minus[i, j])**2
				if(abs(suma) <= 0.00001):
					while_number = 0

			sigma_mniej_minus_minus = sigma_mniej_minus
			sigma_mniej_minus = sigma_mniej
			iteration = iteration + 1


		#koniec petli dla z-

		sigma_r_plus = np.zeros([len(hamilt), len(hamilt)], dtype = 'complex_')
		sigma_l_plus = np.zeros([len(hamilt), len(hamilt)], dtype = 'complex_')
		sigma_r_minus = np.zeros([len(hamilt), len(hamilt)], dtype = 'complex_')
		sigma_l_minus = np.zeros([len(hamilt), len(hamilt)], dtype = 'complex_')

		for i in range(len(new_h0)):
			for j in range(len(new_h0)):
				sigma_l_plus[len(hamilt) - len(new_h0) + i, len(hamilt) - len(new_h0) + j] = sigma_plus[i, j]
		for i in range(len(new_h0)):
			for j in range(len(new_h0)):
				sigma_l_minus[len(hamilt) - len(new_h0) + i, len(hamilt) - len(new_h0) + j] = sigma_mniej[i, j]	
				
		gamma_r = 2 * np.imag(SIGMA_R_plus[T_iter])
		gamma_l = 2 * np.imag(sigma_l_plus)

		G_plus = np.linalg.inv(z * Sd - hamilt - sigma_l_plus - SIGMA_R_plus[T_iter])
		G_minus = np.linalg.inv(z_minus * Sd - hamilt - sigma_l_minus - SIGMA_R_minus[T_iter])

		T = np.matrix.trace(gamma_l @ G_minus @ gamma_r @ G_plus)

		print(E)
		print(T)

		T_vals_re[T_iter] = np.real(T)
		T_vals_im[T_iter] = np.imag(T)
	
		f = f + df
		T_iter = T_iter + 1

plt.plot(T_args, T_vals_re, color='black')
plt.xlabel("energy [ev]")
plt.ylabel("conductance [e^2/h]")
plt.savefig('wykresy/transmisja/transmisja.png')
