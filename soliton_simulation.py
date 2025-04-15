import numpy as np
import matplotlib.pyplot as plt
import math

# Paramètres
Lx = 120  # Taille du domaine en x
Ly = 120  # Taille du domaine en y
Nx = 600  # Nombre de points en x
Ny = 600  # Nombre de points en y
alpha = 1  # Constante dans l'exponentielle initiale
lamb=800*10**-9
v=10
dt=0.01

# Grille spatiale
x = np.linspace(-Lx / 2, Lx / 2, Nx, endpoint=False)
y = np.linspace(-Ly / 2, Ly / 2, Ny, endpoint=False)
x, y = np.meshgrid(x, y)

# Grille de fréquences pour le calcul de la transformée de Fourier
kx = 2 * np.pi * np.fft.fftfreq(Nx, d=Lx / Nx)
ky = 2 * np.pi * np.fft.fftfreq(Ny, d=Ly / Ny)
kx, ky = np.meshgrid(kx, ky)
k_squared = kx**2 + ky**2


#%% condition initial

def IC_sol(x):
    return np.exp(1j*v*x)/np.cosh(x)


#%% Méthode Fresnel équation linéaire

ligne,col=np.shape(x)
ligne =int(ligne/2)

# Ligne centrale des x et k
x_mid = x[ligne]
kx_mid=kx[ligne]

psi = IC_sol(x_mid)

psiL=np.exp(-1j*kx_mid**2*dt/2)
psiN=lambda psi:np.exp(1j*np.abs(psi)**2*dt/2)

# Initialiser la figure pour superposer les courbes
plt.figure(figsize=(10, 6))

for t in range(1000):
    psi*=psiN(psi)
    psi_hat=np.fft.fft(psi)
    psi_hat*=psiL
    psi=np.fft.ifft(psi_hat)
    psi*=psiN(psi)
    
    if t % 100 == 0:
        plt.plot(x_mid,np.abs(psi)**2, 'c')
plt.xlabel("x")
plt.ylabel(r"$|\psi(x)|^2$")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()