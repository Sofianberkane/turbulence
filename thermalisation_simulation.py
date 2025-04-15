import numpy as np
import matplotlib.pyplot as plt
import random
import math
from scipy.optimize import curve_fit

# Paramètres
Lx = 100  # Taille du domaine en x
Ly = 100  # Taille du domaine en y
Nx = 512  # Nombre de points en x
Ny = 256  # Nombre de points en y
deltak = 0.25  # Constante dans l'exponentielle initiale
dt=0.001
I0=np.sqrt(10)
tmax = 4
Nt = int(tmax/dt)

# Grille spatiale
x0 = np.linspace(-Lx / 2, Lx / 2, Nx, endpoint=False)
y0 = np.linspace(-Ly / 2, Ly / 2, Ny, endpoint=False)
x, y = np.meshgrid(x0, y0)

# Grille de fréquences pour le calcul de la transformée de Fourier
kx0 = 2 * np.pi * np.fft.fftfreq(Nx, d=Lx / Nx)
ky0 = 2 * np.pi * np.fft.fftfreq(Ny, d=Ly / Ny)
kx, ky = np.meshgrid(kx0, ky0)
k_squared = kx**2 + ky**2
k=np.sqrt(k_squared)

#%%
def compute_n0(psi):
    psi_hat = np.fft.fft2(psi)
    return np.abs(psi_hat[0, 0]) / np.sum(np.abs(psi_hat))

n0 =[]
time = []
#%% condition initial exp

def IC(k):
    return np.exp(-(k)/deltak**2)+0j

random_phi = np.random.uniform(0, 2 * np.pi, size=kx.shape)

psi_hat=np.sqrt(IC(k_squared))*np.exp(1j*random_phi)


#%% condition initial porte

def IC(nx, ny, a_frac, b_frac):
    x = np.linspace(-1, 1, nx)
    y = np.linspace(-1, 1, ny)
    X, Y = np.meshgrid(x, y)

    # Taille réelle du rectangle dans l'espace normalisé [-1,1]
    a = a_frac
    b = b_frac

    masque = np.where((np.abs(X) <= a / 2) & (np.abs(Y) <= b / 2), 1.0, 0.0)
    return masque

ny,nx = k.shape
masque = IC(nx,ny,1,1)

random_phi = np.random.uniform(0, 2 * np.pi, size=kx.shape)
psi_hat=np.sqrt(masque*k_squared)*np.exp(1j*random_phi)



#%% création de psi avec la CI précédente

psi_keep=np.fft.ifft2(psi_hat)

# plt.imshow(np.abs(psi), extent=(-Lx / 2, Lx / 2, -Ly / 2, Ly / 2), origin="lower")
# plt.colorbar(label="f")
# plt.xlabel("x")
# plt.ylabel("y")
# plt.show()


# CI normalisé

psi_keep=I0*psi_keep/np.sqrt(np.mean(np.abs(psi_keep)**2))

#%% calcul split step fonction random 2D et 1D
psiL=np.exp(-1j*k_squared*dt/2)

def psiN(psi):
    return np.exp(-1j*np.abs(psi)**2*dt/2)


for t in range(Nt):
    psi_keep*=psiN(psi_keep)
    
    psi_hat_keep=np.fft.fft2(psi_keep)
    
    
    psi_hat_keep*=psiL
    
    
    psi_keep=np.fft.ifft2(psi_hat_keep)
    psi_keep*=psiN(psi_keep)
    
    n0.append(compute_n0(psi_keep))
    time.append(t)
    
    psi_hat_keep=np.fft.fft2(psi_keep)
    
    # if t % 100 == 0:
    #     fig, axs = plt.subplots(1, 2,figsize=(15, 8))
    #     # Afficher l'image de abs(psi)**2 dans axs[0]
    #     cax = axs[ 0].imshow(np.abs(psi_keep), extent=(-Lx / 2, Lx / 2, -Ly / 2, Ly / 2), origin="lower")
    #     axs[0].set_xlabel('x', fontsize=18)  # Label pour l'axe X
    #     axs[0].set_ylabel('y', fontsize=18)  # Label pour l'axe Y
        
        
    #     # Tracer abs(psi_mid)**2 dans axs[1]
    #     axs[1].plot(kx0, np.abs(psi_hat_keep[0])**2)
    #     axs[1].set_xlabel("kx", fontsize=18)
    #     axs[1].set_ylabel("Energie de psi", fontsize=18)
    #     # axs[1].set_xscale('log')
    #     # axs[1].set_yscale('log')

#%%
plt.plot(time,n0)
#%%
psi = psi_keep
psi_hat = psi_hat_keep
#%%

def spectre_theo(k, T, mu):
    return T/(k**2+np.abs(mu))

psi_hat_sum = np.zeros_like(psi_hat)
n=0

tmax_plot = 3
Nt_plot = int(tmax_plot/dt)


for t in range(Nt_plot):
    psi*=psiN(psi)
    
    psi_hat=np.fft.fft2(psi)
    
    
    psi_hat*=psiL
    
    
    psi=np.fft.ifft2(psi_hat)
    psi*=psiN(psi)
    
    psi_hat=np.fft.fft2(psi) 
    
    
        
    if t % 100 == 0:
        
        # psi_hat_sum += psi_hat
        # n+=1
        # psi_hat_mean = np.abs(psi_hat_sum/n)
        
        
        params, covariance = curve_fit(spectre_theo, kx0[100:], np.abs(psi_hat[0,100:]))
        T_fit, mu_fit = params
        ST = spectre_theo(kx0[100:], T_fit, mu_fit)
        
        
        
        fig, axs = plt.subplots(2, 2,figsize=(15, 8))
        # Afficher l'image de abs(psi)**2 dans axs[0]
        cax = axs[0,0].imshow(np.abs(psi)**2, extent=(-Lx / 2, Lx / 2, -Ly / 2, Ly / 2), origin="lower")
        axs[0,0].set_xlabel('x')  # Label pour l'axe X
        axs[0,0].set_ylabel('y')  # Label pour l'axe Y
        fig.colorbar(cax, ax=axs[0, 0])  # Ajouter une barre de couleur
        axs[0, 0].set_title("Abs(psi)^2")
        
        
        # Tracer abs(psi_mid)**2 dans axs[1]
        axs[0,1].plot(x0, np.abs(psi[0]))
        axs[0,1].set_title("Abs(psi_mid)^2")
        axs[0,1].set_xlabel("x")
        axs[0,1].set_ylabel("Abs(psi_mid)")
        
        
        # Tracer abs(psihat_mid)**2 dans axs[2] (donc plt)
        axs[0,1].plot(kx0 , np.abs(psi_hat[0,:]),'b', label='exp')
        axs[0,1].plot(kx0[100:] , ST, 'r', label='theo')
        axs[0,1].set_title("Abs(psi_hat_mid)")
        axs[0,1].set_xlabel("kx")
        
        axs[0,1].set_ylabel("Abs(psi_hat_mid)")
        axs[0,1].set_xscale('log')
        axs[0,1].set_yscale('log')
        
        
        
        cax = axs[1, 1].imshow(np.abs(np.fft.ifftshift(psi_hat)), extent=(kx[0][0], kx[-1][-1], ky[0][0], ky[-1][-1]), origin="lower")
        axs[1,1].set_xlabel('kx')  # Label pour l'axe X
        axs[1,1].set_ylabel('ky')  # Label pour l'axe Y
        fig.colorbar(cax, ax=axs[1, 1])  # Ajouter une barre de couleur
        axs[1, 1].set_title("Abs(psi_hat)")
        
        
        # Mettre à jour les titres et labels
        fig.suptitle(f"Time {(t+1000)/1000}")  # Exemple de titre global
        plt.tight_layout()  # Ajuste les positions des subplots pour éviter les chevauchements
        plt.subplots_adjust(hspace=0.6, wspace=0.4)  # Ajuster l'espacement entre les subgraphiques
        
        
        # Affichage final
        plt.show()
        

