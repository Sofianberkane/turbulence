# Simulation de l'Équation de Schrödinger Non Linéaire

Ce dépôt contient deux scripts Python illustrant la résolution numérique de l'équation de Schrödinger non linéaire (NLS) dans deux contextes différents : la propagation d’un **soliton**, et la **thermalisation d’un champ d’onde aléatoire**.

---

## 📁 Contenu

### 1. `soliton_simulation.py`
Ce script simule la propagation d'un **soliton** à l’aide de la méthode du **split-step Fourier** (méthode de Fresnel). On utilise une condition initiale de la forme :

\[
\psi(x, t=0) = \text{sech}(x) \cdot e^{i v x}
\]

On observe une onde localisée qui se propage sans se déformer, typique du comportement des solitons. Plusieurs courbes \( |\psi(x)|^2 \) sont superposées au fil du temps pour illustrer ce mouvement.

### 2. `thermalisation_simulation.py`
Ce script explore la **thermalisation d’un champ d’onde** injecté aléatoirement dans l’espace de Fourier. Il inclut :
- Une initialisation spectrale avec bruit de phase
- L'évolution du système via split-step
- Le suivi de la fraction condensée \( n_0 \)
- L'ajustement du spectre \( |\hat{\psi}(k)| \) par la loi :

\[
n_k^{eq} = \frac{T}{k^2 + |\mu|}
\]

Des visualisations en temps réel du module carré, de son spectre, et de la régression sont proposées.

---

## ⚙️ Dépendances

Les deux scripts nécessitent :

- Python ≥ 3.7
- numpy
- matplotlib
- scipy

Installation rapide :
```bash
pip install numpy matplotlib scipy
