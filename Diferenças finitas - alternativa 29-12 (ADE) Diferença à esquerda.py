# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 10:29:59 2022

@author: jrfri
"""

import numpy as np
import math as m
import matplotlib.pyplot as plt
import scipy.linalg as lin

# constantes

lx = 20 # Tamanho de X
ly = 80 # Tamanho de Y
dx = 1 # Passo na direção X
dy = 1 # Passo na direção Y
nx = int(lx/dx) # Número de nós em X
ny = int(ly/dy) # Número de nós em Y

u0 = 0.05 # Velocidade inicial (Original: 0.05)
D = 0.07407 # Difusividade (Original: 0.07407)
k = u0*dy*dy/D/dx

CC = 1

dimA = (ny-2)*(nx-2)

# Geração da malha:
grid = np.zeros((nx,ny))


# Criação da matriz de solução do problema
A = np.zeros((dimA,dimA))

# Diagonal principal
for i in range(0,dimA):
    A[i][i]= 2+k # Nó do meio

# Diagonais inferior e superior
for i in range(0,dimA-1):
    A[i+1][i]=-k #inferior (Nó da esquerda)
    A[i][i+1]=0 #superior (Nó da direita)
    
# Remoção dos 1's extras: Ny define o passo (Dimensão da matriz D) e Nx define quantos blocos de matrizes teremos
for i in range(1,nx-2):
    A[i*(ny-2)][(i*(ny-2))-1]=0
    A[i*(ny-2)-1][i*(ny-2)]=0
    
# Adição das matrizes identidade
for i in range(ny-2,(ny-2)*(nx-2)): # Começa da posição do fim do primeiro bloco (Nx-2) até anterior ao último bloco
    A[i][i-ny+2]=-1 # Inferior (Nó de trás)
    A[i-ny+2][i]=-1 # Superior (Nó da frente)
    
# Montagem do vetor das condições de contorno
b = np.zeros(dimA)

for i in range(0,ny-2):
    b[i] = CC
    
# Resolução do sistema linear
h = lin.solve(A,b)

# Aplicando as CC's na grid:
    # CC 1: Topo com C = 1
    
for i in range(0,ny):
    grid[0][i]=CC

# Mapeando h de volta pra grid
conth = 0
for i in range(0,nx-2):
    conty = ny-2
    contx = 0
    while conty>0:
        grid[i+1][contx+1] = h[conth]
        conty = conty-1
        conth = conth+1
        contx = contx+1
        
for i in range(0,nx):
    grid[i][ny-1]= grid[i][ny-2]

tau = (A.dot(h))-b # Vetor dos erros locais 
E = lin.solve(A,-tau) # Vetor dos erros globais

norma_tau = lin.norm(tau)
norma_E = lin.norm(E)

plt.imshow(grid, cmap='jet')
plt.colorbar(orientation="horizontal")
plt.show()






