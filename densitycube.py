#code written by Too Hon Lin
import numpy as np
import math
import random
from numpy import matlib

#gendcube2 generates random density cube of two-level system
def gendcube2():
    rho1 = random.uniform(0,1)
    rho2 = random.uniform(0,1)
    rho3 = random.uniform(0,1)
    mat1 = [[rho1,rho2],[rho2,rho3]]
    mat2 = [[rho2,rho3],[rho3,1-rho1]]
    dcube = [mat1,mat2]
    return dcube

#qm2dc maps the quantum matrix to its respective density cube. 
#This function can only accept either 2 or 3 level system
def qm2dc(mat):
    n = len(mat)
    null1 = np.zeros((n,n),dtype=complex)
    null2 = np.zeros((n,n),dtype=complex)
    null3 = np.zeros((n,n),dtype=complex)
    if n == 3:
        dcube = [null1, null2, null3]
    elif n == 2:
        dcube = [null1, null2]
    #hermiticity and normalisation condition
    for i in range(0,n):
        dcube[i][i][i] = mat[i][i]
        for j in range(0,n):
            if (i<j):
                z1 = mat[i][j].imag
                dcube[i][j][j] = math.sqrt(2/3)*z1
                dcube[j][i][j] = dcube[i][j][j]
                dcube[j][j][i] = dcube[i][j][j]
                z2 = mat[i][j].real
                dcube[i][i][j] = math.sqrt(2/3)*z2
                dcube[i][j][i] = dcube[i][i][j]
                dcube[j][i][i] = dcube[i][i][j]
    return dcube

#generate the 4 Pauli cubes
def pauli():
    mat1 = np.mat('1 0;0 1',dtype=complex)
    mat2 = np.mat('0 1;1 0',dtype=complex)
    mat3 = np.mat('0 0-1j;1j 0',dtype=complex)
    mat4 = np.mat('1 0;0 -1',dtype=complex)
    pauli1 = qm2dc(mat1) #first Pauli cube
    pauli2 = math.sqrt(2/3)*qm2dc(mat2) #second Pauli cube
    pauli3 = math.sqrt(2/3)*qm2dc(mat3) #third Pauli cube
    pauli4 = qm2dc(mat4) #4th Pauli cube
    parray = [pauli1, pauli2, pauli3, pauli4] #return array with 4 pauli matrices
    return parray

#inner product
def dot(dcube1,dcube2):
    n = len(dcube1)
    conj = np.conjugate(dcube1)
    inprod = 0
    for i in range(0,n):
        for j in range(0,n):
            for k in range(0,n):
                inprod = inprod + conj[i][j][k]*dcube2[i][j][k]
    return inprod.real

#generate conjugate transpose of the density cube
def getH(dcube):
    n = len(dcube)
    mat1 = np.zeros((n,n),dtype=complex)
    mat2 = np.zeros((n,n),dtype=complex)
    mat3 = np.zeros((n,n),dtype=complex)
    conj = []
    if n == 2:
        conj = [mat1,mat2]
    elif n==3:
        conj = [mat1,mat2,mat3]
    for i in range(0,n):
        for j in range(0,n):
            for k in range(0,n):
                if i==j==k:
                    conj[i][j][k] = dcube[i][j][k]
                elif i==j or i==k or j==k:
                    conj[i][j][k] = dcube[i][k][j]
                else:
                    conj[i][j][k] = dcube[j][i][k].conjugate
    return(conj)

#check if the density cube is hermitian
def ishermit(dcube):
    conj = getH(dcube)
    for i in range(0,n):
        for j in range(0,n):
            for k in range(0,n):
                if conj[i][j][k]!=dcube[i][j][k]:
                    return False
    return True

#ispure takes a state density cube and check if the density cube is pure or mixed
def ispure(dcube):
    inprod = dot(dcube,dcube)
    if inprod == 1:
        return True
    else:
        return False

#isorth takes 2 different state density cubes and check if the density cubes are othogonal
def isorth(dcube1,dcube2):
    inprod = dot(dcube1,dcube2)
    if inprod == 0:
        return True
    else:
        return False
