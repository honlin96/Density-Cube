import numpy as np
import math

def dcube2(mat1,mat2):
    dcube = np.array([mat1,mat2])
    #normalise the cube
    return dcube

def dcube3(mat1,mat2,mat3):
    dcube = np.array([mat1,mat2,mat3])
    #normalise the cube
    return dcube

#generate the 4 Pauli cubes
def pauli():
    mat1 = np.mat('1 0;0 0',dtype=complex)
    mat2 = np.mat('0 0;0 1',dtype=complex)
    mat3 = np.mat('0 1;1 0',dtype=complex)
    mat4 = np.mat('1 0;0 0',dtype=complex)
    mat5 = np.mat('0 0;0 1',dtype=complex)
    mat6 = np.mat('0 1;1 0',dtype=complex)
    mat7 = np.mat('1 0;0 0',dtype=complex)
    mat8 = np.mat('0 0;0 -1',dtype=complex)
    pauli1 = dcube2(mat1,mat2) #first Pauli cube
    pauli2 = math.sqrt(2/3)*dcube2(mat3,mat4) #second Pauli cube
    pauli3 = math.sqrt(2/3)*dcube2(mat5,mat6) #third Pauli cube
    pauli4 = dcube2(mat7,mat8) #4th Pauli cube
    dcube = [pauli1, pauli2, pauli3, pauli4]
    return dcube

#gendcube2 generates random density cube of two-level system
def gendcube2():
    #insert the normalisation condition here
    mat1 = np.matlib.zeros(2)
    mat2 = np.matlib.zeros(2)
    dcube2 = []
    dcube2 = np.array([mat1,mat2])
    return dcube2

#getH returns the conjugate transpose of the density cube, dcube of dimension n.
def getH(dcube,n):
    for i in range(0,n):
        mat = np.asmatrix(dcube[i])
        dcube[i] = mat.getH()
    return dcube

#trace returns the trace of the density cube
def trace(dcube,n):
    global tr
    tr = 0 
    for i in range(0,n):
        mat = np.asmatrix(dcube[i])
        tr = tr + mat[i][i]
    return tr

#ispure takes dcube and n as the input and check if the density cube's is pure or mixed
def ispure(dcube,n):
    #perform inner product on itself
    global tr
    tr = trace(dcube,n)
    if tr == 1:
        return True
    else:
        return False
