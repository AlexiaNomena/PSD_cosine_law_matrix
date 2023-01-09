"""
@author: Nomenjanahary Alexia Raharinirina
@date: November 12, 2022
"""
import numpy as np
import scipy as sp
import sys
import pdb

def CosLM(DX, DY, UX = None, UY = None, fXY = None, c = None):
    """
    @brief Compute the cosine law matrix
    Parameters
    ----------
    DX : shape (M, M), Distance set or Proximity Set associated to a set of points X
    DY : shape (N, N), Distance set or Proximity Set associated to a set of points Y
    UX : shape (M,)  , Distance to theoretical origin point for the set of points X
    UY : shape (N,)    ,Distance to theoretical origin point for the set of points Y   
    fXY  : shape (M, N), Proximity set matrix between the points of X and Y, Compatible with the positions of the points in DX and DY
    c  : Dictionary of parameters with keys "c1", "c2", "c3"
    
    Returns
    -------
    CL_Mat : np.array
            Corresponding Cosine Law Matrix with reference at index 0
    """
    # compute f^0
    F0 = Prox_Mat(DX, DY, UX, UY, fXY)
    # compute cos Mat for W associated with f^0
    CL_Mat0 = (F0[0, :][np.newaxis, :]**2 + F0[:, 0][:, np.newaxis]**2 - F0**2)/2
    
    M = DX.shape[0]
    # compute zeta_f
    CC = np.zeros(CL_Mat0.shape)
    if (UX is not None) or (UY is not None):
        # only pick the components involving f(X,Y) 
        CC[:M, M+1:] = CL_Mat0[:M, M+1:] 
        CC[M+1:,:M] = CL_Mat0[M+1: , :M]
    else:
        CC[:M, M:] = CL_Mat0[:M, M:] 
        CC[M:,:M] = CL_Mat0[M:,:M]
    
    #pdb.set_trace()
    Ri = np.sum(np.abs(CC), axis = 1)
    zeta_f = np.max(Ri)
    if zeta_f == 0:
        sys.exit("Distance/Proximity cannot be zero everywhere")
    # insert the proximity values for the theoretical point z and w_1 = x_1
    c1, c2, c3 = c["c1"], c["c2"], c["c3"]
    if c2 == "default":
        c2 = c1
    if c3 == "default":
        c3 = min(2*c2 + 2, c1 + c2 + 2)
    
    if c3<0:
        sys.exit("c3 is negative but c1, c2, c3 must be positive")
    
    
    D = np.zeros((F0.shape[0]+1, F0.shape[0]+1))
    D[0, 1] = np.sqrt(c1*zeta_f)
    D[0, 2:] =  np.sqrt(F0[0, 1:]**2 + c2*zeta_f)
    D[1:, 0] = D[0, 1:]  
    
    # compute f^vareps
    vareps =  c3*zeta_f
    F_vareps = np.sqrt(F0**2 + vareps) - np.diag(np.sqrt(vareps*np.ones(F0.shape[0])))
    #pdb.set_trace()
    D[1:, 1:] = F_vareps

    # Compute cos Mat for V^z associated to f^vareps with referenc z at index 0       
    CL_Mat = (D[0, :][np.newaxis, :]**2 + D[:, 0][:, np.newaxis]**2 - D**2)/2
    
    return CL_Mat, c1, c2, c3, zeta_f


def Prox_Mat(DX, DY, UX = None, UY = None, fXY = None):
    M = DX.shape[0]
    N = DY.shape[0]
    if (fXY is None) or (np.all(np.isclose(fXY, np.zeros((M, N))))) or (np.any(fXY < 0)):
        sys.exit("fXY non-negative and not zero everywhere is needed \n fXY : Proximity set matrix between the points of X and Y compatible with the positions of the points in DX and DY")

    else:
        if UX is not None:
            # put the distance to the origin after all the points of X
            DX1 = np.concatenate((DX, UX[np.newaxis, :]), axis = 0)
            DX1 = np.concatenate((DX1,np.concatenate((UX.T, np.array([0])), axis = 0)[:, np.newaxis]), axis = 1)
    
        if UY is not None:
            # put the distance to the origin after all the points of X and before all the points of Y
            DY1 = np.concatenate((UY[np.newaxis, :], DY), axis = 0)
            DY1 = np.concatenate((np.concatenate((np.array([0]), UY.T), axis = 0)[:, np.newaxis], DY1), axis = 1)
            
        if (UX is not None) or (UY is not None):
            D = np.zeros((M+N+1, M+N+1))
            D[:M, M+1:] = fXY
            D[M+1:, :M] = D[:M, M+1:].T
            
            if (UX is not None)&(UY is None):
                D[:M+1, :M+1] = DX1
                D[M+1:, M+1:] = DY
        
            elif (UY is not None)&(UX is None):
                D[:M, :M] = DX
                D[M:, M:] = DY1
                
            else:
                D[:M+1, :M+1] = DX1
                D[M:, M:] = DY1
        
        else:
            D = np.zeros((M+N, M+N))
            D[ :M, M:] = fXY
            D[M:,  :M] = D[:M, M:].T
            
            D[:M, :M] = DX
            D[M:, M:] = DY
        
        return D