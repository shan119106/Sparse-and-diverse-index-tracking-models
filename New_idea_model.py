# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 11:47:37 2022

@author: Admin
"""
#importing_data
import random
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import linalg as LA
from scipy.spatial import distance
from sklearn.preprocessing import normalize
from sklearn.cluster import KMeans
import math

#lambda values 
l1 = 0.000000001
l2 = 1.0
l3 = 10000000.0
#input_values_from_Excel
result = pd.read_excel("sp500-monthly-return-10yr.xlsx",usecols = 'B')
#size_of_the_input
n = len(result)-2
#number_of_Assets
m =502
#introducing_the_sparcit_co-efficient
A = [0]*m


#pre-processing_of_portfolio_from_excel
#converting_The_input_into_nparray
y =[]
for i in range(2,n+2):
    y.append(result['S&P 500'][i])
y = np.array(y)
Input = pd.read_excel("sp500-monthly-return-10yr.xlsx",usecols = 'C:SJ')
assets = list(Input)
j =0
X = []
for i in range(n):
    X.append([])
for a in assets:
    for i in range(2,n+2):
        if(pd.isna(Input[a][i])):
            X[i-2].append(0)
        else:
            X[i-2].append(Input[a][i])  
        j=j+1
X = np.array(X)


#objective_function_variables
Q = np.dot(X.transpose(),X)
#bvariable in the second function_of_objective
b = np.dot(y.transpose(),X) 


#converting Q into Positive definite


def nearestPD(A):
    """Find the nearest positive-definite matrix to input

    A Python/Numpy port of John D'Errico's `nearestSPD` MATLAB code [1], which
    credits [2].

    [1] https://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd

    [2] N.J. Higham, "Computing a nearest symmetric positive semidefinite
    matrix" (1988): https://doi.org/10.1016/0024-3795(88)90223-6
    """

    B = (A + A.T) / 2
    _, s, V = LA.svd(B)

    H = np.dot(V.T, np.dot(np.diag(s), V))

    A2 = (B + H) / 2

    A3 = (A2 + A2.T) / 2

    if isPD(A3):
        return A3

    spacing = np.spacing(LA.norm(A))
    # The above is different from [1]. It appears that MATLAB's `chol` Cholesky
    # decomposition will accept matrixes with exactly 0-eigenvalue, whereas
    # Numpy's will not. So where [1] uses `eps(mineig)` (where `eps` is Matlab
    # for `np.spacing`), we use the above definition. CAVEAT: our `spacing`
    # will be much larger than [1]'s `eps(mineig)`, since `mineig` is usually on
    # the order of 1e-16, and `eps(1e-16)` is on the order of 1e-34, whereas
    # `spacing` will, for Gaussian random matrixes of small dimension, be on
    # othe order of 1e-16. In practice, both ways converge, as the unit test
    # below suggests.
    I = np.eye(A.shape[0])
    k = 1
    while not isPD(A3):
        mineig = np.min(np.real(LA.eigvals(A3)))
        A3 += I * (-mineig * k**2 + spacing)
        k += 1

    return A3


def isPD(B):
    """Returns true when input is positive-definite, via Cholesky"""
    try:
        _ = LA.cholesky(B)
        return True
    except LA.LinAlgError:
        return False


if __name__ == '__main__':
    import numpy as np
    for i in range(10):
        for j in range(2, 100):
            A = np.random.randn(j, j)
            B = nearestPD(Q)
            assert (isPD(B))
    print('unit test passed!')
    


EXP = 0
#pre_processing_data_definition
e = 0.5
D = np.identity(m)
De = D -e*np.identity(m)
P = Q -De
L = np.linalg.cholesky(B)
print(L)
ye = np.dot(np.linalg.inv(L.transpose()),b)
zbar = np.random.rand(m)
Summ = np.identity(m)

#spectralclustering

 
k = int(m/3)#number of clusters
D=[]#weight matrix for each asset co-relation
first = X.transpose()
for i in range(m):
    D.append([])
    for j in range(m):
        if(i == j):
            D[i].append(0)
        else:
            D[i].append(distance.sqeuclidean(first[i],first[j]))          
D= np.array(D)
#matrix with sum of weights of each asset
K = D.sum(axis = 1)
#laplacian norm
K = np.sqrt(1/K)
#matrix for laplacian norm
M = np.multiply(K[np.newaxis,:],np.multiply(D,K[:,np.newaxis]))
#EigenValue decomposition
U,Sigma,_ = LA.svd(M,full_matrices = False,lapack_driver='gesvd')
Usubset = U[:,0:k]
#spectral clustering
y_pred_sc = KMeans(n_clusters =k).fit_predict(normalize(Usubset))
#to_seperate_clusters_into_grps
Cluster ={}#cluster_dict
#cluster_dict_making
for i in y_pred_sc:
    if i not in Cluster.keys():
        Cluster[i] =[]
#length_of_kmeans        
cluster_length = len(y_pred_sc)
#clusters_dict_completion
for i in range(cluster_length):
    Cluster[y_pred_sc[i]].append(i)
#definition of A
A = [0]*m
for i in Cluster.keys():
    for j in Cluster[i]:
        if(j%3 == 0):
            A[j] = l1
        elif(j%3 == 1):
            A[j] = l2
        else:
            A[j] = l3
            
A =np.array(A)     
c = A
print('PD checking for A',isPD(np.diag(A)))

      
#PROBLEM_MATHEMATICAL_MODEL
def Zeta(Summ,zbar):
    z = []
    for i in range(m):
        z.append(zbar[i])
    Final = np.identity(m)
    for a in range(m):
        x = 'L1'+ str(a)
        x = [0]*m
        for i in range(m):
            x[i] =[0]*m
            for j in range(m):
                x[i][j] = L[a][i]*L[a][j]
                x[i][j] = x[i][j]*z[i]/De[i][i]
        Summ = Summ + np.array(x)
    Final = Summ
    print('Zeta =',Final)
    result = np.dot(ye,np.dot(Final,ye.transpose())) - np.dot(ye,ye.transpose()) + np.dot(c.transpose(),z)
    return(result)
def Delta_Zeta(z,zbar):
    z1 ={}
    for i in range(m):
        z1[i] = {}
        for j in range(m):
            if(i == j):
                z1[i][j] = math.sqrt(zbar[i])
            else:
                z1[i][j] = 0
    He = [0]*m
    for i in range(m):
        He[i] = np.zeros(m)
        for j in range(m):
            He[i][j] = 0
            for k in range(n):
                He[i][j] += L[i][k] * z1[k][j]
                
    Het = []
    for i in range(m):
        Het.append(0)
    
    for i in range(m):
        Het[i] = [0]*m
        for j in range(n):
            Het[i][j] = He[j][i]
    bzF = []
    for i in range(m):
        bzF.append(b[i]*math.sqrt(zbar[i]))
    H = np.array(He)
    Ht = np.array(Het)
    b1 = np.array(bzF)
    xe = np.dot(np.linalg.inv((np.dot(H,Ht) + De)),b1)
    Final =0
    for i in range(m):
        Final += (De[i][i]*(xe[i]**2)*(z[i]-zbar[i])/zbar[i])
    return(Final)


#----callBack Function----


def Tangent_Cut(model,where):   
    if(where == GRB.Callback.MIPSOL):
        z = model.cbGetSolution(model._z)
        print(z)
        org_model = gp.Model('original_problem')  
        x = x = org_model.addVars(m,vtype = GRB.CONTINUOUS,ub =1,lb =0,name ='x')
        '''Summ = np.identity(m)
        for a in range(m):
            x = 'L1'+ str(a)
            x = [0]*m
            for i in range(m):
                x[i] =[0]*m
                for j in range(m):
                    x[i][j] = L[a][i]*L[a][j]
                    x[i][j] = x[i][j]*z[i]/De[i][i]
            Summ = Summ + np.array(x)
        print(np.linalg.matrix_rank(np.linalg.inv(Summ)))
        print('Optimized =' ,Summ)
        result1 = np.dot(ye,np.dot(np.linalg.inv(Summ),ye.transpose())) - np.dot(ye,ye.transpose()) + gp.quicksum(c[i]*z[i] for i in range(n))
        org_model.setObjective(result1,GRB.MINIMIZE)  
        org_model.setParam(GRB.Param.OutputFlag,0)'''
        expr =0
        for i in range(m):
            for j in range(m):
                    expr += 0.5*Q[i][j]*x[i]*x[j]
        for i in range(m):
            expr -= b[i] * x[i]
            expr += c[i] * z[i]
        for i in range(0,m):
            model.addConstr((z[i] == 0)>>(x[i] == 0))
        org_model.setObjective(eta,GRB.MINIMIZE)
        org_model.optimize()
        print('org_model',org_model.status)
        if(org_model.status != 2):
            print('cut')
            const1 = Zeta(Summ,z)
            const2 = Delta_Zeta(model._z,z)
            model.cbLazy(model._eta >= const1 +const2)
            #result1 = np.dot(ye,np.dot(np.linalg.inv(Summ),ye.transpose())) - np.dot(ye,ye.transpose()) + gp.quicksum(c[i]*z[i] for i in range(n))
        else:
            print("no_cut")
   #np.dot(ye,np.dot(Final,ye.transpose())) - np.dot(ye,ye.transpose())      
        
        
#-----THEOREM-4 Model-----
model = gp.Model("Tangent_cut")
eta = model.addVar(vtype = GRB.CONTINUOUS,name = "N")

z = model.addVars(m,vtype = GRB.BINARY,name ='Z')
model._z = z
x = model.addVars(m,vtype = GRB.CONTINUOUS,ub =1,lb =0,name ='x')

#x = model.addVars(n,vtype = GRB.CONTINUOUS,name ="X")
#----Objective constraint----
expr =0
org_z = 0.5*np.ones(m)
for i in range(m):
    for j in range(m):
            expr += 0.5*Q[i][j]*x[i]*x[j]
for i in range(m):
    expr -= b[i] * x[i]
    expr += c[i] * z[i]
model.addConstr(eta >= expr)
#----indicator_constraint----
for i in range(0,m):
    model.addConstr((z[i] == 0)>>(x[i] == 0))
#----Theorem4 inital tangent cut----
const_part1 = Zeta(Summ,zbar)
const_part2 = Delta_Zeta(z,zbar)
model.addConstr(eta >= const_part1 +const_part2)
#---convexset_constraint----
model.addConstr(gp.quicksum(z[i] for i in range(m)) <= m)
#objective
model.setObjective(eta,GRB.MINIMIZE)
model.Params.DualReductions =0
model.write("model.lp")
model.Params.lazyConstraints = 1
model.setParam(GRB.Param.OutputFlag, 0)
model.optimize(Tangent_Cut)
print("model_status =",model.status)
model_obj = model.getObjective()
print("model_output =",model_obj.getValue())
for i in range(m):
    print(z[i].x,end =',')


      

#OLDMODEL
problem = gp.Model("Index_tracking_problem")
w = problem.addMVar(m,vtype = GRB.CONTINUOUS,lb =0,name="weight")
z = problem.addMVar(m,vtype =GRB.BINARY,lb = 0,name ='z')
for i in range(m):
    EXP += w[i]
    problem.addConstr((z.tolist()[i]==0) >> (w.tolist()[i] ==0))
problem.addConstr(EXP == 1)
a =0
problem.setObjective((w @ B @ w) - b @ w +  z @ np.diag(A) @ z,GRB.MINIMIZE)
problem.optimize()
print(problem.status)
obj = problem.getObjective()
print(obj.getValue)
for i in range(m):
    print("{"+str(i) +":" + str(w[i].x) + '}',end =",")
    


