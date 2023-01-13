# -*- coding: utf-8 -*-
"""
Created on Sun Jul 31 19:35:37 2022

@author: Admin
"""
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import math
#----OuterApproximation------
#-----Variables required
'''n =3
m = 2
l =1
e =1
P = np.array([[25,15,-5],[15,18,0],[-5,0,11]])
D = np.array([[2,0,0],[0,2,0],[0,0,2]])
Q = P + D
c = np.array([1,0,3])
b =np.array([4,2,1])
A = np.array([[1,1,0],[4,0,2]])
G = np.array([[3,0,5]])
w = np.array([0,2])
u = np.array([2])
De = D -e*np.identity(3)
P = Q - De
L = np.linalg.cholesky(P)
ye = np.dot(np.linalg.inv(L),b)
UB = 0
zbar = [0.2,0.1,0.35]'''
matrix_size =5
n =5
m = 2
l =1
e = 0.8
D = np.array([[1.2,0,0,0,0],[0,1.2,0,0,0],[0,0,1.2,0,0],[0,0,0,1.2,0],[0,0,0,0,1.2]])
Q = np.array([[25,15,-5,4,1],[15,18,0,2,-2],[-5,0,11,6,2],[4,2,6,8,2],[1,-2,2,2,5]])
De = D -e*np.identity(n)
P = Q - De
print("eigen_val",np.linalg.eigvals(P))
L = np.linalg.cholesky(P)
c = np.array([1,2,1,1,1])
b =np.array([2,3,5,1,3])
A = np.array([[1,2,6,2,5],[3,4,2,6,1]])
G = np.array([[3,1,4,2,7]])
w = np.array([10,6])
u = np.array([2])
L = np.linalg.cholesky(P)
print("L_rank =",np.linalg.matrix_rank(L))
ye = np.dot(np.linalg.inv(L.transpose()),b)
zbar = [0.1,0.2,0.001,0.5,0.5]
Summ = np.identity(n)
itr = 0
#----zeta functions-----


def Zeta(Summ,zbar):
    z = []
    for i in range(n):
        z.append(zbar[i])
    Final = np.identity(n)
    L1 =[0,0,0,0,0]
    for a in range(n):
        L1[a] = [0,0,0,0,0]
        for i in range(n):
            L1[a][i] =[0,0,0,0,0]
            for j in range(n):
                L1[a][i][j] = L[a][i]*L[a][j]*z[i]/De[i][i]
    Summ = np.identity(n)
    for i in range(n):
        Summ = np.add(Summ ,np.array(L1[i]))
    Final = Summ
    print('Zeta =',Final)
    result = np.dot(ye,np.dot(Final,ye.transpose())) - np.dot(ye,ye.transpose()) + np.dot(c.transpose(),z)
    return(result)
def Delta_Zeta(z,zbar):
    z1 ={}
    for i in range(n):
        z1[i] = {}
        for j in range(n):
            if(i == j):
                z1[i][j] = math.sqrt(zbar[i])
            else:
                z1[i][j] = 0
    He = [0,0,0,0,0]
    for i in range(n):
        He[i] = [0,0,0,0,0]
        for j in range(n):
            He[i][j] = 0
            for k in range(n):
                He[i][j] += L[i][k] * z1[k][j]
                
    Het = []
    for i in range(n):
        Het.append(0)
    
    for i in range(n):
        Het[i] = [0,0,0,0,0]
        for j in range(n):
            Het[i][j] = He[j][i]
    bzF = []
    for i in range(n):
        bzF.append(b[i]*math.sqrt(zbar[i]))
    H = np.array(He)
    Ht = np.array(Het)
    b1 = np.array(bzF)
    xe = np.dot(np.linalg.inv((np.dot(H,Ht) + De)),b1)
    Final =0
    for i in range(n):
        Final += (De[i][i]*(xe[i]**2)*(z[i]-zbar[i])/zbar[i])
    return(Final)


#----callBack Function----


def Tangent_Cut(model,where):   
    if(where == GRB.Callback.MIPSOL):
        z = model.cbGetSolution(model._z)
        print("initial_Z",z)
        
        '''print("initial_z =",z)
        for i in range(n):
            z1[i] = {}
            for j in range(n):
                if(i == j):
                    z1[i][j] = math.sqrt(zbar[i])
                else:
                    z1[i][j] = 0
        org_model = gp.Model('original_problem')
        x = org_model.addVars(n,vtype = GRB.CONTINUOUS,name ="X")'''
        if(z[1] == 0):
            print(1)
            model.cbLazy(model._eta >= 100)
        '''for i in range(n):
            z[i] = org_model.addVar(vtype =GRB.INTEGER,ub = zbar[i],lb = zbar[i])
        L1 =[0,0,0,0,0]
        for a in range(n):
            L1[a] = [0,0,0,0,0]
            for i in range(n):
                L1[a][i] =[0,0,0,0,0]
                for j in range(n):
                    L1[a][i][j] = L[a][i]*L[a][j]*z[i]/De[i][i]
        Summ = np.identity(5)
        for i in range(n):
            Summ = np.add(Summ ,np.array(L1[i]))
        print(Summ)
        first =[0,0,0,0,0]
        for a in range(n):
            first[a] =[0,0,0,0,0]
            for i in range(n):
                for j in range(n):
                    first[a][i] += ye[j]*L1[a][j][i]
        final = 0
        for a in range(n):
            for i in range(n):
                final += first[a][i]*ye[i]                
        result = final + gp.quicksum(c[i]*z[i] for i in range(n)) - gp.quicksum(ye[i]*ye[i] for i in range(n))'''
        '''L1 = [0,0,0,0,0]
        for a in range(n):
            L1[a] = [0,0,0,0,0]
            for i in range(n):
                L1[a][i] = [0,0,0,0,0]
                for j in range(n):
                    L1[a][i][j] = L[a][i]*L[a][j]*z[i]/De[i][i]
        Summ = np.identity(n)
        for i in range(n):
            Summ = np.add(Summ ,np.array(L1[i]))
        print(np.linalg.matrix_rank(np.linalg.inv(Summ)))
        print('Optimized =' ,Summ)
        result1 = np.dot(ye,np.dot(np.linalg.inv(Summ),ye.transpose())) - np.dot(ye,ye.transpose()) + gp.quicksum(c[i]*z[i] for i in range(n))
        org_model.setObjective(result1,GRB.MINIMIZE)  
        org_model.setParam(GRB.Param.OutputFlag,0)
        org_model.optimize()
        expr =0
        for i in range(n):
            for j in range(n):
                    expr += 0.5*Q[i][j]*x[i]*x[j]
        for i in range(n):
            expr -= b[i] * x[i]
            expr += c[i] * z[i]
        #----indicator_constraint----
        for i in range(0,n):
            if(z[i] == 0):
                org_model.addConstr(x[i] == 0)
            else:
                pass
        org_model.setObjective(expr,GRB.MINIMIZE)
        print('org_model',org_model.status)
        if(org_model.status != 2):
            print('cut')
            const1 = Zeta(Summ,z)
            const2 = Delta_Zeta(model._z,z)
            model.cbLazy(model._eta >= const1 +const2)
            #result1 = np.dot(ye,np.dot(np.linalg.inv(Summ),ye.transpose())) - np.dot(ye,ye.transpose()) + gp.quicksum(c[i]*z[i] for i in range(n))
        else:
            print("no_cut")
   #np.dot(ye,np.dot(Final,ye.transpose())) - np.dot(ye,ye.transpose())  '''        
        
        
#-----THEOREM-4 Model-----
model = gp.Model("Tangent_cut")
eta = model.addVar(vtype = GRB.CONTINUOUS,name = "N")
model._eta = eta.getVars()
z = model.addVars(n,vtype = GRB.BINARY,name ='Z')
model._z = z.getVars()
#x = model.addVars(n,vtype = GRB.CONTINUOUS,name ="X")
#----Objective constraint----
expr =0
org_z =[0.5,0.5,0.5,0.5,0.5]
'''for i in range(n):
    for j in range(n):
            expr += 0.5*Q[i][j]*x[i]*x[j]
for i in range(n):
    expr -= b[i] * x[i]
    expr += c[i] * z[i]
model.addConstr(eta >= expr)
#----indicator_constraint----
for i in range(0,n):
    model.addConstr((z[i] == 0)>>(x[i] == 0))'''
#----Theorem4 inital tangent cut----
const_part1 = Zeta(Summ,zbar)
const_part2 = Delta_Zeta(z,zbar)
model.addConstr(eta >= const_part1 +const_part2)
#---convexset_constraint----
model.addConstr(gp.quicksum(z[i] for i in range(n)) <= 3)
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
for i in range(n):
    print(z[i].x)

'''
        for j in range(n):
                L1[i][j] = L1[i][j]*z[i]/De[i][i]
'''