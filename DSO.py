import numpy as np
import gurobipy as gp
from gurobipy import GRB
# python索引下标从0开始
def DSO (Node,Pn,Y,C) :
    Qn=np.zeros(Node)
    Qn[0]=Pn[0]+Pn[3]
    Qn[1]=Pn[1]
    Qn[2]=Pn[2]
    Qn[4]=Pn[4]+Pn[7]
    Qn[6]=Pn[5]+Pn[8]
    Qn[8]=Pn[6]+Pn[9]

    # New op model
    m=gp.Model('DSO')
    # add vars
    Q=m.addVars(Node,Node,name="Q")
    theta = m.addVars(Node, name="theta")
    # define Objective function
    obj = gp.quicksum((Qn[i] - gp.quicksum(Q[i, j] for j in range(Node))) ** 2 for i in range(Node))
    m.setObjective(obj, GRB.MINIMIZE)

    # add s.t.
    m.addConstr(theta[0]==0,name="theta_constraint")

    for i in range(Node):
        for j in range(Node):
            m.addConstr(Q[i, j] == Y[i, j] * (theta[j] - theta[i]), name=f"Q_{i}_{j}" )
            m.addConstr(-C <= Q[i,j], name=f"lb_{i}_{j}")
            m.addConstr(Q[i,j] <= C, name=f"ub_{i}_{j}")
    # op
    m.optimize()

    Q_v = np.array([[Q[i, j].x for j in range(Node)] for i in range(Node)])
    theta_v = np.array([theta[i].x for i in range(Node)])
    dso = np.linalg.norm(Qn - np.sum(Q_v, axis=1), 2)

    return Q_v, theta_v, dso


    #
    # if m.status==GRB.OPTIMAL:
    #     print("DSO op solution")
    #     for i in range(Node):
    #         for j in range(Node):
    #             print(f"Q[{i},{j}] = {Q[i,j].x}")
    #     for i in range(Node):
    #         print(f"theta[{i}]{theta[i].x}")
    #     Q_v = np.array([[Q[i, j].x for j in range(Node)] for i in range(Node)])
    #     theta_v = np.array([theta[i].x for i in range(Node)])
    #     dso = np.linalg.norm(Qn - np.sum(Q_v, axis=1), 2)
    #     m.close()
    #     return Q_v, theta_v, dso
    # else:
    #     print("DSO no solution found")
    #     m.close()
    #     return 0,0,0




