import numpy as np
import gurobipy as gp
from gurobipy import GRB
def agentw(ce,cr,i,Agent,Node,nambda,nu,ap,bp,ar,br,Pmin,Pmax,Rmin,Rmax,Pnm,Pmn,Rnm,Rmn,rho,Q,Pn,upsilon) :
    mw = gp.Model('Agentw')
    Pw = mw.addVars(Agent - 1, name='Pw')
    Rw = mw.addVars(Agent - 1, name='Rw')

    if i == 7:
        obj = (0.5 * ap[0][0] * ( gp.quicksum(Pw[j] for j in range(Agent - 1)) ** 2) + bp[0][0] *  gp.quicksum(
            Pw[j] for j in range(Agent - 1)) +  gp.quicksum(ce[j] * Pw[j] for j in range(Agent - 1))
               + 0.5 * ar[0][0] * ( gp.quicksum(Rw[j] for j in range(Agent - 1)) ** 2) + br[0][0] *  gp.quicksum(
                    Rw[j] for j in range(Agent - 1)) + gp. quicksum(cr[j] * Rw[j] for j in range(Agent - 1))
               +  gp.quicksum((nambda[0][j] * ((Pnm[0][j] - Pmn[0][j]) / 2 - Pw[j]) + nu[0][j] * (
                        (Rnm[0][j] - Rmn[0][j]) / 2 - Rw[j])
                           + ((rho[0][j] / 2) * ((Pnm[0][j] - Pmn[0][j]) / 2 - Pw[j]) ** 2)
                           + ((rho[0][j] / 2) * ((Rnm[0][j] - Rmn[0][j]) / 2 - Rw[j]) ** 2)) for j in range(Agent - 1))
               + upsilon[0][4] * (
                       gp.quicksum(Pw[j] for j in range(Agent - 1)) + Pn[4] -  gp.quicksum(Q[4, j] for j in range(Node)))
               + (1 / 2) * (gp.quicksum(Pw[j] for j in range(Agent - 1)) + Pn[4] - gp.quicksum(
                    Q[4, j] for j in range(Node))) ** 2)
        mw.setObjective(obj, GRB.MINIMIZE)

    elif i == 8:
        obj = (0.5 * ap[0][0] * (gp.quicksum(Pw[j] for j in range(Agent - 1)) ** 2) + bp[0][0] * gp.quicksum(
            Pw[j] for j in range(Agent - 1)) + gp.quicksum(ce[j] * Pw[j] for j in range(Agent - 1))
               + 0.5 * ar[0][0] * (gp.quicksum(Rw[j] for j in range(Agent - 1)) ** 2) + br[0][0] * gp.quicksum(
                    Rw[j] for j in range(Agent - 1)) + gp.quicksum(cr[j] * Rw[j] for j in range(Agent - 1))
               + gp.quicksum((nambda[0][j] * ((Pnm[0][j] - Pmn[0][j]) / 2 - Pw[j]) + nu[0][j] * (
                        (Rnm[0][j] - Rmn[0][j]) / 2 - Rw[j])
                           + ((rho[0][j] / 2) * ((Pnm[0][j] - Pmn[0][j]) / 2 - Pw[j]) ** 2)
                           + ((rho[0][j] / 2) * ((Rnm[0][j] - Rmn[0][j]) / 2 - Rw[j]) ** 2)) for j in range(Agent - 1))
               + upsilon[0][6] * (
                       gp.quicksum(Pw[j] for j in range(Agent - 1)) + Pn[5] - gp.quicksum(Q[6, j] for j in range(Node)))
               + (1 / 2) * (gp.quicksum(Pw[j] for j in range(Agent - 1)) + Pn[5] - gp.quicksum(
                    Q[6, j] for j in range(Node))) ** 2)
        mw.setObjective(obj, GRB.MINIMIZE)

    elif i == 9:
        obj = (0.5 * ap[0][0] * (gp.quicksum(Pw[j] for j in range(Agent - 1)) ** 2) + bp[0][0] * gp.quicksum(
            Pw[j] for j in range(Agent - 1)) + gp.quicksum(ce[j] * Pw[j] for j in range(Agent - 1))
               + 0.5 * ar[0][0] * (gp.quicksum(Rw[j] for j in range(Agent - 1)) ** 2) + br[0][0] * gp.quicksum(
                    Rw[j] for j in range(Agent - 1)) + gp.quicksum(cr[j] * Rw[j] for j in range(Agent - 1))
               + gp.quicksum((nambda[0][j] * ((Pnm[0][j] - Pmn[0][j]) / 2 - Pw[j]) + nu[0][j] * (
                        (Rnm[0][j] - Rmn[0][j]) / 2 - Rw[j])
                           + ((rho[0][j] / 2) * ((Pnm[0][j] - Pmn[0][j]) / 2 - Pw[j]) ** 2)
                           + ((rho[0][j] / 2) * ((Rnm[0][j] - Rmn[0][j]) / 2 - Rw[j]) ** 2)) for j in range(Agent - 1))
               + upsilon[0][8] * (
                       gp.quicksum(Pw[j] for j in range(Agent - 1)) + Pn[6] - gp.quicksum(Q[8, j] for j in range(Node)))
               + (1 / 2) * (gp.quicksum(Pw[j] for j in range(Agent - 1)) + Pn[6] - gp.quicksum(
                    Q[8, j] for j in range(Node))) ** 2)
        mw.setObjective(obj, GRB.MINIMIZE)

    mw.addConstr(gp.quicksum(Pw[j] for j in range(Agent - 1)) >= Pmin, name='Pwmin')
    mw.addConstr(gp.quicksum(Pw[j] for j in range(Agent - 1)) <= Pmax, name='Pwmax')
    mw.addConstr(gp.quicksum(Rw[j] for j in range(Agent - 1)) >= Rmin, name='Rwmin')
    mw.addConstr(gp.quicksum(Rw[j] for j in range(Agent - 1)) <= Rmax, name='Rwmax')
    for j in range(Agent - 1):
        mw.addConstr(Pw[j] >= 0)
        mw.addConstr(Rw[j] <= 0)
    mw.addConstr(gp.quicksum(Pw[j] + Rw[j] for j in range(Agent - 1)) <= Pmax, name='Pw+Rwmax')


    mw.optimize()

    if mw.status == GRB.OPTIMAL:
        print("Agentw op solution")
        Pw = np.array([Pw[j].x for j in range(Node)])
        Rw = np.array([Rw[j].x for j in range(Node)])
        return Pw, Rw
    else:
        print("Agentw no solution found")
