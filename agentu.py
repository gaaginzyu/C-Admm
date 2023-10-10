import numpy as np
import gurobipy as gp
from gurobipy import GRB
def agentu(ce,cr,i,Agent,Node,nambda,nu,ap,bp,ar,br,Pmin,Pmax,Rmin,Rmax,Pnm,Pmn,Rnm,Rmn,rho,Q,Pn,upsilon) :
        mu = gp.Model('Agentu')
        Pu = mu.addVars( Agent - 1, name='Pu')
        Ru = mu.addVars( Agent - 1, name='Ru')
        if i == 3:
            obj = (0.5 * ap[0][0] * (gp.quicksum(Pu[j] for j in range(Agent - 1)) ** 2) + bp[0][0] * gp.quicksum(Pu[j] for j in range(Agent - 1)) + gp.quicksum(ce[j] * Pu[j] for j in range(Agent - 1))
               + 0.5 * ar[0][0] * (gp.quicksum(Ru[j] for j in range(Agent - 1)) ** 2) + br[0][0] * gp.quicksum(Ru[j] for j in range(Agent - 1)) +gp. quicksum(cr[j] * Ru[j] for j in range(Agent - 1))
               + gp.quicksum((nambda[0][j] * ((Pnm[0][j] - Pmn[0][j]) / 2 - Pu[j]) + nu[0][j] * ((Rnm[0][j] - Rmn[0][j]) / 2 - Ru[j])
                + ((rho[0][j] / 2) * ((Pnm[0][j] - Pmn[0][j]) / 2 - Pu[j]) ** 2)
                + ((rho[0][j] / 2) * ((Rnm[0][j] - Rmn[0][j]) / 2 - Ru[j]) ** 2)) for j in range(Agent - 1))
               + upsilon[0][0] * (gp.quicksum(Pu[j] for j in range(Agent - 1)) + Pn[0] - gp.quicksum(Q[0, j] for j in range(Node)))
               + (1 / 2) * (gp.quicksum(Pu[j] for j in range(Agent - 1)) + Pn[0] - gp.quicksum(Q[0, j] for j in range(Node))) ** 2)
            mu.setObjective(obj, GRB.MINIMIZE)
        elif i == 4:
            obj = (0.5 * ap[0][0] * (gp.quicksum(Pu[j] for j in range(Agent - 1)) ** 2) + bp[0][0] * gp.quicksum(
                Pu[j] for j in range(Agent - 1)) + gp.quicksum(ce[j] * Pu[j] for j in range(Agent - 1))
                   + 0.5 * ar[0][0] * (gp.quicksum(Ru[j] for j in range(Agent - 1)) ** 2) + br[0][0] * gp.quicksum(
                        Ru[j] for j in range(Agent - 1)) + gp.quicksum(cr[j] * Ru[j] for j in range(Agent - 1))
                   + gp.quicksum((nambda[0][j] * ((Pnm[0][j] - Pmn[0][j]) / 2 - Pu[j]) + nu[0][j] * (
                            (Rnm[0][j] - Rmn[0][j]) / 2 - Ru[j])
                               + ((rho[0][j] / 2) * ((Pnm[0][j] - Pmn[0][j]) / 2 - Pu[j]) ** 2)
                               + ((rho[0][j] / 2) * ((Rnm[0][j] - Rmn[0][j]) / 2 - Ru[j]) ** 2)) for j in range(Agent - 1))
                   + upsilon[0][4] * (
                           gp.quicksum(Pu[j] for j in range(Agent - 1)) + Pn[7] - gp.quicksum(Q[4, j] for j in range(Node)))
                   + (1 / 2) * (gp.quicksum(Pu[j] for j in range(Agent - 1)) + Pn[7] - gp.quicksum(
                        Q[4, j] for j in range(Node))) ** 2)
            mu.setObjective(obj, GRB.MINIMIZE)
        elif i == 5:
            obj = (0.5 * ap[0][0] * (gp.quicksum(Pu[j] for j in range(Agent - 1)) ** 2) + bp[0][0] * gp.quicksum(
                Pu[j] for j in range(Agent - 1)) +gp. quicksum(ce[j] * Pu[j] for j in range(Agent - 1))
                   + 0.5 * ar[0][0] * (gp.quicksum(Ru[j] for j in range(Agent - 1)) ** 2) + br[0][0] * gp.quicksum(
                        Ru[j] for j in range(Agent - 1)) +gp. quicksum(cr[j] * Ru[j] for j in range(Agent - 1))
                   + gp.quicksum((nambda[0][j] * ((Pnm[0][j] - Pmn[0][j]) / 2 - Pu[j]) + nu[0][j] * (
                            (Rnm[0][j] - Rmn[0][j]) / 2 - Ru[j])
                               + ((rho[0][j] / 2) * ((Pnm[0][j] - Pmn[0][j]) / 2 - Pu[j]) ** 2)
                               + ((rho[0][j] / 2) * ((Rnm[0][j] - Rmn[0][j]) / 2 - Ru[j]) ** 2)) for j in range(Agent - 1))
                   + upsilon[0][6] * (
                           gp. quicksum(Pu[j] for j in range(Agent - 1)) + Pn[8] -gp. quicksum(Q[6, j] for j in range(Node)))
                   + (1 / 2) * (gp.quicksum(Pu[j] for j in range(Agent - 1)) + Pn[8] - gp.quicksum(
                        Q[6, j] for j in range(Node))) ** 2)
            mu.setObjective(obj, GRB.MINIMIZE)
        elif i==6:
            obj = (0.5 * ap[0][0] * (gp.quicksum(Pu[j] for j in range(Agent - 1)) ** 2) + bp[0][0] * gp.quicksum(
                Pu[j] for j in range(Agent - 1)) + gp.quicksum(ce[j] * Pu[j] for j in range(Agent - 1))
                   + 0.5 * ar[0][0] * (gp.quicksum(Ru[j] for j in range(Agent - 1)) ** 2) + br[0][0] *gp. quicksum(
                        Ru[j] for j in range(Agent - 1)) + gp.quicksum(cr[j] * Ru[j] for j in range(Agent - 1))
                   + gp.quicksum((nambda[0][j] * ((Pnm[0][j] - Pmn[0][j]) / 2 - Pu[j]) + nu[0][j] * (
                            (Rnm[0][j] - Rmn[0][j]) / 2 - Ru[j])
                               + ((rho[0][j] / 2) * ((Pnm[0][j] - Pmn[0][j]) / 2 - Pu[j]) ** 2)
                               + ((rho[0][j] / 2) * ((Rnm[0][j] - Rmn[0][j]) / 2 - Ru[j]) ** 2)) for j in range(Agent - 1))
                   + upsilon[0][8] * (
                           gp. quicksum(Pu[j] for j in range(Agent - 1)) + Pn[9] - gp.quicksum(Q[8, j] for j in range(Node)))
                   + (1 / 2) * (gp.quicksum(Pu[j] for j in range(Agent - 1)) + Pn[9] - gp.quicksum(
                        Q[8, j] for j in range(Node))) ** 2)
            mu.setObjective(obj, GRB.MINIMIZE)


        mu.addConstr(gp.quicksum(Pu[j] for j in range(Agent - 1)) >= Pmin, name='Pumin')
        mu.addConstr(gp.quicksum(Pu[j] for j in range(Agent - 1)) <= Pmax, name='Pumax')
        mu.addConstr(gp.quicksum(Ru[j] for j in range(Agent - 1)) >= Rmin, name='Rumin')
        mu.addConstr(gp.quicksum(Ru[j] for j in range(Agent - 1)) <= Rmax, name='Rumax')
        for j in range(Agent - 1):
            mu.addConstr(Pu[j] <= 0)
            mu.addConstr(Ru[j] >= 0)
        mu.addConstr(gp.quicksum(Pu[j] + Ru[j] for j in range(Agent - 1)) <= Pmax, name='Pu+Rumax')


        mu.optimize()

        if mu.status == GRB.OPTIMAL:
            print("Agentu op solution")
            Pu = np.array([Pu[j].x for j in range(Node)])
            Ru = np.array([Ru[j].x for j in range(Node)])
            return Pu, Ru
        else:
            print("Agentu no solution found")
