import numpy as np
import gurobipy as gp
from gurobipy import GRB


def agentu(ce, cr, i, Agent, Node, nambda, nu, ap, bp, ar, br, Pmin, Pmax, Rmin, Rmax, Pnm, Pmn, Rnm, Rmn, rho, Q, Pn,
           upsilon):
    mu = gp.Model('Agentu')
    Pu = mu.addVars(Agent - 1, name='Pu',lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype = GRB.CONTINUOUS)
    Ru = mu.addVars(Agent - 1, name='Ru',lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype = GRB.CONTINUOUS)
    if i == 3:
        obj = (0.5 * ap * (gp.quicksum(Pu[j] for j in range(Agent - 1)) ** 2) +
               bp * gp.quicksum(Pu[j] for j in range(Agent - 1)) +
               gp.quicksum(ce[j] * Pu[j] for j in range(Agent - 1)) +
               0.5 * ar * (gp.quicksum(Ru[j] for j in range(Agent - 1)) ** 2) +
               br * gp.quicksum(Ru[j] for j in range(Agent - 1)) +
               gp.quicksum(cr[j] * Ru[j] for j in range(Agent - 1)) +
               gp.quicksum(nambda[j] * ((Pnm[j] - Pmn[j]) / 2 - Pu[j]) +
                            nu[j] * ((Rnm[j] - Rmn[j]) / 2 - Ru[j]) +
                              (rho[j] / 2) * ((Pnm[j] - Pmn[j]) / 2 - Pu[j]) ** 2  +
                              (rho[j] / 2) * ((Rnm[j] - Rmn[j]) / 2 - Ru[j]) ** 2 for j in range(Agent - 1))+
               upsilon[0] * (gp.quicksum(Pu[j] for j in range(Agent - 1)) + Pn[0] -
                               gp.quicksum(Q[0, j] for j in range(Node))) +
               0.5 * (gp.quicksum(Pu[j] for j in range(Agent - 1)) + Pn[0] -
                      gp.quicksum(Q[0, j] for j in range(Node))) ** 2)
        # mu.setObjective(obj, GRB.MINIMIZE)
    elif i == 4:
        obj = (0.5 * ap * (gp.quicksum(Pu[j] for j in range(Agent - 1)) ** 2) +
               bp * gp.quicksum(Pu[j] for j in range(Agent - 1)) +
               gp.quicksum(ce[j] * Pu[j] for j in range(Agent - 1)) +
               0.5 * ar * (gp.quicksum(Ru[j] for j in range(Agent - 1)) ** 2) +
               br * gp.quicksum(Ru[j] for j in range(Agent - 1)) +
               gp.quicksum(cr[j] * Ru[j] for j in range(Agent - 1))  +
               gp.quicksum(nambda[j] * ((Pnm[j] - Pmn[j]) / 2 - Pu[j]) +
                            nu[j] * ((Rnm[j] - Rmn[j]) / 2 - Ru[j]) +
                              (rho[j] / 2) * ((Pnm[j] - Pmn[j]) / 2 - Pu[j]) ** 2 +
                              ((rho[j] / 2) * ((Rnm[j] - Rmn[j]) / 2 - Ru[j]) ** 2) for j in range(Agent - 1)) +
               upsilon[4] * (gp.quicksum(Pu[j] for j in range(Agent - 1)) + Pn[7] -
                             gp.quicksum(Q[4, j] for j in range(Node))) +
               0.5 * (gp.quicksum(Pu[j] for j in range(Agent - 1)) + Pn[7] -
                      gp.quicksum(Q[4, j] for j in range(Node))) ** 2)
        # mu.setObjective(obj, GRB.MINIMIZE)
    elif i == 5:
        obj = (0.5 * ap * (gp.quicksum(Pu[j] for j in range(Agent - 1)) ** 2) +
               bp * gp.quicksum(Pu[j] for j in range(Agent - 1)) +
               gp.quicksum(ce[j] * Pu[j] for j in range(Agent - 1)) +
               0.5 * ar * (gp.quicksum(Ru[j] for j in range(Agent - 1)) ** 2) +
               br * gp.quicksum(Ru[j] for j in range(Agent - 1)) +
               gp.quicksum(cr[j] * Ru[j] for j in range(Agent - 1)) +
               gp.quicksum(nambda[j] * ((Pnm[j] - Pmn[j]) / 2 - Pu[j]) +
                            nu[j] * ((Rnm[j] - Rmn[j]) / 2 - Ru[j]) +
                              (rho[j] / 2) * ((Pnm[j] - Pmn[j]) / 2 - Pu[j]) ** 2 +
                              (rho[j] / 2) * ((Rnm[j] - Rmn[j]) / 2 - Ru[j]) ** 2 for j in range(Agent - 1)) +
               upsilon[6] * (gp.quicksum(Pu[j] for j in range(Agent - 1)) + Pn[8] -
                             gp.quicksum(Q[6, j] for j in range(Node))) +
               0.5 * (gp.quicksum(Pu[j] for j in range(Agent - 1)) + Pn[8] -
                      gp.quicksum(Q[6, j] for j in range(Node))) ** 2)
        # mu.setObjective(obj, GRB.MINIMIZE)
    elif i == 6:
        obj = (0.5 * ap * (gp.quicksum(Pu[j] for j in range(Agent - 1)) ** 2) +
               bp * gp.quicksum(Pu[j] for j in range(Agent - 1)) +
               gp.quicksum(ce[j] * Pu[j] for j in range(Agent - 1)) +
               0.5 * ar * (gp.quicksum(Ru[j] for j in range(Agent - 1)) ** 2) +
               br * gp.quicksum(Ru[j] for j in range(Agent - 1)) +
               gp.quicksum(cr[j] * Ru[j] for j in range(Agent - 1)) +
               gp.quicksum(nambda[j] * ((Pnm[j] - Pmn[j]) / 2 - Pu[j]) +
                            nu[j] * ((Rnm[j] - Rmn[j]) / 2 - Ru[j]) +
                              (rho[j] / 2) * ((Pnm[j] - Pmn[j]) / 2 - Pu[j]) ** 2 +
                              (rho[j] / 2) * ((Rnm[j] - Rmn[j]) / 2 - Ru[j]) ** 2 for j in range(Agent - 1)) +
               upsilon[8] * (gp.quicksum(Pu[j] for j in range(Agent - 1)) + Pn[9] -
                             gp.quicksum(Q[8, j] for j in range(Node))) +
               0.5 * (gp.quicksum(Pu[j] for j in range(Agent - 1)) + Pn[9] -
                      gp.quicksum(Q[8, j] for j in range(Node))) ** 2)
        # mu.setObjective(obj, GRB.MINIMIZE)

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
        Pu = np.array([Pu[j].x for j in range(Agent - 1)])
        Ru = np.array([Ru[j].x for j in range(Agent - 1)])
        return Pu, Ru
    else:
        print("Agentu no solution found")


