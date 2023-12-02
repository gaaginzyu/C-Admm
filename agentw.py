import numpy as np
import gurobipy as gp
from gurobipy import GRB


def agentw(ce, cr, i, Agent, Node, nambda, nu, ap, bp, ar, br, Pmin, Pmax, Rmin, Rmax, Pnm, Pmn, Rnm, Rmn, rho, Q, Pn,
           upsilon):
    mw = gp.Model('Agentw')
    Pw = mw.addVars(Agent - 1, name='Pw',lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype = GRB.CONTINUOUS)
    Rw = mw.addVars(Agent - 1, name='Rw',lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype = GRB.CONTINUOUS)

    if i == 7:
        obj = (0.5 * ap * (gp.quicksum(Pw[j] for j in range(Agent - 1)) ** 2) +
               bp * gp.quicksum(Pw[j] for j in range(Agent - 1)) +
               gp.quicksum(ce[j] * Pw[j] for j in range(Agent - 1)) +
               0.5 * ar * (gp.quicksum(Rw[j] for j in range(Agent - 1)) ** 2) +
               br * gp.quicksum(Rw[j] for j in range(Agent - 1)) +
               gp.quicksum(cr[j] * Rw[j] for j in range(Agent - 1)) +
               gp.quicksum(nambda[j] * ((Pnm[j] - Pmn[j]) / 2 - Pw[j]) +
                            nu[j] * ((Rnm[j] - Rmn[j]) / 2 - Rw[j]) +
                              (rho[j] / 2) * ((Pnm[j] - Pmn[j]) / 2 - Pw[j]) ** 2  +
                              (rho[j] / 2) * ((Rnm[j] - Rmn[j]) / 2 - Rw[j]) ** 2 for j in range(Agent - 1)) +
               upsilon[4] * (gp.quicksum(Pw[j] for j in range(Agent - 1)) + Pn[4] -
                             gp.quicksum(Q[4, j] for j in range(Node))) +
               0.5 * (gp.quicksum(Pw[j] for j in range(Agent - 1)) + Pn[4] -
                      gp.quicksum(Q[4, j] for j in range(Node))) ** 2)
        # mw.setObjective(obj, GRB.MINIMIZE)
    elif i == 8:
        obj = (0.5 * ap * (gp.quicksum(Pw[j] for j in range(Agent - 1)) ** 2) +
               bp * gp.quicksum(Pw[j] for j in range(Agent - 1)) +
               gp.quicksum(ce[j] * Pw[j] for j in range(Agent - 1)) +
               0.5 * ar * (gp.quicksum(Rw[j] for j in range(Agent - 1)) ** 2) +
               br * gp.quicksum(Rw[j] for j in range(Agent - 1)) +
               gp.quicksum(cr[j] * Rw[j] for j in range(Agent - 1)) +
               gp.quicksum(nambda[j] * ((Pnm[j] - Pmn[j]) / 2 - Pw[j]) +
                              nu[j] * ((Rnm[j] - Rmn[j]) / 2 - Rw[j]) +
                              (rho[j] / 2) * ((Pnm[j] - Pmn[j]) / 2 - Pw[j]) ** 2 +
                              (rho[j] / 2) * ((Rnm[j] - Rmn[j]) / 2 - Rw[j]) ** 2 for j in range(Agent - 1)) +
               upsilon[6] * (gp.quicksum(Pw[j] for j in range(Agent - 1)) + Pn[5] -
                             gp.quicksum(Q[6, j] for j in range(Node))) +
               0.5 * (gp.quicksum(Pw[j] for j in range(Agent - 1)) + Pn[5] -
                      gp.quicksum(Q[6, j] for j in range(Node))) ** 2)
        # mw.setObjective(obj, GRB.MINIMIZE)

    elif i == 9:
        obj = (0.5 * ap * (gp.quicksum(Pw[j] for j in range(Agent - 1)) ** 2) +
               bp * gp.quicksum(Pw[j] for j in range(Agent - 1)) +
               gp.quicksum(ce[j] * Pw[j] for j in range(Agent - 1)) +
               0.5 * ar * (gp.quicksum(Rw[j] for j in range(Agent - 1)) ** 2) +
               br * gp.quicksum(Rw[j] for j in range(Agent - 1)) +
               gp.quicksum(cr[j] * Rw[j] for j in range(Agent - 1)) +
               gp.quicksum(nambda[j] * ((Pnm[j] - Pmn[j]) / 2 - Pw[j]) +
                            nu[j] * ((Rnm[j] - Rmn[j]) / 2 - Rw[j]) +
                              (rho[j] / 2) * ((Pnm[j] - Pmn[j]) / 2 - Pw[j]) ** 2 +
                              (rho[j] / 2) * ((Rnm[j] - Rmn[j]) / 2 - Rw[j]) ** 2 for j in range(Agent - 1)) +
               upsilon[8] * (gp.quicksum(Pw[j] for j in range(Agent - 1)) + Pn[6] -
                             gp.quicksum(Q[8, j] for j in range(Node))) +
               0.5 * (gp.quicksum(Pw[j] for j in range(Agent - 1)) + Pn[6] -
                      gp.quicksum(Q[8, j] for j in range(Node))) ** 2)
        # mw.setObjective(obj, GRB.MINIMIZE)

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
        Pw = np.array([Pw[j].x for j in range(Agent - 1)])
        Rw = np.array([Rw[j].x for j in range(Agent - 1)])
        return Pw, Rw
    else:
        print("Agentw no solution found")
