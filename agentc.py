import numpy as np
import gurobipy as gp
from gurobipy import GRB


def agentc(ce, cr, i, Agent, Node, nambda, nu, ap, bp, ar, br, Pmin, Pmax, Rmin, Rmax, Pnm, Pmn, Rnm, Rmn, rho, Q, Pn,
           upsilon):
    mc = gp.Model('Agentc')
    Pc = mc.addVars(Agent - 1, name='Pc',lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype = GRB.CONTINUOUS)
    Rc = mc.addVars(Agent - 1, name='Rc',lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype = GRB.CONTINUOUS)
    if i == 0:
        obj = (0.5 * ap * (gp.quicksum(Pc[j] for j in range(Agent - 1)) ** 2) +
               bp * gp.quicksum(Pc[j] for j in range(Agent - 1)) +
               gp.quicksum(ce[j] * Pc[j] for j in range(Agent - 1)) +
               0.5 * ar * (gp.quicksum(Rc[j] for j in range(Agent - 1)) ** 2) +
               br * gp.quicksum(Rc[j] for j in range(Agent - 1)) +
               gp.quicksum(cr[j] * Rc[j] for j in range(Agent - 1)) +
               gp.quicksum(nambda[j] * ((Pnm[j] - Pmn[j]) / 2 - Pc[j]) +
                            nu[j] * ((Rnm[j] - Rmn[j]) / 2 - Rc[j])+
                            (rho[j] / 2) * ((Pnm[j] - Pmn[j]) / 2 - Pc[j]) ** 2 +
                            (rho[j] / 2) * ((Rnm[j] - Rmn[j]) / 2 - Rc[j]) ** 2 for j in range(Agent - 1)) +
               upsilon[i] * (gp.quicksum(Pc[j] for j in range(Agent - 1)) + Pn[3] -
                               gp.quicksum(Q[0, j] for j in range(Node))) +
               0.5 * (gp.quicksum(Pc[j] for j in range(Agent - 1)) + Pn[3] -
                        gp.quicksum(Q[0, j] for j in range(Node))) ** 2)
    elif i == 1:
        obj = (0.5 * ap * (gp.quicksum(Pc[j] for j in range(Agent - 1)) ** 2) +
               bp * gp.quicksum(Pc[j] for j in range(Agent - 1)) +
               gp.quicksum(ce[j] * Pc[j] for j in range(Agent - 1)) +
                0.5 * ar * (gp.quicksum(Rc[j] for j in range(Agent - 1)) ** 2) +
               br * gp.quicksum(Rc[j] for j in range(Agent - 1)) +
               gp.quicksum(cr[j] * Rc[j] for j in range(Agent - 1)) +
               gp.quicksum(nambda[j] * ((Pnm[j] - Pmn[j]) / 2 - Pc[j]) +
                             nu[j] * ((Rnm[j] - Rmn[j]) / 2 - Rc[j]) +
                             (rho[j] / 2) * ((Pnm[j] - Pmn[j]) / 2 - Pc[j]) ** 2 +
                             (rho[j] / 2) * ((Rnm[j] - Rmn[j]) / 2 - Rc[j]) ** 2 for j in range(Agent - 1)) +
               upsilon[i] * (gp.quicksum(Pc[j] for j in range(Agent - 1)) -
                               gp.quicksum(Q[1, j] for j in range(Node))) +
               0.5 * (gp.quicksum(Pc[j] for j in range(Agent - 1)) -
                      gp.quicksum(Q[1, j] for j in range(Node))) ** 2)
    elif i == 2:
        obj = (0.5 * ap * (gp.quicksum(Pc[j] for j in range(Agent - 1)) ** 2) +
               bp * gp.quicksum(Pc[j] for j in range(Agent - 1)) +
               gp.quicksum(ce[j] * Pc[j] for j in range(Agent - 1)) +
               0.5 * ar * (gp.quicksum(Rc[j] for j in range(Agent - 1)) ** 2) +
               br * gp.quicksum(Rc[j] for j in range(Agent - 1)) +
               gp.quicksum(cr[j] * Rc[j] for j in range(Agent - 1)) +
               gp.quicksum(nambda[j] * ((Pnm[j] - Pmn[j]) / 2 - Pc[j]) +
                             nu[j] * ((Rnm[j] - Rmn[j]) / 2 - Rc[j]) +
                             (rho[j] / 2) * ((Pnm[j] - Pmn[j]) / 2 - Pc[j]) ** 2 +
                             (rho[j] / 2) * ((Rnm[j] - Rmn[j]) / 2 - Rc[j]) ** 2 for j in range(Agent - 1)) +
               upsilon[i] * (gp.quicksum(Pc[j] for j in range(Agent - 1)) -
                             gp.quicksum(Q[2, j] for j in range(Node))) +
               0.5 * (gp.quicksum(Pc[j] for j in range(Agent - 1)) -
                        gp.quicksum(Q[2, j] for j in range(Node))) ** 2)


    mc.setObjective(obj, GRB.MINIMIZE)

    mc.addConstr(gp.quicksum(Pc[j] for j in range(Agent - 1)) >= Pmin, name='Pcmin')
    mc.addConstr(gp.quicksum(Pc[j] for j in range(Agent - 1)) <= Pmax, name='Pcmax')
    mc.addConstr(gp.quicksum(Rc[j] for j in range(Agent - 1)) >= Rmin, name='Rcmin')
    mc.addConstr(gp.quicksum(Rc[j] for j in range(Agent - 1)) <= Rmax, name='Rcmax')
    for j in range(Agent - 1):
        mc.addConstr(Pc[j] >= 0)
        mc.addConstr(Rc[j] >= 0)

    mc.addConstr(gp.quicksum(Pc[j] + Rc[j] for j in range(Agent - 1)) <= Pmax, name='Pc+RCmax')

    mc.optimize()

    if mc.status == GRB.OPTIMAL:
        Pc = np.array([Pc[j].x for j in range(Agent - 1)])
        Rc = np.array([Rc[j].x for j in range(Agent - 1)])
        return Pc, Rc
    else:
        print("Agentc no solution found")
