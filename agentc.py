import numpy as np
import gurobipy as gp
from gurobipy import GRB
def agentc(ce,cr,i,Agent,Node,nambda,nu,ap,bp,ar,br,Pmin,Pmax,Rmin,Rmax,Pnm,Pmn,Rnm,Rmn,rho,Q,Pn,upsilon) :
        mc =gp.Model('Agentc')
        Pc=mc.addVars(Agent-1,name='Pc')
        Rc=mc.addVars(Agent-1,name='Rc')
        if i == 0:
             obj = (0.5 * ap[0][0] * (gp.quicksum(Pc[j] for j in range(Agent - 1)) ** 2) + bp[0][0] * gp.quicksum(Pc[j] for j in range(Agent - 1)) + gp.quicksum(ce[j] * Pc[j] for j in range(Agent - 1))
                    + 0.5 * ar[0][0] * (gp.quicksum(Rc[j] for j in range(Agent - 1)) ** 2) + br[0][0] * gp.quicksum(Rc[j] for j in range(Agent - 1)) +gp.quicksum( cr[j] * Rc[j] for j in range(Agent - 1))
                    + gp.quicksum((nambda[0][j] * ((Pnm[0][j] - Pmn[0][j]) / 2 - Pc[j])+ nu[0][j] * ((Rnm[0][j] - Rmn[0][j]) / 2 - Rc[j])
                    + ((rho[0][j] / 2) * ((Pnm[0][j] - Pmn[0][j]) / 2 - Pc[j]) ** 2)
                    + ((rho[0][j] / 2) * ((Rnm[0][j] - Rmn[0][j]) / 2 - Rc[j]) ** 2)) for j in range(Agent - 1) )
                    + upsilon[0][i] * (gp.quicksum(Pc[j] for j in range(Agent - 1)) + Pn[3] - gp.quicksum(Q[0, j] for j in range(Node)))
                       + (1 / 2) * (gp.quicksum(Pc[j] for j in range(Agent - 1)) + Pn[3] - gp.quicksum(Q[0, j] for j in range(Node))) ** 2 )
             mc.setObjective(obj, GRB.MINIMIZE)
        elif i==1:

            obj = (0.5 * ap[0][0] * (gp.quicksum(Pc[j] for j in range(Agent - 1)) ** 2) + bp[0][0] * gp.quicksum(Pc[j] for j in range(Agent - 1)) +gp.quicksum(ce[j] * Pc[j] for j in range(Agent - 1))
                   + 0.5 * ar[0][0] * (gp.quicksum(Rc[j] for j in range(Agent - 1)) ** 2) + br[0][0] * gp.quicksum(Rc[j] for j in range(Agent - 1)) +gp.quicksum( cr[j] * Rc[j] for j in range(Agent - 1))
                   + gp.quicksum(nambda[0][j] * ((Pnm[0][j] - Pmn[0][j]) / 2 - Pc[j])
                   + nu[0][j] * ((Rnm[0][j] - Rmn[0][j]) / 2 - Rc[j])
                   + ((rho[0][j] / 2) * ((Pnm[0][j] - Pmn[0][j]) / 2 - Pc[j]) ** 2)
                   + ((rho[0][j] / 2) * ((Rnm[0][j] - Rmn[0][j]) / 2 - Rc[j]) ** 2) for j in range(Agent - 1))
                   + upsilon[0][i] * (gp.quicksum(Pc[j] for j in range(Agent - 1)) - gp.quicksum(Q[1, j] for j in range(Node)))
                      + (1 / 2) * (gp.quicksum(Pc[j] for j in range(Agent - 1)) - gp.quicksum(Q[1, j] for j in range(Node))) ** 2)
            mc.setObjective(obj, GRB.MINIMIZE)
        elif i==2:

            obj = (0.5 * ap[0][0] * (gp.quicksum(Pc[j] for j in range(Agent - 1)) ** 2) + bp[0][0] * gp.quicksum(
                Pc[j] for j in range(Agent - 1)) + gp.quicksum(ce[j] * Pc[j] for j in range(Agent - 1))
                   + 0.5 * ar[0][0] * (gp.quicksum(Rc[j] for j in range(Agent - 1)) ** 2) + br[0][0] * gp.quicksum(
                        Rc[j] for j in range(Agent - 1)) + gp.quicksum(cr[j] * Rc[j] for j in range(Agent - 1))
                   + gp.quicksum(nambda[0][j] * ((Pnm[0][j] - Pmn[0][j]) / 2 - Pc[j])
                              + nu[0][j] * ((Rnm[0][j] - Rmn[0][j]) / 2 - Rc[j])
                              + ((rho[0][j] / 2) * ((Pnm[0][j] - Pmn[0][j]) / 2 - Pc[j]) ** 2)
                              + ((rho[0][j] / 2) * ((Rnm[0][j] - Rmn[0][j]) / 2 - Rc[j]) ** 2) for j in
                              range(Agent - 1))
                   + upsilon[0][i] * (
                               gp.quicksum(Pc[j] for j in range(Agent - 1)) - gp.quicksum(Q[2, j] for j in range(Node)))
                   + (1 / 2) * (
                           gp.quicksum(Pc[j] for j in range(Agent - 1)) - gp.quicksum(Q[2, j] for j in range(Node))) ** 2)
            mc.setObjective(obj, GRB.MINIMIZE)


        mc.addConstr(gp.quicksum(Pc[j] for j in range(Agent-1)) >= Pmin, name='Pcmin')
        mc.addConstr(gp.quicksum(Pc[j] for j in range(Agent-1)) <= Pmax, name='Pcmax')
        mc.addConstr(gp.quicksum(Rc[j] for j in range(Agent-1))>= Rmin, name='Rcmin')
        mc.addConstr(gp.quicksum(Rc[j]  for j in range(Agent-1))<= Rmax, name='Rcmax')
        for j in range(Agent-1):
            mc.addConstr(Pc[j]>= 0 )
            mc.addConstr(Rc[j]>= 0 )
        mc.addConstr(gp.quicksum(Pc[j]+Rc[j] for j in range(Agent-1)) <= Pmax, name='Pc+RCmax')

        mc.optimize()

        if mc.status==GRB.OPTIMAL:
            print("Agentc op solution")
            Pc = np.array([Pc[j].x  for j in range(Node)])
            Rc = np.array([Rc[j].x  for j in range(Node)])
            return Pc, Rc
        else:
            print("Agentc no solution found")



