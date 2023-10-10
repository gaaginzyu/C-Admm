import numpy as np
import scipy.io as scio
from DSO import DSO
from agentc import agentc
from agentu import agentu
from agentw import agentw
# get setup_data
data = scio.loadmat('setup.mat')
# data1 = scio.loadmat('finalresult C=10.mat')
Node=9
Agent=10
dso=1
k=1
resultrp = []
resultrr = []
resultsp = []
resultsr = []
#使用NumPy创建一个长度Node的一维NumPy数组
Pn=np.ones(Agent)
Q=np.zeros(Node)
theta=np.zeros(Node)
D=list()        #创建一个空列表
#值为 1 行向量
sp=data['sp']
sr=data['sr']
rp=data['rp']
rr=data['rr']
# susceptance Y set to 3  9*9矩阵
# Y=np.zeros([Node,Node])
# Y[0,3]=3
# Y[3,0]=3
# Y[3,4]=3
# Y[4,3]=3
# Y[3,8]=3
# Y[8,3]=3
# Y[2,5]=3
# Y[5,2]=3
# Y[5,6]=3
# Y[6,5]=3
# Y[6,7]=3
# Y[7,6]=3
# Y[1,7]=3
# Y[7,1]=3
# Y[4,5]=3
# Y[5,4]=3
# Y[7,8]=3
# Y[8,7]=3
Y = np.array([[0., 0., 0., 3., 0., 0., 0., 0., 0.],
              [0., 0., 0., 0., 0., 0., 0., 3., 0.],
              [0., 0., 0., 0., 0., 3., 0., 0., 0.],
              [3., 0., 0., 0., 3., 0., 0., 0., 3.],
              [0., 0., 0., 3., 0., 3., 0., 0., 0.],
              [0., 0., 3., 0., 3., 0., 3., 0., 0.],
              [0., 0., 0., 0., 0., 3., 0., 3., 0.],
              [0., 3., 0., 0., 0., 0., 3., 0., 3.],
              [0., 0., 0., 3., 0., 0., 0., 3., 0.]])
# line capacity limit C set to 10
C=10
#np的矩阵,访问时都要先写维度upsilon[0][5]。  upsilon=np.zeros(Node)访问时 upsilon[5]
upsilon=np.zeros([1,Node]) #[[0. 0. 0. 0. 0. 0. 0. 0. 0.]]    upsilon=np.zeros(Node) [0. 0. 0. 0. 0. 0. 0. 0. 0.]
Pnm=data['Pnm']
Pmn=data['Pmn']
Rnm=data['Rnm']
Rmn=data['Rmn']
nambda=data['lambda']
nu=data['nu']
rho=data['rho']
ce=data['ce'] #(10,9)
cr=data['cr'] #(10,9)
ap=data['ap']
bp=data['bp']
ar=data['ar']
br=data['br']
Pmin=data['Pmin']
Pmax=data['Pmax']
Rmin=data['Rmin']
Rmax=data['Rmax']
#Pnm是矩阵，只有一行，所以索引第一个[]只能为0，访问具体某个i元素,Pnm[0][i]
#Pnm[0][i]里元素也是二维的列表，也只有一行,访问具体某个j元素也需要加索引[0][j]
while (np.sum(sp)>0.001 or np.sum(sr)>0.001 or np.sum(rp)>0.001 or np.sum(rr)>0.001 or dso>0.01) and k <= 300:
    for i in range(Agent):
        Pn[i]=np.sum(Pnm[0][i])
    Q,theta,dso= DSO(Node,Pn,Y,C)
    D.append(dso)
    print(Q)
    print(theta)


    for i in range(3):
        Pnm[0][i][0][i] = 0 #Pnm[0][i]数组/行向量中 [0][i]个元素置为 0
        Pmn[0][i][0][i] = 0
        Rnm[0][i][0][i] = 0
        Rmn[0][i][0][i]= 0
        nambda[0][i][0][i] = 0
        nu[0][i][0][i] = 0
        rho[0][i][0][i] = 0
        # 返回的Pc Rc 是（9，）的 np 数组  #shapes (9,) (1,10)
        Pc,Rc=agentc(ce[i],cr[i],i,Agent,Node,nambda[0][i],nu[0][i],ap[0][i],bp[0][i],ar[0][i],br[0][i],
                     Pmin[0][i][0][0],Pmax[0][i][0][0],Rmin[0][i][0][0],Rmax[0][i][0][0],
                     Pnm[0][i],Pmn[0][i],Rnm[0][i],Rmn[0][i], rho[0][i],Q,Pn,upsilon)
        Pt =[]
        for j in range(9):
            Pt.append(Pc[j])
        Pt.insert(i, 0.0)
        Rt=[]
        for j in range(9):
            Rt.append(Rc[j])
        Rt.insert(i,0.0)
        sp[0][i]=np.sum((Pt-Pnm[0][i])* (Pt-Pnm[0][i]))
        sr[0][i]=np.sum((Rt-Rnm[0][i])* (Rt-Rnm[0][i]))
        Pnm[0][i]=Pt
        Rnm[0][i]=Rt
        print(Pnm)
        print(Rnm)

    for i in range(3,7):
        Pnm[0][i][0][i] = 0
        Pmn[0][i][0][i] = 0
        Rnm[0][i][0][i] = 0
        Rmn[0][i][0][i] = 0
        nambda[0][i][0][i] = 0
        nu[0][i][0][i] = 0
        rho[0][i][0][i] = 0
        Pu,Ru=agentu(ce[i],cr[i],i,Agent,Node,nambda[0][i],nu[0][i],ap[0][i],bp[0][i],ar[0][i],br[0][i],
                     Pmin[0][i][0][0],Pmax[0][i][0][0],Rmin[0][i][0][0],Rmax[0][i][0][0],
                     Pnm[0][i],Pmn[0][i],Rnm[0][i],Rmn[0][i],rho[0][i],Q,Pn,upsilon)
        Pt =[]
        for j in range(9):
            Pt.append(Pu[j])
        Pt.insert(i, 0.0)
        print(Pt)
        Rt=[]
        for j in range(9):
            Rt.append(Ru[j])
        Rt.insert(i,0.0)
        print(Rt)
        sp[i] = np.sum((Pt - Pnm[0][i]) * (Pt - Pnm[0][i]))
        sr[i] = np.sum((Rt - Rnm[0][i]) * (Rt - Rnm[0][i]))
        Pnm[0][i] = Pt
        Rnm[0][i] = Rt

    for i in range(7,10):
        Pnm[0][i][0][i] =0
        Pmn[0][i][0][i] = 0
        Rnm[0][i][0][i] = 0
        Rmn[0][i][0][i] = 0
        nambda[0][i][0][i] = 0
        nu[0][i][0][i] = 0
        rho[0][i][0][i] = 0
        Pw, Rw = agentw(ce[i], cr[i], i, Agent, Node, nambda[0][i], nu[0][i], ap[0][i], bp[0][i], ar[0][i], br[0][i],
                        Pmin[0][i][0][0], Pmax[0][i][0][0], Rmin[0][i][0][0], Rmax[0][i][0][0],
                        Pnm[0][i], Pmn[0][i], Rnm[0][i], Rmn[0][i],rho[0][i], Q, Pn, upsilon)
        Pt = []
        for j in range(9):
            Pt.append(Pw[j])
        Pt.insert(i, 0.0)
        print(Pt)
        Rt = []
        for j in range(9):
            Rt.append(Rw[j])
        Rt.insert(i, 0.0)
        print(Rt)
        sp[i] = np.sum((Pt - Pnm[0][i]) * (Pt - Pnm[0][i]))
        sr[i] = np.sum((Rt - Rnm[0][i]) * (Rt - Rnm[0][i]))
        Pnm[0][i] = Pt
        Rnm[0][i] = Rt

    for i in range(Agent):
        nambda[0][i][0][i]=0
        nu[0][i][0][i]=0
        Pnm[0][i][0][i]=0
        Pmn[0][i][0][i]=0
        Rnm[0][i][0][i]=0
        Rmn[0][i][0][i]=0
        rho[0][i][0][i]=0

    for i in range(Agent):
        for j in range(Agent):
            Pmn[0][i][0][j] = Pnm[0][j][0][i]

    for i in range(Agent):
        for j in range(Agent):
            Rmn[0][i][0][j] = Rnm[0][j][0][i]

    for i in range(Agent):
        nambda[0][i] = nambda[0][i] - np.array(rho[0][i]) * (np.array(Pnm[0][i]) + np.array(Pmn[0][i])) / 2
        nu[0][i] = nu[0][i] - np.array(rho[0][i]) * (np.array(Rnm[0][i]) + np.array(Rmn[0][i])) / 2

    upsilon[0] = max(upsilon[0] + (Pn[0] + Pn[3] - np.sum(Q[0, :])), 0)
    upsilon[1] = max(upsilon[1] + (Pn[1] - np.sum(Q[1, :])), 0)
    upsilon[2] = max(upsilon[2] + (Pn[2] - np.sum(Q[2, :])), 0)
    upsilon[3] = max(upsilon[3] + (-np.sum(Q[3, :])), 0)
    upsilon[4] = max(upsilon[4] + (Pn[4] + Pn[7] - np.sum(Q[4, :])), 0)
    upsilon[5] = max(upsilon[5] + (-np.sum(Q[5, :])), 0)
    upsilon[6] = max(upsilon[6] + (Pn[6] + Pn[8] - np.sum(Q[6, :])), 0)
    upsilon[7] = max(upsilon[7] + (-np.sum(Q[7, :])), 0)
    upsilon[8] = max(upsilon[8] + (Pn[7] + Pn[9] - np.sum(Q[8, :])), 0)

    for i in range(Agent):
        rp[i] = np.sum((np.array(Pnm[0][i]) + np.array(Pmn[0][i])) * (np.array(Pnm[0][i]) + np.array(Pmn[0][i])))
        rr[i] = np.sum((np.array(Rnm[0][i]) + np.array(Rmn[0][i])) * (np.array(Rnm[0][i]) + np.array(Rmn[0][i])))

    resultrp.append(np.sum(rp))
    resultrr.append(np.sum(rr))
    resultsp.append(np.sum(sp))
    resultsr.append(np.sum(sr))

    k += 1