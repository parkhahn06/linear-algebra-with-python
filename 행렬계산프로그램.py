import sys
import numpy as np
input=sys.stdin.readline
class Matrix:
    def __init__(self):
        self.row, self.col=map(int,input().strip().split())
        data=[list(map(float,input().strip().split())) for _ in range(self.row)]
        self.M=np.array(data, dtype=float)
        self.minrc=min(self.row, self.col)
        self.EPS=1e-10
    def RRA(self, mod, *args):
        if mod=='replacement':
            i, j, c=args
            k=self.M[j, c]/self.M[i, c]
            self.M[j]-=k*self.M[i]
        elif mod=='interchange':
            i, j=args
            self.M[[i, j]]=self.M[[j, i]]
        elif mod=='scaling':
            i, k=args
            self.M[i]*=k
    def REF(self):
        self.pivot_list=[]; r=0
        for c in range(self.col):
            if r==self.row:
                break
            col=self.M[r:, c]
            nz=np.flatnonzero(np.abs(col)>self.EPS)
            if nz.size==0:
                continue
            if nz[0]!=0:
                self.RRA('interchange', r, r+nz[0])
            self.pivot_list.append((r, c))
            for row in range(r+1, self.row):
                if not np.isclose(self.M[row, c], 0.0, atol=self.EPS):
                    self.RRA('replacement', r, row, c)  
            r+=1
        self.rank=len(self.pivot_list)
    def RREF(self):
        for r, c in reversed(self.pivot_list):
            if not np.isclose(self.M[r, c], 1.0, atol=self.EPS):
                self.RRA('scaling', r, 1/self.M[r, c])
            for ro in reversed(range(r)):
                if not np.isclose(self.M[ro, c], 0.0, atol=self.EPS):
                    self.RRA('replacement', r, ro, c)

A=Matrix()
A.REF(); A.RREF()
for x in A.M:
    print(round(x[-1]), end=' ')
            



