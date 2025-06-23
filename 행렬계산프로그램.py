import sys
import numpy as np
from fractions import Fraction
input=sys.stdin.readline
class Matrix:
    def __init__(self):
        self.row, self.col=map(int,input().strip().split())
        data=[list(map(Fraction,input().strip().split())) for _ in range(self.row)]
        self.M=np.array(data, dtype=object)
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
        self.free_variables=[]
        for c in range(self.col):
            if r==self.row:
                break
            col=self.M[r:, c]
            nz=np.flatnonzero(col!=0)
            if nz.size==0:
                if c!=self.col-1:
                    self.free_variables.append(c)
                continue
            if nz[0]!=0:
                self.RRA('interchange', r, r+nz[0])
            self.pivot_list.append((r, c))
            for row in range(r+1, self.row):
                if self.M[row, c]!=0:
                    self.RRA('replacement', r, row, c)  
            r+=1
        self.rank=len(self.pivot_list)
    def RREF(self):
        for r, c in reversed(self.pivot_list):
            if self.M[r, c]!=0:
                self.RRA('scaling', r, 1/self.M[r, c])
            for ro in reversed(range(r)):
                if self.M[ro, c]!=0:
                    self.RRA('replacement', r, ro, c)

A=Matrix()
A.REF(); A.RREF()
inconsistent=False
for x in A.M:
    nz=np.flatnonzero(x!=0)
    if nz.size!=0:
        if nz[0]==A.col-1:
            inconsistent=True
            break
if inconsistent:
    print('Inconsistent System')
else:
    print('Consistent System')
    solutions=[None]*(A.col-1)
    for r, c in A.pivot_list:
        solutions[c]=A.M[r, -1]
    if A.col-1==A.rank:
        print('Uniqueness, Full Rank')
        print('x=('+', '.join(map(str, solutions))+')')
    else:
        print(f'Infinitely Many Solutions Case, Rank={A.rank}')
        print('Parametric Vector Form: p', end='')
        for idx, c in enumerate(A.free_variables):
            print(f' + x{c+1}v{idx+1}', end='')
        print('\nWhere:')
        Acol=A.col-1
        p=[Fraction(0)]*Acol
        for idx, i in enumerate(A.free_variables):
            v=[0]*Acol; v[i]=1
            for r, _ in A.pivot_list:
                v[r]=-A.M[r, i]
                p[r]=A.M[r, -1]
            print(f'v{idx+1}=({', '.join(map(str, v))})')
        print(f'p=({', '.join(map(str, p))})')





