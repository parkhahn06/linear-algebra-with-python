import sys
import numpy as np
from fractions import Fraction
input=sys.stdin.readline
def Mprint(M):
    maxc=0
    M=list([str(s) for s in x] for x in M)
    maxc=max(len(s) for x in M for s in x)
    for x in M:
        x=[s.ljust(maxc) for s in x]
        print('[ '+' '.join(x)+' ]')
class Matrix:
    def __init__(self):
        self.mode=input().strip()
        print('is augmented?:', end=' ')
        self.isaugmented=int(input())
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
        self.free_variables=list(range(self.col-self.isaugmented))
        for c in range(self.col):
            if r==self.row:
                break
            col=self.M[r:, c]
            nz=np.flatnonzero(col!=0)
            if nz.size==0:
                continue
            if nz[0]!=0:
                self.RRA('interchange', r, r+nz[0])
            self.pivot_list.append((r, c))
            for row in range(r+1, self.row):
                if self.M[row, c]!=0:
                    self.RRA('replacement', r, row, c)  
            r+=1
        for _, c in self.pivot_list:
            self.free_variables.remove(c)
        self.rank=len(self.pivot_list)
    def RREF(self):
        self.REF()
        for r, c in reversed(self.pivot_list):
            if self.M[r, c]!=0:
                self.RRA('scaling', r, 1/self.M[r, c])
            for ro in reversed(range(r)):
                if self.M[ro, c]!=0:
                    self.RRA('replacement', r, ro, c)
    def solve(self):
        self.RREF()
        inconsistent=False
        for x in self.M:
            nz=np.flatnonzero(x!=0)
            if nz.size!=0:
                if nz[0]==self.col-1:
                    inconsistent=True
                    break
        if inconsistent:
            print('Inconsistent System')
        else:
            print('Consistent System')
            solutions=[None]*(self.col-1)
            for r, c in self.pivot_list:
                solutions[c]=self.M[r, -1]
            if self.col-1==self.rank:
                print('Uniqueness, Full Rank')
                print('x=('+', '.join(map(str, solutions))+')')
            else:
                print(f'Infinitely Many Solutions Case, Rank={self.rank}')
                print('Parametric Vector Form: p', end='')
                for idx, c in enumerate(self.free_variables):
                    print(f' + x{c+1}v{idx+1}', end='')
                print('\nWhere:')
                Acol=self.col-1
                p=[0]*Acol
                for idx, i in enumerate(self.free_variables):
                    v=[0]*Acol; v[i]=1
                    for r, c in self.pivot_list:
                        v[c]=-self.M[r, i]
                        p[c]=self.M[r, -1]
                    print(f'v{idx+1}=({', '.join(map(str, v))})')
                print(f'p=({', '.join(map(str, p))})')

A=Matrix()
mode=A.mode
if mode=='REF': A.REF(); Mprint(A.M)
elif mode=='RREF': A.RREF(); Mprint(A.M)
elif mode=='solve': A.solve()


