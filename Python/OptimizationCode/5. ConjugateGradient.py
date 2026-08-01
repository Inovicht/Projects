# Jakob Alsaker
# 201808760
# Optimization Programming and Optimization
# 24/10/2022

from numpy import *
from LineSearch import LineSearch

class ConjugateGradient:
    def __init__(self,f,g,x0,nmax=1e3,tol=1e-4,lstol=1e-3):
        self.f = f
        self.g = g
        self.x = x0
        self.n = len(x0)
        self.nmax = nmax
        self.tol = tol
        self.lstol = lstol
        
    def Find_beta(self,c_k,c_k_pre,d_k_pre):
        y_k = c_k-c_k_pre
        Beta_HS = divide(sum(c_k*y_k),sum(d_k_pre*y_k))
        Beta_FR = divide(sum(c_k*c_k),sum(c_k_pre*c_k_pre))
        Beta_PR = divide(sum(c_k*y_k),sum(c_k_pre*c_k_pre))        
        if 0 <= Beta_PR and Beta_PR <= Beta_FR:
            return Beta_PR
        elif Beta_PR > Beta_FR:
            return Beta_FR
        else:
            return 0
    
    def iterating(self):
        k = 0
        i = 0
        d_k_pre = array([0])
        c_k_pre = array([0])
        c_k = array([0])
        x_pre = zeros(self.n)
        while abs(linalg.norm(abs(g(self.x)))-linalg.norm(abs(g(x_pre)))) >= self.tol and k < self.nmax:
            c_k = g(self.x)
            if k % (self.n+1) == 0:
                beta_k = 0
            else:
                beta_k = self.Find_beta(c_k,c_k_pre,d_k_pre)
            i = add(i,3)
            self.d_k = -c_k+beta_k*d_k_pre
            qq = lambda a: f(self.x+a*self.d_k)
            meh = LineSearch(qq,0.1,self.nmax,self.lstol)
            alpha,w = meh.Run()
            i = add(i,w)
            x_pre = self.x
            self.x = self.x+alpha*self.d_k
            d_k_pre = self.d_k
            c_k_pre = c_k
            k = k+1
        return self.x,i


f = lambda x: 1000*x[0]**2+100*x[1]**2+10*x[2]**2+x[3]**2+100*x[0]+10*x[1]+x[2]+0.1*x[3]
g = lambda x: array([2000*x[0]+100,200*x[1]+10,20*x[2]+1,2*x[3]+0.1])
x0 = array([2,2,2,2])

if __name__ == "__main__":
    q = ConjugateGradient(f, g, x0)
    w,k = q.iterating()
    print("Amount of iteration in ConjugateGradient: ",k)
    print("x values found is: ",w)
    print("Value of the function is: ",f(w))
    print("The gradient of the function is: ",g(w))