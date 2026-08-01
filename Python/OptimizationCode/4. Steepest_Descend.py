# Jakob Alsaker
# 201808760
# Optimization Programming and Optimization
# 24/10/2022

from numpy import *
from LineSearch import LineSearch

class SteepestDescend: # Remove Class - make function
    def __init__(self,f,g,x0,n,nmax=5e2,tol=1e-4,lstol=1e-3):
        self.f = f
        self.g = g
        self.x = x0
        self.n = n
        self.nmax = nmax
        self.tol = tol
        self.lstol = lstol
        
    def iterating(self):
        k = 0
        x_pre = zeros(self.n)
        while abs(linalg.norm(abs(g(self.x)))-linalg.norm(abs(g(x_pre)))) >= self.tol or k >= self.nmax:
            self.d_k = -g(self.x)
            qq = lambda a: f(self.x+a*self.d_k)
            meh = LineSearch(qq,0.1,self.nmax,self.lstol)
            alpha,i = meh.Run()
            x_pre = self.x
            self.x = self.x+alpha*self.d_k
            k = k+1
        return self.x,k
        
#f = lambda x: add(power(x[0],2),power(x[1],2))
#g = lambda x: array([multiply(2,x[0]),multiply(2,x[1])])¤x0 = array([0,0])
#n = 2

# f = lambda x: 1000*x[0]**2+100*x[1]**2+10*x[2]**2+x[3]**2+100*x[0]+10*x[1]+x[2]+0.1*x[3]
# g = lambda x: array([2000*x[0]+100,200*x[1]+10,20*x[2]+1,1.1+x[3]*0])
# x0 = array([2,2,2,2])
# n = 4


if __name__ == "__main__":
    q = SteepestDescend(f,g,x0,n)
    w,k = q.iterating()
    print("Amount of iteration in Steepest Descend: ",k)
    print("x values found is: ",w)
    print("Value of the function is: ",f(w))
    print("The gradient of the function is: ",g(w))