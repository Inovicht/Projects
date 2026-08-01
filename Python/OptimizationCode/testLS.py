from LineSearch import LineSearch
from numpy import *

def func1(x): #f = (x-1)^2
    f = subtract(x,1)
    f = power(f,2)
    return f
def func2(x): #f = (x-1.75)^2
    f = subtract(x,1.75)
    f = power(f,2)
    return f
def func3(x): #f = (x-0.01)^2
    f = subtract(x,0.01)
    f = power(f,2)
    return f
def func4(x): # f = (x-5)^2-x
    f = subtract(x,5)
    f = power(f,2)
    f = subtract(f,x)
    return f
def func5(x): # f = -(x-5)^3+5x
    f =subtract(x,5)
    f = power(f,3)
    f = -f
    p = multiply(5,x)
    f = add(f,p)
    return f

delta = 0.7
kmax = 1e5
tol = 1e-3

if __name__ == "__main__":
    func = func2
    I = LineSearch(func,delta,kmax,tol)
    alpha,i = I.Run()
    print("Final Alpha is: %.5f" %alpha)
    print("Final function Value: %.5f" %func(alpha))
    print("Amount of iterations of the cost function: " + str(i))