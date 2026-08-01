from numpy import *
class LineSearch:
    def __init__(self,func,delta=0.01,kmax=10**5,tol=0.01):
        self.func = func
        self.d = delta
        self.kmax = kmax
        self.tol = tol
        
    def Phase0(self,func,i,step=0.01,rho=0.2,beta=0.9): # Does not work
        rho = 0.2
        beta = 0.9
        def grad(func,a=0,nu=0.001):
            grads = divide(func(a)-func(a+nu),a-(a+nu))
            return grads
        a = func(step)
        b = add(func(0),rho*multiply(grad(func),step))
        c = abs(grad(func,step))
        d = beta*abs(grad(func))
        if a > b or c > d:
            k = 0
            while True or k <= self.kmax:
                a = func(step)
                b = add(func(0),0.2*multiply(grad(func),step))
                c = abs(grad(func,step))
                d = 0.9*abs(grad(func))
                i = add(i,8)
                if a <= b and c <= d or k > self.kmax:
                    break
                step = step+0.01
                k = add(k,1)
        return step,i
    
    def Phase1(self,i):
        a_1 = 0
        GR = 1.618
        for n in arange(kmax):
            a_0 = a_1
            a_1 += multiply(self.d,power(GR,n))
            i = add(i,2)
            if self.func(a_0) < self.func(a_1):
                a_0 -= multiply(self.d,power(GR,n-1))
                alpha = array([a_0, a_1])
                del a_0, a_1
                break
        try:
            if alpha is None:
                print("Meh")
        except:
            raise ValueError("Kmax is hit. You have found an unbounded minimum. Change delta")
        return alpha,i
    def Phase2(self,alpha,i):
        tau = 0.618
        a_l = alpha[0]
        a_u = alpha[1]
        while True:
            I = a_u-a_l
            if I < self.tol:
                break
            a_a = a_l+(1-tau)*I
            a_b = a_l+tau*I
            aa = self.func(a_a)
            ab = self.func(a_b)
            if aa < ab:
                a_u = a_b
            elif aa > ab:
                a_l = a_a
            elif aa == ab:
                a_l = a_a
                a_u = a_b
            i = add(i,2)
        alpha = array([a_l,a_u])
        return alpha,i
    def Interpol(self,alpha,i):
        a_u = alpha[1]
        a_l = alpha[0]
        a_i = average(alpha)
        fau = self.func(a_u)
        fal = self.func(a_l)
        fai = self.func(a_i)
        # 3 Eq. 3 unknowns
        a_2 = multiply(divide(1,subtract(a_u,a_i)),subtract(divide(subtract(fau,fal),subtract(a_u,a_l)),divide(subtract(fai,fal),subtract(a_i,a_l))))        
        a_1 = subtract(divide(subtract(fai,fal),subtract(a_u,a_i)),multiply(a_2,add(a_l,a_i)))
        a_0 = subtract(subtract(fal,multiply(a_1,a_l)),multiply(a_2,power(a_l,2)))
        i = add(i,3)
        #Final Alpha
        alpha = -divide(a_1,multiply(2,a_2))
        return alpha,i
                           
    def Run(self):
        i = 1
        self.d,i = self.Phase0(self.func,i,self.d)
        alpha,i = self.Phase1(i)
        alpha,i = self.Phase2(alpha,i)
        alpha_new,i = self.Interpol(alpha,i)
        return alpha_new,i
    

# Assignment specifications
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