from pylab import *

z = 1

def f(x):     return 1.0/x - 1

def fp(x):    return -1.0/x**2

def fpp(x):   return 2.0/x**3

def a(xo):    return 1.0/xo**3

def b(xo):    return -3.0/xo**2

def c(xo):    return 3.0/xo - z

def g(x, xo): return a(xo) * x**2 + b(xo) * x + c(xo)

def r(xo):
    D  = b(xo)**2 - 4*a(xo)*c(xo)
    if(D < 0) : 
        print "Warning, discriminant < 0 so no solution possible."
        return (0,0)
    r1 = (-b(xo) - sqrt(D)) / (2*a(xo))
    r2 = (-b(xo) + sqrt(D)) / (2*a(xo))
    return (r1, r2)

def go():
    global xo
    print r(xo)
    x = linspace(xo-0.7*xo,xo+0.7*xo,100)
    plot(x,f(x))
    plot(x,g(x,xo), '.')
    plot(xo,f(xo),"o")
    (tmp, xo) = r(xo)
    show()
    

xo = 2
go()
go()
## called 11 times and we're in the weeds!
## Also, because we're in such flat land we don't worry about the left-root confusing us by
## being to the right of the initial value.
