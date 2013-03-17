from pylab import *

z = 1

def f(x):     return exp(-x) - z

def fp(x):    return -exp(-x)

def fpp(x):   return exp(-x)

#####################################################

def a(xo):    return fpp(xo)/2.0

def b(xo):    return fp(xo)-2*a(xo)*xo

def c(xo):    return f(xo)-a(xo)*xo**2-b(xo)*xo

def g(x, xo): return a(xo)*x**2 + b(xo)*x + c(xo)

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
    

xo = 1
go()
go()
## called 11 times and we're in the weeds!
## Also, because we're in such flat land we don't worry about the left-root confusing us by
## being to the right of the initial value.
