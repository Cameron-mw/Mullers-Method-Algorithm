︠b291558b-670f-4f69-b99f-48fd84603e9ds︠
def Horner(p,t): #Horners method used for evaluating a polynomial at a given t-value in Muller's and Newton's
    n=len(p)
    b=p[0:n]
    i=1
    while i <= n-1:
        b[i]=p[i]+t*b[i-1]
        i+=1

    n2=len(b)
    q=b[0:n2-1]
    i=1
    while i <= n2-2: #Ignore the last value of b
        q[i]=b[i]+t*q[i-1]
        i+=1
    return b, q[-1] #Return p(x) and p'(t)=q(t)


def Muller(p):
    #Arbitrary Starting Guesses
    r0=-3.534 #r_(n-2)
    r1=2.893 #r_(n-1)
    r2=6.523 #r_n
    i=0
    while (abs(r2-r1) > 0.01) and i < 10: #Want a rough estimate of the root from Muller's using at most 11 iterations
        fr0=Horner(p,r0)[0][-1] #f(r_(n-2))
        fr1=Horner(p,r1)[0][-1] #f(r_(n-1))
        fr2=Horner(p,r2)[0][-1] #f(r_n)

        #Coefficients for the quadratic we want to solve
        a=((r1-r2)*(fr0-fr2)-(r0-r2)*(fr1-fr2))/((r0-r2)*(r1-r2)*(r0-r1))
        b=((r0-r2)^2*(fr1-fr2)-(r1-r2)^2*(fr0-fr2))/((r0-r2)*(r1-r2)*(r0-r1))
        c=fr2

        if abs(b+sqrt(b^2-4*a*c))>abs(b-sqrt(b^2-4*a*c)): #Checking for distance to the root
            r0=r1
            r1=r2
            r2=r2-(2*c)/(b+sqrt(b^2-4*a*c))

        else:
            r0=r1
            r1=r2
            r2=r2-(2*c)/(b-sqrt(b^2-4*a*c))
        i+=1
    return r2

def H_Newton(p,r,E):
    if r.imag()==0: #If the root is real, use regular Newton's
        rn=r-(Horner(p,r)[0][-1]/Horner(p,r)[-1]) #First Iteration
        while abs(rn-r) >= E:
            r=rn
            rn=float(rn-(Horner(p,rn)[0][-1]/Horner(p,rn)[-1])) #Netwon and Horner's method to find root
        return rn
    else: #If the root in complex, use adjusted Newton's (float not accepted for complex numbers)
        rn=r-(Horner(p,r)[0][-1]/Horner(p,r)[-1]) #First Iteration
        while abs(rn-r) >= E:
            r=rn
            rn=(rn-(Horner(p,rn)[0][-1]/Horner(p,rn)[-1])) #Netwon and Horner's method to find root
        return rn

def Quadratic(p): #Finding the roots of the remaining Quadratic
    a=p[0] #Assign each value in p to a quadratic coefficient (p is necessarily length 3)
    b=p[1]
    c=p[2]
    root1=(-2*c)/(b+sqrt(b^2-4*a*c)) #First Root
    root2=(-2*c)/(b-sqrt(b^2-4*a*c)) #Second Root
    return root1, root2

def FindRoots(p,E): #This is the root finding function
    roots=[]
    while len(p) > 3: #For deflated polynomials greater than degree 2
        rn=Muller(p) #Find approx root with Muller's
        root=H_Newton(p,rn,E) #Find precisce root with Newton's
        roots.append(root) #Add root to the roots list
        temp=Horner(p,root)[0] #Copy the coefficient list to a temporary list and evaluate the current root
        temp.remove(Horner(p,root)[0][-1]) #remove the last value from the temporary list
        p=temp #reassign p as this shortened list
        
    roots.append(Quadratic(p)[0]) #Find the first root of the quadratic polynomial 
    roots.append(Quadratic(p)[1]) #Find the second root
    
    for i in xrange(0,len(p)-1): #This is a debugging for-loop for real roots with negligibly small imaginary parts which don't exist in the true root
        if roots[i].imag()==0: #Ignore real roots
            i+=1
        else:
            if abs(roots[i].imag()) <= 10^-5: #10^-5 is an arbitrary choice but it could be lowered if needed, presuming 10^-5 indicates an erroneous complex part of the root
                roots[i]=roots[i].real() #Remove erroneous complex part
                i+=1
            else:
                i+=1 #If imaginary part of the root is correct     
    return roots



poly=[1,2,4,5,2]
FindRoots(poly,10^-8)
︡c3866ea9-b09d-4550-a005-dc3a637e8c71︡{"stdout":"[-0.142387380782455 - 1.66614757361206*I, -0.142387380782455 + 1.66614757361206*I, -0.715225238435090 - 1.55944019550637e-15*I, -1.00000000000000 + 1.55944019550637e-15*I]\n"}︡{"done":true}
︠fd04c3c3-854c-4ecb-8c7f-51c64af1a135s︠
# Functions used to test:
f = []
f.append([1,-4]) # x-4
f.append([1, 1, -6]) # (x-2)(x+3)
f.append([1, 400, 10000, -6000000]) # (x+200)(x+300)(x-100)
f.append([1, -2, 1, -2]) # (x-i)(x+i)(x-2)
f.append([1, 994, -5990, 10000]) # (x-3+i)(x-3-i)(x+1000)
f.append([1, -33, 261, -33, 260]) # (x-i)(x+i)(x-20)(x-13)
f.append([1, 0, 5, 0, 4]) # (x-i)(x+i)(x-2i)(x+2i)
f.append([1, 5, 5, 25, 4, 20]) # (x-i)(x+i)(x-2i)(x+2i)(x+5)
f.append([1, -20.88329, -110.15242564, 25.64032992]) # (x-0.22341)(x+4.55212)(x-25.212)

for i in range(len(f)):
    try:
        FindRoots(f[i],10^(-5))
    except:
        print "The coefficients", f[i], "failed \n"
        continue
︡8dacab27-0eb9-494f-9200-32aed675c67c︡{"stdout":"The coefficients [1, -4] failed \n\n[2, -3]\n[100.000000000000, -200.000000000000, -300.000000000000]\n[2.00000000006306, -2.73305470941835e-11 + 1.00000000006516*I, -3.57314529721549e-11 - 1.00000000006096*I]\n[3.00000000000000 - 0.999999999999999*I, 3.00000000000000 + 0.999999999999999*I, -1000.00000000000]\n[-2.03192455868222e-17 - 1.00000000000000*I, 13.0000000000000, 20.0000000000000, -1.66533453693774e-17 + 1.00000000000000*I]\n[6.80208009174734e-17 - 2.00000000000000*I, -6.81148157585049e-17 + 2.00000000000000*I, 1.36182624096494e-16 + 0.999999999999999*I, -1.36088609255463e-16 - 0.999999999999999*I]\n[5.05475981309420e-17 - 2.00000000000000*I, -4.93018165655520e-17 + 2.00000000000000*I, -1.19275909335025e-16 - 1.00000000000000*I, 3.33066907387547e-17 + 1.00000000000000*I, -5.00000000000000 + 2.22044604925031e-16*I]\n[0.223409999972561, 25.2119999996953, -4.55211999966787]\n"}︡{"done":true}
︠a9ae57d5-a8eb-489e-b9c0-b07aa1fc98c0︠










