# -*- coding: utf-8 -*-
"""
Making array of zeros without numpy: https://stackoverflow.com/questions/4056768/how-to-declare-array-of-zeros-in-python-or-an-array-of-a-certain-size
"""
from pylab import *
from numpy import *
import scipy as sp
import time

C = array([[4,12,-16],[12,37,-43],[-16,-43,98]]) # Test matrix from wikipedia page
random.seed(20)

n = 300
E = random.rand(n,n)*random.random()
#D = 0.5*(D+transpose(D))
E = matmul(E,transpose(E)) # Generating random positive definite matrix

print('Randomly generated matrix E is',E)

test = linalg.cholesky(E)

def choleskymanual(A): 
#    L=[[0]]*len(A) # Numpy independent method
#    for i in range(len(A)):
#        L[i]=L[i]*len(A)
    L = zeros((len(A),len(A)))
    for i in range(len(A)):
        for j in range(i+1):
            lsum=0
            for k in range(j):
                lsum+=sum(L[i,k]*conj(L[j,k]))
            if i>j:
                L[i,j]=(1/L[j,j])*(A[i,j]-lsum)
            else:
                L[j,j]=sqrt(A[j,j]-lsum)
    return L

def cmLDL(A): 
    L = zeros((len(A),len(A)))
    D = zeros((len(A),len(A)))
    for i in range(len(A)):
        for j in range(i+1):
            lsum=0
            lsum2=0 # Require seperate sums as the calculations for elements of D and L use different elements from matrix L
            for k in range(j):
                lsum+=L[j,k]*conj(L[j,k])*D[k,k]
                lsum2+=L[i,k]*conj(L[j,k])*D[k,k]
            D[j,j]=A[j,j]-lsum
#            if i>j:
            L[i,j]=(1/D[j,j])*(A[i,j]-lsum2)
    return L,D


L=choleskymanual(E)
LD,D=cmLDL(E)
#print(L)
Ls=conj(transpose(L))
LDs=conj(transpose(LD))
Anew=matmul(L,Ls)
ADnew=matmul(LD,matmul(D,LDs))

print('L matrix from LL decomposition:',L)
print('L matrix from LDL decomposition:',LD)
print('Diagonal elements from D matrix from LDL decomposition:',diag(D))

#%%
"""
1st Application: Matrix Inversion
"""
start = time.time()
# How to get inverse: [A]][Ainv]=[LL*][Ainv]=I
# Print statements for troubleshooting
def cholinverse(A):
    L=choleskymanual(A)
    Ls=conj(transpose(L))
    I = identity((len(L)))
    U = zeros([len(A),len(A)])
    Ainv = zeros([len(A),len(A)])
    for i in range(len(U)):
        for j in range(i+1):
            if i==j:
                U[i,i]=(I[i,i]/L[i,i])
                Ainv[i,i]=(U[i,i]/Ls[i,i])
#    print(U)
#    print(Ainv)
    for i in range(len(U)):
#        print('i=',i)
#        lsum=0
        for j in range(i+1):
#            print('j=',j)
            lsum=0
            for k in range(i):
                lsum+=L[i,k]*U[k,j]
#                print('k=',k)
#                print(lsum)
#            print(lsum)
#            print(U)
            if i>j:
                U[i,j]=(-lsum)/L[i,i]
    for i in range(len(U)-1,-1,-1):
#        print('i=',i)
#        lsum=0
        for j in range(i,-len(U),-1):
#            print('j=',j)
            lsum=0
            for k in range(i-1,-len(U),-1):
                lsum+=Ls[i,k]*Ainv[k,j]
#                print('k=',k)
#                print(lsum)
#            print(lsum)
#            print(Ainv)
#            if i>j:
            Ainv[i,j]=(U[i,j]-lsum)/Ls[i,i]
            Ainv[j,i]=(U[i,j]-lsum)/Ls[i,i]
    return U,Ainv

AI = linalg.inv(E)
SA  = cholinverse(E)
print('The matrix inversion from numpy.linalg is')
print(AI)
print('The matrix inversion from manual method is')
print(SA[1])

ratioinv=AI/SA[1]
ratioinvf=array([])
for i in range(len(ratioinv)):
    o=mean(ratioinv[i])
    ratioinvf=append(ratioinvf,o)

ratioinvf=mean(ratioinvf)

print('The average ratio between the numpy and manual method at n =',n,' for matrix inversion is',ratioinvf)

#%%
"""
2nd Application: Monte-Carlo Simulation - Corellating unrelated random variables
"""

ranvar = random.normal(0,1,(n,100000)) # Generating n amount of variables for 2 variables from a gaussian distribution
#ranvar[1]=ranvar[1]*10000
#random.seed(40)
#ranvar1 = random.normal(0,1,100000)
#random.seed(60)
#ranvar2 = random.normal(0,1,100000)*10000
#ranvart=zeros(2,dtype=object)
#ranvart[0],ranvart[1]=muarray,ranvar2
#Imat = [[0.6,1],[1,0.4]]

#covmat = matmul(transpose(C),C)
covmat = matmul(transpose(E),E)

L = choleskymanual(covmat)

xcorr = dot(L,ranvar)
#xstd = sqrt(diag(covmat))

plot(abs(xcorr[0]),abs(xcorr[n-2]),'.')
grid()
xlabel('Correlated data set for variable x0')
ylabel('Correlated data set for variable x%s'%(n-2))
title('Correlation between correlated variables x0 and x%s'%(n-2))
figure()
plot(ranvar[0],ranvar[n-1],'.')
grid()
xlabel('Uncorrelated data set for variable x0')
ylabel('Uncorrelated data set for variable x%s'%(n-2))   
title('Correlation between uncorrelated variables x0 and x%s'%(n-2))
#%%
"""
3rd Application: Least-Squares

Because I've already computed the inverse of the matrix via Cholesky, I can simply take Ax=b -> x=Ainvb to get the system of unknown variables x
"""
#B=[0,1,2] # Use for C matrix or other random 3x3 matrix
B = random.rand(n,1)

X = matmul(cholinverse(E)[1],B)
print('The linear equation for a solution')
print(B)
print('via matrix inversion from cholesky gives vector x =')
print(X)

"""
Or I  can use linalg.solve
"""
X = linalg.solve(E,B)
print('The linear equation for a solution')
print(B)
print('via linalg.solve from cholesky gives vector x =')
print(X)


"""
Or I can calculate it more explicitly.

This is a modification of the cholinverse function since they both use back/forward substitution
Ax=b -> A*Ax=A*b -> (LLs)x=A*b -> (L)u=A*b
"""

def cholLS(A,b):
    posdef = matmul(transpose(A),A)
    newb=matmul(transpose(A),b)
    L=choleskymanual(posdef)
    Ls=conj(transpose(L))
    U = zeros([len(b)])
    X = zeros([len(b)])
    for i in range(len(A)):
        U[i]=(newb[i]-dot((L[i]),U))/L[i,i]
    for i in range(len(A)-1,-1,-1):
        X[i]=(U[i]-dot(Ls[i],X))/Ls[i,i]
    return X

print('The linear equation for a solution')
print(B)

Y = cholLS(E,B)

print('via least squares linear regression from cholesky gives vector x =')
print(Y)

ratiols=transpose(X)/Y
ratiolsf=array([])
for i in range(len(ratiols)):
    o=mean(ratiols[i])
    ratiolsf=append(ratiolsf,o)

ratiolsf=mean(ratiolsf)

print('The average ratio between the numpy and manual method at n =',n,' for matrix inversion is',ratioinvf)
print('The average ratio between the numpy and manual method at n =',n,' for least squares is',ratiolsf)

end = time.time()
print('Total time taken to run code = ',end - start)