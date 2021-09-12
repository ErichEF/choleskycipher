# -*- coding: utf-8 -*-
"""
Making array of zeros without numpy: https://stackoverflow.com/questions/4056768/how-to-declare-array-of-zeros-in-python-or-an-array-of-a-certain-size
"""

from numpy import *
from cipherexamples import ciphertest1,ciphertest2,ciphertest3,ciphertest4

random.seed(892373902)

n = 41
E = random.rand(n,n)*random.random()
#D = 0.5*(D+transpose(D))
E = matmul(E,transpose(E)) # Generating "random" positive definite matrix

ciphermatrix = linalg.cholesky(E)

lalphabet=array(['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','.',',','!','?',\
                 '0','1','2','3','4','5','6','7','8','9',''])
calphabet=array(['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','.',',','!','?',\
                 '0','1','2','3','4','5','6','7','8','9',''])

data=loadtxt('ciphertest3.txt',dtype='str')
#data=ciphertest4
splitdata=array(list(data))
spliceddata=array([])

for i in range(len(splitdata)):
    a = list(splitdata[i])
    spliceddata=append(spliceddata,a)
    spliceddata=append(spliceddata,'')

#test=random.randint(0,41)

print(spliceddata)

def encode(x):
    encodedarray=zeros(len(x))
    for k in range(len(lalphabet)):
        lowerletter = lalphabet[k]
        upperletter = calphabet[k]
        lowlettermodifier = ciphermatrix[random.randint(0,41),k]
        caplettermodifier = ciphermatrix[len(ciphermatrix)-1,random.randint(0,41)] # Only pulls from last row so that modifier is never 0
        for i in range(len(x)):
            if x[int(i)] == lowerletter:
                element = ciphermatrix[k,k]
                encodedarray[i]=element+lowlettermodifier
            elif x[int(i)] == upperletter:
                element = ciphermatrix[k,k]+caplettermodifier
                encodedarray[i]=element
            else:
                continue
    return encodedarray

encoded=encode(spliceddata)

savetxt('encodedmatrix.npy',encoded)

print(encoded)