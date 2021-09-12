# -*- coding: utf-8 -*-
"""
Making array of zeros without numpy: https://stackoverflow.com/questions/4056768/how-to-declare-array-of-zeros-in-python-or-an-array-of-a-certain-size
"""

from numpy import *

#for n in range(0,900000000):
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
#data=loadtxt('ciphertest.txt',dtype='str')
#splitdata=array(list(data))
#spliceddata=array([])
#
#for i in range(len(splitdata)):
#    a = list(splitdata[i])
#    spliceddata=append(spliceddata,a)
#
#print(spliceddata)

#test=random.randint(0,41)

encoded=loadtxt('encodedmatrix.npy')

def decode(x):
#    data=array([])
    decodedarray=zeros(len(x),dtype='str')
    for k in range(len(lalphabet)):
        lowlettermodifier = ciphermatrix[random.randint(0,41),k]
        lowerletter = ciphermatrix[k,k]+lowlettermodifier
        caplettermodifier = ciphermatrix[len(ciphermatrix)-1,random.randint(0,41)]
        upperletter = ciphermatrix[k,k]+caplettermodifier
        for i in range(len(x)):
            if x[int(i)] == lowerletter:
                element = lalphabet[k]
                decodedarray[i]=element
            elif x[int(i)] == upperletter:
                element = calphabet[k]
                decodedarray[i]=element
            else:
                continue
    return decodedarray

decoded=decode(encoded)
fulldecode=' '.join(decoded)
print(decoded)
print(fulldecode)