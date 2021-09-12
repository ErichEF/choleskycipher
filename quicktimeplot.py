# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 16:30:58 2019

@author: user
"""

from pylab import *
from numpy import *
import scipy as sp
import time

time = [6.44,8.72,10.8,13.68,17.04,20.81,24.7,29.55,34.75,48.18]
n = linspace(100,200,10)

plot(n,time)
grid()
xlabel('Size of matrices')
ylabel('Time in seconds')
title('Size of matrices versus time taken to compute')