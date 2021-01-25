#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 18:19:48 2021

description : Averages of Arithemetic functions and their plots 
author      : bveitch
version     : 1.0
project     : Averages of arithemtic functions

Supported functions:
    Dirichlet divisor. Averaging the number of divisors of an integer.
    Sum of sigma1 function. Averaging the sum of divisors of an integer.
    Chebyshev. Averaging the von Mangoldt function. Tests Riemann Hypothesis (not really!)
    TO DO:
    Number of integers visible from origin, 'Gauss circle problem'
    Others ...
    
"""
import numpy as np
import math 
import matplotlib.pyplot as plt

def sum_div(n):
    nsqrt=math.floor(math.sqrt(n))
    dsum=-nsqrt*nsqrt
    for i in range(1,nsqrt+1):
        dsum += 2*math.floor(n/i)
    return dsum

def sum_div_pow1(n):
    dsum=0
    for i in range(1,n+1):
        x=math.floor(n/i)
        dsum += x*(x+1)/2
    return dsum

def prime_sieve(n):
    primes=[True]*n
    primes[0]=False
    for p in range(2,len(primes)):
        is_a_p=primes[p-1]
        if(is_a_p):
            for j in range(p*p,n,p):
                primes[j-1]=False
    list_of_primes=[j for j in range(2,n) if primes[j-1]==True]
    return list_of_primes

def chebyshev(n,list_primes):
    csum=0
    for p in list_primes:
        if(p <= math.sqrt(n)):
            count=math.floor(math.log(n)/math.log(p))
            csum += count*math.log(p)
        elif(p<=n):
            csum += math.log(p)
        else:
            break
    return csum

gamma=0.577215664901532860606512090

#Plotting utils:
def plot_sum_divisors(n,fname):
    sum_divisors=[]
    err=[]
    y1=[]
    y05=[]
    for i in range(1,n):
        sumdiv=sum_div(i)
        sum_divisors.append(sumdiv)
        asymp=i*math.log(i)+(2*gamma-1)*i
        err.append(sumdiv-asymp)
        y1.append(math.pow(i,0.32))
        y05.append(math.pow(i,0.5))
    ym=[-y for y in y1]   
    y05m=[-y for y in y05] 
  
    
    fig, ax=plt.subplots()
    ax.set_title("Dirichlet Divisor function, $\Delta$(n), less than {:}".format(n))
    ax.set_xlabel('n')
    ax.plot(err, 'b',label='$\Delta({n})$')
    ax.plot(y1, 'r', label='$n^{0.32}$') 
    ax.plot(ym, 'r') 
    ax.plot(y05,'g', label='$n^{0.5}$')
    ax.plot(y05m, 'g')
    ax.legend(loc='lower left',shadow=True)
    plt.show()
    if(fname is not None):  
        plt.savefig(fname) 

def plot_sum_divisor_pow1(n,fname):
    sum_divisor1=[]
    asymptote=[]
    bigO=[]
    err=[]
    for i in range(2,n):
        sumdiv=sum_div_pow1(i)
        sum_divisor1.append(sumdiv)
        asymp=(math.pi**2/12)*i**2
        bigO.append(i*math.log(i))
        asymptote.append(asymp)
        err.append((sumdiv-asymp))
        
   
    fig, ax=plt.subplots()
    ax.set_title("Sum of divisor function, $\sigma$(n), for n less than {:}".format(n))
    ax.set_xlabel('n')
    ax.plot(err, 'b',label='$\sum_{x < n} \sigma_{1}(n)-({\pi^2/12}) n^2$')
    ax.plot(bigO, 'r',label='$nlog{n}$')
    ax.legend(loc='upper left',shadow=True)
    plt.show()
    if(fname is not None):  
        plt.savefig(fname) 

def plot_chebyshevfn(n,fname):
    list_of_primes=prime_sieve(n)
    cheby=[]
    err=[]
    sqrt=[]
    bigO=[]
    
    for i in range(2,n):
        csum=chebyshev(i,list_of_primes)
        cheby.append(csum)
        asymp=i
        sqrt.append(math.pow(i,0.5))
        bigO.append(math.pow(i,0.5)*math.log(i)**2)
        err.append(csum-asymp)
        
    sqrtm=[-y for y in sqrt]   
    bigOm=[-y for y in bigO]  
    fig, ax=plt.subplots()
    ax.set_title("Chebyshev $\psi$ function for n less than {:}".format(n))
    ax.set_xlabel('n')
    ax.plot(err, 'b',label='$\psi{(n)}-n$')
#These seem crap bound at this range (but assumes RH!!!)
#    ax.plot(bigO,'r', label='$n^{1/2}\log^2{n}$') 
#    ax.plot(bigOm,'r')
    
    ax.plot(sqrt,'g', label='$n^{1/2}$') 
    ax.plot(sqrtm,'g') 
    ax.legend(loc='upper right',shadow=True)
    if(fname is not None):  
        plt.savefig(fname) 

#Actually do the plotting:
fn='SumDivisors' 

if(fn=='Dirichlet'):
    plot_sum_divisors(100000,'Dirichlet.jpg')
elif(fn=='SumDivisors'):
    plot_sum_divisor_pow1(100000,'Sigma_Err.jpg')
elif(fn=='Chebyshev'):
    plot_chebyshevfn(10000,'Chebyshev.jpg')
else:
    ValueError("Invalid function type")