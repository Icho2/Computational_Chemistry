#Problem Set 1

#Question 1
#a
primes = []
for i in range(2, 101):
    flag = True
    for j in range(2,i):
        if i % j == 0:
            flag = False
    if flag == True:
        primes.append(i)
##print(primes)

#b
form_primes = []
for i in range(1,11):
    odd_num = (2**i) - 1
    for num in primes:
        if odd_num == num:
            form_primes.append(odd_num)
##print(form_primes)

#Question 2
#a
import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt
c = 1
x_mean = 0
x_delta = 1
def F(x):
    return c*np.exp((-(x-x_mean)**2) / (2*(x_delta**2)))

x_vals = np.linspace(-5, 5, 100)
y_vals = F(x_vals)
#plt.plot(x_vals, y_vals)
plt.title('Gaussian Distribution')
plt.xlabel('x')
plt.ylabel('F(x)')
##plt.show()

#b)
#here I am changing the mean 
x_mean = 2
x_vals = np.linspace(-5, 5, 100)
y_vals = F(x_vals)
#plt.plot(x_vals, y_vals)
plt.title('Gaussian Distribution')
plt.xlabel('x')
plt.ylabel('F(x)')
##plt.show()

#here I change the delta x
x_mean = 0
x_delta = 2
x_vals = np.linspace(-5, 5, 100)
y_vals = F(x_vals)
#plt.plot(x_vals, y_vals)
plt.title('Gaussian Distribution')
plt.xlabel('x')
plt.ylabel('F(x)')
##plt.show()

##print("The mean of x parameter causes the distribution to translate to the left or right. Altering the delta of x causes the distribution to widen or narrow.")

#c)
#wrote this down on my ipad

#d)
from scipy import integrate

c = 1/(x_delta * np.sqrt(2 * np.pi))
x_mean = 0
x_delta = 1
lower_limit = -2.5 * x_delta
upper_limit =  2.5 * x_delta
result,error = integrate.quad(F, lower_limit, upper_limit)
#print("Decimal percentage of distribution: ",result)

#Question 3
#a

