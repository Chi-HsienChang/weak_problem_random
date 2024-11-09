import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import numpy as np

x=[]
y=[]



with open('fe.txt', 'r') as file:
    for line in file:
        words = line.strip().split()
        print(words)
        x.append(int(words[0]))
        y.append(int(words[1]))

fig, ax = plt.subplots()

ax.set_xscale('log')
ax.set_yscale('log')

ax.scatter(x, y, color='blue', marker='^', label='Experiment')

coefficients = np.polyfit(np.log(x), np.log(y), 1)
polynomial = np.poly1d(coefficients)

x_fit = np.logspace(np.log10(x[0]), np.log10(x[-1]), 100)
ax.plot(x_fit, np.exp(polynomial(np.log(x_fit))), color='black', label=r'Trend Line: $O(l^{{ {:.2f} }}log(l))$'.format(polynomial[1]))


ax.set_xlabel('problem size')
ax.set_ylabel('population size')
plt.title('')
plt.legend()

PATH = './example_plot/'
for num in x:
    PATH = PATH + str(num) + '_'
PATH += '.png'

plt.savefig(PATH, dpi=300, bbox_inches='tight')
