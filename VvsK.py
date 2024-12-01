import numpy as np
import matplotlib.pyplot as plt

global mu0
mu0 = 3/2**(2/3)
print(mu0**2)
# Define the function Vqq
def V_func(K_bar, mu_bar, sigma):
    K = K_bar*mu0**0.5
    mu = mu_bar*mu0
    # sigma = sigma*mu0**1.5
    return ((K**4)/8 - (mu*K**2)/4 - sigma*K/2)/mu0**2

# Generate K values and calculate V
K = np.linspace(-2, 2.5, 100)
V = V_func(K, 2, 1)

# Plot the function
plt.plot(K, V, label=r'$V(K)$', color='blue')
plt.axhline(0, color='black', linewidth=0.8, linestyle='solid')  # x-axis
plt.axvline(0, color='black', linewidth=0.8, linestyle='solid')  # y-axis

# Add labels and title
plt.xlabel(r'$K$', fontsize=12)
plt.ylabel(r'$V(K)$', fontsize=12)
plt.title('Potential Function Plot', fontsize=14)

# Activate grid and show legendq
# plt.grid(True)
# plt.ylim(-0.2,0.2)
plt.legend()
plt.show()