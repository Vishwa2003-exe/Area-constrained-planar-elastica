import numpy as np
import matplotlib.pyplot as plt

global mu0
mu0 = 3/2**(2/3)
print(mu0**2)

def V_func(K_bar, mu_bar, sigma):
    K = K_bar*mu0**0.5
    mu = mu_bar*mu0
    # sigma = sigma*mu0**1.5
    return ((K**4)/8 - (mu*K**2)/4 - sigma*K/2)/mu0**2

def Kdash_func(K_bar, mu_bar, sigma, E):
    res = 2*(E-V_func(K_bar, mu_bar, sigma))*mu0**2
    res = (res**0.5) / mu0
    return res, -1*res

# Generate K values and calculate V
K = np.linspace(-1.2, 2, 10000)

# Plot the function
colors = ['blue', 'green']
n = [2,5]
# Plotting
for idx, E in enumerate([0.253, 0.32]):
    Ebar = E / mu0**2
    p, m = Kdash_func(K, 0.5 * mu0, 1, Ebar)
    color = colors[idx]
    plt.plot(K, p, color=color, label=f'n = {n[idx]} (Positive Kdash)')
    plt.plot(K, m, color=colormu )

# Add axes lines
plt.axhline(0, color='black', linewidth=0.8, linestyle='solid')  # x-axis
plt.axvline(0, color='black', linewidth=0.8, linestyle='solid')  # y-axis

# Add labels and title
plt.xlabel(r'$K$', fontsize=12)
plt.ylabel(r"${K'}$", fontsize=12)
plt.title("$K'$ vs $K$", fontsize=14)

# Activate grid and show legend
plt.grid(True)
plt.legend()
plt.show()