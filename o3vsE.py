import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy.special import ellipk


mu0 = 3 / 2 ** (2 / 3)

def Kminmax(E, mu):
    sigma = 1
    coefficients = [1, 0, -2 * mu, -4 * sigma, -8 * E]
    roots = np.roots(coefficients)
    real_roots = [root.real for root in roots if abs(root.imag) < 1e-6]
    real_roots = np.sort(real_roots)

    if len(real_roots) > 2:
        kmax = real_roots[2]
        kmin = real_roots[1]
    else:
        kmax = max(real_roots)
        kmin = min(real_roots)

    return kmin, kmax

def khi(k):
    return ellipk(k ** 2)

def PI(a1square, k):
    n = a1square
    def integrand(theta):
        return 1 / ((1 - n * np.sin(theta) ** 2) * np.sqrt(1 - k ** 2 * np.sin(theta) ** 2))

    result, _ = quad(integrand, 0, np.pi/2)
    return result

def theta(E, mu):
    b, a = Kminmax(E, mu)
    a1sq = 4 / (a + b) + a * b - 0.25 * (a + b) ** 2
    b1 = -0.5 * (a + b)
    A = np.sqrt((a - b1) ** 2 + a1sq)
    B = np.sqrt((b - b1) ** 2 + a1sq)
    g = 1 / np.sqrt(A * B)
    alpha = (A - B) / (A + B)
    alpha2 = (b * A - a * B) / (a * B + b * A)
    k = np.sqrt(((a - b) ** 2 - ((A - B) ** 2)) / (4 * A * B))
    res = alpha2 * khi(k) + ((alpha - alpha2) / (1 - alpha ** 2)) * PI(((alpha ** 2) / (alpha ** 2 - 1)), k)
    res = res * 8 * g * (a * B + b * A) / (A - B)
    res = np.real(res)
    return res

def VK3(mu):
    sigma = 1
    m = (108 + 108 * np.sqrt(1 - (mu / mu0) ** 3 + 0j)) ** (1 / 3)
    kp = m / 12 + mu / m
    K = 2 * kp
    res = ((K ** 4) / 8 - (mu * K ** 2) / 4 - sigma * K / 2)
    return res

Ebar = np.linspace(0, 1, 1000)
mubars = [-2.5, -1, 0.5, 1.01, 1.5, 2]

# Create a plot
plt.figure(figsize=(8, 6))

# Loop through the mu values and plot the data
for mubar in mubars:
    mu = mubar * mu0
    E = np.real(Ebar * mu0 ** 2 + VK3(mu))
    y = [theta(E[i], mu) / (2 * np.pi) for i in range(len(E))]

    plt.plot(Ebar, y, label=r'$\mu = {:.2f}$'.format(mubar))

# Plot the horizontal lines y = 1/n for n values = -4, -3, -2, 2, 3, 4
n_values = [-4, -3, -2, 2, 3, 4]
for n in n_values:
    plt.axhline(y=1/n, color='black', linestyle=':', linewidth=1)
    text = "n = "+str(n)

    # Add labels positioned over the lines
    plt.text(Ebar[-1], 1/n, text, color='black', fontsize=10, verticalalignment='bottom', horizontalalignment='right')

# Add labels, title, and legend
plt.xlabel(r'$\overline{E}$', fontsize=14)
plt.ylabel(r'$\theta_i / (2 \pi)$', fontsize=14)
plt.title('Theta vs Energy for Different Values of $\mu$', fontsize=16)
plt.legend(title=r'$\mu$ Values', fontsize=12)

# Show the plot
plt.tight_layout()
plt.show()


