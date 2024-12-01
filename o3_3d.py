import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy.special import ellipk


mu0 = 3/2**(2/3)

def Kminmax(E, mu):
    sigma = 1
    coefficients = [1, 0, -2 * mu, -4 * sigma, -8 * E]
    roots = np.roots(coefficients)
    # print(np.imag(roots))
    real_roots = [root.real for root in roots if abs(root.imag) < 1e-6]
    real_roots = np.sort(real_roots)


    if len(real_roots) > 2:
        kmax = real_roots[2]
        kmin = real_roots[1]
    else:
        kmax = max(real_roots)
        kmin = min(real_roots)
    # real_roots = np.sort(real_roots)

    return kmin, kmax

def khi(k):
    res = ellipk(k ** 2)
    return res


def PI(a1square, k):
    n = a1square
    def integrand(theta):
        return 1 / ((1 - n * np.sin(theta) ** 2) * np.sqrt(1 - k ** 2 * np.sin(theta) ** 2))

    result, _ = quad(integrand, 0, np.pi/2)
    return result


def theta(E, mu):
    b,a = Kminmax(E, mu)
    a1sq = 4/(a+b) + a*b - 0.25*(a+b)**2
    b1 = -0.5*(a+b)
    A = np.sqrt((a-b1)**2 + a1sq)
    B = np.sqrt((b-b1)**2 + a1sq)
    g = 1/np.sqrt(A*B)
    alpha = (A-B)/(A+B)
    alpha2 = (b*A - a*B) / (a*B + b*A)
    k = np.sqrt(((a-b)**2 - ((A-B)**2)) / (4*A*B))
    # print(k)
    res = alpha2*khi(k) + ((alpha - alpha2) / (1 - alpha**2)) * PI(((alpha**2) / (alpha**2 - 1)), k)
    res = res*8*g*(a*B + b*A)/(A-B)
    res = np.real(res)
    return res

def VK3(mu):
    sigma = 1
    m = (108 + 108 * np.sqrt(1 - (mu / mu0) ** 3 + 0j)) ** (1 / 3)
    kp = m / 12 + mu / m
    K = 2 * kp
    res = ((K ** 4) / 8 - (mu * K ** 2) / 4 - sigma * K / 2)
    return res

def filter(array):
    for i in range(len(array)):
        if np.isnan(array[i]):  # Check if the value is NaN
            array[i] = array[i - 1]  # Replace with the previous value
    return array



# Assuming VK3, theta, mu0, and filter are already defined

range1_points = 5000  # Number of points in [0, 0.4]
range2_points = 10000  # Number of points in [0.4, 0.7]
range3_points = 5000 # Number of points in [0.7, 1]

# Create the ranges
Ebar_range1 = np.linspace(0, 0.3, range1_points, endpoint=False)
Ebar_range2 = np.linspace(0.3, 0.5, range2_points, endpoint=False)
Ebar_range3 = np.linspace(0.5, 1, range3_points)

# Concatenate the ranges
Ebar = np.concatenate((Ebar_range1, Ebar_range2, Ebar_range3))

mubars = np.linspace(-0.5, 1, 10)

# Create meshgrid for 3D plot
E_mesh, mu_mesh = np.meshgrid(Ebar, mubars)

# Initialize Z (O) values
O_mesh = np.zeros_like(E_mesh)

# Calculate O for each combination of Ebar and mubars
for i, mubar in enumerate(mubars):
    mu = mubar * mu0
    E = np.real(Ebar * mu0 ** 2 + VK3(mu))
    O = []
    for e in E:
        val = theta(e, mu) / (2 * np.pi)
        O.append(val)

    # Filter NaN values and update O_mesh
    # O = filter(O)
    O_mesh[i, :] = np.array(O)



fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

surface = ax.plot_surface(
    mu_mesh, E_mesh, O_mesh,
    cmap='viridis',
    edgecolor='none',
    alpha=0.9
)

cbar = fig.colorbar(surface, ax=ax, shrink=0.5, aspect=10)
cbar.set_label('O value', rotation=270, labelpad=15)

ax.view_init(elev=25, azim=135)

ax.set_xlabel(r'$\overline{\mu}$', fontsize=14)
ax.set_ylabel(r'$\overline{E}$', fontsize=14)
ax.set_zlabel(r'$\theta_i / (2 \pi)$', labelpad = 10)
ax.set_title('3D Plot of $\theta_i / (2 \pi)$ vs $\overline{E}$ and $\overline{\mu}$', fontsize=14)

ax.set_zlim(-1, 1)

plt.show()



