import numpy as np
import matplotlib.pyplot as plt

# Define mu0
mu0 = 3 / 2 ** (2 / 3)

# Define theaimax function
def theaimax(mu, n):
    # m is calculated from the given formula
    m = (108 + 108 * np.sqrt(1 - (mu / mu0) ** 3 + 0j)) ** (1 / 3)
    kp = m / 12 + mu / m
    km = m / 12 - mu / m  # Adjusting as per usual formulas

    # Calculate K based on n
    if n == 1:
        K = -1 * kp + 1j * np.sqrt(3) * km
    elif n == 3:
        K = 2 * kp
    else:
        raise ValueError("Incorrect value of n. Please choose n = 1 or n = 3.")

    # Compute the result (res)
    res = np.sqrt(3 - mu / K ** 2)  # This will result in a complex number if applicable
    res = np.sqrt(2) / res
    return res

# Generate the plot
mubar1 = np.linspace(1.01, 2, 100)
thetabar1 = theaimax(mubar1 * mu0, 1)

mubar2 = np.linspace(-3, 2, 100)
thetabar2 = theaimax(mubar2 * mu0, 3)

# Plot for n = 1
plt.plot(mubar1, -np.real(thetabar1), label=r'$\Theta_{0i}^{\max}$ for n=1')

# Plot for n = 3
plt.plot(mubar2, np.real(thetabar2), label=r'$\Theta_{0i}^{\max}$ for n=3')

# Add axis labels and title
plt.axhline(0, color='black', linewidth=0.8, linestyle='solid')  # x-axis
plt.axvline(0, color='black', linewidth=0.8, linestyle='solid')  # y-axis

# Add labels and title
plt.xlabel(r'$\mu / \mu_0$', fontsize=12)
plt.ylabel(r'$\Theta_{0i}^{\max}$', fontsize=12)
plt.title('Potential Function Plot', fontsize=14)

# Add a legend to differentiate the curves
plt.legend()

# Show the plot
plt.show()
