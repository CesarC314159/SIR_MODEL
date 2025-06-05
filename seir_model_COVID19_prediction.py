import numpy as np
from scipy.integrate import odeint
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# Sample aggregated data
aggregated_data = [6547, 26969, 21670, 101291, 314106, 146799, 76933, 98787, 194342, 329811, 375642, 180655, 152136, 172573, 88544, 53347, 295274, 647111, 310518, 78647, 46783, 532580, 1328421, 227990, 43192, 82125, 235497, 323730, 320381, 216503, 87891, 45454, 57749, 102697, 107409, 59314, 41786, 31372, 26125, 28405, 46518, 91930, 57063, 27072, 21915, 43089, 70581, 46427, 20540, 11626, 14760, 38162, 94260, 62034, 32376, 11352, 7775, 12154, 16453, 12985, 10342, 9733, 10798]

# Time points (one per one-month interval)
t = np.arange(len(aggregated_data))

# SEIR model differential equations
def seir_model(y, t, beta, sigma, gamma):
    S, E, I, R = y
    N = S + E + I + R
    dSdt = -beta * S * I / N
    dEdt = beta * S * I / N - sigma * E
    dIdt = sigma * E - gamma * I
    dRdt = gamma * I
    return [dSdt, dEdt, dIdt, dRdt]

# Cost function to minimize
def residuals(params, t, data):
    beta, sigma, gamma = params
    E0 = data[0]
    y0 = [21e6 - E0, E0, 0, 0]
    sol = odeint(seir_model, y0, t, args=(beta, sigma, gamma))
    I = sol[:, 2]
    return np.sum((I - data)**2)

# Initial guess for parameters
initial_guess = [0.1, 0.1, 0.1]

# Optimization
result = minimize(residuals, initial_guess, args=(t, aggregated_data), method='L-BFGS-B', bounds=[(0, 1), (0, 1), (0, 1)])
beta_opt, sigma_opt, gamma_opt = result.x

# Solve SEIR with optimized parameters
E0 = aggregated_data[0]
y0 = [21e6 - E0, E0, 0, 0]
sol = odeint(seir_model, y0, t, args=(beta_opt, sigma_opt, gamma_opt))
S, E, I, R = sol.T

# Plot results
plt.figure(figsize=(20, 10))
plt.plot(t, aggregated_data, 'o', label='Actual Data')
plt.plot(t, I, '-', label='SEIR Model Prediction')
plt.xlabel('Time (1-month intervals)')
plt.ylabel('Number of Infected Individuals')
plt.title('SEIR Model Fit to COVID-19 Data')
plt.legend()
plt.grid(True)
plt.show()

# Print optimized parameters
print(f"Optimized Parameters: Beta: {beta_opt:.4f}, Sigma: {sigma_opt:.4f}, Gamma: {gamma_opt:.4f}")
