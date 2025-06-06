import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# Load the hare and lynx data
data = pd.read_csv('hare_lynx.csv')

# Extract the data
years = data['Year'].values
hares = data['Hare'].values
lynxes = data['Lynx'].values

# Lotka-Volterra model differential equations
def lotka_volterra(y, t, alpha, beta, delta, gamma):
    x, y = y
    dxdt = alpha * x - beta * x * y
    dydt = delta * x * y - gamma * y
    return [dxdt, dydt]

# Cost function to minimize
def residuals(params, t, data):
    alpha, beta, delta, gamma = params
    y0 = data[0]
    sol = odeint(lotka_volterra, y0, t, args=(alpha, beta, delta, gamma))
    return np.sum((sol[:, 0] - data[:, 0])**2 + (sol[:, 1] - data[:, 1])**2)

# Initial guess for parameters
initial_guess = [0.1, 0.1, 0.1, 0.1]

# Combine hare and lynx data for optimization
combined_data = np.column_stack((hares, lynxes))

# Optimization
result = minimize(residuals, initial_guess, args=(years, combined_data), method='L-BFGS-B', bounds=[(0, 1), (0, 1), (0, 1), (0, 1)])
alpha_opt, beta_opt, delta_opt, gamma_opt = result.x

# Solve Lotka-Volterra with optimized parameters
y0 = combined_data[0]
sol = odeint(lotka_volterra, y0, years, args=(alpha_opt, beta_opt, delta_opt, gamma_opt))
hares_pred, lynxes_pred = sol.T

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(years, hares, 'o-', label='Hare Data')
plt.plot(years, lynxes, 'o-', label='Lynx Data')
plt.plot(years, hares_pred, '-', label='Hare Model Prediction')
plt.plot(years, lynxes_pred, '-', label='Lynx Model Prediction')
plt.xlabel('Year')
plt.ylabel('Population')
plt.title('Lotka-Volterra Model Fit to Hare and Lynx Data')
plt.legend()
plt.grid(True)
plt.show()
