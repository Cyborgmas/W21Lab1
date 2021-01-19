import spinmob as sm
import numpy as n
import scipy.optimize as opt
import os

csv_dir = os.path.dirname(__file__) + "/data/"

cst = sm.data.load(csv_dir + "constants.csv")
runs = sm.data.load(csv_dir + "run_data.csv")

P0 = cst[0][0] * 1000 # kPa -> Pa
P0_err = cst[1][0] * 1000
m = cst[2][0]
m_err = cst[3][0]
d = cst[4][0] / 100      # cm -> m
d_err = cst[5][0] / 100
V = cst[6][0]
V_err = cst[7][0]

def sinFct(t, A, B, OMEGA, gamma, C):
    return (n.exp(-gamma*t) * (A*n.cos(OMEGA*t) + B*n.sin(OMEGA*t))) + C

parameters = []
cov = []

for i in range(10):
    x = runs[i * 4]
    y = runs[(i * 4) + 1]
    data = opt.curve_fit(sinFct, x, y)
    parameters.append(data[0])
    cov.append(data[1])

avg_params =  sum(parameters) / 10

avg_cov = sum(cov) / 10

err_sq = [avg_cov[i][i] / 10 for i in range(5)]

x = n.linspace(0, 10, 250)
y = sinFct(x, avg_params[0], avg_params[1], avg_params[2], avg_params[3], avg_params[4])

x_data = [runs[i*4] for i in range(10)]
y_data = [runs[i*4 + 1] for i in range(10)]

x_data.append(x)
y_data.append(y)

data_names = ["Run #" + str(i) for i in range(10)]
data_names.append("Fitted Data")

#sm.plot.xy.data(x_data, y_data, label=data_names)
#sm.plot.xy.data(x, y, label="Fitted Data")

omega = avg_params[2]
period = 2*n.pi / omega

omega_err = n.sqrt(err_sq[2])
period_err = omega_err * (period / omega)

print("Period = %.4f +- %.4f" % (period, period_err))

A = n.pi * (d / 2)**2
A_err = d_err * A * 2 / d 

g = 9.81

P = P0 + (m * g / A)
P_err = n.sqrt(P0_err**2 + (m_err * g / A)**2 + (A_err*m*g/A**2)**2)

gamma_ = 4*n.pi**2 * m * V / (A**2 * P * period**2)
gamma_err = gamma_ * n.sqrt( (m_err / m)**2 + (V_err / V)**2 + (P_err / P)**2 + (2 * A_err / A)**2 + (2 * period_err / period)**2 )

print("Gamma = %.3f +- %.3f" % (gamma_, gamma_err))



