import math
import random as rd
import numpy as np
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from expressions import *

# markers - linestyle - color

def tv1():
    print("--- Teilversuch 1 ---")

    # degrees
    alpha_angles = [ 8.06, 15.42, 23.00 ]
    alpha_uncertainties = [ 0.12, 0.12, 0.12 ]
    beta_angles = [ 7.24, 13.75, 20.44 ]
    beta_uncertainties = [ 0.09, 0.09, 0.09 ]

    '''
    for n in range(3):
        print(math.sin(alpha_angles[n] / 180 * math.pi) / (n + 1))
    
    print()

    for n in range(3):
        print(math.sin(beta_angles[n] / 180 * math.pi) / (n + 1))

    exit()
    '''

    theta_alpha_vars = [ Var(alpha_angles[i] / 180 * math.pi, alpha_uncertainties[i] / 180 * math.pi, f"alpha_{i}") for i in range(3) ]
    theta_beta_vars =  [ Var(beta_angles[i] / 180 * math.pi, beta_uncertainties[i] / 180 * math.pi, f"beta_{i}") for i in range(3) ]

    print(f"K_alpha-Linie:")
    tv1_calculate_lambdas(theta_alpha_vars, 71.08e-12)

    print(f"K_beta-Linie:")
    tv1_calculate_lambdas(theta_beta_vars, 63.09e-12)
    
def tv1_calculate_lambdas(theta_vars, literature_value):
    d_2 = Const(2 * 282.01e-12)
    variables = [ ]
    total_var = None
    values = [ ]
    for n in range(3):
        l = Div(Mult(d_2, Sin(theta_vars[n])), Const(n + 1))
        value = l.eval()
        uncertainty = gaussian(l, [ theta_vars[n] ])
        print(f"lambda = {value:8.4g}, \u0394lambda = {uncertainty:8.04g}")
        values.append(value)

    '''
        newvar = Var(value, uncertainty, f"lambda_{n + 1}")
        if total_var is None:
            total_var = newvar
        else:
            total_var = Add(total_var, newvar)
        variables.append(newvar)

    total_var = Div(total_var, Const(3))
    total = total_var.eval()
    uncertainty = gaussian(total_var, variables)
    '''

    total = np.mean(values)
    uncertainty = np.std(values)
    print(f"Durchschnitt: lambda = {total:8.4g}, \u0394lambda = {uncertainty:8.4g}")
    dev = abs(total - literature_value) / uncertainty
    print(f"Abweichung vom Literaturwert: {dev:.1f} Standardabweichungen")

def tv4():
    print("--- Teilversuch 4 ---")

    print("Auswertung vor Ort:")

    print("Alpha-Linien:")
    e = [ 3.81, 6.02, 7.74, 17.12, 21.83, 9.34 ]
    Z = [ 22, 26, 29, 42, 47, 79 ]
    eval_tv4(e, Z, 2)

    print("Beta-Linien:")
    e = [ 17.06, 17.06, 17.12, 19.35, 24.77, 11.10 ]
    eval_tv4(e, Z, 3)

# n_1 = 2 => alpha, n_1 = 3 => beta
def eval_tv4(e, Z, n_1):
    R_y = 13.6e-3 # in keV
    for i in range(len(e)):
        sqrt_ratio = (e[i] / R_y) ** 0.5
        sigma = Z[i] - sqrt_ratio * (1 - n_1 ** -2) ** -0.5

        print(f"Z = {Z[i]}, n_1 = {n_1}: sqrt_ratio = {sqrt_ratio:8.2f}, sigma = {sigma:8.2f}")

tv1()
tv4()
