import numpy as np


# These are the numbers you enter
# ----------------------------------
avg_pressure = 103386.45        # [Pa]
avg_wavelength = 632.8          # [nm]
avg_temp = 23                   # [C]
avg_humidity = 0.58             # [fraction]
avg_co2 = 408.05                # [ppm]

std_pressure = 40               # [Pa]
std_wavelength = 10             # [nm]
std_temp = 0.5                  # [C]
std_humidity = 0.01             # [fraction]
std_co2 = 5                     # [ppm]
# ----------------------------------


# These values are from ciddor
# ----------------------------------
a_0 = 1.58123e-6
a_1 = -2.9331e-8
a_2 = 1.1043e-10

b_0 = 5.707e-6
b_1 = -2.051e-8

c_1 = -2.376e-6
c_0 = 1.9898e-4

d = 1.83e-11

e = -7.65e-9

A = 1.2378847e-5
B = -0.019121316
C = 33.93711047
D = -6343.1645
    
k_0 = 238.0185
k_1 = 5792105
k_2 = 57.362
k_3 = 167917
    
w_0 = 295.235
w_1 = 2.6422
w_2 = -0.03238
w_3 = 0.004028

R = 8.31451

α = 1.00062
β = 3.14e-8
γ = 5.6e-7
# ----------------------------------


# Generates inputs from a normal dist for markov sampling
def generate_inputs():
    pressure = np.random.normal(avg_pressure, std_pressure)
    wavelength = np.random.normal(avg_wavelength, std_wavelength)
    temperature = np.random.normal(avg_temp, std_temp)
    humidity = np.random.normal(avg_humidity, std_humidity)
    co2 = np.random.normal(avg_co2, std_co2)
    return pressure, temperature, wavelength, co2, humidity

# determines compressibility
def compressibility(T, p, x_w):
    t = T - 273.15 
    return 1 - ((p / T) * (a_0 + (a_1 * T) + (a_2 * t**2) + ((b_0 + (b_1 * t)) * x_w) + ((c_0 + (c_1 * t)) * x_w**2))) + ((p / T)**2 * (d + (e * x_w**2)))

# computes index of refraction
def index(P, T, W, co2, H):
    t = T
    λ = W
    σ = 1 / λ
    T = t + 273.15 # Lord Kelvin
    p = P
    h = H

    # saturation vapor pressure of water vapor in air at T
    svp = np.exp((A * T * T) + (B * T) + (C) + (D / T))

    # enhancement factor of water vapor in air
    f = α + (β * p) + (γ * t * t)

    # molar fraction of water vapor in moist air
    x_w= ((f * h * svp) / p)

    Z = compressibility(T, P, x_w)
    Z_a = compressibility(288.15, 101325, 0)    # Dry Air
    Z_w = compressibility(293.15, 1333, 1)      # Pure Water Vapor

    M_w = 0.018015                                          # Molar mass of dry air (kg / mol)
    M_a = 0.001 * (28.9635 + (1.2011e-5 * (co2 - 400)))     # Molar mass of water vapor (kg/ mol)


    ρ_axs = ((101325 * M_a) / (Z_a * R * 288.15))           # density of standard air
    ρ_ws = ((1333 * M_a) / (Z_w * R * 293.15)) * (M_w / M_a)# density of of standard water vapor

    ρ_w = (p * M_w * x_w) / (Z * R * T)                     # density of water vaport component
    ρ_a = (p * M_a * (1 - x_w)) / (Z * R * T)               # density of the dry component of air


    # refractive index of water vapor at standard temperature and pressure.
    n_ws = 1 + (0.00000001 * (1.022 * (w_0 + (w_1 * σ* σ) + (w_2 * σ* σ* σ* σ) + (w_3 * σ* σ* σ* σ* σ* σ))))

    # refractive index of standard air at 15 C, 1 atm, 0% humidity, 450 ppm CO2
    n_as = 1 + ((0.00000001) * ((k_1 / (k_0 - (σ* σ))) + (k_3 / (k_2 - (σ* σ)))))

    # refractive index of standard air at 15 C, 1 atm, 0% humidity, co2 ppm CO2
    n_axs = 1 + ((n_as - 1) * (1 + (0.000000534 * (co2 - 450))))


    # black magic
    N_prop = 1 + ((ρ_a / ρ_axs) * (n_axs - 1)) + ((ρ_w / ρ_ws) * (n_ws - 1))
    return N_prop


def markov():  # Takes 10 thousand samples and returns average index of refraction and standard deviation
    values = []
    for i in range(10000):
        P, T, W, co2, H = generate_inputs()
        values.append(index(P, T, W, co2, H))
    values = np.asarray(values)
    return (np.average(values), np.std(values))

[n, δn] = markov()

print('The index of refraction of air in our laboratory is {:.7f} ± {:.7f}'.format(round(n, 8), round(δn, 7)))