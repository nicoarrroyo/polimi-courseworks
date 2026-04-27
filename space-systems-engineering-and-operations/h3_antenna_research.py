import numpy as np
import matplotlib.pyplot as plt

# --- 1. SMOS Mission Global Parameters ---
R_earth     = 6371              # [km] Earth radius
k           = 1.38 * 10**(-23)  # [J/K] Boltzmann Constant
BER         = 10**(-7)          # Required ESA Downlink BER 
L_i         = 1.85              # [dB] Implementation Losses (TODO-EST)
L_add       = 7.0               # [dB] Additional/Atmospheric Losses (TODO-EST)
T_Rgs       = 250.0             # [K] Ground Station System Noise Temp (TODO-EST)

# --- 2. Downlink Scenarios Setup ---
# Scenario 0: X-band Science Data to Svalbard GS
# Scenario 1: X-band Science Data to ESAC GS
# Scenario 2: S-band TT&C HKTM to Kiruna GS

# Frequencies [GHz] 
# X-band allocation (8025-8400 MHz) , S-band downlink (2200-2290 MHz) 
f_d         = np.array([8.2, 8.2, 2.245])     

# Data Rates [Mbps]
# X-band Science: 16.8 Mbps , S-band max HKTM: 722 kbps 
R_b         = np.array([16.8, 16.8, 0.722])   

# Ground Station Dish Diameters [m]
# Svalbard: 11m, ESAC: 3.5m, Kiruna: 15m
D_gs        = np.array([11.0, 3.5, 15.0])    

# Spacecraft Antenna Properties
# X-band: Moderate-gain isoflux  (TODO-EST 6 dBi)
# S-band: LGA peak gain ~4 dBi 
antenna_gain = np.array([6.0, 6.0, 4.0]) # [dBi]

# Spacecraft Transmit Power [W]
# X-band: TWTA (TODO-EST 15W), S-band: SSPA  (TODO-EST 5W)
P_Tss = np.array([15.0, 15.0, 5.0])      

# --- 3. Budget Calculations ---
P_Tss_dBW   = 10 * np.log10(P_Tss)     # [W] -> [dBW]
EIRP        = antenna_gain + P_Tss_dBW # [dBW] Effective Isotropic Radiated Power

G_gs        = 17.8 + 20 * np.log10(D_gs) + 20 * np.log10(f_d)   # [dB] Ground Station Gain
noise_gs    = 10 * np.log10(k) + 10 * np.log10(T_Rgs)           # [dBW/Hz] GS Noise
GRN0R       = G_gs - noise_gs                                   # [dB/K] Gain-to-noise density ratio

# Evaluate altitudes around the SMOS nominal 763 km orbit 
alt_range   = np.linspace(500, 1000, 500)                    # [km] 
slant_range = np.sqrt((R_earth + alt_range)**2 - R_earth**2) # [km] Approximation for overhead
lambda_     = (3 * 10**8) / (f_d * 10**9)                    # [m] Wavelength

# Free Space Loss [dB] 
# (Shape is 500 altitudes x 3 scenarios)
L_fs_dB     = 20 * np.log10(4 * np.pi * slant_range[:, None] * 10**3 / lambda_) 

# Initialize arrays for Carrier-to-Noise and Energy-per-bit
CN0R        = np.zeros((500, 3))
EBN0        = np.zeros((500, 3))

for i in range(3):
    CN0R[:, i] = EIRP[i] + GRN0R[i] - L_fs_dB[:, i] - L_add # [dBHz] Carrier to Noise Density Ratio 
    EBN0[:, i] = CN0R[:, i] - 10 * np.log10(R_b[i] * 10**6) # [dB] Energy per bit-to-noise density

# Find required EBN0 for the target BER
EBN0_est    = np.linspace(0, 15, 100)                   # [dB]
P_e         = 0.5 * np.exp(-10**(EBN0_est / 10))        # Simplified Prob of Error curve
EBN0_req    = np.interp(BER, P_e[::-1], EBN0_est[::-1]) # [dB] Required E_b/N_0

link_margin = EBN0 - EBN0_req - L_i                     # [dB]

# --- 4. Plotting ---
plt.figure(figsize=(10, 6))

plt.plot(alt_range, link_margin[:, 0], 'b-', label='X-band -> Svalbard (11m)')
plt.plot(alt_range, link_margin[:, 1], 'r-', label='X-band -> ESAC (3.5m)')
plt.plot(alt_range, link_margin[:, 2], 'g-', label='S-band -> Kiruna (15m)')

# Plot SMOS actual altitude
plt.axvline(x=763, color='k', linestyle=':', label='SMOS Nominal Altitude (763km)')
plt.axhline(y=0, color='k')

plt.ylabel('Link Margin [dB]')
plt.xlabel('Altitude [km]')
plt.title('SMOS TTMTC Downlink Link Margins vs Altitude')
plt.legend(loc='lower left')
plt.grid(True)
plt.show()
