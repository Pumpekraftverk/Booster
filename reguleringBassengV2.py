"""
Reguleringssystem for "basseng" til Fjone kraftverk
"""

# %% Importerer pakker

import matplotlib.pyplot as plt
import numpy as np

# %% Viktig parametere til simulering

A = 9                 # Tversnittareal på basseng [m**2]
T = 5                 # Tidskonstant [s]
t_delay = 1           # Forsinkelse [s]

# %% Verdier for boosterpumpe 

u_max = 20            # Maks flow [m**3/s]
u_min = 0             # Minste flow [m**3/s]

# %% Verdier for hovedpumpe

Fut = 18              # Flow i hovedpumpe [m**3/s]
amplitude_for = 0.5   # Ampltitude for forstyrrelse/ målestøy
periodetid_for = 50   # Periodetid for forstyrrelse/ målestøy

# %% Høydeverdier for bassenget

h_min = 0             # Minste høyde [m]
h_max = 14            # Maks høyde [m]
h_lav = 10            # Laveste høyde for å unngå kavitasjon [m]
h_sp = 12             # Settpunkt [m]
h_k = 12              # Initialbetingelse til høyde [m]

# %% Tidsinnstillinger

dt = 0.1                               # Tidssteg [s]
t_start = 0                            # Starttid [s]
t_stop = 5*60                          # Sluttid [s]
N_sim = int((t_stop - t_start)/dt) + 1 # Antall iterasjoner til simuleringsløkka

# %% Tidsforsinkelse

u_delayed_init = Fut # Startflow for booster [m**3/s]
N_delay = int(round(t_delay/dt)) + 1 # Str på delay array
delay_array_k = np.zeros(N_delay) + u_delayed_init # Initaliserer delay array

# %% Initialiserer arrays for plotting

u_array = np.zeros(N_sim)    # Booster
p_array = np.zeros(N_sim)    # Hovedpumpe
h_array = np.zeros(N_sim)    # Høyde til bassenget
h_sp_array = np.zeros(N_sim) # Settpunkt til bassenget
t_array = np.zeros(N_sim)    # Tid
ek_array = np.zeros(N_sim)   # Avvik mellom settpunkt og høyde på bassenget

# %% PI-tuning med skogestad-metode for rask nivåregulering

def skogestad_metode(deltaH, deltaF, A):
    Tc = A*(deltaH/deltaF)
    Ki = 1/A
    
    Kp = 1/(Ki*Tc)
    Ti = 2*Tc
   
    return Kp, Ti

#%% Parametere til Skogestad-metoden

deltaH = 2            # Maks endring i høydenivå til basseng [m]
deltaF = 2            # Maks forskjell i flow på booster og hovedpumpe [m**3/s]

Kp, Ti = skogestad_metode(deltaH, deltaF, A)

#%% Pådrag fra forover- og tilbakekobling til booster

u_man = 0           # Manuelt pådrag [m**3/s]
u_i = 0             # Initaliserer I-leddet til PI-regulator

def regulator_kontroller(Kp, Ti, u_i, u_man, u_forover):

    # P-leddet
    u_p = Kp*e_k

    # I-leddet
    u_i_k = (Kp/Ti)*dt*e_k
    u_i += u_i_k

    # Anti wind-up
    u_i = np.clip(u_i, u_min - u_forover, u_max - u_forover)

    # Totalt pådrag
    u_k = u_forover + u_man + u_p + u_i
 
    # Setter grenseverdier
    u_k = np.clip(u_k, u_min, u_max)

    return u_k, u_i

#%% Tidsforsinkelse med tidskonstant for å simulere ikke-ideelle forhold i systemet

def system_prosess(u_k, delay_array_k, Fb):
    
    # Forsinkelse
    u_delayed_k = delay_array_k[-1]
    delay_array_k[1:] = delay_array_k[0:-1]
    delay_array_k[0] = u_k
    delay_array_kp1 = delay_array_k
    
    # Eulers forovermetode for tidskonstantdynamikk med tidsforsinkelse
    dFbdt = (u_delayed_k - Fb)/T
    Fb_kp1 = dFbdt*dt + Finn
    
    return Fb_kp1, delay_array_kp1

#%% Forstyrrelser under pumping fra hovedpumpe

def forstyrrelse(amplitude, periodetid, tid):
    A = amplitude
    Tp = periodetid
    
    return A*np.sin((2*np.pi*tid)/Tp)

#%% Simuleringsløkke

Finn = u_delayed_init    # Initialbetingelse til booster [m**3/s]

for k in range(0, N_sim):
    
    # Tid
    t_k = k*dt
    
    # Settpunkt
    h_sp_k = h_sp
    
    # Error
    e_k = h_sp_k - h_k
    
    # Forstyrrelser i hovedpumpe
    Fut_k = Fut + forstyrrelse(amplitude_for, periodetid_for, t_k)
    
    # Måling av flow til hovedpumpe som er foroverkoblet
    u_forover = Fut 
    
    # Pådrag fra regulator til boosterpumpe
    u_k, u_i = regulator_kontroller(Kp, Ti, u_i, u_man, u_forover)
    
    # Lagrer verdier i arrayer til plotting
    ek_array[k] = e_k
    t_array[k] = t_k
    h_array[k] = h_k
    u_array[k] = Finn
    p_array[k] = Fut_k
    h_sp_array[k] = h_sp_k
    
    # Pådrag fra boosterpumpe
    Finn_kp1, delay_array_kp1 = system_prosess(u_k, delay_array_k, Finn)
    
    # Eulers metode for å finne høydeendring
    dhdt_k = (Finn - Fut_k)/A
    h_kp1 = dhdt_k*dt + h_k
    
    # Time index shift
    h_k = h_kp1
    Finn = Finn_kp1
    delay_array_k = delay_array_kp1
    
    # Setter grenseverdier på høyde til bassenget
    h_k = np.clip(h_k, h_min, h_max)
    
#%% Printer snitt og makserror

t0 = 50
t1 = 250

mean_h = np.mean(h_array[int(t0/dt):int(t1/dt)])
mean_h_sp = np.mean(h_sp_array[int(t0/dt):int(t1/dt)])
mean_e = abs(mean_h_sp - mean_h)

#print('Avviket mellom', t0, 'sek og', t1, 'sek er', round(mean_e, 3), '[m]')

mean_e_array = abs(h_sp_array-h_array)
print('Gjennomsnittet til avviket er:', abs(round(np.mean(mean_e_array),3)), '[m]')
print('Maks avvik er:', abs(round((np.max(mean_e_array)),3)), '[m]')

# %% Plotting

plt.close('all')
plt.figure(1)

plt.subplot(2, 1, 1)
plt.plot(t_array, h_sp_array, 'r', label='h_sp')
plt.plot(t_array, h_array, 'b', label='h')
plt.plot(t_array, np.zeros(N_sim)+h_lav, 'c',
          label='h_l')
plt.legend()
plt.grid()
plt.xlim(t_start, t_stop)
plt.xlabel('t [s]')
plt.ylabel('høyde [m]')

plt.subplot(2, 1, 2)
plt.plot(t_array, u_array, 'm', label='Booster')
plt.plot(t_array, p_array, 'g', label='Pumpe')
plt.legend()
plt.grid()
plt.xlim(t_start, t_stop)
plt.xlabel('t [s]')
plt.ylabel('Vannføring [m**3/s]')

plt.savefig('sim_basseng2.png')
plt.show()