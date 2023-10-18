import numpy as np
import matplotlib.pyplot as plt
from moduli import ruote_fattorigrippaggio as rfg
from moduli.ruote_interpolazione import interp as ip
from datetime import datetime
from tabulate import tabulate
import os
import shutil
import pandas as pd

################################################            VARIABILI GLOBALI           ################################################  

start_time = datetime.now()
out_file = open("Grippaggio.txt", "w")                     # scrivo il file txt

####    Importazione
VG = pd.read_csv('variabili_globali.csv', index_col=0)
CNC = pd.read_csv('caratteristiche_nc.csv', index_col=0)
VC = pd.read_csv('variabili_cinematiche.csv', index_col=0)
PC = pd.read_csv('correzione.csv', index_col=0)
IP = pd.read_csv('ipotesi_progetto.csv', index_col=0)
IF = pd.read_csv('influence_factorsP.csv', index_col=0)

####    VARIABILI GLOBALI
theta_p = (VG.loc['r1','valori'])                                               # [rad] angolo di pressione 
beta = 0 
DentiP, DentiR = (VG.loc['r2','valori']), (VG.loc['r3','valori'])               # denti ruota e pignone
u_real = DentiR/DentiP                                                          # gear ratio 
m_uni =  (VG.loc['r5','valori']) 

alpha_wt = alpha_wt = (PC.loc['r21','pignone'])                                     # [rad] angolo di pressione trasversale di funzionamento
alpha_wt_deg = np.degrees(alpha_wt) 
rho_rel =  (PC.loc['r23','pignone'])                                                # [mm] curvatura relativa


Ft = (PC.loc['r18','pignone'])                                                      # [N] forza circonferenziale
b = (PC.loc['r12','pignone'])                                                       # [mm] larghezza di fascia

####    DIMENSIONI&VELOCITÀ
vel_rot_p = (VC.loc['r1','pignone'])                            # [1/s] velocità di rotazione pignone
raggio_corretto_pignone = (PC.loc['r7','pignone'])*0.001        # [m]
v_p = vel_rot_p*raggio_corretto_pignone                         # [m/s] velocità lineare al diametro primitivo/di lavoro pignone
diametro_p_p = 2*(CNC.loc['r1','pignone'])                      # [mm] diametro primitiva pignone
alphat1 = (PC.loc['r22','pignone'])                             # [rad] angolo in testa del pignone

vel_rot_r = (VC.loc['r1','ruota'])                              # [1/s] velocità di rotazione ruota
raggio_corretto_ruota = (PC.loc['r7','ruota'])*0.001            # [m]
v_r = vel_rot_r*raggio_corretto_ruota                           # [m/s] velocità lineare al diametro primitivo/di lavoro ruota
diametro_p_r = 2*(CNC.loc['r1','ruota'])                        # [mm] diametro primitiva ruota
alphat2 = (PC.loc['r22','ruota'])                               # [rad] angolo in testa del pignone

####    IPOTESI DI PROGETTO
####    Materiale ruote 16NiCr11
E_modulus = (IP.loc['r1','valori'])                             # [N/mm2] modulo elastico
v_poisson = (IP.loc['r2','valori'])                             # coefficiente di Poisson
rho_val = (IP.loc['r10','valori'])                              # [kg/m3] densità
lambda_val = (IP.loc['r11','valori'])                           # [N/sK] capacità termica
c_val = (IP.loc['r12','valori'])                                # [Nm/kgK] calore specifico
Ra = (IP.loc['r13','valori'])

####    Lubrificatione
eta_oil = (IP.loc['r20','valori'])                          # [mPas] viscosità dinamica                   
C = (IP.loc['r21','valori'])                                # polialfaolfine/esteri                  
                                            

####     Coefficienti calcolati
KA = (IF.loc['r1','valori'])
KV = (IF.loc['r2','valori'])
KH_beta = (IF.loc['r3','valori'])
KH_alpha = (IF.loc['r5','valori'])
KMP =  1










################################################            GRIPPAGGIO V.VULLO GEARS VOL3           ################################################  
####    Grafico e calcolo della velocità di slip
inizio = (PC.loc['r25','ruota'])
fine = (PC.loc['r24','pignone'])
g1 = np.arange(-inizio, 0, 0.01)
vg1 = np.array([np.abs(rfg.vg(v_p, gg, 0,  diametro_p_p)) for gg in g1])
vg_1 = np.max(vg1)

g2 = np.arange(0, fine, 0.01)
vg2 = np.array([rfg.vg(v_r, gg, 0,  diametro_p_r) for gg in g2])
vg_2 = np.max(vg2)

plt.figure(num=1)
plt.plot(g1, vg1, "b", label="pignone")
plt.plot(g2, vg2, "r", label ="ruota")
plt.xlabel('Segmento dei contatti [mm]')
plt.ylabel(r'$v_{g1}, v_{g2}$, slip velocity [m/s]')
plt.grid()
plt.legend()
plt.savefig('slipvelocity.png', dpi=300)

tab_vel = [["pitch line velocity", round(v_p,2)],
           ["slip velocity pignone", round(vg_1,2)],
           ["slip velocity ruota", round(vg_2,2)]]

st0 = "velocità"
st00 = tabulate(tab_vel)
out_file.write(st0 + '\n'), print(st0)
out_file.write(st00 + '\n'), print(st00)





################################################            Coefficienti di contatto          ################################################  
BM1 = np.sqrt(lambda_val*rho_val*c_val*0.001)
BM2 = BM1





################################################            flash temperature theta_fl          ################################################
################################################            Thermoelastic factor & Reduced modulus          ################################################      
Er = (E_modulus)/(1-v_poisson**2)

XM = ((Er**(1/4))/BM1)*np.sqrt(1000)





################################################            Approacch factor          ################################################      
XJ = 1





################################################           Load sharing factor          ################################################         
gamma_A = -u_real*(np.tan(alphat2)/np.tan(alpha_wt) - 1) 
gamma_B = (np.tan(alphat1)/np.tan(alpha_wt)) - 1 - (2*np.pi/(DentiP*np.tan(alpha_wt)))
gamma_D = -u_real*(np.tan(alphat2)/np.tan(alpha_wt)) - 1 + (2*np.pi/(DentiP*np.tan(alpha_wt)))
gamma_E = (np.tan(alphat1)/np.tan(alpha_wt) - 1)

X_gammaAB = (7-2)/15 + (1/3)*(g1-gamma_A)/(gamma_B-gamma_A)
X_gamma = 1
X_gammaDE = (7-2)/15 + (1/3)*(gamma_E-g2)/(gamma_E-gamma_D)





################################################           Geometrical factor          ################################################         
####    Angle factor
x_betaalpha = ip(20, 0.978, 22, 1.007, alpha_wt_deg)               # interpolazione dalla tabella 7.1

XG1 = rfg.calcola_XG(x_betaalpha, u_real, gamma_A, gamma_E)
XG2 = rfg.calcola_XG(x_betaalpha, u_real, gamma_E, gamma_A)





################################################          coefficiente d'attrito          ################################################         
####    transverse unit load
wbt = KA*KV*KH_beta*KH_alpha*KMP*(Ft/b)

####    cumulative velocity vector
vsigmac = 2*v_p*np.sin(alpha_wt)

####    local radius of curvature
rho1 = ((1+gamma_A)/(1 + u_real))*(diametro_p_p/2)*np.sin(alpha_wt)
rho2 = ((u_real-gamma_E)/(1 + u_real))*(diametro_p_r/2)*np.sin(alpha_wt)

####    local relative radius of curvature
rho_rel = (rho1*rho2)/(rho1 + rho2)

####    lubricant factor
XL = C * eta_oil**(-0.05)

####    roughness factor
XR = (Ra)**(0.25)

####    Friction coefficient
MU_m = rfg.calcola_mum(XL, XR, vsigmac, rho_rel, wbt)





################################################          theta_fl          ################################################         
theta_fl1 = rfg.calcola_thetafl(MU_m, XM, XJ, XG1, X_gamma, wbt, v_p, diametro_p_p)
theta_fl2 = rfg.calcola_thetafl(MU_m, XM, XJ, XG2, X_gamma, wbt, v_p, diametro_p_r)

tab_thetafl = [["Coefficienti di contatto", round(BM1,2)],
                ["Thermoelastic factor", round(XM,2)],
                 ["Approacch factor", XJ], 
                 ["Load sharing factor", X_gamma],
                 ["Angle factor", x_betaalpha],
                 ["Posizione adimensionale gammaA", round(gamma_A, 2) ], 
                 ["Posizione adimensionale gammaE", round(gamma_E, 2) ],
                 ["Geometrical factor1", round(XG1,2)],
                 ["Geometrical factor2", round(XG2,2)],
                  ["transverse unit load", round(wbt,2)],
                  ["cumulative velocity vector", round(vsigmac,2)],
                  ["local relative radius of curvature", round(rho_rel,2)],
                  ["lubricant factor", round(XL,4)],
                  ["roughness factor", round(XR,4)],
                  ["Friction coefficient", round(MU_m,4)],
                  ["flash temperature1", round(theta_fl1,2)],
                    ["flash temperature2", round(theta_fl2,2)],
                 ]

st2 = "flash temperature"
st3 = tabulate(tab_thetafl)
out_file.write(st2 + '\n'), print(st2)
out_file.write(st3 + '\n'), print(st3)





################################################          grafico temperatura flash - velocità          ################################################         
vp1 = np.arange(0, v_p, 0.01)
theta_fl1_plot = [rfg.calcola_thetafl(MU_m, XM, XJ, XG1, X_gamma, wbt, vv, diametro_p_p) for vv in vp1]
theta_fl2_plot = [rfg.calcola_thetafl(MU_m, XM, XJ, XG2, X_gamma, wbt, vv, diametro_p_r) for vv in vp1]

plt.figure(num=2)
plt.plot(vp1, theta_fl1_plot , "b", label="pignone")
plt.plot(vp1,theta_fl2_plot , "r", label ="ruota")
plt.xlabel('Pitch line velocity [m/s]')
plt.ylabel(r'$\Theta_{fl}$, flash temperature [C]')
plt.grid()
plt.legend()
plt.savefig('flashtemeprature_v.png', dpi=300)

################################################          grafico temperatura flash - carico          ################################################         
wbt1 = np.arange(0, wbt, 0.01)
theta_fl1_plot1 = [rfg.calcola_thetafl(MU_m, XM, XJ, XG1, X_gamma, ww, v_p, diametro_p_p) for ww in wbt1]
theta_fl2_plot1 = [rfg.calcola_thetafl(MU_m, XM, XJ, XG2, X_gamma, ww, v_p, diametro_p_r) for ww in wbt1]

plt.figure(num=3)
plt.plot(wbt1, theta_fl1_plot1 , "b", label="pignone")
plt.plot(wbt1,theta_fl2_plot1 , "r", label ="ruota")
plt.xlabel('transverse unit load [N/mm]')
plt.ylabel(r'$\Theta_{fl}$, flash temperature [C]')
plt.grid()
plt.legend()
plt.savefig('flashtemeprature_w.png', dpi=300)





################################################          verifica della temperatura flash          ################################################         
####    Banda di pressione hertziana
bh = rfg.calcola_bh(Ft, Er, b, theta_p, rho_rel)

alpha_wn = np.arcsin((np.sin(alpha_wt)*np.cos(beta)))
wbn = wbt/(np.cos(alpha_wn)*np.cos(beta))

####  Numeri di Péclet
Pe1 = rfg.calcola_Pe(vg_1, bh*0.001, rho_val, c_val, lambda_val, np.pi/2)
Pe2 = rfg.calcola_Pe(vg_2, bh*0.001, rho_val, c_val, lambda_val, np.pi/2)

st4 = f"NON è possibile utilizzare la relazione generale 7.20: {round(Pe1,2)}, {round(Pe2,2)}"
st5 = f"È possibile utilizzare la relazione generale 7.20: {round(Pe1,2)}, {round(Pe2,2)}"
if Pe1<5 and Pe2<5:
    out_file.write(st4 + '\n'), print(st4)
else: 
    out_file.write(st5 + '\n'), print(st5)





################################################          interfacial bulk temperature          ################################################         
####    theta_fl_media
theta_fl_media = (theta_fl1 + theta_fl2)/(gamma_E-gamma_A)
st61 = f"La temperatura flash media: {round(theta_fl_media,2)}"
out_file.write(st61 + '\n'), print(st61)

####  lubrication system factor
XS = 0.2                                                    # a bagno d'olio

####  multiple mating pinion factor
XMP = 1

####    overral bulk temeprature
theta_oil = 85                                              # fig 8.4 vullo
theta_M = theta_oil + 0.47*XS*XMP*theta_fl_media

st6 = f"La temperatura all'interfaccia è: {round(theta_M,2)}"
out_file.write(st6 + '\n'), print(st6)

theta_Msp = 0.0021*v_p**2 - 0.1188*v_p + 77.088
st7 = f"Dalla formula sperimentale si ottiene invece: {round(theta_Msp,2)}"
out_file.write(st7 + '\n'), print(st7)





################################################          interfacial contactl temperature          ################################################         
theta_fl = max(theta_fl1, theta_fl2)

theta_B = theta_M + theta_fl

st8 = f"L'olio lubrificante deve resistere alla massima temperatura possibile di: {round(theta_B,2)}"
out_file.write(st8 + '\n'), print(st8)


################################################          RUMORE          ################################################ 
L = rfg.calcola_L(u_real, beta, KV, 1.3, 15)
st9 = f"Il rumore prodotto dal riduttore alla distanza di 1m è pari a: {round(L)} dB"
out_file.write(st9 + '\n'), print(st9)

beta_val = np.arange(0, np.radians(20), 0.01)
l_plot = [rfg.calcola_L(u_real, bb, KV, 1.3, 15) for bb in beta_val]

plt.figure(num=4)
plt.plot(beta_val, l_plot)
plt.xlabel("angolo d'elica [rad]")
plt.ylabel('Livello di rumorosità [dB]')
plt.grid()
plt.savefig('rumorosità.png', dpi=300)



out_file.close()
################################################            ORGANIZZAZIONE FILE        ################################################
# Ottieni il percorso della directory in cui si trova lo script
percorso_script = os.path.dirname(os.path.abspath(__file__))

# Costruisci il percorso completo della cartella di origine
cartella_origine = os.path.join(percorso_script)

# Crea una cartella per le figure, tabelle e log all'interno della directory dello script (senza errori se esistono già)
os.makedirs(os.path.join(cartella_origine, 'figure'), exist_ok=True)
os.makedirs(os.path.join(cartella_origine, 'tabelle'), exist_ok=True)
os.makedirs(os.path.join(cartella_origine, 'log'), exist_ok=True)

# Elenca tutti i file nella cartella di origine
files = os.listdir(cartella_origine)

# Sposta i file nelle cartelle corrispondenti
for file in files:
    if file.endswith('.png') or file.endswith('.jpg'):
        shutil.move(os.path.join(cartella_origine, file), os.path.join(cartella_origine, 'figure', file))
    elif file.endswith('.xlsx'):
        shutil.move(os.path.join(cartella_origine, file), os.path.join(cartella_origine, 'tabelle', file))
    elif file.endswith('.txt'):
        shutil.move(os.path.join(cartella_origine, file), os.path.join(cartella_origine, 'log', file))



end_time = datetime.now()
print('Elapsed Time: {}'.format(end_time - start_time))
plt.show()