import numpy as np
import matplotlib.pyplot as plt
from moduli import ruote_fattorigeometrici as rfge
from moduli.ruote_interpolazione import interp as interp
from tabulate import tabulate
import pandas as pd
import openpyxl
import sqlite3
import os
import shutil
from datetime import datetime

################################################            CORPO            ################################################         
start_time = datetime.now()

out_file = open("Geometria.txt", "w")                       # scrivo il file log

####    INIZIALIZZAZIONE DELLE VARIABILI GLOBALI

materiale_path = "ipotesipreliminari/materiale.xlsx"
workbook = openpyxl.load_workbook(materiale_path)           # Inizializzo i valori caratteristici necessari del materiale
sheet = workbook.active                                     # Questo seleziona il foglio attivo
HV = sheet["C4"].value                                      # Durezza Vickers

# bending
Af = sheet["C16"].value
Bf = sheet["C17"].value
sigma_flim = HV*Af + Bf                              # [N/mm2] nominal stress number for bending


theta = np.radians(20)                                      # Angolo di pressione
beta = np.radians(20)                                       # Angolo d'elica  





####    Definizione del fattore moltiplicativo
x = (1 + np.mean(np.array([8,1,3])))/10
stu1 = f"Il tuo fattore moltiplicativo è: {x}"
out_file.write(stu1 + '\n'), print(stu1)





####    Motore in ingresso
P_in = x*30*1e3                     # W
P_re = 0.98*0.99*0.98*P_in          # W
omega_in = x*13000                  # rpm  
omega_in_s = 2*np.pi*omega_in/60    # rad/s
tau_calc = x*0.5 + 0.2

omega_input = input("Se il riduttore è multistadio, inserisci la velocità di rotazione corrispondente in [rpm]: ")
omega_in = float(omega_input)
omega_in_s = 2*np.pi*omega_in/60

tab_car = [ ["Pin [W]", P_in],
     ["Pre [W]", round(P_re,2)], 
     ["Velrot [rpm]", omega_in],
     ["Velrot [rad/s]", round(omega_in_s, 2)],
     ["tau", tau_calc],]

stu2 = "Le caratteristiche del motore in ingresso saranno le seguenti:"
stu3 = tabulate(tab_car, headers=[])
out_file.write(stu2 + '\n'), print(stu2)
out_file.write(stu3 + '\n'), print(stu3)

st00 = f"Il rapporto di trasmissione richiesto è: {tau_calc}"
out_file.write(st00 + '\n'), print(st00)

tau_input = input('Se il riduttore è multistadio, inserisci il valore del rapporto di riduzione che scegli: ')
tau = float(tau_input)





####    Numero minimo di denti & Fattore di ricoprimento
k_prime_pignone = 1 
k_prime_ruota = 1.25

prima_denti_pignone = rfge.z_min_ideale(k_prime_pignone, tau, theta) 

prima_denti_ruota = rfge.z_min_ideale(k_prime_ruota, tau, theta) 

prima_fattore_ricoprimento = rfge.epsilon_ideale(prima_denti_ruota, k_prime_ruota, prima_denti_pignone, k_prime_pignone, theta)

st1 = f"In prima approssimazione, il minimo numero di denti di pignone e ruota è: {int(prima_denti_pignone)} {int(prima_denti_ruota)}"
st2 = f"Con un fattore di ricoprimento associato di: {round(prima_fattore_ricoprimento, 2)}"
out_file.write(st1 + '\n'), print(st1)
out_file.write(st2 + '\n'), print(st2)





####    Plot del numero minimo di denti
tau_valori = np.arange(0, 1, 0.01)      # Ruote cilindriche esterne
z_min_valori_pignone = [2 * k_prime_pignone * (1 + np.sqrt(1 + tt * (1 + tt) * np.sin(theta)**2)) / ((2 + tt) * np.sin(theta)**2) for tt in tau_valori ]
z_min_valori_ruota = [2 * k_prime_ruota * (1 + np.sqrt(1 + tt * (1 + tt) * np.sin(theta)**2)) / ((2 + tt) * np.sin(theta)**2) for tt in tau_valori ]

fig, axs = plt.subplots(1, 2)
### Pignone
axs[0].plot(z_min_valori_pignone, tau_valori)
axs[0].plot([prima_denti_pignone, prima_denti_pignone], [0, 0.45], 'r')
axs[0].plot([0,prima_denti_pignone], [0.45, 0.45], 'r')
axs[0].set_xlabel('Numero minimo di denti')
axs[0].set_ylabel('Rapporto di trasmissione' r'$~\tau>0$')
axs[0].set_title('Per il pignone')
axs[0].grid(True)
### Ruota
axs[1].plot(z_min_valori_ruota, tau_valori)
axs[1].plot([prima_denti_ruota, prima_denti_ruota], [0, 0.45], 'r')
axs[1].plot([0,prima_denti_ruota], [0.45, 0.45], 'r')
axs[1].set_xlabel('Numero minimo di denti')
axs[1].set_ylabel("Rapporto di trasmissione" r"$~\tau>0$")
axs[1].set_title('Per la ruota')
axs[1].grid(True)

fig.suptitle('Numero minimo di denti', fontsize=16)
plt.tight_layout()  
plt.savefig('numerominimodidenti.png', dpi=300)




####   Minimi denti intagliabili
reale_pignone = rfge.z_min_intagliabili(k_prime_pignone, theta)
reale_ruota = rfge.z_min_intagliabili(k_prime_ruota, theta)

st3 = f"Il minimo numero di denti intagliabili per pignone e ruota è: {int(reale_pignone)}, {int(reale_ruota)}"
out_file.write(st3 + '\n'), print(st3)



####   Numeri pratici
pratico_pignone = reale_pignone*(5/6)
pratico_ruota = reale_ruota*(5/6)

epsilon_pratico = ((pratico_ruota + pratico_pignone)/(2*np.pi))*np.tan(theta)

st4 = f"Il minimo numero di denti PRATICO di pignone e ruota è: {int(pratico_pignone)}, {int(pratico_ruota)}"
st5 = f"Con un fattore di ricoprimento associato di: {round(epsilon_pratico, 2)}"
out_file.write(st4 + '\n'), print(st4)
out_file.write(st5 + '\n'), print(st5)





####   Reale numero di denti 
DentiP_input = input('Inserisci il numero dei denti del pignone: ') 
DentiP = int(DentiP_input)
DentiR_guess = DentiP + 15
toll = 1e-5
DentiR = int(rfge.z_min_reale(DentiP, DentiR_guess, tau, toll))

tau_re = DentiP/DentiR

epsilon_reale = ((DentiR + DentiP)/(2*np.pi))*np.tan(theta)

st6 = f"Si scelgono tuttavia un numero di denti di pignone e ruota pari a: {DentiP}, {DentiR}"
st7 = f"In modo da garantire un fattore di ricoprimento pari a: {round(epsilon_reale, 2)}"
st8 = f"Ed un rapporto di riduzione velocità pari a: {round(tau_re, 3)}"
out_file.write(st6 + '\n'), print(st6)
out_file.write(st7 + '\n'), print(st7)
out_file.write(st8 + '\n'), print(st8)





####   Velocità e momenti
vel_ruota_s = (DentiP/DentiR)*omega_in_s        #rad/s

momento_torc_pign = P_re/omega_in_s             #Nm
momento_torc_ruota = P_re/vel_ruota_s           #Nm

tab_car = [ ["M_t [Nm]", round(momento_torc_pign, 2), round(momento_torc_ruota, 2)],
     ["Velocità [rad/s]", round(omega_in_s,2), round(vel_ruota_s, 2)]]

st9 = "A queste ruote si associano le seguenti caratteristiche"
st10 = tabulate(tab_car, headers=["Pignone", "Ruota"])
out_file.write(st9 + '\n'), print(st9)
out_file.write(st10 + '\n'), print(st10)




####   Calcolo del modulo secondo Lewis
k_pignone = interp(22, 0.651, 24, 0.629, 23)
k_ruota = interp(50, 0.461, 60, 0.430, 51)
print(round(k_pignone,3), round(k_ruota,3))
lambda_val = 10

modulo_pignone = rfge.calcola_modulo(k_pignone, momento_torc_pign*1e3, lambda_val, (sigma_flim/2)/5)
modulo_ruota = rfge.calcola_modulo(k_ruota, momento_torc_ruota*1e3, lambda_val, (sigma_flim/2)/5)

st11 = f"I moduli associati secondo Lewis sono perciò: {round(modulo_pignone,3)} {round(modulo_ruota,3)}"
out_file.write(st11 + '\n'), print(st11)




####   Gandezze caratteristiche
modulo = input('Inserisci il valore del modulo scelto: ')
m = float(modulo)
a = m
u = 1.25*m

####   Raggi primitivi
r_pignone = rfge.raggi_primitive(m, DentiP)
r_ruota = rfge.raggi_primitive(m, DentiR)

####   Raggio di testa
Rtp = r_pignone + a
Rtr = r_ruota + a

####   Raggio di piede
Rpp = r_pignone - u
Rpr = r_ruota - u

####   Raggio fondamentale
rho_pignone = r_pignone*np.cos(theta)
rho_ruota = r_ruota*np.cos(theta)

####   Passo
passo_pignone = m*np.pi
passo_ruota = m*np.pi

####   Interasse
interasse = r_pignone+r_ruota

####   Spessore del dente sulla primitiva
spp = passo_pignone/2
spr = passo_ruota/2

tab_grand = [ ["Primitive [mm]", r_pignone, r_ruota],
     ["Testa [mm]", Rtp, Rtr], 
     ["Piede [mm]", round(Rpp,2), round(Rpr,2)],
     ["Fondamentali [mm]", round(rho_pignone,2), round(rho_ruota,2)],
     ["Passo [mm]", round(passo_pignone,2), round(passo_ruota,2)],
     ["Interasse [mm]", round(interasse,2), "="],
     ["Sp su Prim [mm]", round(spp,2), round(spr,2)],]

st9 = "A queste ruote si associano le seguenti caratteristiche"
st12 = tabulate(tab_grand, headers=["Pignone", "Ruota"])
out_file.write(st9 + '\n'), print(st9)
out_file.write(st12 + '\n'), print(st12)




####   Verifiche
####   Spessore del dente
ev_incognitap = rfge.ev_gamma(spp, r_pignone, theta)
ev_incognitar = rfge.ev_gamma(spr, r_ruota, theta)


####   Angolo di pressione in testa
theta_tp = np.arccos(rho_pignone/(Rtp))
theta_tr = np.arccos(rho_ruota/(Rtr))


spessore_testap = 2*(Rtp)*(ev_incognitap - rfge.ev(theta_tp))
spessore_testar = 2*(Rtr)*(ev_incognitar - rfge.ev(theta_tr))


spessore_lim = 0.2*m

tab_spesst = [["ev_gamma [rad]", round(ev_incognitap, 4), round(ev_incognitar, 4)], 
              ["Angolo in testa [rad]", round(theta_tp, 4), round(theta_tr, 4)],
              ["Spessore testa [mm]", round(spessore_testap, 2), round(spessore_testar, 2)],
              ["Spessore limite [mm]", spessore_lim, spessore_lim]]

st13 = "La verifica sullo spessore è superata!"
st14 = "Ritenta, spessore in testa troppo basso!"
st15 = tabulate(tab_spesst, headers=["Pignone", "Ruota"])

if round(spessore_testap, 2)>spessore_lim and round(spessore_testar, 2)>spessore_lim:
    out_file.write(st13 + '\n'), print(st13)
else:
    out_file.write(st14 + '\n'), print(st14)

out_file.write(st15 + '\n'), print(st15)





####   Interferenza
phi = np.pi/2 - theta

accesso_teorico = r_pignone*np.cos(phi) # T1C

recesso_teorico = r_ruota*np.cos(phi) # CT2

accesso_reale = np.sqrt(Rtr**2 - rho_ruota**2) - recesso_teorico

recesso_reale = np.sqrt(Rtp**2 - rho_pignone**2) - accesso_teorico 

tab1 = [ ["Teorico", round(accesso_teorico,3), round(recesso_teorico,3)],
     ["Reale", round(accesso_reale,3), round(recesso_reale,3)]]

st16 = "La verifica sull'interferenza è superata"
st17 = "Ritenta, c'è interferenza!"
st18 = tabulate(tab1, headers=["Accesso", "Recesso"])

if round(accesso_reale,3) < round(accesso_teorico,3) and round(recesso_reale,3) < round(recesso_teorico,3):
    out_file.write(st16 + '\n'), print(st16)
else:
    out_file.write(st17 + '\n') , print(st17)

out_file.write(st18 + '\n'), print(st18)




####   Fattore di ricoprimento
contatti = accesso_reale+recesso_reale

epsilon_verifica = contatti/(passo_pignone*np.cos(theta))

st19 = f"Il fattore di ricoprimento è {round(epsilon_verifica, 2)}"
st20 = f"Il fattore di ricoprimento non è abbastanza! {round(epsilon_verifica, 2)}"

if round(epsilon_verifica, 2)>1.2:
    out_file.write(st19 + '\n'), print(st19)
else:
    out_file.write(st20 + '\n'), print(st20)





####   Strisciamenti specifici
inizio = -accesso_reale
fine = recesso_reale
step = 0.01
delta = np.arange(inizio, fine, step)

ks_pignone_plot1 = [rfge.ks_pignone(tau_re, dd, r_pignone, theta) for dd in delta]
ks_ruota_plot1 = [rfge.ks_ruota(tau_re, dd, r_ruota, theta) for dd in delta]

plt.figure(num=2)
plt.plot(delta, ks_pignone_plot1, "b", label = "pignone")
plt.plot(delta, ks_ruota_plot1, "r", label = "ruota")
plt.xlabel('Segmento dei contatti [mm]')
plt.ylabel(r'$k_s$')
plt.title('Grafico degli strisciamenti specifici (x=0  x\'=0)')
plt.legend()
plt.grid()
plt.savefig('strisciamenti_nc.png', dpi=300)

ksmax_pignone = max(np.abs(np.min(ks_pignone_plot1)), np.abs(np.max(ks_pignone_plot1)))
ksmax_ruota = max(np.abs(np.min(ks_ruota_plot1)), np.abs(np.max(ks_ruota_plot1)))

delta_k = np.abs(ksmax_pignone - ksmax_ruota) 

tab2 = [ ["ksmax", round(ksmax_pignone, 3), round(ksmax_ruota,3)],
        ["delta k", round(delta_k, 3) ]]

st21 = "La prima verifica sugli strisciamenti è superata"
st22 = "La seconda verifica sugli strisciamenti è superata"
st23 = "Nessuna delle verifiche è superata, è necessaria la correzione!"
st24 = tabulate(tab2, headers=["Pignone", "Ruota"])

if ksmax_pignone < 1.4 and ksmax_ruota < 1.4:
    out_file.write(st21 + '\n'), print(st21)
elif delta_k<0.2:
    out_file.write(st22 + '\n'), print(st22)
else:
    out_file.write(st23 + '\n'), print(st23)

out_file.write(st24 + '\n'), print(st24)





################################################            DATABASE        ################################################
out_file.close()                                                       # chiudo e salvo il file log


####    Variabili globali
variabili_globali = {
    'variabili':["angolo di pressione [rad]", "denti pignone", "denti ruota", "rapporto di trasmissione", "modulo"],
    'valori':[ theta, DentiP, DentiR, tau_re, m]
}
index_labels1 = ['r1', 'r2', 'r3', 'r4', 'r5']

VG = pd.DataFrame(variabili_globali, index=index_labels1)

# Connessione al database SQLite
conn = sqlite3.connect("variabili_globali.db")

# Salva il DataFrame nel database
VG.to_sql('variabili_globali', conn, if_exists='replace', index_label='index')

# Chiudi la connessione
conn.close()

# Tabella su excel sempre facilmente consultabile
VG.to_excel('variabili_globali.xlsx', index=True)





####    Caratteristiche non corrette
caratteristiche_nc = {
    'variabili':["primitive", "fondamentali", "testa", "piede", "passo", "interasse", "sp su primitiva"],
    'pignone':[ r_pignone, round(rho_pignone,2), Rtp, round(Rpp,2), round(passo_pignone,2), round(interasse,2), round(spp,2)],
    'ruota':[r_ruota, round(rho_ruota,2), Rtr, round(Rpr,2), round(passo_pignone,2), round(interasse,2), round(spp,2)]
}

index_labels2 = ['r1', 'r2', 'r3', 'r4', 'r5', 'r6', 'r7']

CNC = pd.DataFrame(caratteristiche_nc, index=index_labels2)

# Connessione al database SQLite
conn1 = sqlite3.connect("caratteristiche_nc.db")

# Salva il DataFrame nel database
CNC.to_sql('caratteristiche_nc', conn1, if_exists='replace', index_label='index')

# Chiudi la connessione
conn1.close()

CNC.to_excel('caratteristiche_nc.xlsx', index=True)





####    Variabili cinematiche
variabili_cinematiche = {
    'variabili':["velocità di rotazione (rad/s)", "momento torcente"],
    'pignone':[round(omega_in_s,2), round(momento_torc_pign, 2)],
    'ruota':[round(vel_ruota_s, 2),  round(momento_torc_ruota, 2)]
    
}

index_labels3 = ['r1', 'r2']

VC = pd.DataFrame(variabili_cinematiche, index=index_labels3)

# Connessione al database SQLite
conn2 = sqlite3.connect("variabili_cinematiche.db")

# Salva il DataFrame nel database
VC.to_sql('variabili_cinematiche', conn2, if_exists='replace', index_label='index')

# Chiudi la connessione
conn2.close()

VC.to_excel('variabili_cinematiche.xlsx', index=True)









################################################            ORGANIZZAZIONE FILE        ################################################
# Ottieni il percorso della directory in cui si trova lo script
percorso_script = os.path.dirname(os.path.abspath(__file__))

# Costruisci il percorso completo della cartella di origine
cartella_origine = os.path.join(percorso_script)

# Crea una cartella per le figure, tabelle e log all'interno della directory dello script (senza errori se esistono già)
os.makedirs(os.path.join(cartella_origine, 'figure'), exist_ok=True)
os.makedirs(os.path.join(cartella_origine, 'tabelle'), exist_ok=True)
os.makedirs(os.path.join(cartella_origine, 'log'), exist_ok=True)
os.makedirs(os.path.join(cartella_origine, 'database'), exist_ok=True)


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
    elif file.endswith('.db'):
        shutil.move(os.path.join(cartella_origine, file), os.path.join(cartella_origine, 'database', file))




end_time = datetime.now()
print('Elapsed Time: {}'.format(end_time - start_time))
plt.show()
 

