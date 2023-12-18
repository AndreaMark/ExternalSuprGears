import numpy as np
import matplotlib.pyplot as plt
from moduli import ruote_fattorigeometrici as rfge
from moduli import ruote_fattoripitting as rfp
from moduli import ruote_correzione as rc
import pandas as pd
from tabulate import tabulate
import sqlite3
import os
import shutil
from datetime import datetime

################################################            CORPO            ################################################  
start_time = datetime.now()

out_file = open("Correzione.txt", "w")                                # scrivo il file log

######################################################    Importazione Valori dai database
####    VALORI GLOBALI
vg_path = "database/variabili_globali.db"
# Connessione
conn = sqlite3.connect(vg_path)
cur = conn.cursor()
cur.execute("SELECT valori FROM variabili_globali")
vg = [row[0] for row in cur.fetchall()]
theta_p = vg[0]                                       # Angolo di pressione
DentiP = vg[1]
DentiR = vg[2]                                        # Numero di denti di pignone e ruota
tau_re = DentiP/DentiR                                # Reale rapporto di trasmissione
m = vg[4]                                             # modulo e addendum
a = m
beta = 0  
conn.close()





####    CARATTERISTICHE GEOMETRICHE NON CORRETTE
cnc_path = "database/caratteristiche_nc.db"
conn = sqlite3.connect(cnc_path)
#   PIGNONE
cur1 = conn.cursor()
cur1.execute("SELECT pignone FROM caratteristiche_nc")
cncP = [row[0] for row in cur1.fetchall()]
#   RUOTA
cur2 = conn.cursor()
cur2.execute("SELECT ruota FROM caratteristiche_nc")
cncR = [row[0] for row in cur2.fetchall()]
####    Pignone
r_pignone = cncP[0]     # Primitiva
rho_pignone = cncP[1]   # Fondamentale
####    Ruota
r_ruota = cncR[0]       # Primitiva
rho_ruota = cncR[1]     # Fondamentale
conn.close()





####    VARIABILI CINEMATICHE
vc_path = "database/variabili_cinematiche.db"
conn = sqlite3.connect(vc_path)
#   PIGNONE
cur1 = conn.cursor()
cur1.execute("SELECT pignone FROM variabili_cinematiche")
vcP = [row[0] for row in cur1.fetchall()]
#   RUOTA
cur2 = conn.cursor()
cur2.execute("SELECT ruota FROM variabili_cinematiche")
vcR = [row[0] for row in cur2.fetchall()]

Momento_pignone = vcP[1]                        # Nm
Momento_ruota = vcR[1]                          # Nm
conn.close()





####    Raggi di testa di ruota e pignone NON CORRETTI
Rtp = r_pignone + a
Rtr = r_ruota + a

####    Calcolo del segmento dei contatti NON CORRETTO
phi = np.pi/2 - theta_p
accesso_reale = np.sqrt(Rtr**2 - rho_ruota**2) - r_ruota*np.cos(phi)
recesso_reale = np.sqrt(Rtp**2 - rho_pignone**2) - r_pignone*np.cos(phi)
delta = np.arange(-accesso_reale, recesso_reale, 0.01)





####    Correzione dentatura
correzione = rc.correzione_dentature(tau_re, delta, theta_p, 0.3, 0.2, rho_pignone, rho_ruota, r_pignone, r_ruota, DentiP, DentiR, m)
st97 = f"Il numero di iterazioni è: {correzione[10]}"
out_file.write(st97 + '\n'), print(st97)
####    Variabili post correzione
x_pignone_approssimato = round(correzione[0], 3)
x_ruota_approssimato = round(correzione[1], 3)
theta_l = correzione[2]
accesso_reale_new = correzione[3]
recesso_reale_new = correzione[4]

X = 2 * (x_pignone_approssimato+x_ruota_approssimato)/(DentiP + DentiR)
Y = (np.cos(theta_p)/np.cos(theta_l) - 1)
K = ((DentiP + DentiR)/2)*(X-Y)
clearence = (1/4 - K)*m                     #Vullo1 p256...

st98 = f"Fattori vullo {X}, {Y}, {K}"
out_file.write(st98 + '\n'), print(st98)

st99 = f"clearance, Km0 {clearence}, {K*m}"
out_file.write(st99 + '\n'), print(st99)



raggio_testa_pignone_corretto = 0.5*m*(DentiP + 2*(1+x_pignone_approssimato-K))
raggio_testa_ruota_corretto = 0.5*m*(DentiR + 2*(1+x_ruota_approssimato-K))
modulo_lavoro = correzione[5]
# working circular pith
passo_primitiva_lavoro = np.pi*m*(np.cos(theta_p)/np.cos(theta_l))
spessore_primitiva_lavoro = passo_primitiva_lavoro/2


raggio_primitiva_pignone_corretto = correzione[6]
raggio_primitiva_ruota_corretto = correzione[7]
raggio_piede_pignone_corretto = 0.5*m*(DentiP - 2*(5/4 - x_pignone_approssimato))
raggio_piede_ruota_corretto = 0.5*m*(DentiR - 2*(5/4 - x_ruota_approssimato))
y = correzione[8]
y_prime = correzione[9]
interasse = round(raggio_primitiva_pignone_corretto + raggio_primitiva_ruota_corretto, 2)

####    SEGMENTO DEI CONTATTI CORRETTO E PLOT
delta_new = np.arange(-accesso_reale_new, recesso_reale_new, 0.01)

ks_pignone_plot2 = [rfge.ks_pignone(tau_re, dd, r_pignone, theta_l) for dd in delta_new]
ks_ruota_plot2 = [rfge.ks_ruota(tau_re, dd, r_ruota, theta_l) for dd in delta_new]


titolo = 'Grafico degli strisciamenti specifici (x={}  x\'={})'.format(x_pignone_approssimato, x_ruota_approssimato)
plt.figure(num=3)
plt.plot(delta_new, ks_pignone_plot2, "b", label = "pignone")
plt.plot(delta_new, ks_ruota_plot2, "r", label = "ruota")
plt.xlabel('Segmento dei contatti [mm]')
plt.ylabel(r'$k_s$')
plt.title(titolo)
plt.legend()
plt.grid()
plt.savefig('strisciamenti_c.png', dpi=300)

ks_max_pignone_corretto = max(np.abs(np.min(ks_pignone_plot2)), np.abs(np.max(ks_pignone_plot2)))
ks_max_ruota_corretto = max(np.abs(np.min(ks_ruota_plot2)), np.abs(np.max(ks_ruota_plot2)))
delta_ks_corretto = np.abs(ks_max_pignone_corretto - ks_max_ruota_corretto)

tab2 = [ ["ksmax", round(ks_max_pignone_corretto, 3), round(ks_max_ruota_corretto,3)],
        ["delta k", round(delta_ks_corretto, 3), round(delta_ks_corretto, 3) ]]

st1 = "Le verifiche sullo strisciamento sono superate!"
st2 = tabulate(tab2, headers=["Pignone", "Ruota"])
out_file.write(st1 + '\n'), print(st1)
out_file.write(st2 + '\n'), print(st2)





####    Verifiche corrette
####    Spessore della testa del dente corretto
evgamma_corr_Pignone = rfge.ev_gamma(spessore_primitiva_lavoro, raggio_primitiva_pignone_corretto, theta_l)     #rfge.ev(theta_p) + (2/DentiP)*(np.pi/4 + x_pignone_approssimato*np.tan(theta_p)) 
evgamma_corr_Ruota = rfge.ev_gamma(spessore_primitiva_lavoro, raggio_primitiva_ruota_corretto, theta_l)         #rfge.ev(theta_p) + (2/DentiR)*(np.pi/4 + x_ruota_approssimato*np.tan(theta_p))

####    Angolo di pressione in testa
theta_t_corr_pignone = np.arccos(rho_pignone/raggio_testa_pignone_corretto)
theta_t_corr_ruota = np.arccos(rho_ruota/raggio_testa_ruota_corretto)

spessore_testa_corretto_pignone = 2*raggio_testa_pignone_corretto*(evgamma_corr_Pignone - rfge.ev(theta_t_corr_pignone))
spessore_testa_corretto_ruota = 2*raggio_testa_ruota_corretto*(evgamma_corr_Ruota - rfge.ev(theta_t_corr_ruota))

spessore_lim = 0.2*m

tab_spesst = [["ev_gamma [rad]", round(evgamma_corr_Pignone, 4), round(evgamma_corr_Ruota, 4)], 
              ["Angolo in testa [rad]", round(theta_t_corr_pignone, 4), round(theta_t_corr_ruota, 4)],
              ["Spessore testa [mm]", round(spessore_testa_corretto_pignone, 2), round(spessore_testa_corretto_ruota, 2)],
              ["Spessore limite [mm]", spessore_lim, spessore_lim]]

st3 = "La verifica sullo spessore è superata!"
st4 = "Ritenta, spessore in testa troppo basso!"
st5 = tabulate(tab_spesst, headers=["Pignone", "Ruota"])

if spessore_testa_corretto_pignone>spessore_lim and spessore_testa_corretto_ruota>spessore_lim :
    out_file.write(st3 + '\n'), print(st3)
else:
    out_file.write(st4 + '\n'), exit(st4)

out_file.write(st5 + '\n'), print(st5)





####    Interferenza
phi_new = np.pi/2 - theta_l

t1c = raggio_primitiva_pignone_corretto*np.cos(phi_new)     # Accesso Teorico

ct2 = raggio_primitiva_ruota_corretto*np.cos(phi_new)       # Recesso Teorico



tab1 = [ ["Teorico", round(t1c,3), round(ct2,3)],
     ["Reale", round(accesso_reale_new,3), round(recesso_reale_new,3)]]

st6 = f"La verifica sull'interferenza è superata {round((t1c + ct2),2)}, {round((accesso_reale_new + recesso_reale_new),2)}"
st7 = "C'è qualcosa che non va!"
st8 = tabulate(tab1, headers=["Accesso", "Recesso"])

if accesso_reale_new < t1c and recesso_reale_new < ct2:
    out_file.write(st6 + '\n'), print(st6)
else:
    out_file.write(st7 + '\n'), exit(st7)

out_file.write(st8 + '\n'), print(st8)





####    Fattore di ricoprimento
# Basic Fanelli
contatti_new = accesso_reale_new + recesso_reale_new              
passo = m*np.pi
epsilon_verifica_new = contatti_new/(passo*np.cos(theta_p))


# Vullo Vol1 eq 3.44 p94
# fatt1 = 1/(2*np.pi*np.cos(theta_p))
# fatt2 = np.sqrt((np.sin(theta_p)**2 * DentiR**2 ) + 4*k2*(DentiR + k2))
# fatt3 = np.sqrt((np.sin(theta_p)**2 * DentiP**2 ) + 4*k1*(DentiP + k1))
# fatt4 = (DentiP + DentiR)*np.sin(theta_p)

# epsilon_verifica_new = fatt1 * (fatt2 + fatt3 - fatt4)


st9 = f"Il fattore di ricoprimento è {round(epsilon_verifica_new, 1)}"
st10 = "Il fattore di ricoprimento non è abbastanza!"

if round(epsilon_verifica_new, 1)>=1.2:
    out_file.write(st9 + '\n'), print(st9)
else:
    out_file.write(st10 + '\n'), print(st9)





####    Grandezze caratteristiche
####    Spessori primitive di taglio
s_taglio_pignone = 2*r_pignone*(evgamma_corr_Pignone - rfge.ev(theta_p))
s_taglio_ruota = 2*r_ruota*(evgamma_corr_Ruota - rfge.ev(theta_p))

####    Spessori primitive di lavoro
#s_lavoro_pignone = 2*raggio_primitiva_pignone_corretto*(evgamma_corr_Pignone - rfge.ev(theta_l))
#s_lavoro_ruota = 2*raggio_primitiva_pignone_corretto*(evgamma_corr_Ruota - rfge.ev(theta_l))





####    Altre grandezze necessarie
addendum_pignone = round(raggio_testa_pignone_corretto - raggio_primitiva_pignone_corretto, 2)
addendum_ruota = round(raggio_testa_ruota_corretto - raggio_primitiva_ruota_corretto, 2)

dedendum_pignone = round(raggio_primitiva_pignone_corretto - raggio_piede_pignone_corretto, 2)
dedendum_ruota = round(raggio_primitiva_ruota_corretto - raggio_piede_ruota_corretto, 2)

# deltad = (abs(r_pignone-raggio_primitiva_pignone_corretto) + abs(r_ruota-raggio_primitiva_ruota_corretto))
# gioco_radiale = 2*deltad*np.sin(theta_l)

gioco_radiale = interasse - raggio_testa_pignone_corretto-raggio_piede_ruota_corretto

alpha_t = np.arctan(np.tan(theta_p)/np.cos(beta))                                           # angolo di pressione trasversale di riferimento UNI 8862 tab4

alpha_wt = np.arccos(((2*r_pignone)*np.cos(alpha_t))/(2*raggio_primitiva_pignone_corretto))                     # angolo di pressione trasversale di funzionamento UNI8862 tab4

rho_red = rfp.curvatura_relativa((2*raggio_piede_pignone_corretto), (2*raggio_piede_ruota_corretto), alpha_wt)   # curvatura relativa





####    Tabella riassuntiva
tab_grandezze = [ ["Coefficienti di correzione", x_pignone_approssimato, x_ruota_approssimato],
                    ["Angolo di lavoro [rad]", round(theta_l, 5), '='], 
                    ["Angolo di lavoro [grad]", round(np.degrees(theta_l), 2), '='],
                    ["Modulo [mm]", round(modulo_lavoro,3), '='],
                    ["Addendum [mm]", addendum_pignone, addendum_ruota],
                    ["Dedendum [mm]", dedendum_pignone, dedendum_ruota],
                    ["Raggio primitiva di lavoro [mm]", round(raggio_primitiva_pignone_corretto, 2), round(raggio_primitiva_ruota_corretto, 2)],
                    ["Raggio di testa [mm]", round(raggio_testa_pignone_corretto,2), round(raggio_testa_ruota_corretto,2)],
                    ["Raggio di piede [mm]", round(raggio_piede_pignone_corretto,2), round(raggio_piede_ruota_corretto,2)],
                    ["Raggio Fondamentale [mm]", round(rho_pignone,2), round(rho_ruota,2)],
                    ["Distanze tra primitive di lavoro e di taglio", round(y, 3), round(y_prime, 3)], 
                    ["Larghezza di fascia [mm]", 25, 25],
                    ["Spessore su primitiva di taglio [mm]", round(s_taglio_pignone,2), round(s_taglio_ruota,2)],
                    ["Spessore su primitiva di lavoro [mm]", round(spessore_primitiva_lavoro,2), round(spessore_primitiva_lavoro,2)],
                    ["Interasse [mm]", interasse, "="],
                    ["Gioco radiale [mm]", round(gioco_radiale, 2), "="],
                    ["Angolo di pressione trasversale di riferimento [rad]", round(alpha_t, 5), round(alpha_t, 5)],
                    ["Angolo di pressione trasversale di funzionamento [rad]",round(alpha_wt, 5), round(alpha_wt, 5) ],
                    ["Curvatura relativa", round(rho_red, 2), round(rho_red, 2)], 
                    ["Accesso", round(accesso_reale_new,2), "\\"],
                    ["Recesso", "\\", round(recesso_reale_new,2)],]

st11 = "A queste ruote si associano le seguenti grandezze caratteristiche CORRETTE"
st12 = tabulate(tab_grandezze, headers=["Pignone", "Ruota"])
out_file.write(st11 + '\n'), print(st11)
out_file.write(st12 + '\n'), print(st12)





####    Forze agenti
raggio_pignone_forza = raggio_primitiva_pignone_corretto*0.001                      # m
raggio_ruota_forza = raggio_primitiva_ruota_corretto*0.001                        # m

F_tangenziale_pignone = Momento_pignone/raggio_pignone_forza
F_tangenziale_ruota = Momento_ruota/raggio_ruota_forza

F_assiale_pignone = F_tangenziale_pignone*(np.tan(beta)/np.cos(theta_l))
F_assiale_ruota = F_tangenziale_ruota*(np.tan(beta)/np.cos(theta_l))

F_radiale_pignone = F_tangenziale_pignone*np.tan(theta_l)
F_radiale_ruota = F_tangenziale_ruota*np.tan(theta_l)

tab_forze = [["Forze circonferenziali [N]", round(F_tangenziale_pignone, 2), round(F_tangenziale_ruota, 2)],
             ["Forze radiali [N]", round(F_radiale_pignone, 2), round(F_radiale_ruota, 2) ], 
             ["Forze assiali [N]", round(F_assiale_pignone,2), round(F_assiale_ruota,2)]]

st13 = "Le forze agenti su queste ruote sono:"
st14 = tabulate(tab_forze, headers=["Pignone", "Ruota"])
out_file.write(st13 + '\n'), print(st13)
out_file.write(st14 + '\n'), print(st14)





################################################            DATABASE        ################################################
out_file.close()                                                        # chiudo e salvo il file log

post_correzione = {
    'variabili':["Coefficienti di correzione",
                 "Angolo di lavoro [rad]",
                 "Angolo di lavoro [grad]",
                 "Modulo [mm]",
                 "Addendum [mm]",
                 "Dedendum [mm]",
                 "Raggio primitiva di lavoro [mm]",
                 "Raggio di testa [mm]",
                 "Raggio di piede [mm]",
                 "Raggio Fondamentale [mm]",
                 "Distanze tra primitive di lavoro e di taglio",
                 "Larghezza di fascia [mm]",
                 "Spessore su primitiva di taglio [mm]",
                 "Spessore su primitiva di lavoro [mm]",
                 "Interasse [mm]",
                 "Gioco radiale [mm]", 
                 "Fattore di ricoprimento",
                 "Forze circonferenziali [N]",
                 "Forze radiali [N]",
                 "Angolo di pressione trasversale di riferimento [rad]",
                 "Angolo di pressione trasversale di funzionamento [rad]",
                 "Angolo in testa [rad]", 
                 "Curvatura relativa [mm]", 
                 "Accesso",
                 "Recesso",                             
                ],
    'pignone':[x_pignone_approssimato,
               round(theta_l, 5),
               round(np.degrees(theta_l), 2),           
               round(modulo_lavoro,3),
               addendum_pignone,
               dedendum_pignone,
               round(raggio_primitiva_pignone_corretto, 2),
               round(raggio_testa_pignone_corretto,2),
               round(raggio_piede_pignone_corretto,2),
               round(rho_pignone,2),
               round(y, 3),
               25,
               round(s_taglio_pignone,2),
               round(spessore_primitiva_lavoro,2),
               interasse,
               round(gioco_radiale, 2),
               round(epsilon_verifica_new, 1),
               round(F_tangenziale_pignone, 2),
               round(F_radiale_pignone, 2),
               round(alpha_t, 5),
               round(alpha_wt, 5),
               round(theta_t_corr_pignone,5), 
               round(rho_red, 2), 
               round(accesso_reale_new,2),
               0,
               ],
    'ruota':[x_ruota_approssimato,
             round(theta_l,5),
             round(np.degrees(theta_l), 2),            
             round(modulo_lavoro,3),
             addendum_ruota,
             dedendum_ruota,
             round(raggio_primitiva_ruota_corretto, 2),
             round(raggio_testa_ruota_corretto,2),
             round(raggio_piede_ruota_corretto,2),
             round(rho_ruota,2),
             round(y_prime, 3),
             25,
             round(s_taglio_ruota,2),
             round(spessore_primitiva_lavoro,2),
             interasse,
             round(gioco_radiale, 2),
             round(epsilon_verifica_new, 1),
             round(F_tangenziale_ruota, 2),
             round(F_radiale_ruota, 2),
             round(alpha_t, 5),
             round(alpha_wt, 5),
             round(theta_t_corr_ruota,5),
             round(rho_red, 2), 
             0,
             round(recesso_reale_new,2)],             
}
index_labels1 = ['r1', 'r2', 'r3', 'r4', 'r5', 
                 'r6', 'r7', 'r8', 'r9', 'r10', 
                 'r11', 'r12', 'r13', 'r14', 'r15', 
                 'r16', 'r17', 'r18', 'r19', 'r20',
                 'r21', 'r22', 'r23', 'r24', 'r25']

PC = pd.DataFrame(post_correzione, index=index_labels1)

# Connessione al database SQLite
conn = sqlite3.connect("post_correzione.db")

# Salva il DataFrame nel database
PC.to_sql('post_correzione', conn, if_exists='replace', index_label='index')

# Chiudi la connessione
conn.close()

# Tabella su excel sempre facilmente consultabile
PC.to_excel('correzione.xlsx', index=True)









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