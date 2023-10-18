import numpy as np
import matplotlib.pyplot as plt
from moduli import ruote_fattorigeometrici as rfge
from moduli import ruote_fattoripitting as rfp
from moduli import ruote_correzione as rc
import pandas as pd
from tabulate import tabulate
import os
import shutil
from datetime import datetime

################################################            CORPO            ################################################  
start_time = datetime.now()

out_file = open("Correzione.txt", "w")                                # scrivo il file txt

####    Importazione
VG = pd.read_csv('variabili_globali.csv', index_col=0)
CNC = pd.read_csv('caratteristiche_nc.csv', index_col=0)
VC = pd.read_csv('variabili_cinematiche.csv', index_col=0)





####    Inizializzazione variabili globali
theta_p = (VG.loc['r1','valori'])                                       # Angolo di pressione
DentiP = (VG.loc['r2','valori'])
DentiR = (VG.loc['r3','valori'])                                        # Numero di denti di pignone e ruota
tau_re = DentiP/DentiR                                                  # Reale rapporto di trasmissione
m = (VG.loc['r5','valori'])                                             # modulo e addendum
a = m
beta = 0                                                                # angolo d'elica





####    Parametri iniziali
####    Raggi primitivi delle ruote NON corrette
r_pignone = (CNC.loc['r1','pignone'])
r_ruota = (CNC.loc['r1','ruota'])

####    Raggi fondamentali, questi raggi non subiranno correzione
rho_pignone = (CNC.loc['r2','pignone'])
rho_ruota = (CNC.loc['r2','ruota'])

####    Raggi di testa di ruota e pignone NON CORRETTI
Rtp = r_pignone + a
Rtr = r_ruota + a

####    Calcolo del segmento dei contatti NON CORRETTO
phi = np.pi/2 - theta_p
accesso_reale = np.sqrt(Rtr**2 - rho_ruota**2) - r_ruota*np.cos(phi)
recesso_reale = np.sqrt(Rtp**2 - rho_pignone**2) - r_pignone*np.cos(phi)
delta = np.arange(-accesso_reale, recesso_reale, 0.01)





####    Correzione dentatura
correzione = rc.correzione_dentature(tau_re, delta, theta_p, 0.31, 0.20, rho_pignone, rho_ruota, r_pignone, r_ruota, DentiP, DentiR, m)





####    SEGMENTO DEI CONTATTI CORRETTO E PLOT
delta_new = np.arange(-correzione[3], correzione[4], 0.01)

ks_pignone_plot2 = [rfge.ks_pignone(tau_re, dd, r_pignone, correzione[2]) for dd in delta_new]
ks_ruota_plot2 = [rfge.ks_ruota(tau_re, dd, r_ruota, correzione[2]) for dd in delta_new]

x_pignone_approssimato = round(correzione[0], 3)
x_ruota_approssimato = round(correzione[1], 3)
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
evgamma_corr_Pignone = rfge.ev(theta_p) + (2/DentiP)*(np.pi/4 + x_pignone_approssimato*np.tan(theta_p)) 
evgamma_corr_Ruota = rfge.ev(theta_p) + (2/DentiR)*(np.pi/4 + x_ruota_approssimato*np.tan(theta_p))

####    Angolo di pressione in testa
theta_t_corr_pignone = np.arccos(rho_pignone/(correzione[5]))
theta_t_corr_ruota = np.arccos(rho_ruota/(correzione[6]))

spessore_testa_corretto_pignone = 2*(correzione[5])*(evgamma_corr_Pignone - rfge.ev(theta_t_corr_pignone))
spessore_testa_corretto_ruota = 2*(correzione[6])*(evgamma_corr_Ruota - rfge.ev(theta_t_corr_ruota))

spessore_lim = 0.2*correzione[7]

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
phi_new = np.pi/2 - correzione[2]

t1c = r_pignone*np.cos(phi_new)     # Accesso Teorico

ct2 = r_ruota*np.cos(phi_new)       # Recesso Teorico

accesso_reale_new = correzione[3]

recesso_reale_new = correzione[4]

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
contatti_new = accesso_reale_new + recesso_reale_new

passo_new = correzione[7]*np.pi

epsilon_verifica_new = contatti_new/(passo_new*np.cos(correzione[2]))

st9 = f"Il fattore di ricoprimento è {round(epsilon_verifica_new, 1)}"
st10 = "Il fattore di ricoprimento non è abbastanza!"

if round(epsilon_verifica_new, 1)>=1.2:
    out_file.write(st9 + '\n'), print(st9)
else:
    out_file.write(st10 + '\n'), exit(st10)





####    Grandezze caratteristiche
####    Spessori primitive di taglio
s_taglio_pignone = 2*r_pignone*(evgamma_corr_Pignone - rfge.ev(theta_p))
s_taglio_ruota = 2*r_ruota*(evgamma_corr_Ruota - rfge.ev(theta_p))

####    Spessori primitive di lavoro
s_lavoro_pignone = 2*correzione[8]*(evgamma_corr_Pignone - rfge.ev(correzione[2]))
s_lavoro_ruota = 2*correzione[9]*(evgamma_corr_Ruota - rfge.ev(correzione[2]))





####    Altre grandezze necessarie
addendum_pignone = round(correzione[5]-correzione[8], 2)
addendum_ruota = round(correzione[6]-correzione[9], 2)
dedendum_pignone = round(correzione[8]-correzione[10], 2)
dedendum_ruota = round(correzione[9]-correzione[11], 2)
interasse = round(correzione[8] + correzione[9], 2)
deltad = (abs(r_pignone-correzione[8]) + abs(r_ruota-correzione[9]))
gioco_radiale = 2*deltad*np.sin(correzione[2])

alpha_t = np.arctan(np.tan(theta_p)/np.cos(beta))                                           # angolo di pressione trasversale di riferimento UNI 8862 tab4

alpha_wt = np.arccos(((2*r_pignone)*np.cos(alpha_t))/(2*correzione[8]))                     # angolo di pressione trasversale di funzionamento UNI8862 tab4

rho_red = rfp.curvatura_relativa((2*correzione[10]), (2*round(correzione[11])), alpha_wt)   # curvatura relativa





####    Tabella riassuntiva
tab_grandezze = [ ["Coefficienti di correzione", x_pignone_approssimato, x_ruota_approssimato],
                    ["Angolo di lavoro [rad]", round(correzione[2], 5), round(correzione[2],5)], 
                    ["Angolo di lavoro [grad]", round(np.degrees(correzione[2]), 2), round(np.degrees(correzione[2]), 2)],
                    ["Modulo [mm]", round(correzione[7],3), round(correzione[7],3)],
                    ["Addendum [mm]", addendum_pignone, addendum_ruota],
                    ["Dedendum [mm]", dedendum_pignone, dedendum_ruota],
                    ["Raggio primitiva di lavoro [mm]", round(correzione[8], 2), round(correzione[9], 2)],
                    ["Raggio di testa [mm]", round(correzione[5],2), round(correzione[6],2)],
                    ["Raggio di piede [mm]", round(correzione[10],2), round(correzione[11],2)],
                    ["Raggio Fondamentale [mm]", round(rho_pignone,2), round(rho_ruota,2)],
                    ["Distanze tra primitive di lavoro e di taglio", round(correzione[12], 3), round(correzione[13], 3)], 
                    ["Larghezza di fascia [mm]", 25, 25],
                    ["Spessore su primitiva di taglio [mm]", round(s_taglio_pignone,2), round(s_taglio_ruota,2)],
                    ["Spessore su primitiva di lavoro [mm]", round(s_lavoro_pignone,2), round(s_lavoro_ruota,2)],
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
Momento_pignone = (VC.loc['r2','pignone'])                      # Nm
Momento_ruota = (VC.loc['r2','ruota'])                          # Nm
raggio_pignone_forza = correzione[8]*0.001                      # m
raggio_ruota_forza = correzione[9]*0.001                        # m

F_tangenziale_pignone = Momento_pignone/raggio_pignone_forza
F_tangenziale_ruota = Momento_ruota/raggio_ruota_forza

F_assiale_pignone = F_tangenziale_pignone*(np.tan(0)/np.cos(theta_t_corr_pignone))
F_assiale_ruota = F_tangenziale_ruota*(np.tan(0)/np.cos(theta_t_corr_ruota))

F_radiale_pignone = F_tangenziale_pignone*np.tan(correzione[2])
F_radiale_ruota = F_tangenziale_ruota*np.tan(correzione[2])

tab_forze = [["Forze circonferenziali [N]", round(F_tangenziale_pignone, 2), round(F_tangenziale_ruota, 2)],
             ["Forze radiali [N]", round(F_radiale_pignone, 2), round(F_radiale_ruota, 2) ], 
             ["Forze assiali [N]", round(F_assiale_pignone,2), round(F_assiale_ruota,2)]]

st13 = "Le forze agenti su queste ruote sono:"
st14 = tabulate(tab_forze, headers=["Pignone", "Ruota"])
out_file.write(st13 + '\n'), print(st13)
out_file.write(st14 + '\n'), print(st14)





################################################            DATABASE        ################################################
out_file.close()                                                        # chiudo e salvo il file .txt

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
               round(correzione[2], 5),
               round(np.degrees(correzione[2]), 2),           
               round(correzione[7],3),
               addendum_pignone,
               dedendum_pignone,
               round(correzione[8], 2),
               round(correzione[5],2),
               round(correzione[10],2),
               round(rho_pignone,2),
               round(correzione[12], 3),
               25,
               round(s_taglio_pignone,2),
               round(s_lavoro_pignone,2),
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
             round(correzione[2],5),
             round(np.degrees(correzione[2]), 2),            
             round(correzione[7],3),
             addendum_ruota,
             dedendum_ruota,
             round(correzione[9], 2),
             round(correzione[6],2),
             round(correzione[11],2),
             round(rho_ruota,2),
             round(correzione[13], 3),
             25,
             round(s_taglio_ruota,2),
             round(s_lavoro_ruota,2),
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

PC.to_csv('correzione.csv', index=True)
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