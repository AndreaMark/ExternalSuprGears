# SpurExternalGears
Questa cartella contiene tutto il necessario a dimensionare e verificare delle ruote dentate cilindriche a denti dritti, esterni. 

All'interno della cartella "moduli" sono contenuti tutte le funzioni necessarie al calcolo dei fattori geometrici e di resistenza. 
  - ruote_fattorigeometrici
    Contiene tutte le funzioni necessarie al primo dimensionamento della ruota dentata: numero di denti, grandezze caratteristiche, strisciamenti...
  - ruote_correzione
    Contiene una funzione iterativa necessario al calcolo dei coefficienti di correzione e delle caratteristiche corrette della ruota dentata.
  - ruote_fattoriinfluenza
    Contiene tutte quelle funzioni necessarie al calolo dei fattori di influenza secondo norma ISO 6336-1
  - ruote_fattoriflessione
    Contiene la maggiore parte delle funzioni necessarie al calcolo dei coefficienti correttivi di flessione secondo norma ISO 6336-3
  - ruote_fattoripitting
    Contiene la maggior parte delle funzione necessarie al calcolo dei coefficienti correttivi di pressione da contatto secondo norma ISO 6336-2
  - ruote_fattorigrippaggio
    Contiene la maggior parte delle funzioni necessario al calcolo dei coefficiente e delle temperature di esercizio all'interno di un riduttore secondo norma ISO 6336-20 (Vullo-Gears-Vol2)
  - ruote_interpolazione
    Contiene una funzione utile al calcolo dell'iterpolazione lineare di c oefficienti in tabella.

L'ordine in cui runnare gli script è chiaramente il seguente: 
1 - 2 - 3.1 - 3.2 - 4

Ogni script creerà un database dal quale lo script precedente attingerà valori, una tabella al fine di visualizzare i valori, un file "log" dove registererà i print e delle immagini. 
Ogni file viene salvato nella sua relativa cartella tranne i database, che per semplicità esecutiva vengono salvati all'interno della stessa cartella dello script. 

In questo modo è necessario modificare all'interno dello script solo i valori riguardanti la scelta del materiale all'interno del file "1", inserire alcuni fattori di input, e poi compilare i successivi script senza dover reinserire i dati per ogni script. 

(Il file "4" produce solo un log e delle immagini)


Andrea Marchegiani 
18/10/2023
