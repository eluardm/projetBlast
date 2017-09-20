#!/usr/bin/python3
# -*-coding:utf-8 -*

from mdBlast import *
from time import time, gmtime, strftime
from datetime import datetime
import sys
import argparse


tps1 = time()

############################################################################
#                                                                          #
#                           GESTION DES PARAMETRES                         #
#                                                                          #
############################################################################

# Dictionnaire de paramètres (Valeurs par défaut)
parametres = {
    "requete": None,
    "bdd": "../data/testBDD.fasta",
    "matrice": "../data/blosum62.txt",
    "version": 1,
    "w": 3,
    "seuilScore": 0,
    "seuilEvalue": 10,
    "nbHSP": 250,
    "out": "../resultats/results_blast.txt"}

'''
Gestion des paramètres:
8 paramètres :
-h: aide
-r: Chemin vers la séquence requête au format fasta.
-b: Chemin vers la BDD au format fasta ou saved
    (précédement construite par une premiere utilisattion avec un fichier fasta).
    (Par défaut utilisation de la table de hachage pré-construite fournie)
-m: Chemin vers la matrice de score. (par défaut, matrice blosum62 fournie)
-v: Version de Blast (1 ou 2). (2 non fonctionelle)
-w: Taille des mots. (Par défaut, 3)
-s: Score minimal des HSP en sortie. (Par défaut, 0)
-e: Seuil de e-value. (Par défaut, 10)
-n: Nombre de HSP en sortie. (Par défaut, 250)
-o: fichier de sortie.
'''

if(len(sys.argv) > 16):
    sys.exit("erreur: trop de parametres (voir aide -h)")

parser = argparse.ArgumentParser()
parser.add_argument("-r", help="Chemin vers la séquence requête au format fasta.")
parser.add_argument("-b", help="Chemin vers la BDD au format fasta. \n(précédement construite par une premiere utilisattion avec un fichier fasta)\n(Par défaut, utilisation de la table de hachage pré-construite fournie)")
parser.add_argument("-m", help="Chemin vers la matrice de score. (par défaut, matrice blosum62 fournie.")
parser.add_argument("-v", type=int, help="Version de Blast (1 ou 2). (2 non fonctionelle)")
parser.add_argument("-w", type=int, help="Taille des mots. (Par défaut, 3)")
parser.add_argument("-s", type=int, help="Score minimal des HSP en sortie. (Par défaut, 0)")
parser.add_argument("-e", type=float, help="Seuil de e-value. (Par défaut, 10)")
parser.add_argument("-n", type=int, help="Nombre de HSP en sortie. (Par défaut, 250)")
parser.add_argument("-o", help="Fichier de sortie")

args = parser.parse_args()
parametres = verifParametres(args, parametres)


############################################################################
#                                                                          #
#                            PROGRAMME PRINCIPAL                           #
#                                                                          #
############################################################################

# Importer la séquence requête----------------------------------------------------
print("Importation séquence requête... ")
seqReq = importSeqFasta(parametres["requete"])
# seqReq => (id, séquence)
print("OK\n")

# Importer la matrice de score----------------------------------------------------
print("Importation de la matrice de scores... ")
blosum62 = MatScore(parametres["matrice"])
print("OK\n")

# Création de la table de hachage de mots de W lettres pour la séquence requête---
if parametres["w"] > len(seqReq[1])/2:
    sys.exit("erreur: Taille des mots trop grande.")

print("Construction du dictionaire de mots similaires pour la séquence requête... ")
dicoWRequete = dicoMots(seqReq[1], parametres["w"], blosum62)
print("OK\n")

# Importer la base de données----------------------------------------------------
# Si la base est fournie au format fasta, elle est sauvegardée au format saved
# afin de pouvoir être importée plus rapidement les fois suivantes.
print("Importation de la base de donnée... ")
if(parametres["bdd"].split('.')[-1] == "fasta"):
    bdd = BDD(parametres["bdd"])
    print("OK\n")
    # Sauvegarde de la base de donnée sur disque au format saved.
    print("Sauvegarde de la base de donnée... ")
    res = dumpBDDSaved(parametres["bdd"], bdd)
    print('OK ('+res+')\n')
else:
    bdd = loadBDDSaved(parametres["bdd"])
    print("OK\n")


# Recherche des hits et extension----------------------------------------------------------
print("BLAST en cours... \n")
# Classe au fur et à mesure les HSP par ordre croissant
lScoreCroissant = []

dicoHSP = {}
for motReq in dicoWRequete:
    for posMotReq in dicoWRequete[motReq]:
        if motReq in bdd.dico_3wBDD.keys():
            for infoMotBDD in bdd.dico_3wBDD[motReq]:
                # infoMotBDD => (idSequenceBDD, position)
                # Si le hit est présent dans un HSP existant, alors il est ignoré.
                # Approximation pour gagner du temps de calcul.                
                if(not (infoMotBDD[0] in dicoHSP and (dicoHSP[infoMotBDD[0]].q_start <= infoMotBDD[1] <= dicoHSP[infoMotBDD[0]].q_end) and (dicoHSP[infoMotBDD[0]].r_start <= posMotReq <= dicoHSP[infoMotBDD[0]].r_end))):

# Extension des hits---------------------------------------------------------
                    hsp_tmp = HSP(motReq, posMotReq, infoMotBDD)
                    hsp_tmp.extension(seqReq[1], bdd.dico_seqBDD[infoMotBDD[0]], blosum62)
                    hsp_tmp.e_value(bdd.sizeBDD)
# Sélection des HSP à conserver----------------------------------------------
                    # Ignorer les HSP de score trop faible ou de e-value trop élevée.
                    if(hsp_tmp.score_tot < parametres["seuilScore"] or hsp_tmp.eval > parametres["seuilEvalue"]):
                        continue
                    # Ignorer les HSP si il existe déjà un HSP de score supérieur sur la même séquence.
                    if((hsp_tmp.hit.infoMotBDD[0] not in dicoHSP) or (hsp_tmp.hit.infoMotBDD[0] in dicoHSP and hsp_tmp.score_tot > dicoHSP[hsp_tmp.hit.infoMotBDD[0]].score_tot)):
                        hsp_tmp.pourcentageID_SIM(seqReq[1], bdd.dico_seqBDD[infoMotBDD[0]], blosum62)
                        dicoHSP[hsp_tmp.hit.infoMotBDD[0]] = hsp_tmp
                        # Classement des HSP par ordre croissant de score.
                        if(hsp_tmp.hit.infoMotBDD[0] in lScoreCroissant):
                            lScoreCroissant.remove(hsp_tmp.hit.infoMotBDD[0])
                        i = 0
                        for i, hsp in enumerate(lScoreCroissant):
                            if(dicoHSP[hsp].score_tot > hsp_tmp.score_tot):
                                break
                        lScoreCroissant.insert(i, hsp_tmp.hit.infoMotBDD[0])

# Ecriture sortie---------------------------------------------------
lScoreCroissant.reverse()
i = 0

with open(parametres["out"], "w") as fOut:
    while i < min(parametres["nbHSP"], len(lScoreCroissant)):
        fOut.write("### HSP "+str(i+1)+" ###\n")
        fOut.write(dicoHSP[lScoreCroissant[i]].ecriture(seqReq, bdd.dico_seqBDD[dicoHSP[lScoreCroissant[i]].hit.infoMotBDD[0]]))
        i += 1

print("Nombre de HSP: "+str(len(lScoreCroissant)))
print("\nDurée d'exécution : ")
print(strftime('%H:%M:%S', gmtime(time()-tps1)))
