# -*-coding:utf-8 -*

# IMPORT MODULES------------------------------------------------------------
import re
from Bio import SeqIO
from math import exp
import os.path
import pickle
import sys


# FONCTIONS-----------------------------------------------------------------
def verifParametres(args, parametres):
    """
    Verification des paramètres. Retourne le dictionaire de paramètres
    """

    # Verification du paramètre -r (Séquence requete fasta)
    extension = None
    if args.r is not None:
        if os.path.isfile(args.r):
            if (args.r.split('/')[-1]).split('.')[-1].lower() == "fasta":
                parametres["requete"] = args.r
            else:
                sys.exit("erreur: Séquence requete doit être au format fasta.")
        else:
            sys.exit("erreur: Fichier séquence requete inexistant.")
    else:
        sys.exit("erreur: Séquence requete doit être spécifié.")

    # Verifier le paramètre -b (Base de donnée au format fasta)
    extension = None
    if args.b is not None:
        if os.path.isfile(args.b):
            if (args.b.split('/')[-1]).split('.')[-1].lower() in ["fasta", "saved"]:
                parametres["bdd"] = args.b
            else:
                sys.exit("erreur: Base de données doit être au format fasta ou saved.")
        else:
            sys.exit("erreur: Fichier base de données inexistant.")
    elif os.path.isfile(parametres["bdd"]):
        pass
    else:
        sys.exit("erreur: Fichier base de données inexistant.")

    # Verifier le paramètre -m (Matrice de score)
    if args.m is not None:
        if os.path.isfile(args.m):
                parametres["matrice"] = args.m
        else:
            sys.exit("erreur: Fichier matrice de score inexistant.")
    elif os.path.isfile("./blosum62.txt"):
        pass
    else:
        sys.exit("erreur: Fichier matrice de score inexistant.")

    # Verifier le paramètre -v (Version de Blast)
    if args.v is not None:
        if args.v in ["1", "2"]:
            parametres["version"] = int(args.v)
        else:
            sys.exit("erreur: version incorrecte")

    # Verifier le paramètre -w (Taille des mots)
    if args.w is not None:
        if(args.w >= 1):
            parametres["w"] = int(args.w)
        else:
            sys.exit("erreur: Taille de mots doit être superieur à 0.")

    # Verifier le paramètre -s (Score minimal des HSP en sortie)
    if args.s is not None:
            parametres["seuilScore"] = int(args.s)

    # Verifier le paramètre -e (Seuil de e-value)
    if args.e is not None:
        if float(args.e) > 0.0:
            parametres["seuilEvalue"] = float(args.e)
        else:
            sys.exit("erreur: Seuil evalue doit être superieur à 0.")

    # Verifier le paramètre -n (Nombre de HSP en sortie)
    if args.n is not None:
        if int(args.n) > 0:
            parametres["nbHSP"] = int(args.n)
        else:
            sys.exit("erreur: Nombre de HSP en sortie doit être superieur à 0.")

    # Verifier le paramètre -o (Fichier de sortie)
    if args.o is not None:
        try:
            with open(args.o,"w"): pass
        except IOError:
            sys.exit("erreur: Le fichier de sortie n'pas pu être ouvert")
        parametres["out"] = args.o

    return parametres


def loadBDDSaved(file):
    """
    Charge une base de donnée sous forme d'un dictionaire binaire avec pickle.
    """

    with open(file, 'rb') as fInput:
        # Incompatible Python 2.7 si le fichier est écrit avec Python3.
        bdd = pickle.load(fInput)
    return bdd


def dumpBDDSaved(fileName, bdd):
    """
    Enregistre une base de donnée sous forme d'un dictionaire binaire
    à partir d'un fichier fasta avec pickle.
    """

    regex = re.compile("(.+)fasta$")
    res = regex.search(fileName.lower())
    with open(res.group(1)+'saved', 'wb') as fOutput:
        # Paramètre protocole=2, ne corrige pas le problème de compatibilité.
        pickle.dump(bdd, fOutput)
    return res.group(1)+'saved'


def mot_similaire(mot1, mat):
    """
    Construit la liste des mot similaire à 'mot1' selon la matrice de score donnée.
    """

    T = 13
    motSim = []
    i = 0
    j = 0
    k = 0
    for i in range(len(mat.indx_AA)-1):
        for j in range(len(mat.indx_AA)-1):
            for k in range(len(mat.indx_AA)-1):
                mot2 = mat.indx_AA[i] + mat.indx_AA[j] + mat.indx_AA[k]
                score = mat.score(mot1[0], mat.indx_AA[i]) + mat.score(mot1[1], mat.indx_AA[j]) + mat.score(mot1[2], mat.indx_AA[k])
                if(score >= T):
                    motSim.append(mot2)
    return motSim


def importSeqFasta(file):
    """
    Import une sequence fasta et retourne l'id et la séquence.
    """

    regex = re.compile("^.{2}\|([a-zA-Z0-9]+)\|")
    record = SeqIO.read(file, "fasta")
    res = regex.search(record.id)
    return (res.group(1), str(record.seq))


def dicoMots(seq, w, mat):
    """
    Construit la table de hachage contenant, pour chaque mot
    de taille w de la séquence seq, la ou les position(s) du mot
    dans la séquence.
    """

    dico = {}
    for i in range(len(seq) + 1 - w):
        if seq[i:i+w] not in dico.keys():
            tmp_mot = mot_similaire(seq[i:i+w], mat)
            for mot in tmp_mot:
                dico[mot] = [i]
        else:
            dico[seq[i:i+w]].append(i)
    return dico


def scoreExtension(seq, seqbdd, mat, posSeq=0, posSeqbdd=0, sens=1):
    """
    Calcul le score (à partir de la matrice de score mat) de l'allongement
    d'un hit à la position posSeq (séquence seq) et posSeqbdd (séquence seqbdd)
    vers la droite (sens = 1) ou la gauche (sens = -1).
    Retourne les position d'arret de l'alignement sur les 2 seq et le score.
    """

    s = [0]
    sum_score = 0
    seuil = 15
    debut = 0
    if(sens == -1):
        debut = -1
        l = -(min(posSeq, posSeqbdd)+1)
        if(posSeq == 0 or posSeqbdd == 0):
            return (posSeq, posSeqbdd, 0)
    elif(sens == 1):
        # tester que pos < len
        l = min((len(seq)-posSeq), (len(seqbdd)-posSeqbdd))
    else:
        return None
    # print(posSeq)
    # print(posSeqbdd)
    for i in range(debut, l, sens):
        score_AA = mat.score(seq[posSeq+i], seqbdd[posSeqbdd+i])
        sum_score += score_AA
        # print(score_AA)
        s.append(s[i]+score_AA)
        maxi = max(s)
        if((maxi - s[i+1]) >= seuil):
            break
    return (posSeq+i, posSeqbdd+i, sum_score)


# CLASSES-----------------------------------------------------------------
class Hit:
    """
    Information sur un hit associé à un HSP.
    ammorce : mot de w lettre.
    posMotReq : position de l'ammorce sur la séquence requête.
    infoMotBDD : tuple  - ID de la séquence de la BDD
                        - position de l'ammorce sur la séquence de la BDD
    """
    def __init__(self, mot, posMotReq, infoMotBDD):
        self.ammorce = mot
        self.posMotReq = posMotReq
        self.infoMotBDD = infoMotBDD


class MatScore:
    """
    Information sur la matrice de score.
    matrice : matrice des valeurs de scores.
    indx_AA : liste des AA -> indice des AA en ligne et colonne.
    """

    def __init__(self, FileName):
        self.matrice = []
        self.indx_AA = []
        with open(FileName, "r") as fileMatrice:
            self.indx_AA = fileMatrice.readline().split()
            for ligne in fileMatrice:
                self.matrice.append([int(x) for x in ligne.split()[1:]])

    def __repr__(self):
        out = "  "+"  ".join(self.indx_AA)+"\n"
        for i, ligne in enumerate(self.matrice):
            out += str(self.indx_AA[i])+"  "+"  ".join([str(x) for x in ligne])+"\n"
        return out

    def score(self, aa1, aa2):
        """
        Retourne le score de substitution de 2 acides aminés aa1 et aa2.
        """
        return self.matrice[self.indx_AA.index(aa1)][self.indx_AA.index(aa2)]


class HSP:
    """
    Information sur un HSP.
    r_end: position du début de l'HSP sur la séquence requête.
    q_end: position de fin de l'HSP sur la séquence requête.
    r_start: position du début de l'HSP sur la séquence de la BDD.
    q_start: position de fin de l'HSP sur la séquence de la BDD.
    score_tot: score total du HSP sur toute sa longueur.
    hit: Information du hit associé au HSP.
    identite: Pourcentage d'identité entre la séquence requête et la
              séquence de la base de donnée associé au HSP.
    similarite: Pourcentage de similarité (selon la matrice de score) entre la séquence requête et la
                séquence de la base de donnée associé au HSP.
    length: longueur du HSP.
    eval: E-value du HSP.
    """

    def __init__(self, motReq, posMotReq, infoMotBDD):
        self.r_end = -1
        self.q_end = -1
        self.r_start = -1
        self.q_start = -1
        self.score_tot = 0
        self.hit = Hit(motReq, posMotReq, infoMotBDD)
        self.identite = -1
        self.similarite = -1
        self.length = -1
        self.eval = -1

    def __repr__(self):
        # strHit = "Hit : '{}', position requête {}, id séquence bdd {}, position séquence bdd {}\n".format(self.hit.ammorce, self.hit.posMotReq, self.hit.infoMotBDD[0], self.hit.infoMotBDD[1])
        strHSP = "Postion séquence requête : {} - {}\nPosition séquence {}: {} - {}\nScore : {}\nId. : {:4.2f}%\nSim. : {:4.2f}%\nEvalue : {:4.2f}\n".format(self.r_start, self.r_end, self.hit.infoMotBDD[0], self.q_start, self.q_end, self.score_tot, self.identite, self.similarite, self.eval)
        return strHSP

    def extension(self, seq, seqbdd, mat):
        """ Extension du hit vers la droite puis vers la gauche.
            Determine le score total du HSP et sa longeur.
        """
        s = 0
        s2 = 0
        # Vers la droite
        (self.r_end, self.q_end, s) = scoreExtension(seq, seqbdd, mat, self.hit.posMotReq, self.hit.infoMotBDD[1], 1)
        # print("droite")
        # print(self.r_end,self.q_end, s)

        # Vers la gauche
        (self.r_start, self.q_start, s2) = scoreExtension(seq, seqbdd, mat, self.hit.posMotReq, self.hit.infoMotBDD[1], -1)
        # print("gauche")
        # print(self.r_start,self.q_start, s2)
        self.length = self.r_end - self.r_start
        self.score_tot = s+s2

    def pourcentageID_SIM(self, seq, seqbdd, mat):
        """
        Calcul le pourcentage d'identité et de similarité
        (substitution à score positif ou nul selon la matrice de score) ayant un score de
        entre les 2 séquence le long de l'HSP.
        """

        for i in range(0,  self.length + 1):
            if (seq[i] == seqbdd[i]):
                self.identite += 1
            if(mat.score(seq[i], seqbdd[i]) >= 0):
                self.similarite += 1
        self.identite = self.identite*100.0/self.length
        self.similarite = self.similarite*100.0/self.length

    def e_value(self, M):
        """
        Calcul approximatif de la e-value de l'HSP
        """

        # M = taille de la base de donnée
        K = 0.128
        L = 0.311
        self.eval = K*M*self.length*exp(-L*self.score_tot)


    def ecriture(self, seq, seqbdd):
        """
        Format l'écriture des informations du HSP dans un fichier.
        """

        strHSP = "Score : {}\nId. : {:4.2f}%\nSim. : {:4.2f}%\nEvalue : {}\n\n".format(self.score_tot, self.identite, self.similarite, self.eval)
        seqr = seq[1]
        i = 0
        for i in range(int(self.length/60)):
            strHSP += "{:6}{:5} {} {:5}\n".format(seq[0], self.r_start+i*60+1, seqr[:60], self.r_start+i*60+60)
            strHSP += "{:6}{:5} {} {:5}\n".format(self.hit.infoMotBDD[0], self.q_start+i*60+1, seqbdd[:60], self.q_start+i*60+60)
            seqr = seqr[60:]
            seqbdd = seqbdd[60:]

        strHSP += "{:6}{:5} {} {:5}\n".format(seq[0], self.r_start+int(self.length/60)*60, seqr[:self.length % 60], self.r_end)
        strHSP += "{:6}{:5} {} {:5}\n\n\n".format(self.hit.infoMotBDD[0], self.q_start+int(self.length/60)*60, seqbdd[:self.length % 60], self.q_end)

        return strHSP


class BDD:
    """
    Information sur la base de donnée.
    sizeBDD: taille de la base de donnée, nombre d'acides aminés.
    fileFasta: chemin vers le fichier fasta contenant les séquences de la base de donnée.
    dico_seqBDD: dictionnaire avec les id de séquences en clés et les sequences en valeurs.
    dico_3wBDD: dictionnaire de mot de W lettres (clés)
                associé à une liste de tuple : Id sequénce, poistion dans la séquence.
    """

    def __init__(self, fileFasta, W=3):
        self.sizeBDD = 0
        self.fileFasta = fileFasta
        self.dico_seqBDD = self.readFasta()
        self.dico_3wBDD = self.create3wBDD(W)

    def __repr__(self):
        out = ""
        for mot in self.dico_3wBDD:
            out += mot+"\n"
            out += str(self.dico_3wBDD[mot])+"\n\n"
        return out

    def readFasta(self):
        """
        Lecture du fichier fasta et création du dictionaire des séquences.
        """

        with open(self.fileFasta, "r") as fasta:
            dico_seq = {}
            regex = re.compile("^.{2}\|([a-zA-Z0-9]+)\|")
            for record in SeqIO.parse(fasta, "fasta"):
                res = regex.search(record.id)
                dico_seq[res.group(1)] = str(record.seq)
        return dico_seq

    def create3wBDD(self, W):
        """
        Création du dictionnaire de mot de W lettres
        et calcul de la taille de la base de données.
        """

        bdd = {}
        for ID in self.dico_seqBDD:
            for i in range(len(self.dico_seqBDD[ID])+1-W):
                if self.dico_seqBDD[ID][i:i+W] not in bdd.keys():
                    bdd[self.dico_seqBDD[ID][i:i+W]] = [(ID, i)]
                    self.sizeBDD += 1
                else:
                    bdd[self.dico_seqBDD[ID][i:i+W]].append((ID, i))
                    self.sizeBDD += 1
        self.sizeBDD *= W
        return bdd
