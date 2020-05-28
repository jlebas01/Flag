/*=============================================================*/
/*= E.Incerti (incerti@univ-eiffel.fr)                        =*/
/*= M2/IMAC3 - ANIMATION / MOTEUR PHYSIQUE                    =*/
/*= Syst�me Masses/Ressort en dimension 1 : module liaison    =*/
/*=============================================================*/

#ifndef  _LINK_
#define _LINK_

#include <PMat.h>

/* identificateurs de type                                                 */

/* Module produisant une force constante -- '�quivalent' du point fixe     */
#define _FRC_CONST     0
/* Ressort lin�aire (de Hook) - force proportionnelle � l'�longation       */
#define _RESSORT       1
/* Frein lin�aire - force proportionnelle � la vitesse relative des points */
#define _FREIN         2
/* Ressort+Frein int�gr�s                                                  */
#define _RESSORT_FREIN 3
/* Liaison condtionnelle de rupture -- condition sur l'�long. ou la force  */
#define _CONDIT_RF     4
/* But�e conditionelle : ressort frein inactif au-del� d'un seuil de dist. */
#define _RF_BUTEE      5

#define _FRC_EXT      6

typedef struct _link {
    int type;      /* type pour usages divers         */
    /*------------------------------------------------------*/
    double k;      /* raideur pour les ressorts          */
    double l;      /* longueur � vide pour les ressorts  */
    double z;      /* viscosit� pour les freins          */
    double s;      /* seuil pour les liens conditionnels */
    double frcx;    /* force � distribuer sur M1 et M2    */
    double frcy;
    double frcz;
    double ax;
    double ay;
    double az;
    double fx;
    double fy;
    double fz;
    double t;
    PMat
            *M1, *
            M2; /* les 2 PMat de branchement          */
    bool on_off;                /* flag d'activit�     */
    void (*setup)(struct _link *); /* l'algo de calcul    */
    /*------------------------------------------------------*/
    void (*draw)(struct _link *);  /* fonction de dessin  */
}
        Link;

/* les  constructeurs */
/*! Cr�ation d'un module Force Constante (exp. Gravit�)                                  !*/
void FrcConst(Link *L, double frc, const char *axe);

/*! Cr�ation d'un ressort lin�aire (de Hook)                                             !*/
void Ressort(Link *L, double k);

/*! Cr�ation d'un frein cin�tique lin�aire                                               !*/
void Frein(Link *L, double z);

/*! Cr�ation d'un ressort+frein : les 2 pr�c�dents combin�s                              !*/
void RessortFrein(Link *L, double k, double z);

/*! Cr�ation d'une but�e visco-�lastique seuill�e, non lin�aire                          !*/
void RF_Butee(Link *L, double k, double z, double s);

/*! Cr�ation d'une liaison de rupture avec condition sur l'�longation                    !*/
void RF_CondPos(Link *L, double k, double z, double s);

/*! Connexion d'une Liaison entre 2 points Mat                                           !*/
void Connect(PMat *M1, Link *L, PMat *M2);

void FrcExt(Link *L, double ax, double ay, double az, double fx, double fy, double fz);
#endif
