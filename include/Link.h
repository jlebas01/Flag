/*=============================================================*/
/*= E.Incerti (incerti@univ-eiffel.fr)                        =*/
/*= M2/IMAC3 - ANIMATION / MOTEUR PHYSIQUE                    =*/
/*= Système Masses/Ressort en dimension 1 : module liaison    =*/
/*=============================================================*/

#ifndef  _LINK_
#define _LINK_

#include <PMat.h>

/* identificateurs de type                                                 */

/* Module produisant une force constante -- 'équivalent' du point fixe     */
#define _FRC_CONST     0
/* Ressort linéaire (de Hook) - force proportionnelle à l'élongation       */
#define _RESSORT       1
/* Frein linéaire - force proportionnelle à la vitesse relative des points */
#define _FREIN         2
/* Ressort+Frein intégrés                                                  */
#define _RESSORT_FREIN 3
/* Liaison condtionnelle de rupture -- condition sur l'élong. ou la force  */
#define _CONDIT_RF     4
/* Butée conditionelle : ressort frein inactif au-delà d'un seuil de dist. */
#define _RF_BUTEE      5

#define _FRC_EXT      6

typedef struct _link {
    int type;      /* type pour usages divers         */
    /*------------------------------------------------------*/
    double k;      /* raideur pour les ressorts          */
    double l;      /* longueur à vide pour les ressorts  */
    double z;      /* viscosité pour les freins          */
    double s;      /* seuil pour les liens conditionnels */
    double frcx;    /* force à distribuer sur M1 et M2    */
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
    bool on_off;                /* flag d'activité     */
    void (*setup)(struct _link *); /* l'algo de calcul    */
    /*------------------------------------------------------*/
    void (*draw)(struct _link *);  /* fonction de dessin  */
}
        Link;

/* les  constructeurs */
/*! Création d'un module Force Constante (exp. Gravité)                                  !*/
void FrcConst(Link *L, double frc, const char *axe);

/*! Création d'un ressort linéaire (de Hook)                                             !*/
void Ressort(Link *L, double k);

/*! Création d'un frein cinétique linéaire                                               !*/
void Frein(Link *L, double z);

/*! Création d'un ressort+frein : les 2 précédents combinés                              !*/
void RessortFrein(Link *L, double k, double z);

/*! Création d'une butée visco-élastique seuillée, non linéaire                          !*/
void RF_Butee(Link *L, double k, double z, double s);

/*! Création d'une liaison de rupture avec condition sur l'élongation                    !*/
void RF_CondPos(Link *L, double k, double z, double s);

/*! Connexion d'une Liaison entre 2 points Mat                                           !*/
void Connect(PMat *M1, Link *L, PMat *M2);

void FrcExt(Link *L, double ax, double ay, double az, double fx, double fy, double fz);
#endif
