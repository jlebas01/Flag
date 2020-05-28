/*=============================================================*/
/*= E.Incerti (incerti@univ-eiffel.fr)                        =*/
/*= M2/IMAC3 - ANIMATION / MOTEUR PHYSIQUE                    =*/
/*= Système Masses/Ressort en dimension 1 : module matériel   =*/
/*=============================================================*/

#ifndef  _PMAT_
#define _PMAT_

#include <g3x.h>

#define _POINTFIXE 0
#define _PARTICULE 1

typedef struct _pmat {
    int type;     /* type pour usages divers            */
    /*-----------------------------------------------------*/
    double m;        /* paramètre de masse                 */
    double pos;      /* position courante                  */
    double vitx;      /* vitesse  courante                  */
    double vity;      /* vitesse  courante                  */
    double vitz;      /* vitesse  courante                  */
    double frcx;      /* buffer de force                    */
    double frcy;      /* buffer de force                    */
    double frcz;      /* buffer de force                    */
    void (*setup)(struct _pmat *, double h); /* intégrateur */
    /*-----------------------------------------------------*/
    G3Xcolor col;   /* couleur RGBA (4 float)              */
    double x;     /* 2è coord. pour positionnement 2D    */
    double z;
    double lastX;
    double lastY;
    double lastZ;
    double coeff;
    void (*draw)(struct _pmat *);     /* fonction de dessin */
} PMat;

/*! Création d'une particule mobile !*/
/* avec l'intégrateur LeapFrog */
void MassLF(PMat *M, double pos, double x, double z, double m, double coeff);

/* variante, avec l'intégrateur Euler Explicite */
void MassEE(PMat *M, double pos, double x, double z, double m, double coeff);

/*! Création d'une masse fixe !*/
void Fixe(PMat *M, double pos, double x, double z);

#endif
