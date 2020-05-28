/*=============================================================*/
/*= E.Incerti (incerti@univ-eiffel.fr)                        =*/
/*= M2/IMAC3 - ANIMATION / MOTEUR PHYSIQUE                    =*/
/*= Système Masses/Ressort en dimension 1 : module matériel   =*/
/*=============================================================*/

#include <PMat.h>


/*! Algorithme de la masse libre 2D (2° loi de Newton)    !*/
/*! intégration : leapfrog                                !*/
static void leapfrog(PMat *M, double h) {
    M->vitx += M->coeff * (h * M->frcx / M->m);   /* 1° intégration : vitesse   V(n+1)=V(n)+h*F(n)/m - EXplicite */
    M->vity += M->coeff *(h * M->frcy / M->m);   /* 1° intégration : vitesse   V(n+1)=V(n)+h*F(n)/m - EXplicite */
    M->vitz += M->coeff *(h * M->frcz / M->m);   /* 1° intégration : vitesse   V(n+1)=V(n)+h*F(n)/m - EXplicite */
    M->lastX = M->x;
    M->lastY = M->pos;
    M->lastZ = M->z;
    M->x += h * M->vitx;        /* 2° intégration : position  X(n+1)=X(n)+h*V(n+1) - IMplicite */
    M->pos += h * M->vity;        /* 2° intégration : position  X(n+1)=X(n)+h*V(n+1) - IMplicite */
    M->z += h * M->vitz;        /* 2° intégration : position  X(n+1)=X(n)+h*V(n+1) - IMplicite */
    M->frcx = 0.;              /* on vide le buffer de force */
    M->frcy = 0.;
    /*M->frcz = 0.;*/
    M->frcz = 0.0;
}

/*! intégration : Euler Explicite                         !*/
static void eulerexp(PMat *M, double h) {
    M->lastX = M->x;
    M->lastY = M->pos;
    M->lastZ = M->z;
    M->x += h * M->vitx;        /* 1° intégration : position  X(n+1)=X(n)+h*V(n)   - EXplicite */
    M->pos += h * M->vity;        /* 1° intégration : position  X(n+1)=X(n)+h*V(n)   - EXplicite */
    M->z += h * M->vitz;        /* 1° intégration : position  X(n+1)=X(n)+h*V(n)   - EXplicite */
    M->vitx += h * M->frcx / M->m;   /* 2° intégration : vitesse   V(n+1)=V(n)+h*F(n)/m - EXplicite */
    M->vity += h * M->frcy / M->m;   /* 2° intégration : vitesse   V(n+1)=V(n)+h*F(n)/m - EXplicite */
    M->vitz += h * M->frcz / M->m;   /* 2° intégration : vitesse   V(n+1)=V(n)+h*F(n)/m - EXplicite */
    M->frcx = 0.;              /* on vide le buffer de force */
    M->frcy = 0.;
    M->frcz = 0.;
}

/*! Algorithme du point fixe (position constante)         !*/
static void pointfixe(PMat *M, double h) { /* ne fait rien, à part vider le buffer de force       */
    M->frcx = 0.;
    M->frcy = 0.;
    M->frcz = 0.;
    /* ne sert que si on veut "figer" une particule mobile */
    M->vitx = 0.;
    M->vity = 0.;
    M->vitz = 0.;
    M->lastX = 0.;
    M->lastY = 0.;
    M->lastZ = 0.;
    M->coeff = 1.0;
}


static void drawcirc(PMat *M) { /* rayon indexé sur la masse */
    g3x_Material(M->col, 0.25, 0.7, 0.5, 0.5, 0.5f);
    if (M->type == _POINTFIXE) {
        glPushMatrix();
        glTranslatef(M->x, M->pos, M->z);
        glutSolidSphere(0.03 * M->m, 40, 40);
        glPopMatrix();
    } else {
        glPushMatrix();
        glTranslatef(M->x, M->pos, M->z);
        glutSolidSphere(0.03 * M->m, 40, 40);
        glPopMatrix();
    }
}

static void drawdot(PMat *M) { /* rayon indexé sur la masse */
    g3x_Material(M->col, 0.25, 0.7, 0.5, 0.5, 1.);
    if (M->type == _POINTFIXE) {
        glPushMatrix();
        glTranslatef(M->x, M->pos, M->z);
        glutSolidSphere(0.03 * M->m, 40, 40);
        glPopMatrix();
    } else {
        glPushMatrix();
        glTranslatef(M->x, M->pos, M->z);
        glutSolidSphere(0.03 * M->m, 40, 40);
        glPopMatrix();
    }
}


/*! Création d'une masse libre !*/
/*  Création d'une particule libre : attribution des paramètres de position et masse (vitesse nulle)  */
/*  avec l'intégrateur LeapFrog */
extern void MassLF(PMat *M, double pos, double x, double z, double m, double coeff) {
    M->type = _PARTICULE;
    /* paramètres pour le moteur physique */
    M->pos = pos;
    M->x = x;
    M->z = z;
    M->vitx = 0.;
    M->vity = 0.;
    M->vitz = 0.;
    M->frcx = 0.;
    M->frcy = 0.;
    M->frcz = 0.;
    M->coeff = coeff;
    M->m = m;
    M->setup = &leapfrog;
    /* paramètres graphiques */
    M->lastX = 0.0;
    M->lastY = 0.0;
    M->lastZ = 0.0;
    memcpy(M->col, G3Xb, sizeof(G3Xcolor));
    M->draw = &drawdot;
}

/* variante, avec l'intégrateur Euler Explicite */
extern void MassEE(PMat *M, double pos, double x, double z, double m, double coeff) {
    M->type = _PARTICULE;
    /* paramètres pour le moteur physique */
    M->pos = pos;
    M->x = x;
    M->z = z;
    M->vitx = 0.;
    M->vity = 0.;
    M->vitz = 0.;
    M->frcx = 0.;
    M->frcy = 0.;
    M->frcz = 0.;
    M->coeff = coeff;
    M->m = m;
    M->setup = &eulerexp;
    /* paramètres graphiques */
    M->lastX = 0.0;
    M->lastY = 0.0;
    M->lastZ = 0.0;
    memcpy(M->col, G3Xg, sizeof(G3Xcolor));
    M->draw = &drawcirc;
}

/*! Création d'une masse fixe !*/
extern void Fixe(PMat *M, double pos, double x, double z) {
    M->type = _POINTFIXE;
    /* paramètres pour le moteur physique */
    M->pos = pos;
    M->x = x;
    M->z = z;
    M->vitx = 0.;
    M->vity = 0.;
    M->vitz = 0.;
    M->frcx = 0.;
    M->frcy = 0.;
    M->frcz = 0.;
    M->coeff = 1.0;
    M->m = 1.; /* juste pour le dessin */
    M->setup = &pointfixe;
    /* paramètres graphiques */
    M->lastX = 0.0;
    M->lastY = 0.0;
    M->lastZ = 0.0;
    memcpy(M->col, G3Xr, sizeof(G3Xcolor));
    M->draw = &drawcirc;
}

