/*=============================================================*/
/*= E.Incerti (incerti@univ-eiffel.fr)                        =*/
/*= M2/IMAC3 - ANIMATION / MOTEUR PHYSIQUE                    =*/
/*= Système Masses/Ressort en dimension 1 : module liaison    =*/
/*=============================================================*/

#include <g3x.h>

#include <PMat.h>
#include <Link.h>

/*================================================*/
/*=      ALGORITHMES DES MODULES DE LIAISON      =*/
/*================================================*/

/*! Algorithme Force Constante             !*/
/*  1 seul point mat. M1 ou M2, au choix... */
static void AlgoFrcConst(Link *L) {
    PMat *M = (L->M1 ? L->M1 : L->M2);
    M->frcx += L->frcx;
    M->frcy += L->frcy;
    M->frcz += L->frcz;
}

static void AlgoFrcExt(Link *L) {
    PMat *M = (L->M1 ? L->M1 : L->M2);

   /* int aMin = -200;
    int aMax = 200;
    int fMin = 100;
    int fMax = 400;

    double ax = (rand() % ((aMax + 1) - aMin) + aMin) / 100.0f;
    double ay  = (rand() % ((aMax + 1) - aMin) + aMin) / 100.0f;
    double az  = (rand() % ((aMax + 1) - aMin) + aMin) / 100.0f;

    double Fx  = (rand() % ((fMax + 1) - fMin) + fMin) / 100.0f;
    double Fy  = (rand() % ((fMax + 1) - fMin) + fMin) / 100.0f;
    double Fz  = (rand() % ((fMax + 1) - fMin) + fMin) / 100.0f;*/

    double fx = L->ax * cos(L->t * L->fx);
    double fy = L->ay * cos(L->t * L->fy);
    double fz = L->az * cos(L->t * L->fz);

    /*double fx = ax * cos(L->t * Fx);
    double fy = ay * cos(L->t * Fy);
    double fz = az * cos(L->t * Fz);*/

    L->t++;
    M->frcx += fx;
    M->frcy += fy;
    M->frcz += fz;
}

/*! Algo. ressort de Hook (linéaire)                            !*/
/* produit une force : F=k*(d-l0).AB                             */
/* où (d-l0) est l'allongement et AB le vecteur inter particules */
static void AlgoRessort(Link *L) {
    double d = sqrt(
            pow((L->M2->x - L->M1->x), 2.0) + pow((L->M2->pos - L->M1->pos), 2.0) + pow((L->M2->z - L->M1->z), 2.0));
    double M1M2x = (L->M2->x - L->M1->x) / d;
    double M1M2y = (L->M2->pos - L->M1->pos) / d;
    double M1M2z = (L->M2->z - L->M1->z) / d;
    double fx = L->k * (d - L->l) * M1M2x;
    double fy = L->k * (d - L->l) * M1M2y;
    double fz = L->k * (d - L->l) * M1M2z;
    L->M1->frcx += fx * L->M1->coeff;
    L->M1->frcy += fy * L->M1->coeff;
    L->M1->frcz += fz * L->M1->coeff;
    L->M2->frcx -= fx * L->M2->coeff;
    L->M2->frcy -= fy * L->M2->coeff;
    L->M2->frcz -= fz * L->M2->coeff;
}

/*! Algo. frein cinétique                           !*/
/* produit une force : F=z*(V2-V1)                   */
/* où (V2-V1) est la vitesse relative des particules */
static void AlgoFrein(Link *L) {
    double fx = L->z * (L->M2->vitx - L->M1->vitx);
    double fy = L->z * (L->M2->vity - L->M1->vity);
    double fz = L->z * (L->M2->vitz - L->M1->vitz);

    L->M1->frcx += fx;
    L->M1->frcy += fy;
    L->M1->frcz += fz;
    L->M2->frcx -= fx;
    L->M2->frcy -= fy;
    L->M2->frcz -= fz;
}

/*! Algo. ressort+frein         !*/
/*combine les 2 algos précédents */
static void AlgoRessortFrein(Link *L) {
    /*double d = L->M2->pos - L->M1->pos;*/
    double d = sqrt(
            pow((L->M2->x - L->M1->x), 2.0) + pow((L->M2->pos - L->M1->pos), 2.0) + pow((L->M2->z - L->M1->z), 2.0));
    double M1M2x = (L->M2->x - L->M1->x) / d;
    double M1M2y = (L->M2->pos - L->M1->pos) / d;
    double M1M2z = (L->M2->z - L->M1->z) / d;
    /*double f = L->k * (d - L->l) + L->z * (L->M2->vit - L->M1->vit);*/
    double fx = L->k * (d - L->l) * M1M2x + L->z * (L->M2->vitx - L->M1->vitx);
    double fy = L->k * (d - L->l) * M1M2y + L->z * (L->M2->vity - L->M1->vity);
    double fz = L->k * (d - L->l) * M1M2z + L->z * (L->M2->vitz - L->M1->vitz);

    L->M1->frcx += fx * L->M1->coeff;
    L->M1->frcy += fy * L->M1->coeff;
    L->M1->frcz += fz * L->M1->coeff;
    L->M2->frcx -= fx * L->M2->coeff;
    L->M2->frcy -= fy * L->M2->coeff;
    L->M2->frcz -= fz * L->M2->coeff;
}

/*! Algo. butée visco-élastique                    !*/
/*! active uniquement si dist. < seuil             !*/
/* comme précédemment mais avec un seuil d'activité */
static void AlgoRF_Butee(Link *L) {
    /*double d = L->M2->pos - L->M1->pos;*/
    double d = sqrt(
            pow((L->M2->x - L->M1->x), 2.0) + pow((L->M2->pos - L->M1->pos), 2.0) + pow((L->M2->z - L->M1->z), 2.0));
    if (d > L->s) return;

    double M1M2x = (L->M2->x - L->M1->x) / d;
    double M1M2y = (L->M2->pos - L->M1->pos) / d;
    double M1M2z = (L->M2->z - L->M1->z) / d;

    double fx = L->k * (d - L->l) * M1M2x + L->z * (L->M2->vitx - L->M1->vitx);
    double fy = L->k * (d - L->l) * M1M2y + L->z * (L->M2->vity - L->M1->vity);
    double fz = L->k * (d - L->l) * M1M2z + L->z * (L->M2->vitz - L->M1->vitz);

    L->M1->frcx += fx * L->M1->coeff;
    L->M1->frcy += fy * L->M1->coeff;
    L->M1->frcz += fz * L->M1->coeff;
    L->M2->frcx -= fx * L->M2->coeff;
    L->M2->frcy -= fy * L->M2->coeff;
    L->M2->frcz -= fz * L->M2->coeff;
}


/*! Algo. lien visco-élastique "cassable"          !*/
/*! si d>seuil, la liaison devient inactive        !*/
static void AlgoRF_CondPos(Link *L) { /* si la liaison est déjà cassée : rien */
    if (!L->on_off) return;

    /*double d = L->M2->pos - L->M1->pos;*/
    double d = sqrt(
            pow((L->M2->x - L->M1->x), 2.0) + pow((L->M2->pos - L->M1->pos), 2.0) + pow((L->M2->z - L->M1->z), 2.0));
    /* si l'allongement est trop fort : la liaison casse */
    if (d > L->s) {
        L->on_off = 0;
        return;
    }

    double M1M2x = (L->M2->x - L->M1->x) / d;
    double M1M2y = (L->M2->pos - L->M1->pos) / d;
    double M1M2z = (L->M2->z - L->M1->z) / d;

    double fx = L->k * (d - L->l) * M1M2x + L->z * (L->M2->vitx - L->M1->vitx);
    double fy = L->k * (d - L->l) * M1M2y + L->z * (L->M2->vity - L->M1->vity);
    double fz = L->k * (d - L->l) * M1M2z + L->z * (L->M2->vitz - L->M1->vitz);

    L->M1->frcx += fx * L->M1->coeff;
    L->M1->frcy += fy * L->M1->coeff;
    L->M1->frcz += fz * L->M1->coeff;
    L->M2->frcx -= fx * L->M2->coeff;
    L->M2->frcy -= fy * L->M2->coeff;
    L->M2->frcz -= fz * L->M2->coeff;
}

/* quelques couleurs (cf. <g2x_Colors.h>) */
static G3Xcolor Lcols[8] = {G3Xr, G3Xo, G3Xy, G3Xg, G3Xb, G3Xc, G3Xm, G3Xk};

static void DrawM1M2(Link *L) {
    glLineWidth(3.0f);
    g3x_Material(Lcols[L->type], 0.25, 0.7, 0.5, 0.5, 1.);
    glBegin(GL_LINES);
    glVertex3f(L->M1->x, L->M1->pos, L->M1->z);
    glVertex3f(L->M2->x, L->M2->pos, L->M2->z);
    glEnd();
}


static void DrawFrc(Link *L) {
    /*PMat *M = (L->M1 != NULL ? L->M1 : L->M2);
    glLineWidth(1.0f);
    g3x_Material(Lcols[L->type], 0.25, 0.7, 0.5, 0.5, 1.);
    glBegin(GL_LINES);
    glVertex3f(M->x, M->pos, M->z);
    glVertex3f(M->x * 0.001 * L->frcx, M->pos + 0.001 * L->frcy, M->z * 0.001 * L->frcz);
    glEnd();*/
    /*g2x_Line(M->x,M->pos,M->x,M->pos+0.001*L->frc,Lcols[L->type],1);*/
}

/*================================================*/
/*= FONCTIONS DE CREATION DESdoubl MODULES DE LIAISON =*/
/*================================================*/

/*! Création d'un module Force Constante (exp. Gravité) !*/
extern void FrcConst(Link *L, double force_const, const char *axe) {
    L->type = _FRC_CONST;

    /* paramètres pour le moteur physique */
    if (NULL != strchr(axe, 'x')) {
        L->frcx = force_const;
    }
    if (NULL != strchr(axe, 'y')) {
        L->frcy = force_const;
    }
    if (NULL != strchr(axe, 'z')) {
        L->frcz = force_const;
    }
    L->k = 0.;
    L->z = 0.;
    L->s = 0.;
    L->l = 0.;
    L->ax = 0.0f;
    L->ay = 0.0f;
    L->az = 0.0f;
    L->fx = 0.0f;
    L->fy = 0.0f;
    L->fz = 0.0f;
    L->t = 0.0f;
    L->on_off = true;
    L->setup = &AlgoFrcConst;
    /* paramètres graphiques              */
    L->draw = &DrawFrc;
}

extern void FrcExt(Link *L, double ax, double ay, double az, double fx, double fy, double fz) {
    L->type = _FRC_EXT;

    L->k = 0.;
    L->z = 0.;
    L->s = 0.;
    L->l = 0.;
    L->ax = ax;
    L->ay = ay;
    L->az = az;
    L->fx = fx;
    L->fy = fy;
    L->fz = fz;
    L->t = 0.0f;
    L->on_off = true;
    L->setup = &AlgoFrcExt;
    /* paramètres graphiques              */
    L->draw = &DrawFrc;
}

/*! Création d'un ressort linéaire (de Hook) !*/
extern void Ressort(Link *L, double k) {
    L->type = _RESSORT;
    /* paramètres pour le moteur physique */
    L->k = k;
    L->z = 0.;
    L->s = 0.;
    L->l = 0.0f;
    L->ax = 0.0f;
    L->ay = 0.0f;
    L->az = 0.0f;
    L->fx = 0.0f;
    L->fy = 0.0f;
    L->fz = 0.0f;
    L->t = 0.0f;
    L->on_off = true;
    L->setup = &AlgoRessort;
    /* paramètres graphiques              */
    L->draw = &DrawM1M2;
}

/*! Création d'un frein cinétique linéaire !*/
extern void Frein(Link *L, double z) {
    L->type = _FREIN;
    /* paramètres pour le moteur physique */
    L->k = 0.;
    L->z = z;
    L->s = 0.;
    L->ax = 0.0f;
    L->ay = 0.0f;
    L->az = 0.0f;
    L->fx = 0.0f;
    L->fy = 0.0f;
    L->fz = 0.0f;
    L->t = 0.0f;
    L->on_off = true;
    L->setup = &AlgoFrein;
    /* paramètres graphiques              */
    L->draw = &DrawM1M2;
}

/*! Création d'un ressort+frein !*/
extern void RessortFrein(Link *L, double k, double z) {
    L->type = _RESSORT_FREIN;
    /* paramètres pour le moteur physique */
    L->k = k;
    L->z = z;
    L->s = 0.;
    L->l = 0.;
    L->ax = 0.0f;
    L->ay = 0.0f;
    L->az = 0.0f;
    L->fx = 0.0f;
    L->fy = 0.0f;
    L->fz = 0.0f;
    L->t = 0.0f;
    L->on_off = true;
    L->setup = &AlgoRessortFrein;
    /* paramètres graphiques              */
    L->draw = &DrawM1M2;
}

/*! Création d'une butée visco-élastique seuillée !*/
extern void RF_Butee(Link *L, double k, double z, double s) {
    L->type = _RF_BUTEE;
    /* paramètres pour le moteur physique */
    L->k = k;
    L->z = z;
    L->s = s; /* seuil de distance pour détachement   */
    L->l = 0.;
    L->ax = 0.0f;
    L->ay = 0.0f;
    L->az = 0.0f;
    L->fx = 0.0f;
    L->fy = 0.0f;
    L->fz = 0.0f;
    L->t = 0.0f;
    L->on_off = true;
    L->setup = &AlgoRF_Butee;
    /* paramètres graphiques              */
    L->draw = &DrawM1M2;
}


/*! Création d'une liaison de rupture avec condition sur l'élongation !*/
extern void RF_CondPos(Link *L, double k, double z, double s) {
    L->type = _CONDIT_RF;
    /* paramètres pour le moteur physique */
    L->k = k;
    L->z = z;
    L->s = s; /* seuil d'elongation pour rupture      */
    L->l = 0.;
    L->ax = 0.0f;
    L->ay = 0.0f;
    L->az = 0.0f;
    L->fx = 0.0f;
    L->fy = 0.0f;
    L->fz = 0.0f;
    L->t = 0.0f;
    L->on_off = true;
    L->setup = &AlgoRF_CondPos;
    /* paramètres graphiques              */
    L->draw = &DrawM1M2;
}


/*! Connexion d'une Liaison entre 2 points Mat !*/
/* pas utilisé dans sytèmes de particules       */
/* ça sert surtout pour les topologies fixes    */
extern void Connect(PMat *M1, Link *L, PMat *M2) {
    L->M1 = M1;
    L->M2 = M2;
    if (M1 == NULL || M2 == NULL) return;
    double distance = sqrt(pow((M2->x - M1->x), 2.0) + pow((M2->pos - M1->pos), 2.0) + pow((M2->z - M1->z), 2.0));
    L->l = distance;
}
