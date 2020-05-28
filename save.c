/*=============================================================*/
/*= E.Incerti (incerti@univ-eiffel.fr)                        =*/
/*= M2/IMAC3 - ANIMATION / MOTEUR PHYSIQUE                    =*/
/*= Système Masses/Ressort en dimension 1 : Oscillateur       =*/
/*= Modèle 0 : juste une masse, un point fixe et un ressort   =*/
/*= Rupture d'équilibre par impulsion de force ou déplacement =*/
/*=============================================================*/
#include <g3x.h>

#include <PMat.h>
#include <Link.h>

#define NB_FIXE 20
#define NB_MOVE 20

/*! Paramètres de la fenetre graphique !*/
static int pixwidth = 640, pixheight = 640;
static double xmin = -6., ymin = -6., xmax = +6., ymax = +6.;

/*! VARIABLES DE SIMULATION !*/
static double h;          /* pas d'échantillonnage 1./Fech      */
static int Fe;         /* fréquence d'échantillonnage        */
static int Fa = 1.0;       /* on affiche 1 echantillon sur Faff  */
static double m, k, z;      /* paramètres physiques               */
static char PouF = 'F'; /* type de perturbation initiale      */
static double tempo = 0.02; /* temporisation de la simul          */



/*! systeme de particules  !*/
/*static PMat M[NB_MOVE], S[NB_FIXE];*/
static PMat M[NB_FIXE][NB_MOVE];
static PMat M2[NB_FIXE][NB_MOVE];
static Link *L = NULL;
static Link *tmp = NULL;
static int NB_LINK = 0;

/*=================================================================*/
/*= sic.                                                          =*/
/*=================================================================*/

void alloc_one_more_link(Link **L) {
    tmp = (Link *) realloc(*L, ++NB_LINK * sizeof(Link));
    if (NULL == tmp) {
        exit(EXIT_FAILURE);
    }
    *L = tmp;
    /* free(tmp);*/
}

void reset(void) {
    int i = 0;
    int j = 0;
    double px = 0.02f, py = 0.0f, pz = 0.0f;
    for (i = 0; i < NB_FIXE; i++) {
        M[i][0].x = 0.0f;
        M[i][0].pos = 0.0f;
        M[i][0].z = -pz;

        M2[i][0].x = px;
        M2[i][0].pos = 0.0f;
        M2[i][0].z = -pz;
        pz += 0.1f;
    }
    pz = 0.0f;
    for (i = 0; i < NB_FIXE; i++) {
        py = 0.1f;
        for (j = 1; j < NB_MOVE; j++) {
            M[i][j].x = 0.0f;
            M[i][j].pos = py;
            M[i][j].z = -pz;

            M2[i][j].x = px;
            M2[i][j].pos = py;
            M2[i][j].z = -pz;
            py += 0.1f;
        }
        pz += 0.1f;
    }
    for (i = 0; i < NB_FIXE; i++) {
        for (j = 0; j < NB_MOVE; j++) {

            M[i][j].vitx = M[i][j].frcx = 0.;
            M[i][j].vity = M[i][j].frcy = 0.;
            M[i][j].vitz = M[i][j].frcz = 0.;

            M2[i][j].vitx = M2[i][j].frcx = 0.;
            M2[i][j].vity = M2[i][j].frcy = 0.;
            M2[i][j].vitz = M2[i][j].frcz = 0.;
            switch (PouF) {
                case 'P' :
                case 'p' : /* contrainte initiale : particule déplacée */
                    /*M.pos = .5;*/
                    break;
                case 'F' :
                    /*M[i][j].frcz += -(M[i][j].m * 9.80665f) ;*/
                case 'f' : /* Autre forme de contrainte initiale : impulsion de force */
                    M[i][j].frcx = 0.7f * SQR(Fe);   /* notez le calibrage par Fe² */
                    M[i][j].frcy = 0.7f * SQR(Fe);   /* notez le calibrage par Fe² */
                    M[i][j].frcz = 0.7f * SQR(Fe);   /* notez le calibrage par Fe² */

                    M2[i][j].frcx = 0.7f * SQR(Fe);   /* notez le calibrage par Fe² */
                    M2[i][j].frcy = 0.7f * SQR(Fe);   /* notez le calibrage par Fe² */
                    M2[i][j].frcz = 0.7f * SQR(Fe);   /* notez le calibrage par Fe² */
                    break;
            }
        }
    }
}

/*=-------------------------=*/
/*=  Constructon du modèle  =*/
/*=-------------------------=*/
void Modeleur(void) {
    Fe = 100;           /* paramètre du simulateur Fe=1/h                  */
    h = 1. / Fe;
    m = 1.0f;            /* la plupart du temps, toutes les masses sont à 1 */
    k = 0.04 * SQR(Fe);  /* raideurs   -> à multiplier par mb*Fe²                    */
    z = 0.007f * Fe;      /* viscosités -> à multiplier par mb*Fe                     */
    /*z=0.0024;
    k=0.00965;*/

    /*! les objets : une masse mobile (M) reliée à un pt. fixe (S) par un ressort-frein (L) !*/
    int i = 0;
    int j = 0;
    double px = 0.02f, py = 0.0f, pz = 0.0f;
    for (i = 0; i < NB_FIXE; i++) {
        Fixe(&M[i][0], 0.0f, py, -pz);
        Fixe(&M2[i][0], px, py, -pz);
        pz += 0.1f;
    }
    pz = 0.0f;
    double CRD = 1.0f; /*Coefficient de rendement décroissant pour la masse*/
    double max = 100.0f;
    double min = 10.0f;
    double reliquat = ((max - min) / (NB_MOVE - 1) / 100.0f);

    for (i = 0; i < NB_FIXE; i++) {
        py = 0.1f;
        CRD = 1.0f;
        for (j = 1; j < NB_MOVE; j++) {
            /*Fixe(&M[i][j], 0.0f, px, -pz);*/
            /* la même chose, mais avec l'intégrateur Euler Explicite */
            /*MassEE(&M[i][j], 0.0f, px, -pz, m);*/
            MassLF(&M[i][j], 0.0f, py, -pz, m * CRD);
            MassLF(&M2[i][j], px, py, -pz, m * CRD);
            CRD -= reliquat;
            py += 0.1f;
        }
        pz += 0.1f;
    }

    /* connection du ressort entre les 2 masses */


    for (i = 0; i < NB_FIXE; i++) {
        for (j = 0; j < NB_MOVE; j++) {
            /* lien horizontal +1 */
            if (j + 1 < NB_MOVE) {
                alloc_one_more_link(&L);
                RessortFrein(&L[NB_LINK - 1], k, z);
                Connect(&M[i][j], &L[NB_LINK - 1], &M[i][j + 1]);

                alloc_one_more_link(&L);
                RessortFrein(&L[NB_LINK - 1], k, z);
                Connect(&M2[i][j], &L[NB_LINK - 1], &M2[i][j + 1]);

                alloc_one_more_link(&L);
                RessortFrein(&L[NB_LINK - 1], k, z);
                Connect(&M[i][j], &L[NB_LINK - 1], &M2[i][j + 1]);

                alloc_one_more_link(&L);
                RessortFrein(&L[NB_LINK - 1], k, z);
                Connect(&M2[i][j], &L[NB_LINK - 1], &M[i][j + 1]);
            }
            /* lien horizontal +2 */
            if (j + 2 < NB_MOVE) {
                alloc_one_more_link(&L);
                RessortFrein(&L[NB_LINK - 1], k, z);
                Connect(&M[i][j], &L[NB_LINK - 1], &M[i][j + 2]);

                alloc_one_more_link(&L);
                RessortFrein(&L[NB_LINK - 1], k, z);
                Connect(&M2[i][j], &L[NB_LINK - 1], &M2[i][j + 2]);
            }
            /* lien vertical +1 */
            if (i + 1 < NB_FIXE) {
                alloc_one_more_link(&L);
                RessortFrein(&L[NB_LINK - 1], k, z);
                Connect(&M[i][j], &L[NB_LINK - 1], &M[i + 1][j]);

                alloc_one_more_link(&L);
                RessortFrein(&L[NB_LINK - 1], k, z);
                Connect(&M2[i][j], &L[NB_LINK - 1], &M2[i + 1][j]);
            }
            /* lien vertical +2 */
            if (i + 2 < NB_FIXE) {
                alloc_one_more_link(&L);
                RessortFrein(&L[NB_LINK - 1], k, z);
                Connect(&M[i][j], &L[NB_LINK - 1], &M[i + 2][j]);

                alloc_one_more_link(&L);
                RessortFrein(&L[NB_LINK - 1], k, z);
                Connect(&M2[i][j], &L[NB_LINK - 1], &M2[i + 2][j]);
            }
            /* lien diagonal nord-ouest -> sud-est*/
            if ((i + 1 < NB_FIXE) && (j + 1 < NB_MOVE)) {
                alloc_one_more_link(&L);
                RessortFrein(&L[NB_LINK - 1], k, z);
                Connect(&M[i][j], &L[NB_LINK - 1], &M[i + 1][j + 1]);

                alloc_one_more_link(&L);
                RessortFrein(&L[NB_LINK - 1], k, z);
                Connect(&M2[i][j], &L[NB_LINK - 1], &M2[i + 1][j + 1]);
            }
            /* lien diagonal nord-est -> sud-ouest*/
            if ((i + 1 < NB_FIXE) && (j - 1 >= 0)) {
                alloc_one_more_link(&L);
                RessortFrein(&L[NB_LINK - 1], k, z);
                Connect(&M[i][j], &L[NB_LINK - 1], &M[i + 1][j - 1]);

                alloc_one_more_link(&L);
                RessortFrein(&L[NB_LINK - 1], k, z);
                Connect(&M2[i][j], &L[NB_LINK - 1], &M2[i + 1][j - 1]);
                /*printf("[%d][%d]\n",i+1,j-1);*/
            }
            alloc_one_more_link(&L);
            RessortFrein(&L[NB_LINK - 1], k, z);
            Connect(&M[i][j], &L[NB_LINK - 1], &M2[i][j]);
        }
    }


    printf("nb liens %d \n", NB_LINK);

    reset();

}

/*=------------------------------=*/
/*=        Initialisation        =*/
/*=------------------------------=*/
void init(void) {
    /* construction du modèle   */
    Modeleur();

    /* les paramètres réglables */
    g3x_CreateScrollv_i("Fa", &Fa, 1, 20, 1, "Fa");
    g3x_CreateScrollv_d("tmp", &tempo, 0., .1, 1., "tempo");
    g3x_CreateScrollh_d("k", &k, k * 0.01, k * 5., 1, "k");
    g3x_CreateScrollh_d("z", &z, z * 0.01, z * 5., 1, "z");
    g3x_CreatePopUp("reset", reset, "reset");
}


/*=------------------------------=*/
/*=          L'Affichage         =*/
/*=------------------------------=*/
void dessin(void) { /* fréquence d'affichage réglable */
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    glEnable(GL_BLEND);
    /*glLoadIdentity();*/
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    int i = 0;
    int j = 0;
    g3x_Material(G3Xr, 0.25, 0.7, 0.5, 0.5, 1.);
    for (i = 0; i < NB_FIXE; i++) {
        for (j = 0; j < NB_MOVE; j++) {
            M[i][j].draw(&M[i][j]);
            M2[i][j].draw(&M2[i][j]);
        }
    }
    for (i = 0; i < NB_LINK; i++) {
        L[i].draw(&L[i]);
    }
    g3x_SetRefreshFreq(Fa);
    for (i = 0; i < NB_LINK; i++) {
        L[i].k = k;
        L[i].z = z; /* mise à jour des paramètres => scrollbar */
    }
}


/*=------------------------------=*/
/*=        Le Simulateur         =*/
/*=------------------------------=*/
void Moteur_Physique(void) {
    /* d'abord les PMat */
    int i = 0;
    int j = 0;
    for (i = 0; i < NB_FIXE; i++) {
        for (j = 0; j < NB_MOVE; j++) {
            M[i][j].setup(&M[i][j], h);
            M2[i][j].setup(&M2[i][j], h);
        }
    }
    /* ensuite les Link */
    for (i = 0; i < NB_LINK; i++) {
        L[i].setup(&L[i]);
    }

    g3x_tempo(tempo); /* temporisation, si ça va trop vite */
}


/*=-------------------------=*/
/*= Fonction de libération  =*/
/*= appelée à la fin.       =*/
/*=-------------------------=*/
void quit(void) {
    /* rien ici */
}

/*=-------------------------=*/
/*=                         =*/
/*=-------------------------=*/
int main(int argc, char *argv[]) {

    /* initialisation de la fenêtre graphique et paramétrage Gl */
    g3x_InitWindow(*argv, pixwidth, pixheight);
    /* param. géométrique de la caméra. cf. gluLookAt(...) */
    g3x_SetPerspective(40., 100., 1.);
    /* position, orientation de la caméra */
    g3x_SetCameraSpheric(0.25 * PI, +0.25 * PI, 6., (G3Xpoint) {0., 0., 0.});

    /*g3x_CreateScrollv_i("n",&n,3,100,1,"...");*/

    /* définition des fonctions */
    g3x_SetInitFunction(init);
    g3x_SetAnimFunction(Moteur_Physique);
    g3x_SetDrawFunction(dessin);
    g3x_SetExitFunction(quit);

    /*g3x_SetScrollWidth(2);
    g3x_CreateScrollv_i("n",&n,3,100,1.,"résolution");
    g3x_CreateScrollh_i("n",&n,3,100,1.,"résolution");*/
    /* Exit();*/
    /* boucle d'exécution principale */

    return g3x_MainStart();
}
