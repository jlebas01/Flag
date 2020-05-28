/*=============================================================*/
/*= E.Incerti (incerti@univ-eiffel.fr)                        =*/
/*= M2/IMAC3 - ANIMATION / MOTEUR PHYSIQUE                    =*/
/*= Système Masses/Ressort en dimension 1 : Corde             =*/
/*= Modèle 2 : avec force extérieure ("gravité")              =*/
/*=============================================================*/

#include <g2x.h>

#include <PMat.h>
#include <Link.h>

/*! Paramètres de la fenetre graphique !*/
static int      pixwidth=640,pixheight=640;
static double   xmin=-6.,ymin=-6.,xmax=+6.,ymax=+6.;

/*! VARIABLES DE SIMULATION !*/
static int      Fe;        /* fréquence d'échantillonnage        */
static double   h;         /* pas d'échantillonnage 1./Fech      */
static int      Fa=1;      /* on affiche 1 echantillon sur Faff  */
static double   m,k,z,g;   /* paramètres physiques               */
static double   tempo=0.02;/* temporisation de la simul          */ 

/*! systeme "Masses-Ressorts" : les particules et les liaisons !*/
static int      nbm;
static PMat    *TabM=NULL;
static int      nbl;
static Link    *TabL=NULL;

/*=================================================================*/
/*= sic.                                                          =*/
/*=================================================================*/
void reset(void)
{	
  int i;
  for (i=0;i<nbm;i++) TabM[i].pos=TabM[i].vit=TabM[i].frc=0.;
}


/*=-------------------------=*/
/*=  Constructon du modèle  =*/
/*=-------------------------=*/
bool Modeleur(void)
{
  /*! le Modèle : un tableau de particules, un tableau de liaisons !*/
  nbm = 11;              /* 9 particules et 2 points fixes     */
  if (!(TabM=(PMat*)calloc(nbm,sizeof(PMat)))) return false;
  nbl = (nbm-1)+(nbm-2); /* 10 ressorts-freins + 9 "gravités"  */
  if (!(TabL=(Link*)calloc(nbl,sizeof(Link)))) return false;

  Fe= 100;           /* paramètre du simulateur Fe=1/h                  */ 
  h = 1./Fe;
  m = 1.;            /* la plupart du temps, toutes les masses sont à 1 */
  /* avec les paramètres ci-dessous : limite de stabilité               */
  k = 0.866*SQR(Fe); /* raideurs   -> à multiplier par mb*Fe²           */
  z = 0.08*Fe;       /* viscosirés -> à multiplier par mb*Fe            */
  g = -5.*Fe;        /* la gravité : elle aussi à calibrer avec Fe      */
  

  /*! les particules !*/
  int i;
  PMat* M=TabM;
  Fixe(M++,0.,-5.);
  for (i=1;i<nbm-1;i++) MassLF(M++,0.,-5.+i,m); 
  Fixe(M++,0.,+5.);

  /*! les liaisons !*/  
  Link* L=TabL;
  for (i=0;i<nbm-1;i++) RessortFrein(L++,k,z);
  for (i=0;i<nbm-2;i++) FrcConst    (L++,g  ); 
  
  /*! les connections => topologie fixe !*/
  L=TabL;
  /* les ressorts inter-particules */
  M=TabM;
  while (L<TabL+nbm-1) 
  {
    Connect(M,L,M+1);
    M++;
    L++;
  }
  /* les "gravités" individuelles  */
  M=TabM+1; /* 1° particule mobile */
  while (L<TabL+nbl) 
  {
    Connect(M,L,NULL);
    M++;
    L++;
  }

  reset();
  return true;
}


/*=------------------------------=*/
/*=        Initialisation        =*/
/*=------------------------------=*/
void init(void)
{
  /* construction du modèle   */
	if (!Modeleur()) exit(1);
  
  /* les paramètres réglables */
  g2x_CreateScrollv_i("Fa",&Fa,1,20,1,"Fa");
  g2x_CreateScrollv_d("tmp",&tempo,0.,.1,1.,"tempo");
	g2x_CreateScrollh_d("k",&k,k*0.01,k*5.,1,"k");
	g2x_CreateScrollh_d("z",&z,z*0.01,z*5.,1,"z");
	g2x_CreatePopUp("reset",reset,"reset");
}

/*=------------------------------=*/
/*=          L'Affichage         =*/
/*=------------------------------=*/
void dessin(void)
{
	/* fréquence d'affichage réglable */
	g2x_SetRefreshFreq(Fa);
  PMat *M=TabM;
  while (M<TabM+nbm)  { M->draw(M); ++M; }
  Link *L=TabL;  
  while (L<TabL+nbl)  { L->draw(L); L->k=k; L->z=z; ++L; } /* mise à jour des paramètres => scrollbar */
}


/*=------------------------------=*/
/*=        Le Simulateur         =*/
/*=------------------------------=*/
void Moteur_Physique(void)
{  
  PMat *M=TabM;
  while (M<TabM+nbm) { M->setup(M,h); ++M; } 
  Link *L=TabL;
  while (L<TabL+nbl) { L->setup(L)  ; ++L; } 
  g2x_tempo(tempo); /* temporisation, si ça va trop vite */
}


/*=-------------------------=*/
/*= Fonction de libération  =*/
/*= appelée à la fin.       =*/
/*=-------------------------=*/
void quit(void)
{
  if (TabM!=NULL) free(TabM); 
  if (TabL!=NULL) free(TabL); 
}

/*=-------------------------=*/
/*=                         =*/
/*=-------------------------=*/
int main(int argc, char* argv[])
{
  g2x_InitWindow(*argv,pixwidth,pixheight);  
  g2x_SetWindowCoord(xmin,ymin,xmax,ymax);

  g2x_SetInitFunction(init);
  g2x_SetAnimFunction(Moteur_Physique);
  g2x_SetDrawFunction(dessin);          
  g2x_SetExitFunction(quit);      

  return g2x_MainStart();
}
