/*=============================================================*/
/*= E.Incerti (incerti@univ-eiffel.fr)                        =*/
/*= M2/IMAC3 - ANIMATION / MOTEUR PHYSIQUE                    =*/
/*= Syst�me Masses/Ressort en dimension 1 : Corde             =*/
/*= Mod�le 1 : sans force ext�rieure (pas de "gravit�")       =*/
/*= Rupture d'�quilibre par impulsion de force ou d�placement =*/  
/*=============================================================*/

#include <g2x.h>

#include <PMat.h>
#include <Link.h>

/*! Param�tres de la fenetre graphique !*/
static int    pixwidth=640,pixheight=640;
static double xmin=-6.,ymin=-6.,xmax=+6.,ymax=+6.;

/*! VARIABLES DE SIMULATION !*/
static double   h;          /* pas d'�chantillonnage 1./Fech      */
static int      Fe;         /* fr�quence d'�chantillonnage        */
static int      Fa=1;       /* on affiche 1 echantillon sur Faff  */
static double   m,k,z;      /* param�tres physiques               */
static char     PouF = 'f'; /* type de perturbation initiale      */
static double   tempo=0.02; /* temporisation de la simul          */ 

/*! systeme "Masses-Ressorts" : les particules et les liaisons !*/
static int      nbm;
static PMat    *TabM=NULL;
static int      nbl;
static Link    *TabL=NULL;

/*=================================================================*/
/*= sic.                                                          =*/
/*=================================================================*/
void reset()
{	
  int i;
  for (i=0;i<nbm;i++) TabM[i].pos=TabM[i].vit=TabM[i].frc=0.;
  switch (PouF)
  {
    case 'P' : case 'p' : /* contrainte initiale : particule d�plac�e */
      TabM[1].pos=1.;
      break;
    case 'F' : case 'f' : /* Autre forme de contrainte initiale : impulsion de force */
      TabM[1].frc=1.*SQR(Fe); /* notez le calibrage par Fe� */
      break;
  }
}
/*=-------------------------=*/
/*=  Constructon du mod�le  =*/
/*=-------------------------=*/
bool Modeleur(void)
{
  /*! le Mod�le : un tableau de particules, un tableau de liaisons !*/
  nbm = 11;
  if (!(TabM=(PMat*)calloc(nbm,sizeof(PMat)))) return false;
  nbl = nbm-1;
  if (!(TabL=(Link*)calloc(nbl,sizeof(Link)))) return false;

  
  Fe= 100;           /* param�tre du simulateur Fe=1/h                  */ 
  h = 1./Fe;
  m = 1.;            /* la plupart du temps, toutes les masses sont � 1 */
  /* avec les param�tres ci-dessous : limite de stabilit�               */
  k = 0.865*SQR(Fe); /* raideurs   -> � multiplier par mb*Fe�           */
  z = 0.08*Fe;       /* viscosir�s -> � multiplier par mb*Fe            */

  /*! les particules !*/
  int i;
  PMat* M=TabM;
  Fixe(M++,0.,-5.);
  for (i=1;i<nbm-1;i++) MassLF(M++,0.,-5.+i,m); 
  Fixe(M++,0.,+5.);
    
  /*! les liaisons !*/  
  Link* L=TabL;
  for (i=0;i<nbl;i++) RessortFrein(L++,k,z);
  
  /*! les connections => topologie fixe !*/
  L=TabL;
  M=TabM;
  while (L<TabL+nbl) 
  {
    Connect(M,L,M+1);
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
  /* construction du mod�le   */
	if (!Modeleur()) exit(1);
  
  /* les param�tres r�glables */
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
	/* fr�quence d'affichage r�glable */
	g2x_SetRefreshFreq(Fa);
  PMat *M=TabM;
  while (M<TabM+nbm)  { M->draw(M); ++M; }
  Link *L=TabL;
  while (L<TabL+nbl)  { L->draw(L); L->k=k; L->z=z; ++L; } /* mise � jour des param�tres => scrollbar */
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
  g2x_tempo(tempo); /* temporisation, si �a va trop vite */
}


/*=-------------------------=*/
/*= Fonction de lib�ration  =*/
/*= appel�e � la fin.       =*/
/*=-------------------------=*/
void quit(void)
{
  if (TabM) free(TabM); 
  if (TabL) free(TabL); 
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
