/*=================================================================*/
/*= E.Incerti - incerti@u-pem.fr                                   =*/
/*= Université Paris-Est-Marne-la-Vallée                          =*/
/*= Balistic launch : analytic and step-by-step solutions         =*/
/*=================================================================*/
#include <g2x.h>

/* window dimensions (pixels)    */
uint   pixwidth=1400,pixheight=700;
/* drawing zone (floating-point) */
double xmin=-5.,ymin=-5.,xmax=+205.,ymax=+105.;

G2Xpoint  P0;          /* initial position and speed */  
G2Xvector V0;
double g,m,za,kc,zc;   /* 'Physical' parameters      */
double dt;             /* analytic time step         */


/*=================================================================*/
/*= Initialize                                                    =*/
/*=================================================================*/
void Init(void)
{
  g = 10.; /* gravity                                           */
  m = 1.0;  /* ball mass                                         */
  za= 0.01;  /* air friction (elementary kinetic damping)         */
  kc= 0.9;  /* collision treatment : inverse kinetics            */        
  zc= 0.9;  /*                                                   */
  P0=(G2Xpoint) {  0., 25.};   /* initial position and speed     */
  V0=(G2Xvector){ 10., 40.};
  dt=0.05;     /* analytic   time step */

  /* scrollbars and buttons */
	g2x_CreateScrollv_d("za",&za,0.   ,0.1,1.,"za");
	g2x_CreateScrollh_d("dt",&dt,0.001,0.1,1.,"dt");
}



/*======================================================================*/
/*= Computes an draws analytic solution for initial conditions (p0,vO) =*/
/*======================================================================*/
void Analytic(G2Xpoint* pi, G2Xvector* vi)
{
	G2Xpoint  p,q;	
  G2Xvector v;
	double    t = 0;
  double    w=m/za,e;
	
  switch (G2Xzero(za))
  {
    case true : /* without "air" friction */
      q=*pi;  
      do
      {
        t+=dt;
			  /* setup position */
			  p.x = pi->x + vi->x*t;
			  p.y = pi->y + vi->y*t- 0.5*g*t*t;		
			  /* setup speed    */
			  v.x = vi->x;
			  v.y = vi->y - g*t;		
        /* draw           */
        g2x_Line(q.x,q.y,p.x,p.y,G2Xwa,1);
        g2x_Plot(p.x,p.y,G2Xwb,3);
        q=p;
     } while (p.y+dt*v.y>0.);
    break;  
    default : /* with "air" friction */
      q=*pi;  
      do
      {
        t += dt;
        e  = exp(-t/w);
			  /* setup position */
			  p.x = pi->x + w*(      vi->x)*(1.-e);
			  p.y = pi->y + w*(g*w + vi->y)*(1.-e) - (g*w)*t;		
			  /* setup speed    */
			  v.x = (    vi->x)*e;
			  v.y = (g*w+vi->y)*e - g*w;	
        /* draw           */
        g2x_Line(q.x,q.y,p.x,p.y,G2Xwa,1);
        g2x_Plot(p.x,p.y,G2Xwb,3);
        q=p;
      } while (p.y+dt*v.y>0.);
  }
  /* setup initial position */
  *pi = q;
  /* setup initial speed    */
  vi->x =  zc*v.x;
  vi->y = -kc*v.y;
}


/*=================================================================*/
/*=                                                               =*/
/*=================================================================*/
void Draw(void)
{
	G2Xpoint  p=P0;	
	G2Xvector v=V0;		
	g2x_Axes();
  do
  {
    Analytic(&p,&v);
  }
  while (!G2Xzero(v.x) && p.x<xmax);
}


/*===========================================================================*/
/*= Cleaning function                                                       =*/
/*===========================================================================*/
void Quit(void)
{
  /* nothing to do here */
}

/*===========================================================================*/
/*=                                                                         =*/
/*===========================================================================*/
int main(int argc, char* argv[])
{	  
  /* window statement */
  g2x_InitWindow(argv[0],pixwidth,pixheight);
  g2x_SetWindowCoord(xmin,ymin,xmax,ymax);
  /* dialog functions */
  g2x_SetInitFunction(Init);
  g2x_SetExitFunction(Quit);
	g2x_SetDrawFunction(Draw);  

  /* simulation loop -- passed to glutMainLoop() */
  return g2x_MainStart();
}
