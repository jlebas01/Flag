https://github.com/jlebas01/Flag.git/*=================================================================*/
/*= E.Incerti - incerti@upem.fr                                   =*/
/*= Universit� Paris-Est-Marne-la-Vall�e                          =*/
/*= Balistic launch : analytic and step-by-step solutions         =*/
/*=================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include <g2x.h>

/* pixel window dimensions          */
uint   pixwidth=1400,pixheight=700;
/* floating-point window dimensions */
double xmin=-5.,ymin=-5.,xmax=+205.,ymax=+105.;

/* 'Physical' parameters                                                       */
double m,za,kc,zc,r;
/* initial position and speed */
G2Xpoint  P0;
G2Xvector V0;
/* current positions Euler, Implicit, Leapfrog */
G2Xpoint  Pe,Pi,Pl;
/* speed, acceleration  Euler, Implicit, Leapfrog, gravity */
G2Xvector Ve,Vi,Vl,Ae,Ai,Al,G;

/* analytic & simulation time step */
double dt,h;

/* FLAGS to show/hide simulation process */
bool FLAG_EXPLICIT_EULER=true,
     FLAG_IMPLICIT_EULER=true,
     FLAG_LEAPFROG=true;
     
/*=================================================================*/
/*= sic.                                                          =*/
/*=================================================================*/
void reset(void)
{	
	Pe=Pi=Pl=P0;
	Ve=Vi=Vl=V0;
	Ae.x=Ai.x=Al.x=0.;
	Ae.y=Ai.y=Al.y=0.;
}

/*=================================================================*/
/*= Initialize                                                    =*/
/*=================================================================*/
void Init(void)
{
  G =(G2Xvector){  0.,-10.}; /* gravity force                    */
  m = 1.0;  /* ball mass                                         */
  za= 0.01; /* air friction (elementary kinetic damping)         */
  kc= 0.9;  /* collision treatment : inverse kinetics            */        
  zc= 0.9;  /*                                                   */
  r = 2.0;  /* ball ray                                          */
  P0=(G2Xpoint) {  0., 25.};   /* initial position and speed     */
  V0=(G2Xvector){ 10., 40.};
  dt=0.05;     /* analytic   time step */
  h =0.05;     /* simulation time step */
  reset();

  /* scrollbars and buttons */
	g2x_CreateScrollv_d("za",&za,0.,0.1,1.,"za");
	g2x_CreateScrollh_d("dt",&dt,0.001,0.1,1.,"dt");
	g2x_CreateScrollh_d(" h",&h ,0.001,0.1,1.,"h");
	g2x_CreatePopUp("reset",reset,"reset");
	g2x_CreateSwitch("EE",&FLAG_EXPLICIT_EULER,"affiche/masque Euler Explicit");
	g2x_CreateSwitch("IE",&FLAG_IMPLICIT_EULER,"affiche/masque Euler Implicit");
	g2x_CreateSwitch("LF",&FLAG_LEAPFROG      ,"affiche/masque LeapFrog      ");
}

/*======================================================================*/
/*= Computes an draws analytic solution for initial conditions (p0,vO) =*/
/*======================================================================*/
void Analytic(G2Xpoint* pi, G2Xvector* vi)
{
	G2Xpoint  p,q;	
  G2Xvector v;
	double    t = 0;
  double    e,w=m/za;
	
  q=*pi;  
	/* two sets of "continuous time" formulation : with or without "air friction" */
  switch (G2Xzero(za)) /* G2Xzero(x) <=> (fabs(x)<EPSILON?true:false) */
  {
    case true : /* without "air" friction */
      do
      {
        t+=dt;
			  /* setup current position */
			  p.x = pi->x + vi->x*t + 0.5*G.x*t*t;
			  p.y = pi->y + vi->y*t + 0.5*G.y*t*t;		
        /* draw                   */
        g2x_Line(q.x,q.y,p.x,p.y,G2Xwa,1);
        g2x_Plot(p.x,p.y,G2Xwb,3);
        q=p;
      } while (p.y>r);
      /* setup collision speed    */
      v.x = vi->x + G.x*t;
      v.y = vi->y + G.y*t;
    break;  
    default :  /* with "air" friction */
      do
      {
        t += dt;
        e  = exp(-t/w);
			  /* setup current position */
			  p.x = pi->x + (G.x*w)*t + w*(vi->x - G.x*w)*(1.-e);		
			  p.y = pi->y + (G.y*w)*t + w*(vi->y - G.y*w)*(1.-e);		
        /* draw                   */
        g2x_Line(q.x,q.y,p.x,p.y,G2Xwa,1);
        g2x_Plot(p.x,p.y,G2Xwb,3);
        q=p;
      } while (p.y>r);
      /* setup collision speed    */
      v.x = (vi->x - G.x*w)*e + G.x*w;	
      v.y = (vi->y - G.y*w)*e + G.y*w;	
  }
  /* setup initial position */
  *pi = q;
  /* setup initial speed    */
  vi->x =  zc*v.x;
  vi->y = -kc*v.y;
}

/*=================================================================*/
/*= Analytic solution is drawn and saved as backgroung image      =*/
/*=================================================================*/
void Background(void)
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



/*=================================================================*/
/*= step-by-step function called by the simulation loop           =*/
/*=================================================================*/
/*= Explicit method =*/
void Simul_Explicit(void)
{  
	/* 1 - set up ball position */
	Pe.x = Pe.x + h*Ve.x;
	Pe.y = Pe.y + h*Ve.y;
	/* 2 - set up ball velocity */
	Ve.x = Ve.x + h*Ae.x; 
	Ve.y = Ve.y + h*Ae.y; 
  /* 3 - computes the resulting force applied to the ball       */
	Ae.x = G.x-(za/m)*Ve.x;
	Ae.y = G.y-(za/m)*Ve.y;
	if (Pe.y<=r) 	/* collision detection */
	{
		Ve.x *=  zc;  
		Ve.y *= -kc; /* collision treatment : inverse kinetics */
		Pe.y  = r;
	}
}

/*= Implicit method =*/
void Simul_Implicit(void)
{  
	/* 1 - set up ball velocity */
	double w = m/(m+za*h);
	Vi.x = (Vi.x + G.x*h)*w; 
	Vi.y = (Vi.y + G.y*h)*w; 
	/* 2 - set up ball position */
	Pi.x = Pi.x + h*Vi.x;
	Pi.y = Pi.y + h*Vi.y;
  /* 3 - computes the resulting force applied to the ball       */
	Ai.x = G.x-(za/m)*Vi.x;
	Ai.y = G.y-(za/m)*Vi.y;
	if (Pi.y<=r) 	/* collision detection */
	{
		Vi.x *=  zc;  
		Vi.y *= -kc; /* collision treatment : inverse kinetics */
		Pi.y  = r;
	}
}


/*= Hybrid method =*/
void Simul_LeapFrog(void)
{  
  /* 1 - computes the resulting force applied to the ball       */
	Al.x = G.x-(za/m)*Vl.x;
	Al.y = G.y-(za/m)*Vl.y;
	/* 2 - set up ball velocity */
	Vl.x = Vl.x + h*Al.x; 
	Vl.y = Vl.y + h*Al.y; 
	/* 3 - set up ball position */
	Pl.x = Pl.x + h*Vl.x;
	Pl.y = Pl.y + h*Vl.y;
	if (Pl.y<=r) 	/* collision detection */
	{
		Vl.x *=  zc;  
		Vl.y *= -kc; /* collision treatment : inverse kinetics */
		Pl.y  = r;
	}
}

void Anim(void)
{
  if (FLAG_EXPLICIT_EULER) Simul_Explicit();
  if (FLAG_IMPLICIT_EULER) Simul_Implicit();
  if (FLAG_LEAPFROG)       Simul_LeapFrog();

	/* the ball leaves the window, back to initial conditions */
	if (Pe.y<=0.   || Pi.y<=0.   || Pl.y<=0. || 
	    Pe.x>=xmax || Pi.x>=xmax || Pl.x>=xmax) reset();
  else g2x_tempo(0.1*h);
}

/*===========================================================================*/
/*=                                                                         =*/
/*===========================================================================*/
void Draw(void)
{
	Background(); /* call the background image */
	if (FLAG_EXPLICIT_EULER) { g2x_Circle(Pe.x,Pe.y,r,G2Xb ,3); /* Explicit Euler */ 
                             g2x_StaticPrint(50,pixheight-20,G2Xb,'l',"Explicit");
                           } 
  if (FLAG_IMPLICIT_EULER) { g2x_Circle(Pi.x,Pi.y,r,G2Xg ,3); /* Imlicit Euler  */
	                           g2x_StaticPrint(50,pixheight-30,G2Xg,'l',"Implicit");
                           }
	if (FLAG_LEAPFROG)       { g2x_Circle(Pl.x,Pl.y,r,G2Xr ,3); /* Leapfrog       */
	                           g2x_StaticPrint(50,pixheight-40,G2Xr,'l',"LeapFrog");
                           }
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
  g2x_InitWindow("Balistic",pixwidth,pixheight);
  g2x_SetWindowCoord(xmin,ymin,xmax,ymax);
  /* dialog functions */
  g2x_SetInitFunction(Init);
	g2x_SetDrawFunction(Draw);  
	g2x_SetAnimFunction(Anim);  
  g2x_SetExitFunction(Quit);

  /* simulation loop -- passed to glutMainLoop() */
  return g2x_MainStart();
}
