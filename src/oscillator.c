/* ============================================================= */
/*  E.Incerti (incerti@univ-mlv.fr)                              */
/*  M2/IMAC3 - ANIMATION & MOTEURS PHYSIQUES - Octobre 2010      */
/*                                                               */
/*  Oscillateur mécanique élémentaire :                          */
/* test comparatif des méthodes Euler (Ex/Im)plicite et Leafrog  */
/* - Utilisation des paramètres réduit Z=hz/m et K=h²k/m         */
/* ============================================================= */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <g2x.h>

#define EXPLICIT 0
#define IMPLICIT 1
#define LEAPFROG 2

typedef struct
{
	double    A,B;
	double    x,r;
	double    f;
	G2Xcolor  col;
	bool      on;
	bool      OSC,DIV;
	char      name[64];
	double    y;
} method;
/* les 3 méthodes : Euler Explicite, Euler Implicite, LeapFrog */
method EE,EI,LF;

bool   SA=true; /* flag pour solution analytique */

/*uint   pixwidth=1340,pixheight=400;*/
uint   pixwidth=1800,pixheight=600;
double xmin=0.,ymin=-1.,xmax=+1.,ymax=+1.;
double yamp=1.,xamp=1.;

int    FECH=100;
double h,u; 
double m,k,z,K,Z;      
double x0,v0,x1;
double fc;

void setupmeth(method* M)
{
	double D;
	D=SQR(M->A)+4*M->B;
	M->OSC =(D<0.?true:false);
	if (M->OSC == IMPLICIT) 
	{
		M->DIV =(-M->B<=1?false:true);
		M->f   =(  M->DIV?0.:(1./(2.*PI*h))*acos(0.5*M->A/sqrt(-M->B)));	
	}	
	else 
    M->DIV=(((-M->A+M->B<=1)&&(M->A+M->B<=1))?false:true);
}

void setup_EI(void)
{
	EI.A=(2.+Z)/(1+K+Z);           
	EI.B=-1/(1+K+Z);
	setupmeth(&EI);
}

void setup_EE(void)
{
	EE.A=2.-Z;           
	EE.B=-1+Z-K; 
	setupmeth(&EE);
}

void setup_LF(void)
{
	LF.A=2.-Z-K;               
	LF.B=Z-1;
	setupmeth(&LF);
}

void drawmeth(method M)
{
	g2x_Print(xmin+0.01*xmax,ymin+M.y,M.col ,'l',M.name);	
	if (M.DIV) 
	{
		g2x_Print(xmin+0.05*xmax,ymin+M.y,M.col ,'l',"A=%lf B=%lf -> DIV",M.A,M.B);	
		return;
	}
	if (M.OSC) g2x_Print(xmin+0.05*xmax,ymin+M.y,M.col ,'l',"A=%lf B=%lf -> f=%lf",M.A,M.B,M.f);	
	else       g2x_Print(xmin+0.05*xmax,ymin+M.y,M.col ,'l',"A=%lf B=%lf         ",M.A,M.B);	
	
	double y,x,r,t;
	M.x=x0;
	M.r=x0;
	for (t=u;t<xmax;t+=u)
	{
		y =  M.x;
		M.x= M.A*y+M.B*M.r;
		M.r= y;
		x  = (M.x>0.)?MIN(M.x*yamp,0.999*ymax):MAX(M.x*yamp,0.999*ymin);
		r  = (M.r>0.)?MIN(M.r*yamp,0.999*ymax):MAX(M.r*yamp,0.999*ymin);
		g2x_Line(t-u,r,t,x,M.col,2);
	}	
}

/* SOLUTION ANALYTIQUE */
double D,Dmin,Dmax,r1,r2;
double c1,c2;
double (*sol)(double);

double sol_1(double t) { return c1*exp(r1*t) + c2*exp(r2*t); } /*! régime amorti    !*/
double sol_2(double t) { return (c1*t+c2)*exp(r1*t);         } /*! régime critique  !*/
double sol_3(double t) { return c1*exp(r1*t)*cos(r2*t+c2);   } /*! régime oscillant !*/

/*! mise à jour !*/
void setup_analytic(void)
{
	double D=(z*z-4.*k*m); /* discriminant */	
	/* régime critique */
	if (G2Xzero(D))
	{
	 	r1=r2=-z/(2.*m);
		c1=v0-r1*x0;
		c2=x0;
		sol=&sol_2;
		fc=-1.;
		return;
	} 
	/* régime amorti */
	if (D>0.)
	{
		r1=(-z-sqrt(D))/(2.*m);
		r2=(-z+sqrt(D))/(2.*m);
		c1=(v0-r2*x0)/(r1-r2);
		c2=x0-c1;
		sol=&sol_1;
		fc=0.;
		return;
	}
	/* régime oscillant */
	r1=-z/(2.*m);
	r2=sqrt(-D)/(2.*m);
	c2=atan((r1*x0-v0)/(r2*x0));
	c1=x0/cos(c2);
	sol=&sol_3;
	fc=sqrt(k/m-(z*z)/(4*m*m));
}

void draw_analytic(void)
{
	double x=x0,r=x0,t;
	for (t=u;t<xmax;t+=u)
	{
		r=x;
		x=(*sol)(t*xamp);
		g2x_Line(t-u,r,t,x,G2Xwb ,1);
	}	
}

void dessin(void)
{ 
  h = 1./FECH;
  u =(xmax-xmin)*h/xamp; /* pas graphique */

  /* paramètres réduits       */
	K = h*h*k/m;
	Z =   h*z/m;
  /* mise à jour des méthodes */
	setup_EE();
	setup_EI();
	setup_LF();
	setup_analytic();
	
	if (SA   ) draw_analytic();
	if (EE.on) drawmeth(EE);
	if (EI.on) drawmeth(EI);
	if (LF.on) drawmeth(LF);

	g2x_Print(xmin+0.07*xmax,ymax-0.1,G2Xk ,'l'," h=%lf   ",h );
	g2x_Print(xmin+0.01*xmax,ymax-0.2,G2Xk ,'l'," k=%.4lf ",k );
	g2x_Print(xmin+0.07*xmax,ymax-0.2,G2Xk ,'l'," K=%.4lf ",K );
	g2x_Print(xmin+0.01*xmax,ymax-0.3,G2Xk ,'l'," z=%.4lf ",z ); 
	g2x_Print(xmin+0.07*xmax,ymax-0.3,G2Xk ,'l'," Z=%.4lf ",Z ); 
	
	if (fc>0.) g2x_Print(xmin+0.01*xmax,ymax-0.1,G2Xyc ,'l',"fc=%lf",fc/(2.*PI));

}

void init(void)
{
	x0=x1=0.50;
	v0=0.;
	m =1.;
	k =1000.;
	z =1.;
  FECH = 1000;
	h =1./FECH;

  /* paramètres réduits       */
	K=h*h*k/m;
	Z=  h*z/m;	

	sprintf(EE.name,"Euler Expl."); EE.y=0.1; memcpy(EE.col,G2Xb  ,sizeof(G2Xcolor));
	sprintf(EI.name,"Euler Impl."); EI.y=0.2; memcpy(EI.col,G2Xg  ,sizeof(G2Xcolor));
	sprintf(LF.name,"Leapfrog   "); LF.y=0.3; memcpy(LF.col,G2Xr  ,sizeof(G2Xcolor));

	/* les scrollbars : 3 horiz. / 2 vertic. */
	g2x_CreateScrollh_i("Fe"  ,&FECH, 1   ,1000  , 1.0,"sampling rate  ");	
	g2x_CreateScrollh_d("k"   ,&k   , 0.  , 1.e+4, 1.0,"spring stiffness");
	g2x_CreateScrollh_d("z"   ,&z   , 0.  , 1.e+1, 1.0,"damper viscosity");
	g2x_CreateScrollv_d("xmax",&xamp, 0.05, 50.  , 1.0,"xmax");
	g2x_CreateScrollv_d("yamp",&yamp, 0.01, 1.0  , 1.0,"yamp");

  /* les boutons switch */
	g2x_CreateSwitch("SA",&SA   ,"Sol. Analytique");
	g2x_CreateSwitch("EE",&EE.on,"Euler Explicite");
  g2x_CreateSwitch("EI",&EI.on,"Euler Implicite");
  g2x_CreateSwitch("LF",&LF.on,"Leapfrog       ");

}

/*=============================*/
/*=                           =*/
/*=============================*/
int main(int argc, char* argv[])
{ 
	/* Paramètres de la fenetre graphique */
  g2x_InitWindow(*argv,pixwidth,pixheight);
  g2x_SetWindowCoord(xmin,ymin,xmax,ymax);

	g2x_SetInitFunction(init);
	g2x_SetDrawFunction(dessin);

  return g2x_MainStart();
}
