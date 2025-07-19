//CONSISTENT TRANSPORT CLSVOF ALGORITHM USING PLIC SCHEME FOR VOLUME FRACTIONS AND ENO2 SCHEME FOR LEVEL SETS.
//BOUNDARY CONDITIONS ARE AS FOLLOWS.
//LEFT, RIGHT, AND BOTTOM BOUNDARY - NO SLIP AND NO PENETRATION
//TOP BOUNDARY - OUTFLOW
class CLSVOF:public MBASE
{
	int **tag;	//tags of interfacial cells
	double **Ft;	//intermediate volume fractions
	double **Phit;	//intermediate level set functions
	double **rhot;	//intermediate density field
	double **A_xt,**A_yt;	//intermediate advection terms for X and Y momentum equations
	double mass_act;	//actual mass
	void updt_ghost(double **Phi);	//update the ghost cells of cell centered field
	double ENO2(int flag,int i,int j,double **Phi,double V);	//ENO2 scheme
	double UP1(int flag,int i,int j,double **Phi,double V);	//1st order upwind scheme
	double QUICK(int flag,int i,int j,double **Phi,double V);	//QUICK scheme
	void tag_X(double **Fa);	//tagging algorithm in X direction
	void tag_Y(double **Fa);	//tagging algorithm in Y direction
	void recon_int(int i,int j,double F,double **Phi,VECTOR *N,double *s);	//reconstruct the interface
	double vol_frac_flux(int t,int flag,int i,int j,double F,double **Phi,double V);	//calculate the advection volume fraction flux
	void adv_X(int t);	//solve the advection equation(implicit discretization in X direction)
	void adv_Y(int t);	//solve the advection equation(implicit discretization in Y direction)
	void reinit(int t);	//LS reinitialization algorithm
	public:
			CLSVOF(); ~CLSVOF();
			void ini(double xc,double yc);
			void den_ini();	//initialize density field at n+1 time step
			void solve(int n);	//CLSVOF advection algorithm
			void prop_updt();	//update density and viscosity field for next time step
			void mass_err();	//calculate mass error
			void lsvf_write(int t);	//tecplot file output
			void ls_complete(int t);	//ls file output including the ghost cells
};
CLSVOF::CLSVOF()
{
	tag=new int*[J+2];
	Phit=new double*[J+2];
	Ft=new double*[J+2];
	rhot=new double*[J+1];
	A_xt=new double*[J+1];
	A_yt=new double*[J+1];
	for(int i=0;i<J+2;i++)
	{
		tag[i]=new int[I+2];
		Phit[i]=new double[I+2];
		Ft[i]=new double[I+2];
		if(i<J+1)
		{
			rhot[i]=new double[I+1];
			A_xt[i]=new double[I+1];
			A_yt[i]=new double[I+1];
		}
	}
	cout<<"CLSVOF: MEMORY ALLOCATED"<<endl;
}
CLSVOF::~CLSVOF()
{
	for(int i=0;i<J+2;i++)
	{
		delete[] tag[i];
		delete[] Phit[i];
		delete[] Ft[i];
		if(i<J+1)
		{
			delete[] rhot[i];
			delete[] A_xt[i];
			delete[] A_yt[i];
		}
	}
	delete[] tag;
	delete[] Phit;
	delete[] Ft;
	delete[] rhot;
	delete[] A_xt;
	delete[] A_yt;
	cout<<"CLSVOF: MEMORY RELEASED"<<endl;
}
void CLSVOF::ini(double xc,double yc)
{
	INI ms(Xm,Ym,CX,CY,F,Phi,xc,yc);
	ms.VF();	//initial exact volume fractions and level sets are calculated here
	updt_ghost(F);	//update ghost cells of F
	ms.LS(tag);
	mass_act=0.0;
	for(int j=1;j<=J;j++)
		for(int i=1;i<=I;i++)
			mass_act+=F[j][i];
	for(int j=0;j<=J+1;j++)	//initialize density, viscosity, and advection field (including ghost nodes)
	{
		for(int i=0;i<=I+1;i++)
		{
			rho_n[j][i]=rho_1*F[j][i]+rho_0*(1.0-F[j][i]);
			mu[j][i]=mu_1*F[j][i]+mu_0*(1.0-F[j][i]);
			if((i>=1)&&(i<=I)&&(j>=1)&&(j<=J))	//only inner domain
			{
				A_x[j][i]=u[j][i];
				A_y[j][i]=v[j][i];
			}
		}
	}
}
void CLSVOF::den_ini()
{
	for(int j=0;j<=J+1;j++)	//calculate density field at n+1 step (including ghost nodes)
		for(int i=0;i<=I+1;i++)
			rho_np1[j][i]=rho_1*F[j][i]+rho_0*(1.0-F[j][i]);
}
void CLSVOF::prop_updt()
{
	for(int j=0;j<=J+1;j++)	//including ghost nodes
	{
		for(int i=0;i<=I+1;i++)
		{
			rho_n[j][i]=rho_np1[j][i];
			mu[j][i]=mu_1*F[j][i]+mu_0*(1.0-F[j][i]);
			if((i>=1)&&(i<=I)&&(j>=1)&&(j<=J))	//only inner domain
			{
				A_x[j][i]=u[j][i];
				A_y[j][i]=v[j][i];
			}
		}
	}
}
void CLSVOF::updt_ghost(double **Phia)
{
	for(int j=1;j<=J;j++)	//left and right ghost nodes (Neumann bc)
	{
		Phia[j][0]=Phia[j][1];
		Phia[j][I+1]=Phia[j][I];
	}
	for(int i=1;i<=I;i++)	//bottom and top ghost nodes (Neumann bc)
	{
		Phia[0][i]=Phia[1][i];
		Phia[J+1][i]=Phia[J][i];
	}
}
double CLSVOF::ENO2(int flag,int i,int j,double **Phia,double V)
{
	if(abs(V)<=EPS) return 0.0;
	double plus,minus;	//plus and minus flux
	if(flag==0)	//X direction flux
	{
		plus=Phia[j][i+1]-0.5*MINMOD((Phia[j][i+1]-Phia[j][i]),(Phia[j][i+2]-Phia[j][i+1]));
		minus=Phia[j][i]+0.5*MINMOD((Phia[j][i]-Phia[j][i-1]),(Phia[j][i+1]-Phia[j][i]));
	}
	else if(flag==1)	//Y direction flux
	{
		plus=Phia[j+1][i]-0.5*MINMOD((Phia[j+1][i]-Phia[j][i]),(Phia[j+2][i]-Phia[j+1][i]));
		minus=Phia[j][i]+0.5*MINMOD((Phia[j][i]-Phia[j-1][i]),(Phia[j+1][i]-Phia[j][i]));
	}
	if(V>0.0) return minus;
	else return plus;
}
double CLSVOF::UP1(int flag,int i,int j,double **Phia,double V)
{
	if(abs(V)<=EPS) return 0.0;
	if(flag==0)	//X direction flux
	{
		if(V>0.0) return Phia[j][i];
		else return Phia[j][i+1];
	}
	else if(flag==1)	//Y direction flux
	{
		if(V>0.0) return Phia[j][i];
		else return Phia[j+1][i];
	}
	else { cout<<"CLSVOF: ERROR IN UPWIND SCHEME!"<<endl; return 0; }
}
double CLSVOF::QUICK(int flag,int i,int j,double **Phia,double V)
{
	if(abs(V)<=EPS) return 0.0;
	if(flag==0)	//X direction flux
	{
		if(V>0.0) return (0.75*Phia[j][i]+0.375*Phia[j][i+1]-0.125*Phia[j][i-1]);
		else return (0.75*Phia[j][i+1]+0.375*Phia[j][i]-0.125*Phia[j][i+2]);
	}
	else if(flag==1)	//Y direction flux
	{
		if(V>0.0) return (0.75*Phia[j][i]+0.375*Phia[j+1][i]-0.125*Phia[j-1][i]);
		else return (0.75*Phia[j+1][i]+0.375*Phia[j][i]-0.125*Phia[j+2][i]);
	}
	else { cout<<"CLSVOF: ERROR IN QUICK SCHEME!"<<endl; return 0; }
}
void CLSVOF::tag_X(double **Fa)
{
	for(int j=0;j<=J+1;j++)	//reinitialize the cell tags (including ghost cells)
		for(int i=0;i<=I+1;i++)
			tag[j][i]=0;
	for(int j=1;j<=J;j++)	//inner domain
	{
		for(int i=1;i<=I;i++)
		{
			if((Fa[j][i]>TRUNC_l)&&(Fa[j][i]<(1.0-TRUNC_u)))	//locate interfacial cell
			{
				tag[j][i]=1;	//tag interfacial cell
				tag[j][i+1]=1; tag[j][i-1]=1;	//tag neighbouring cells in X direction
			}
		}
	}
}
void CLSVOF::tag_Y(double **Fa)
{
	for(int j=0;j<=J+1;j++)	//reinitialize the cell tags (including ghost cells)
		for(int i=0;i<=I+1;i++)
			tag[j][i]=0;
	for(int j=1;j<=J;j++)	//inner domain
	{
		for(int i=1;i<=I;i++)
		{
			if((Fa[j][i]>TRUNC_l)&&(Fa[j][i]<(1.0-TRUNC_u)))	//locate interfacial cell
			{
				tag[j][i]=1;	//tag interfacial cell
				tag[j+1][i]=1; tag[j-1][i]=1;	//tag neighbouring cells in Y direction
			}
		}
	}
}
void CLSVOF::reinit(int t)
{
	double h=MAX2(dx,dy);
	double temp,a,b,gamma=0.5*h;	//gamma is the distance parameter
//--------------------INITIALIZATION SCHEME (BASED ON ADVECTED VF FIELD)--------------------------------------------------
	for(int j=0;j<=J+1;j++)	//reinitialize the tag values and calculate LS from volume fractions (including ghost nodes)
	{
		for(int i=0;i<=I+1;i++)
		{
			tag[j][i]=0;
			Phi[j][i]=(2.0*F[j][i]-1.0)*gamma;
		}
	}
	for(int j=1;j<=J;j++)	//tagging in horizontal direction (inner domain)
	{
		for(int i=1;i<=I;i++)
		{
			if(SGN(Phi[j][i]*Phi[j][i+1])!=1)	//LS changes sign
			{
				if(abs(Phi[j][i])<=abs(Phi[j][i+1])) tag[j][i]=1;	//tag cell having minimum LS value
				else if(i!=I) tag[j][i+1]=1;
			}
		}
	}
	for(int i=1;i<=I;i++)	//tagging in vertical direction (inner domain)
	{
		for(int j=1;j<=J;j++)
		{
			if(SGN(Phi[j][i]*Phi[j+1][i])!=1)	//LS changes sign
			{
				if(abs(Phi[j][i])<=abs(Phi[j+1][i])) tag[j][i]=1;	//tag cell having minimum LS value
				else if(j!=J) tag[j+1][i]=1;
			}
		}
	}
	for(int j=0;j<=J+1;j++)	//reset LS in untagged cells (including ghost cells)
		for(int i=0;i<=I+1;i++)
			if(tag[j][i]!=1) Phi[j][i]=SGN(Phi[j][i])*100.0;
	for(int j=1;j<=J;j++)	//initialize level sets of the ghost cells(left and right boundaries)
	{
		Phi[j][0]=Phi[j][1]; tag[j][0]=tag[j][1];
		Phi[j][I+1]=Phi[j][I]; tag[j][I+1]=tag[j][I];
	}
	for(int i=0;i<=I+1;i++)	//initialize level sets of the ghost cells(bottom and top boundaries)
	{
		Phi[0][i]=Phi[1][i]; tag[0][i]=tag[1][i];
		Phi[J+1][i]=Phi[J][i]; tag[J+1][i]=tag[J][i];
	}
//---------------SOLUTION OF DISCRETE EQUATIONS(including the ghost cells)-----------------------------
	for(int sweep=1,i_ini,j_ini,di,dj;sweep<=4;sweep++)	//Gauss-Siedel sweeps
	{
		switch(sweep)	//direction of each Gauss-Siedel sweep
		{
			case 1: j_ini=0; i_ini=0;
					dj=1; di=1;
					break;
			case 2: j_ini=0; i_ini=I+1;
					dj=1; di=-1;
					break;
			case 3: j_ini=J+1; i_ini=I+1;
					dj=-1; di=-1;
					break;
			case 4: j_ini=J+1; i_ini=0;
					dj=-1; di=1;
					break;
			default: break;
		}
		for(int j=j_ini;((j>=0)&&(j<=J+1));j+=dj)	//sweep the domain in the required direction (SMART LOOPS!)
		{
			for(int i=i_ini;((i>=0)&&(i<=I+1));i+=di)
			{
				if(tag[j][i]==1) continue;	//interface cells are not updated
				if(i==0) a=Phi[j][i+1];	//left boundary
				else if(i==(I+1)) a=Phi[j][i-1];	//right boundary
				else	//inner domain
				{
					if(SGN(Phi[j][i])==1.0) a=MIN2(Phi[j][i+1],Phi[j][i-1]);
					else a=MAX2(Phi[j][i+1],Phi[j][i-1]);
				}
				if(j==0) b=Phi[j+1][i];	//bottom boundary
				else if(j==(J+1)) b=Phi[j-1][i];	//top boundary
				else	//inner domain
				{
					if(SGN(Phi[j][i])==1.0) b=MIN2(Phi[j+1][i],Phi[j-1][i]);
					else b=MAX2(Phi[j+1][i],Phi[j-1][i]);
				}
				if(SGN(Phi[j][i])==1.0)	//positive viscosity solution
				{
					if((abs(a-b)-h)>=-SMALL) temp=MIN2(a,b)+h;
					else temp=0.5*(a+b+sqrt(2.0*pow(h,2.0)-pow((a-b),2.0)));
					Phi[j][i]=MIN2(Phi[j][i],temp);
				}
				else	//negative viscosity solution
				{
					if((abs(a-b)-h)>=-SMALL) temp=MAX2(a,b)-h;
					else temp=0.5*(a+b-sqrt(2.0*pow(h,2.0)-pow((a-b),2.0)));
					Phi[j][i]=MAX2(Phi[j][i],temp);
				}
			}
		}
	}
}
void CLSVOF::recon_int(int i,int j,double Fa,double **Phia,VECTOR *N,double *s)
{
	double sx,sy,sm,s_max,s_min,s_c;	//geometric parameters
	double F_c;	//minimum volume fraction in the cell
	*N=VECTOR((0.5*(Phia[j][i+1]-Phia[j][i-1])/dx),(0.5*(Phia[j+1][i]-Phia[j-1][i])/dy));
	*N=(*N).unit();
	sx=abs((*N).x)*dx; sy=abs((*N).y)*dy; sm=sx+sy;
	F_c=MIN2(Fa,(1.0-Fa));
	s_max=MAX2(sx,sy); s_min=MIN2(sx,sy);
	if(F_c<(0.5*s_min/s_max))
		s_c=sqrt(2.0*F_c*sx*sy);
	else
		s_c=F_c*s_max+0.5*s_min;
	if(Fa>0.5)
		*s=sm-s_c;
	else
		*s=s_c;
}
double CLSVOF::vol_frac_flux(int t,int flag,int i,int j,double Fa,double **Phia,double V)
{
	if(abs(V)<=EPS) return 0.0;	//zero velocity
	else if((abs(Fa)<=TRUNC_l)||((1.0-Fa)<=TRUNC_u)) return (Fa*V);	//single phase cell
	VECTOR N; double s;	//equation of interface
	double Dx,Dy,Vol_adv;
	int flag1=0;
	double dsx,dsy,sm,ds_max,ds_min,s_c;	//geometric parameters
	double F_c,dv0;	//minimum volume fraction in the cell
	recon_int(i,j,Fa,Phia,&N,&s);	//reconstruct interface
	if(flag==0)	//volume fraction flux in X direction
	{
		if((N.x*V)<0.0) {Dx=dx-abs(V)*dt; flag1=1;}
		else Dx=abs(V)*dt;
		Dy=dy;
	}
	else	//volume fraction flux in Y direction
	{
		if((N.y*V)<0.0) {Dy=dy-abs(V)*dt; flag1=1;}
		else Dy=abs(V)*dt;
		Dx=dx;
	}
	dsx=abs(N.x)*Dx; dsy=abs(N.y)*Dy;
	sm=dsx+dsy;
	ds_max=MAX2(dsx,dsy); ds_min=MIN2(dsx,dsy);
	s_c=MIN2(s,(MAX2((sm-s),0.0)));
	if(s_c<ds_min) F_c=0.5*pow(s_c,2.0)/(dsx*dsy);
	else F_c=(s_c-0.5*ds_min)/ds_max;
	if(s>(0.5*sm)) dv0=Dx*Dy*(1.0-F_c);
	else dv0=Dx*Dy*F_c;
	if(flag1==1) Vol_adv=Fa*dx*dy-dv0;
	else Vol_adv=dv0;
	if(flag==0)	//volume fraction flux in X direction
		return (SGN(V)*Vol_adv/(dt*dy));
	else	//volume fraction flux in Y direction
		return (SGN(V)*Vol_adv/(dt*dx));
}
void CLSVOF::adv_X(int t)
{
	double F_adv[2],LS_adv[2],rho_adv[2],u_f[2],v_f[2];	//advection fluxes for volume fractions, LS, density, and advected velocity for X and Y momentum equation
	double dtdx=dt/dx,dtdy=dt/dy;	//constant ratios for the advection equation
	int up;	//donor cell index
	tag_X(F);	//tag cells for advection in X direction
	for(int j=1;j<=J;j++)	//implicit discretization in X direction
	{
		for(int i=1;i<=I;i++)
		{
			if(i==1)
			{
				F_adv[0]=0.0;	//fluid does not cross the boundary
				LS_adv[0]=UP1(0,1,j,Phi,u_EW[j][0]);	//upwind scheme for left boundary flux
				rho_adv[0]=0.0;
				u_f[0]=0.0;	//no penetration
				v_f[0]=0.0;	//no slip
			}
			if(i==I)
			{
				F_adv[1]=0.0;	//fluid does not cross the boundary
				LS_adv[1]=UP1(0,I,j,Phi,u_EW[j][I]);	//upwind scheme for right boundary flux
				rho_adv[1]=0.0;
				u_f[1]=0.0;	//no penetration
				v_f[1]=0.0;	//no slip
			}
			else	//inner domain
			{
				up=(u_EW[j][i]<0.0)?(i+1):i;
				F_adv[1]=vol_frac_flux(t,0,up,j,F[j][up],Phi,u_EW[j][i]);
				LS_adv[1]=ENO2(0,i,j,Phi,u_EW[j][i]);	//ENO2 scheme for inner domain
				rho_adv[1]=(rho_1-rho_0)*F_adv[1]+rho_0*u_EW[j][i];
				if((tag[j][i]==0)&&(tag[j][i+1]==0)&&(i>1)&&(i<I-1)&&(j>1)&&(j<J-1))	//use QUICK scheme
				{
					u_f[1]=QUICK(0,i,j,A_x,u_EW[j][i]);
					v_f[1]=QUICK(0,i,j,A_y,u_EW[j][i]);
				}
				else	//use UP1 scheme
				{
					u_f[1]=UP1(0,i,j,A_x,u_EW[j][i]);
					v_f[1]=UP1(0,i,j,A_y,u_EW[j][i]);
				}
			}
			Ft[j][i]=(F[j][i]-dtdx*(F_adv[1]-F_adv[0]))/(1.0-dtdx*(u_EW[j][i]-u_EW[j][i-1]));
			Phit[j][i]=(Phi[j][i]-dtdx*(u_EW[j][i]*LS_adv[1]-u_EW[j][i-1]*LS_adv[0]))/(1.0-dtdx*(u_EW[j][i]-u_EW[j][i-1]));
			rhot[j][i]=(rho_n[j][i]-dtdx*(rho_adv[1]-rho_adv[0]))/(1.0-dtdx*(u_EW[j][i]-u_EW[j][i-1]));
			A_xt[j][i]=(rho_n[j][i]*A_x[j][i]-dtdx*(u_f[1]*rho_adv[1]-u_f[0]*rho_adv[0]))/(1.0-dtdx*(u_EW[j][i]-u_EW[j][i-1]));
			A_yt[j][i]=(rho_n[j][i]*A_y[j][i]-dtdx*(v_f[1]*rho_adv[1]-v_f[0]*rho_adv[0]))/(1.0-dtdx*(u_EW[j][i]-u_EW[j][i-1]));
			F_adv[0]=F_adv[1];	//updation for the next cell
			LS_adv[0]=LS_adv[1];
			rho_adv[0]=rho_adv[1];
			u_f[0]=u_f[1]; v_f[0]=v_f[1];
			A_xt[j][i]/=rhot[j][i];	//extract intermediate velocity field
			A_yt[j][i]/=rhot[j][i];
			if(abs(Ft[j][i])<=TRUNC_l) { Ft[j][i]=0.0; rhot[j][i]=rho_0; }	//clipping of volume fractions and density
			else if(abs(1.0-Ft[j][i])<=TRUNC_u) { Ft[j][i]=1.0; rhot[j][i]=rho_1; }
		}
	}
	updt_ghost(Ft);	//update ghost cells of Ft
	updt_ghost(Phit);	//update ghost cells of Phit
	tag_Y(Ft);	//tag cells for advection in Y direction
	for(int i=1;i<=I;i++)	//explicit discretization in Y direction
	{
		for(int j=1;j<=J;j++)
		{
			if(j==1)
			{
				F_adv[0]=0.0;	//fluid does not cross the boundary
				LS_adv[0]=UP1(1,i,1,Phit,v_NS[0][i]);	//upwind scheme for bottom boundary flux
				rho_adv[0]=0.0;
				u_f[0]=0.0;	//no slip
				v_f[0]=0.0;	//no penetration
			}
			if(j==J)
			{
				F_adv[1]=0.0;	//fluid does not cross the boundary
				LS_adv[1]=UP1(1,i,J,Phit,v_NS[J][i]);	//upwind scheme for top boundary flux
				u_f[1]=A_xt[J][i];	//outflow
				v_f[1]=A_yt[J][i];	//outflow
				rho_adv[1]=rho_0*v_NS[J][i];	//outflow
			}
			else
			{
				up=(v_NS[j][i]<0.0)?(j+1):j;
				F_adv[1]=vol_frac_flux(t,1,i,up,Ft[up][i],Phit,v_NS[j][i]);
				LS_adv[1]=ENO2(1,i,j,Phit,v_NS[j][i]);	//ENO2 scheme for inner domain
				rho_adv[1]=(rho_1-rho_0)*F_adv[1]+rho_0*v_NS[j][i];
				if((tag[j][i]==0)&&(tag[j+1][i]==0)&&(i>1)&&(i<I-1)&&(j>1)&&(j<J-1))	//use QUICK scheme
				{
					u_f[1]=QUICK(1,i,j,A_xt,v_NS[j][i]);
					v_f[1]=QUICK(1,i,j,A_yt,v_NS[j][i]);
				}
				else	//use UP1 scheme
				{
					u_f[1]=UP1(1,i,j,A_xt,v_NS[j][i]);
					v_f[1]=UP1(1,i,j,A_yt,v_NS[j][i]);
				}
			}
			F[j][i]=Ft[j][i]*(1.0+dtdy*(v_NS[j][i]-v_NS[j-1][i]))-dtdy*(F_adv[1]-F_adv[0]);
			Phi[j][i]=Phit[j][i]*(1.0+dtdy*(v_NS[j][i]-v_NS[j-1][i]))-dtdy*(v_NS[j][i]*LS_adv[1]-v_NS[j-1][i]*LS_adv[0]);
			rho_np1[j][i]=rhot[j][i]*(1.0+dtdy*(v_NS[j][i]-v_NS[j-1][i]))-dtdy*(rho_adv[1]-rho_adv[0]);
			A_x[j][i]=rhot[j][i]*A_xt[j][i]*(1.0+dtdy*(v_NS[j][i]-v_NS[j-1][i]))-dtdy*(u_f[1]*rho_adv[1]-u_f[0]*rho_adv[0]);
			A_y[j][i]=rhot[j][i]*A_yt[j][i]*(1.0+dtdy*(v_NS[j][i]-v_NS[j-1][i]))-dtdy*(v_f[1]*rho_adv[1]-v_f[0]*rho_adv[0]);
			F_adv[0]=F_adv[1];	//updation for the next cell
			LS_adv[0]=LS_adv[1];
			rho_adv[0]=rho_adv[1];
			u_f[0]=u_f[1]; v_f[0]=v_f[1];
			A_x[j][i]/=rho_np1[j][i];	//extract velocity field
			A_y[j][i]/=rho_np1[j][i];
			if(abs(F[j][i])<=TRUNC_l) { F[j][i]=0.0; rho_np1[j][i]=rho_0; }	//clipping of volume fractions and density
			else if(abs(1.0-F[j][i])<=TRUNC_u) { F[j][i]=1.0; rho_np1[j][i]=rho_1; }
		}
	}
	updt_ghost(F);	//update ghost cells of F
	updt_ghost(rho_np1);	//update ghost cells of rho
}
void CLSVOF::adv_Y(int t)
{
	double F_adv[2],LS_adv[2],rho_adv[2],u_f[2],v_f[2];	//advection fluxes for volume fractions, LS, density, and advected velocity for X and Y momentum equation
	double dtdx=dt/dx,dtdy=dt/dy;	//constant ratios for the advection equation
	int up;	//donor cell index
	tag_Y(F);	//tag cells for advection in Y direction
	for(int i=1;i<=I;i++)	//implicit discretization in Y direction
	{
		for(int j=1;j<=J;j++)
		{
			if(j==1)
			{
				F_adv[0]=0.0;	//fluid does not cross the boundary
				LS_adv[0]=UP1(1,i,1,Phi,v_NS[0][i]);	//upwind scheme for bottom boundary flux
				rho_adv[0]=0.0;
				u_f[0]=0.0;	//no slip
				v_f[0]=0.0;	//no penetration
			}
			if(j==J)
			{
				F_adv[1]=0.0;	//fluid does not cross the boundary
				LS_adv[1]=UP1(1,i,J,Phi,v_NS[J][i]);	//upwind scheme for top boundary flux
				u_f[1]=A_x[J][i];	//outflow
				v_f[1]=A_y[J][i];	//outflow
				rho_adv[1]=rho_0*v_NS[J][i];	//outflow
			}
			else
			{
				up=(v_NS[j][i]<0.0)?(j+1):j;
				F_adv[1]=vol_frac_flux(t,1,i,up,F[up][i],Phi,v_NS[j][i]);
				LS_adv[1]=ENO2(1,i,j,Phi,v_NS[j][i]);	//ENO2 scheme for inner domain
				rho_adv[1]=(rho_1-rho_0)*F_adv[1]+rho_0*v_NS[j][i];
				if((tag[j][i]==0)&&(tag[j-1][i]==0)&&(i>2)&&(i<I-2)&&(j>2)&&(j<J-2))	//use QUICK scheme
				{
					u_f[1]=QUICK(1,i,j,A_x,v_NS[j][i]);
					v_f[1]=QUICK(1,i,j,A_y,v_NS[j][i]);
				}
				else	//use UP1 scheme
				{
					u_f[1]=UP1(1,i,j,A_x,v_NS[j][i]);
					v_f[1]=UP1(1,i,j,A_y,v_NS[j][i]);
				}
			}
			Ft[j][i]=(F[j][i]-dtdy*(F_adv[1]-F_adv[0]))/(1.0-dtdy*(v_NS[j][i]-v_NS[j-1][i]));
			Phit[j][i]=(Phi[j][i]-dtdy*(v_NS[j][i]*LS_adv[1]-v_NS[j-1][i]*LS_adv[0]))/(1.0-dtdy*(v_NS[j][i]-v_NS[j-1][i]));
			rhot[j][i]=(rho_n[j][i]-dtdy*(rho_adv[1]-rho_adv[0]))/(1.0-dtdy*(v_NS[j][i]-v_NS[j-1][i]));
			A_xt[j][i]=(rho_n[j][i]*A_x[j][i]-dtdy*(u_f[1]*rho_adv[1]-u_f[0]*rho_adv[0]))/(1.0-dtdy*(v_NS[j][i]-v_NS[j-1][i]));
			A_yt[j][i]=(rho_n[j][i]*A_y[j][i]-dtdy*(v_f[1]*rho_adv[1]-v_f[0]*rho_adv[0]))/(1.0-dtdy*(v_NS[j][i]-v_NS[j-1][i]));
			F_adv[0]=F_adv[1];	//updation for the next cell
			LS_adv[0]=LS_adv[1];
			rho_adv[0]=rho_adv[1];
			u_f[0]=u_f[1]; v_f[0]=v_f[1];
			A_xt[j][i]/=rhot[j][i];	//extract intermediate velocity field
			A_yt[j][i]/=rhot[j][i];
			if(abs(Ft[j][i])<=TRUNC_l) { Ft[j][i]=0.0; rhot[j][i]=rho_0; }	//clipping of volume fractions and density
			else if(abs(1.0-Ft[j][i])<=TRUNC_u) { Ft[j][i]=1.0; rhot[j][i]=rho_1; }
		}
	}
	updt_ghost(Ft);	//update ghost cells of Ft
	updt_ghost(Phit);	//update ghost cells of Phit
	tag_X(Ft);	//tag cells for advection in X direction
	for(int j=1;j<=J;j++)	//explicit discretization in X direction
	{
		for(int i=1;i<=I;i++)
		{
			if(i==1)
			{
				F_adv[0]=0.0;	//fluid does not cross the boundary
				LS_adv[0]=UP1(0,1,j,Phit,u_EW[j][0]);	//upwind scheme for left boundary flux
				rho_adv[0]=0.0;
				u_f[0]=0.0;	//no penetration
				v_f[0]=0.0;	//no slip
			}
			if(i==I)
			{
				F_adv[1]=0.0;	//fluid does not cross the boundary
				LS_adv[1]=UP1(0,I,j,Phit,u_EW[j][I]);	//upwind scheme for right boundary flux
				rho_adv[1]=0.0;
				u_f[1]=0.0;	//no penetration
				v_f[1]=0.0;	//no slip
			}
			else	//inner domain
			{
				up=(u_EW[j][i]<0.0)?(i+1):i;
				F_adv[1]=vol_frac_flux(t,0,up,j,Ft[j][up],Phit,u_EW[j][i]);
				LS_adv[1]=ENO2(0,i,j,Phit,u_EW[j][i]);	//ENO2 scheme for inner domain
				rho_adv[1]=(rho_1-rho_0)*F_adv[1]+rho_0*u_EW[j][i];
				if((tag[j][i]==0)&&(tag[j][i-1]==0)&&(i>2)&&(i<I-2)&&(j>2)&&(j<J-2))	//use QUICK scheme
				{
					u_f[1]=QUICK(0,i,j,A_xt,u_EW[j][i]);
					v_f[1]=QUICK(0,i,j,A_yt,u_EW[j][i]);
				}
				else	//use UP1 scheme
				{
					u_f[1]=UP1(0,i,j,A_xt,u_EW[j][i]);
					v_f[1]=UP1(0,i,j,A_yt,u_EW[j][i]);
				}
			}
			F[j][i]=Ft[j][i]*(1.0+dtdx*(u_EW[j][i]-u_EW[j][i-1]))-dtdx*(F_adv[1]-F_adv[0]);
			Phi[j][i]=Phit[j][i]*(1.0+dtdx*(u_EW[j][i]-u_EW[j][i-1]))-dtdx*(u_EW[j][i]*LS_adv[1]-u_EW[j][i-1]*LS_adv[0]);
			rho_np1[j][i]=rhot[j][i]*(1.0+dtdx*(u_EW[j][i]-u_EW[j][i-1]))-dtdx*(rho_adv[1]-rho_adv[0]);
			A_x[j][i]=rhot[j][i]*A_xt[j][i]*(1.0+dtdx*(u_EW[j][i]-u_EW[j][i-1]))-dtdx*(u_f[1]*rho_adv[1]-u_f[0]*rho_adv[0]);
			A_y[j][i]=rhot[j][i]*A_yt[j][i]*(1.0+dtdx*(u_EW[j][i]-u_EW[j][i-1]))-dtdx*(v_f[1]*rho_adv[1]-v_f[0]*rho_adv[0]);
			F_adv[0]=F_adv[1];	//updation for the next cell
			LS_adv[0]=LS_adv[1];
			rho_adv[0]=rho_adv[1];
			u_f[0]=u_f[1]; v_f[0]=v_f[1];
			A_x[j][i]/=rho_np1[j][i];	//extract velocity field
			A_y[j][i]/=rho_np1[j][i];
			if(abs(F[j][i])<=TRUNC_l) { F[j][i]=0.0; rho_np1[j][i]=rho_0; }	//clipping of volume fractions and density
			else if(abs(1.0-F[j][i])<=TRUNC_u) { F[j][i]=1.0; rho_np1[j][i]=rho_1; }
		}
	}
	updt_ghost(F);	//update ghost cells of F
	updt_ghost(rho_np1);	//update ghost cells of rho
}
void CLSVOF::solve(int n)
{
	if((n%2)==0) adv_Y(n);	//Strang splitting
	else adv_X(n);
	reinit(n);
}
void CLSVOF::lsvf_write(int t)
{
	string fname="ls_vol_frac_"+to_string(t)+".dat";
	ofstream p_out(fname);
	p_out<<"TITLE = \"LEVEL SETS AND VOLUME FRACTIONS\""<<endl;
	p_out<<"FILETYPE = SOLUTION"<<endl;
	p_out<<"VARIABLES = \"F\",\"Phi\""<<endl;
	p_out<<"ZONE T=\""<<t*dt<<"\", I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK, VARLOCATION=([1,2]=CELLCENTERED), SOLUTIONTIME="<<t*dt<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<F[j][i];
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<Phi[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"CLSVOF: LEVEL SETS AND VOLUME FRACTIONS FILE OUTPUT SUCCESSFUL AT n = "<<t<<endl;
}
void CLSVOF::mass_err()
{
	double mass=0.0;
	for(int j=1;j<=J;j++)
		for(int i=1;i<=I;i++)
			mass+=F[j][i];
	cout<<"CLSVOF: MASS ERROR = "<<(mass-mass_act)<<endl;
}
void CLSVOF::ls_complete(int t)
{
	string fname="ls_comlpete_"+to_string(t)+".dat";
	ofstream p_out(fname);
	p_out<<"TITLE = \"LEVEL SETS INCLUDING GHOST CELLS\""<<endl;
	p_out<<"FILETYPE = FULL"<<endl;
	p_out<<"VARIABLES = \"X\",\"Y\",\"Phi\""<<endl;
	p_out<<"ZONE T=\""<<t*dt<<"\", I="<<I+2<<", J="<<J+2<<", DATAPACKING=BLOCK, SOLUTIONTIME="<<t*dt<<endl;
	double x=0.0-0.5*dx,y=0.0-0.5*dy;
	for(int j=0;j<=J+1;j++)	//X coordinates
	{
		for(int i=0;i<=I+1;i++)
		{
			p_out<<" "<<x;
			x+=dx;
		}
		p_out<<endl;
		x=0.0-0.5*dx;
	}
	p_out<<endl;
	for(int j=0;j<=J+1;j++)	//Y coordinates
	{
		for(int i=0;i<=I+1;i++)
			p_out<<" "<<y;
		p_out<<endl;
		y+=dy;
	}
	p_out<<endl;
	for(int j=0;j<=J+1;j++)
	{
		for(int i=0;i<=I+1;i++)
			p_out<<" "<<Phi[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"CLSVOF: COMPLETE LEVEL SET FILE OUTPUT SUCCESSFUL AT n = "<<t<<endl;
}
