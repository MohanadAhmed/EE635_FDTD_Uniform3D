#include <stdio.h>
#include <algorithm>
#include "Uniform3DGrid.h"
#define sqr(x) pow(x, 2.0)

#define LOOP_TRIPLE(is, ie, js, je, ks,ke) \
	for(int i = is; i < ie; i++){	\
		for(int j = js; j < je; j++) {\
			for(int k = ks; k < ke; k++){
			
#define LOOP_CLOSE()	}}}

#define IDX(i,j,k) ((i)*(jb*kb) + (j)*(kb) + (k))
#define IDXP(i,j,k, p) (IDX(i,j,k)*npoles + p)
#define FIELDS_NUM 6

void DumpSigmaArray(std::vector<Uniform3DGrid::Value>& arr, const char* filename);

void Uniform3DGrid::Allocate()
{
	ie = pml_back + round(xx/ddx) + pml_front; ic = pml_back + round(xx/(2*ddx));
	je = pml_left + round(yy/ddy) + pml_right; jc = pml_left + round(yy/(2*ddy));
	ke = pml_bottom + round(zz/ddz) + pml_top; kc = pml_bottom + round(zz/(2*ddz));
	
	ib = ie + 1; jb = je + 1; kb = ke + 1;
	
	Value ds = sqr(1/ddx) + sqr(1/ddy) + sqr(1/ddz);
	dt = dt_factor / (c0 * sqrt(ds));
	
	printf("Center = %d , %d, %d\n", ic, jc, kc);
	
	PrintSimData(nullptr);
	const int array_size = ib*jb*kb;
	printf("Array Size = %d\n", array_size);
	
	{	/***** Allocate for Field Arrays ***********/
		Ex.resize(array_size, 0);
		Ey.resize(array_size, 0);
		Ez.resize(array_size, 0);
		
		Dx.resize(array_size, 0);
		Dy.resize(array_size, 0);
		Dz.resize(array_size, 0);
		
		Hx.resize(array_size, 0);
		Hy.resize(array_size, 0);
		Hz.resize(array_size, 0);
		
		Bx.resize(array_size, 0);
		By.resize(array_size, 0);
		Bz.resize(array_size, 0);
	}
	
	{	/***** Allocate for Field Arrays ***********/
		Dxy.resize(array_size, 0);
		Dxz.resize(array_size, 0);
		Dyx.resize(array_size, 0);
		Dyz.resize(array_size, 0);
		Dzx.resize(array_size, 0);
		Dzy.resize(array_size, 0);
		
		Bxy.resize(array_size, 0);
		Bxz.resize(array_size, 0);
		Byx.resize(array_size, 0);
		Byz.resize(array_size, 0);
		Bzx.resize(array_size, 0);
		Bzy.resize(array_size, 0);
	}
	
	{	/***** Allocate for Coefficients ***********/
		caDxy.resize(array_size, 1);
		caDxz.resize(array_size, 1);
		caDyx.resize(array_size, 1);
		caDyz.resize(array_size, 1);
		caDzx.resize(array_size, 1);
		caDzy.resize(array_size, 1);
		
		cbDxy.resize(array_size, (dt/ddy));
		cbDxz.resize(array_size, (dt/ddz));
		cbDyx.resize(array_size, (dt/ddx));
		cbDyz.resize(array_size, (dt/ddz));
		cbDzx.resize(array_size, (dt/ddx));
		cbDzy.resize(array_size, (dt/ddy));
		
		daBxy.resize(array_size, 1);
		daBxz.resize(array_size, 1);
		daByx.resize(array_size, 1);
		daByz.resize(array_size, 1);
		daBzx.resize(array_size, 1);
		daBzy.resize(array_size, 1);
		
		dbBxy.resize(array_size, (dt/ddy));
		dbBxz.resize(array_size, (dt/ddz));
		dbByx.resize(array_size, (dt/ddx));
		dbByz.resize(array_size, (dt/ddz));
		dbBzx.resize(array_size, (dt/ddx));
		dbBzy.resize(array_size, (dt/ddy));
	}
	
	{	/***** Initialize PML Regions of Coefficients ***********/
		const double pmln = 4;
		//const Value R = 1e-6;
		Value SigmaMax = 3.0/dt;
		
		std::vector<Value> s_Dx_back(pml_back,0);
		std::vector<Value> s_Dx_front(pml_front,0);
		std::vector<Value> s_Dy_left(pml_left,0);
		std::vector<Value> s_Dy_right(pml_right,0);
		std::vector<Value> s_Dz_bottom(pml_bottom,0);
		std::vector<Value> s_Dz_top(pml_top,0);
		                   
		std::vector<Value> s_Bx_back(pml_back,0);
		std::vector<Value> s_Bx_front(pml_front,0);
		std::vector<Value> s_By_left(pml_left,0);
		std::vector<Value> s_By_right(pml_right,0);
		std::vector<Value> s_Bz_bottom(pml_bottom,0);
		std::vector<Value> s_Bz_top(pml_top,0);
		
		
		for(int i = 0; i < pml_back;i++){
			// Value SigmaMax = - ((pmln + 1) * c0 * log(R))/(2*pml_back*ddx);
			s_Dx_back[i] = SigmaMax * pow( (1.0*i+1.0)/pml_back, pmln);
			s_Bx_back[i] = SigmaMax * pow( (1.0*i+0.5)/pml_back, pmln);
		}
		
		for(int i = 0; i < pml_front;i++){
			// Value SigmaMax = - ((pmln + 1) * c0 * log(R))/(2*pml_front*ddx);
			s_Dx_front[i] = SigmaMax * pow( (1.0*i)/pml_front, pmln);
			s_Bx_front[i] = SigmaMax * pow( (1.0*i+0.5)/pml_front, pmln);
		}
		for(int i = 0; i < pml_left;i++){
			// Value SigmaMax = - ((pmln + 1) * c0 * log(R))/(2*pml_left*ddy);
			s_Dy_left[i] = SigmaMax * pow( (1.0*i+1.0)/pml_left, pmln);
			s_By_left[i] = SigmaMax * pow( (1.0*i+0.5)/pml_left, pmln);
		}
		for(int i = 0; i < pml_right;i++){
			// Value SigmaMax = - ((pmln + 1) * c0 * log(R))/(2*pml_right*ddy);
			s_Dy_right[i] = SigmaMax * pow( (1.0*i)/pml_right, pmln);
			s_By_right[i] = SigmaMax * pow( (1.0*i+0.5)/pml_right, pmln);
		}
		for(int i = 0; i < pml_bottom;i++){
			// Value SigmaMax = - ((pmln + 1) * c0 * log(R))/(2*pml_bottom*ddz);
			s_Dz_bottom[i] = SigmaMax * pow( (1.0*i+1.0)/pml_bottom, pmln);
			s_Bz_bottom[i] = SigmaMax * pow( (1.0*i+0.5)/pml_bottom, pmln);
		}
		for(int i = 0; i < pml_top;i++){
			// Value SigmaMax = - ((pmln + 1) * c0 * log(R))/(2*pml_top*ddz);
			s_Dz_top[i] = SigmaMax * pow( (1.0*i)/pml_top, pmln);
			s_Bz_top[i] = SigmaMax * pow( (1.0*i+0.5)/pml_top, pmln);
		}
		
		DumpSigmaArray(s_Dx_back, "params/sigma_Dx_back");
		DumpSigmaArray(s_Dx_front, "params/sigma_Dx_front");
		DumpSigmaArray(s_Dy_left, "params/sigma_Dy_left");
		DumpSigmaArray(s_Dy_right, "params/sigma_Dy_right");
		DumpSigmaArray(s_Dz_bottom, "params/sigma_Dz_bottom");
		DumpSigmaArray(s_Dz_top, "params/sigma_Dz_top");
		
		DumpSigmaArray(s_Bx_back, "params/sigma_Bx_back");
		DumpSigmaArray(s_Bx_front, "params/sigma_Bx_front");
		DumpSigmaArray(s_By_left, "params/sigma_By_left");
		DumpSigmaArray(s_By_right, "params/sigma_By_right");
		DumpSigmaArray(s_Bz_bottom, "params/sigma_Bz_bottom");
		DumpSigmaArray(s_Bz_top, "params/sigma_Bz_top");
		
		/*Back PML:X*/
		for(int i = pml_back; i > 0; i--){
			for(int  j = 0; j < je;j++){
				for(int k = 0; k < ke;k++){
					int z = pml_back - i;
					caDyx[IDX(i,j,k)] = (2 - dt*s_Dx_back[z])/(2 + dt*s_Dx_back[z]);
					cbDyx[IDX(i,j,k)] = (dt/ddx)/(1 + dt*s_Dx_back[z]*0.5);
					daBzx[IDX(i,j,k)] = (2 - dt*s_Bx_back[z])/(2 + dt*s_Bx_back[z]);
					dbBzx[IDX(i,j,k)] = (dt/ddx)/(1 + dt*s_Bx_back[z]*0.5);
					
					caDzx[IDX(i,j,k)] = (2 - dt*s_Dx_back[z])/(2 + dt*s_Dx_back[z]);
					cbDzx[IDX(i,j,k)] = (dt/ddx)/(1 + dt*s_Dx_back[z]*0.5);
					daByx[IDX(i,j,k)] = (2 - dt*s_Bx_back[z])/(2 + dt*s_Bx_back[z]);
					dbByx[IDX(i,j,k)] = (dt/ddx)/(1 + dt*s_Bx_back[z]*0.5);
				}
			}
		}
		printf("BackDone\n");
		/*Front PML:X*/
		for(int i = ib - pml_front; i < ib; i++){
			for(int  j = 0; j < je;j++){
				for(int k = 0; k < ke;k++){
					int z = i - (ib - pml_front);
					caDyx[IDX(i,j,k)] = (2 - dt*s_Dx_front[z])/(2 + dt*s_Dx_front[z]);
					cbDyx[IDX(i,j,k)] = (dt/ddx)/(1 + dt*s_Dx_front[z]*0.5);;
					daBzx[IDX(i,j,k)] = (2 - dt*s_Bx_front[z])/(2 + dt*s_Bx_front[z]);
					dbBzx[IDX(i,j,k)] = (dt/ddx)/(1 + dt*s_Bx_front[z]*0.5);
					
					caDzx[IDX(i,j,k)] = (2 - dt*s_Dx_front[z])/(2 + dt*s_Dx_front[z]);
					cbDzx[IDX(i,j,k)] = (dt/ddx)/(1 + dt*s_Dx_front[z]*0.5);
					daByx[IDX(i,j,k)] = (2 - dt*s_Bx_front[z])/(2 + dt*s_Bx_front[z]);
					dbByx[IDX(i,j,k)] = (dt/ddx)/(1 + dt*s_Bx_front[z]*0.5);
				}
			}
		}
		printf("FrontDone\n");
		/*Left PML:Y*/
		for(int i = 0; i < ie; i++){
			for(int  j = pml_left; j > 0;j--){
				for(int k = 0; k < ke;k++){
					int z = pml_left - j;
					caDxy[IDX(i,j,k)] = (2 - dt*s_Dy_left[z])/(2 + dt*s_Dy_left[z]);
					cbDxy[IDX(i,j,k)] = (dt/ddx)/(1 + dt*s_Dy_left[z]*0.5);
					daBzy[IDX(i,j,k)] = (2 - dt*s_By_left[z])/(2 + dt*s_By_left[z]);
					dbBzy[IDX(i,j,k)] = (dt/ddx)/(1 + dt*s_By_left[z]*0.5);
					    
					caDzy[IDX(i,j,k)] = (2 - dt*s_Dy_left[z])/(2 + dt*s_Dy_left[z]);
					cbDzy[IDX(i,j,k)] = (dt/ddx)/(1 + dt*s_Dy_left[z]*0.5);
					daBxy[IDX(i,j,k)] = (2 - dt*s_By_left[z])/(2 + dt*s_By_left[z]);
					dbBxy[IDX(i,j,k)] = (dt/ddx)/(1 + dt*s_By_left[z]*0.5);
				}
			}
		}
		printf("LeftDone\n");
		/*Right PML:Y*/
		for(int i = 0; i < ie; i++){
			for(int  j = jb - pml_right; j < je;j++){
				for(int k = 0; k < ke;k++){
					int z = j - (jb - pml_right);
					caDxy[IDX(i,j,k)] = (2 - dt*s_Dy_right[z])/(2 + dt*s_Dy_right[z]);
					cbDxy[IDX(i,j,k)] = (dt/ddy)/(1 + dt*s_Dy_right[z]*0.5);
					daBzy[IDX(i,j,k)] = (2 - dt*s_By_right[z])/(2 + dt*s_By_right[z]);
					dbBzy[IDX(i,j,k)] = (dt/ddy)/(1 + dt*s_By_right[z]*0.5);
					    
					caDzy[IDX(i,j,k)] = (2 - dt*s_Dy_right[z])/(2 + dt*s_Dy_right[z]);
					cbDzy[IDX(i,j,k)] = (dt/ddx)/(1 + dt*s_Dy_right[z]*0.5);
					daBxy[IDX(i,j,k)] = (2 - dt*s_By_right[z])/(2 + dt*s_By_right[z]);
					dbBxy[IDX(i,j,k)] = (dt/ddy)/(1 + dt*s_By_right[z]*0.5);
				}
			}
		}
		printf("RightDone\n");
		/*Bottom PML:Z*/
		for(int i = 0; i < ie; i++){
			for(int j = 0; j < je; j++){
				for(int  k = pml_bottom; k > 0;k--){
					int z = pml_bottom - k;
					/*Faulty*/
					caDxz[IDX(i,j,k)] = (2 - dt*s_Dz_bottom[z])/(2 + dt*s_Dz_bottom[z]);
					cbDxz[IDX(i,j,k)] = (dt/ddz)/(1 + dt*s_Dz_bottom[z]*0.5);
					daByz[IDX(i,j,k)] = (2 - dt*s_Bz_bottom[z])/(2 + dt*s_Bz_bottom[z]);
					dbByz[IDX(i,j,k)] = (dt/ddz)/(1 + dt*s_Bz_bottom[z]*0.5);
					    
					caDyz[IDX(i,j,k)] = (2 - dt*s_Dz_bottom[z])/(2 + dt*s_Dz_bottom[z]);
					cbDyz[IDX(i,j,k)] = (dt/ddz)/(1 + dt*s_Dz_bottom[z]*0.5);
					daBxz[IDX(i,j,k)] = (2 - dt*s_Bz_bottom[z])/(2 + dt*s_Bz_bottom[z]);
					dbBxz[IDX(i,j,k)] = (dt/ddz)/(1 + dt*s_Bz_bottom[z]*0.5);
				}
			}
		}
		printf("BottomDone\n");
		/*Top PML:Z*/
		for(int i = 0; i < ie; i++){
			for(int j = 0; j < je; j++){
				for(int  k = (kb - pml_top); k < ke;k++){
					int z = k - (kb - pml_top);
					caDxz[IDX(i,j,k)] = (2 - dt*s_Dz_top[z])/(2 + dt*s_Dz_top[z]);
					cbDxz[IDX(i,j,k)] = (dt/ddz)/(1 + dt*s_Dz_top[z]*0.5);
					daByz[IDX(i,j,k)] = (2 - dt*s_Bz_top[z])/(2 + dt*s_Bz_top[z]);
					dbByz[IDX(i,j,k)] = (dt/ddz)/(1 + dt*s_Bz_top[z]*0.5);
					    
					caDyz[IDX(i,j,k)] = (2 - dt*s_Dz_top[z])/(2 + dt*s_Dz_top[z]);
					cbDyz[IDX(i,j,k)] = (dt/ddz)/(1 + dt*s_Dz_top[z]*0.5);;
					daBxz[IDX(i,j,k)] = (2 - dt*s_Bz_top[z])/(2 + dt*s_Bz_top[z]);
					dbBxz[IDX(i,j,k)] = (dt/ddz)/(1 + dt*s_Bz_top[z]*0.5);
				}
			}
		}
		printf("TopDone\n");
		std::vector<Value> orient(array_size, 0);
		{	/********* Init Orientation Array ************/
			LOOP_TRIPLE(ie/4, ie/2, je/4, je/2, ke/4, ke/2)
				orient[IDX(i,j,k)] = 1;
			LOOP_CLOSE();
			LOOP_TRIPLE(ie/4, ie/2, je/4, je/2, ke/2, 3*ke/4)
				orient[IDX(i,j,k)] = 2;
			LOOP_CLOSE();
			LOOP_TRIPLE(ie/4, ie/2, je/2, 3*je/4, ke/4, ke/2)
				orient[IDX(i,j,k)] = 3;
			LOOP_CLOSE();
			LOOP_TRIPLE(ie/4, ie/2, je/2, 3*je/4, ke/2, 3*ke/4)
				orient[IDX(i,j,k)] = 4;
			LOOP_CLOSE();
			
			LOOP_TRIPLE(ie/2, 3*ie/4, je/4, je/2, ke/4, ke/2)
				orient[IDX(i,j,k)] = 5;
			LOOP_CLOSE();
			LOOP_TRIPLE(ie/2, 3*ie/4, je/4, je/2, ke/2, 3*ke/4)
				orient[IDX(i,j,k)] = 6;
			LOOP_CLOSE();
			LOOP_TRIPLE(ie/2, 3*ie/4, je/2, 3*je/4, ke/4, ke/2)
				orient[IDX(i,j,k)] = 7;
			LOOP_CLOSE();
			LOOP_TRIPLE(ie/2, 3*ie/4, je/2, 3*je/4, ke/2, 3*ke/4)
				orient[IDX(i,j,k)] = 8;
			LOOP_CLOSE();
		}
		StoreSimParams("params/simData.txt");
		
		DumpParameterArray(orient, "orient");
		
		DumpParameterArray(caDxy,"caDxy" );
		DumpParameterArray(caDxz,"caDxz" );
		DumpParameterArray(caDyx,"caDyx" );
		DumpParameterArray(caDyz,"caDyz" );
		DumpParameterArray(caDzx,"caDzx" );
		DumpParameterArray(caDzy,"caDzy" );
		
		DumpParameterArray(cbDxy,"cbDxy" );
		DumpParameterArray(cbDxz,"cbDxz" );
		DumpParameterArray(cbDyx,"cbDyx" );
		DumpParameterArray(cbDyz,"cbDyz" );
		DumpParameterArray(cbDzx,"cbDzx" );
		DumpParameterArray(cbDzy,"cbDzy" );
		
		DumpParameterArray(daBxy,"daBxy" );
		DumpParameterArray(daBxz,"daBxz" );
		DumpParameterArray(daByx,"daByx" );
		DumpParameterArray(daByz,"daByz" );
		DumpParameterArray(daBzx,"daBzx" );
		DumpParameterArray(daBzy,"daBzy" );
		
		DumpParameterArray(dbBxy,"dbBxy" );
		DumpParameterArray(dbBxz,"dbBxz" );
		DumpParameterArray(dbByx,"dbByx" );
		DumpParameterArray(dbByz,"dbByz" );
		DumpParameterArray(dbBzx,"dbBzx" );
		DumpParameterArray(dbBzy,"dbBzy" );
	}
	
	{	/***** Allocate for General Algorithm Variable ***********/
		Px.resize(array_size*npoles, 0);
		Px1.resize(array_size*npoles, 0);
		Px2.resize(array_size*npoles, 0);
		
		Py.resize(array_size*npoles, 0);
		Py1.resize(array_size*npoles, 0);
		Py2.resize(array_size*npoles, 0);
		
		Pz.resize(array_size*npoles, 0);
		Pz1.resize(array_size*npoles, 0);
		Pz2.resize(array_size*npoles, 0);
		
		epsx.resize(array_size, eps0);
		epsy.resize(array_size, eps0);
		epsz.resize(array_size, eps0);
		
		C1x.resize(array_size*npoles, 0);
		C1y.resize(array_size*npoles, 0);
		C1z.resize(array_size*npoles, 0);
		
		C2x.resize(array_size*npoles, 0);
		C2y.resize(array_size*npoles, 0);
		C2z.resize(array_size*npoles, 0);
		
		C3x.resize(array_size*npoles, 0);
		C3y.resize(array_size*npoles, 0);
		C3z.resize(array_size*npoles, 0);
		
		material_map.resize(array_size, 0);
	}
	
	{	/***** Allocate for General Algorithm Variable ***********/
		profile.resize(jb*kb);
	}
}

void Uniform3DGrid::Simulate(Source src)
{
	Value time;
	probes_data.resize((nmax*FIELDS_NUM*probes.size()));
#ifdef CAVITY_TEST
	LOOP_TRIPLE(ie/4, ie/2, je/4, je/2, ke/4, ke/2)
		Hz[IDX(i,j,k)] = 10;
	LOOP_CLOSE();
	LOOP_TRIPLE(ie/4, ie/2, je/4, je/2, ke/4, ke/2)
		Hx[IDX(i,j,k)] = 10;
	LOOP_CLOSE();
#endif
	
	printf("Probes:------------------------\n");
	for(auto prb: probes){
		printf("%50s \t %3d %3d %3d\n", prb.name, prb.x, prb.y, prb.z);
	}
	
	DumpParameterArray(profile, "profile");
	
	for(nstep = 0; nstep < nmax;nstep++)
	#pragma omp parallel num_threads(2)
	{	
/***********************UpdateSplitElectricFlux*******************************/
		#pragma omp for
		for(int i = 0; i < ie;i++){//X
			for(int  j = 1; j < je;j++){
				for(int k = 1; k < ke;k++){
					Dxy[IDX(i,j,k)] = caDxy[IDX(i,j,k)]*Dxy[IDX(i, j, k)] + 
					cbDxy[IDX(i,j,k)]*(Hz[IDX(i,j,k)] - Hz[IDX(i,j - 1,k)]);
					Dxz[IDX(i,j,k)] = caDxz[IDX(i,j,k)]*Dxz[IDX(i, j, k)] -
					cbDxz[IDX(i,j,k)]*(Hy[IDX(i,j,k)] - Hy[IDX(i,j, k - 1)]);
					Dx[IDX(i,j,k)] = Dxy[IDX(i,j,k)] + Dxz[IDX(i,j,k)];
				}
			}
		}	
		#pragma omp for
		for(int i = 1; i < ie;i++){//Y
			for(int  j = 0; j < je;j++){
				for(int k = 1; k < ke;k++){
					Dyz[IDX(i,j,k)] = caDyz[IDX(i,j,k)]*Dyz[IDX(i, j, k)] + 
					cbDyz[IDX(i,j,k)]*(Hx[IDX(i,j,k)] - Hx[IDX(i, j, k - 1)]);
					Dyx[IDX(i,j,k)] = caDyx[IDX(i,j,k)]*Dyx[IDX(i, j, k)] -
					cbDyx[IDX(i,j,k)]*(Hz[IDX(i,j,k)] - Hz[IDX(i - 1, j, k)]);
					Dy[IDX(i,j,k)] = Dyz[IDX(i,j,k)] + Dyx[IDX(i,j,k)];
				}
			}
		}
		#pragma omp for
		for(int i = 1; i < ie;i++){//Z
			for(int  j = 1; j < je;j++){
				for(int k = 0; k < ke;k++){
					Dzx[IDX(i,j,k)] = caDzx[IDX(i,j,k)]*Dzx[IDX(i, j, k)] + 
					cbDzx[IDX(i,j,k)]*(Hy[IDX(i,j,k)] - Hy[IDX(i - 1, j ,k)]);
					Dzy[IDX(i,j,k)] = caDzy[IDX(i,j,k)]*Dzy[IDX(i, j, k)] -
					cbDzy[IDX(i,j,k)]*(Hx[IDX(i,j,k)] - Hx[IDX(i, j - 1, k)]);
					Dz[IDX(i,j,k)] = Dzx[IDX(i,j,k)] + Dzy[IDX(i,j,k)];
				}
			}
		}
		/*END***UpdateSplitElectricFlux*/
		
	
		#pragma omp for
		for(int i = 0; i < ie;i++){//X
			for(int  j = 1; j < je;j++){
				for(int k = 1; k < ke;k++){
					for(int p = 0; p < npoles; p++){
						Px[IDXP(i,j,k,p)] = C1x[IDXP(i,j,k,p)]*Px1[IDXP(i,j,k,p)] + 
											C2x[IDXP(i,j,k,p)]*Px2[IDXP(i,j,k,p)] + 
											C3x[IDXP(i,j,k,p)]*Ex[IDX(i,j,k)];
					}
				}
			}
		}
		#pragma omp for
		for(int i = 1; i < ie;i++){//Y
			for(int  j = 0; j < je;j++){
				for(int k = 1; k < ke;k++){
					for(int p = 0; p < npoles; p++){
						Py[IDXP(i,j,k,p)] = C1y[IDXP(i,j,k,p)]*Py1[IDXP(i,j,k,p)] + 
											C2y[IDXP(i,j,k,p)]*Py2[IDXP(i,j,k,p)] + 
											C3y[IDXP(i,j,k,p)]*Ey[IDX(i,j,k)];
					}
				}
			}
		}
		#pragma omp for
		for(int i = 1; i < ie;i++){//Z
			for(int  j = 1; j < je;j++){
				for(int k = 0; k < ke;k++){
					for(int p = 0; p < npoles; p++){
						Pz[IDXP(i,j,k,p)] = C1z[IDXP(i,j,k,p)]*Pz1[IDXP(i,j,k,p)] + 
											C2z[IDXP(i,j,k,p)]*Pz2[IDXP(i,j,k,p)] + 
											C3z[IDXP(i,j,k,p)]*Ez[IDX(i,j,k)];
					}
				}
			}
		}
	
		
		/*UpdateElectricField*/
		#pragma omp for
		for(int i = 0; i < ie;i++){
			for(int  j = 0; j < je;j++){
				for(int k = 0; k < ke;k++){
					Value PxT = 0;Value PyT = 0; Value PzT = 0;
					for(int  p = 0; p < npoles; p++){
						PxT += Px[IDXP(i,j,k,p)];
						PyT += Py[IDXP(i,j,k,p)];
						PzT += Pz[IDXP(i,j,k,p)];
					}
					Ex[IDX(i, j, k)] = (Dx[IDX(i, j, k)] - PxT)/epsx[IDX(i,j,k)];
					Ey[IDX(i, j, k)] = (Dy[IDX(i, j, k)] - PyT)/epsy[IDX(i,j,k)];
					Ez[IDX(i, j, k)] = (Dz[IDX(i, j, k)] - PzT)/epsz[IDX(i,j,k)];
				}
			}
		}
		//Ez[IDX(ic, jc, kc)] = src(nstep, time);
		/*UpdateSplitMagneticFlux*/
		#pragma omp for
		for(int i = 1; i < ie;i++){//X
			for(int  j = 0; j < je;j++){
				for(int k = 0; k < ke;k++){
					Bxz[IDX(i, j, k)] = daBxz[IDX(i,j,k)]*Bxz[IDX(i, j, k)] + 
					dbBxz[IDX(i,j,k)]*(Ey[IDX(i,j,k+1)] - Ey[IDX(i, j ,k)]);
					Bxy[IDX(i, j, k)] = daBxy[IDX(i,j,k)]*Bxy[IDX(i, j, k)] - 	//AKHHHHHHHHH THIS IS BADDDDD!! daBxy dbBxy!!
					dbBxy[IDX(i,j,k)]*(Ez[IDX(i,j+1,k)] - Ez[IDX(i, j ,k)]);
					Bx[IDX(i,j,k)] = Bxz[IDX(i, j, k)] + Bxy[IDX(i, j, k)];
				}
			}
		}
		#pragma omp for
		for(int i = 0; i < ie;i++){//Y
			for(int  j = 1; j < je;j++){
				for(int k = 0; k < ke;k++){
					Byx[IDX(i, j, k)] = daByx[IDX(i,j,k)]*Byx[IDX(i, j, k)] +
					dbByx[IDX(i,j,k)]*(Ez[IDX(i+1,j,k)] - Ez[IDX(i, j ,k)]);
					Byz[IDX(i, j, k)] = daByz[IDX(i,j,k)]*Byz[IDX(i, j, k)] - 
					dbByz[IDX(i,j,k)]*(Ex[IDX(i,j,k+1)] - Ex[IDX(i, j ,k)]);
					By[IDX(i,j,k)] = Byx[IDX(i,j,k)] + Byz[IDX(i,j,k)];
				}
			}
		}	
		#pragma omp for
		for(int i = 0; i < ie;i++){//Z
			for(int  j = 0; j < je;j++){
				for(int k = 1; k < ke;k++){
					Bzy[IDX(i, j, k)] = daBzy[IDX(i,j,k)]*Bzy[IDX(i, j, k)] + 
					dbBzy[IDX(i,j,k)]*(Ex[IDX(i,j+1,k)] - Ex[IDX(i, j ,k)]);
					Bzx[IDX(i, j, k)] = daBzx[IDX(i,j,k)]*Bzx[IDX(i, j, k)] -
					dbBzx[IDX(i,j,k)]*(Ey[IDX(i+1,j,k)] - Ey[IDX(i, j ,k)]);
					Bz[IDX(i,j,k)] = Bzy[IDX(i,j,k)] + Bzx[IDX(i,j,k)];
				}
			}
		}

		/*UpdateMagneticField*/
		#pragma omp for
		for(int i = 0; i < ie;i++){
			for(int  j = 0; j < je;j++){
				for(int k = 0; k < ke;k++){
					Hx[IDX(i, j, k)] = Bx[IDX(i, j, k)]/mu0;
					Hy[IDX(i, j, k)] = By[IDX(i, j, k)]/mu0;
					Hz[IDX(i, j, k)] = Bz[IDX(i, j, k)]/mu0;
				}
			}
		}
		
		#pragma omp master
		{
			/* Swap Polarization Arrays */
			// printf("%x\n", Px.data());
			// printf("%x\n", Px1.data());
			// printf("%x\n---\n", Px2.data());
			
			std::swap(Px, Px1);
			std::swap(Px, Px2);
			std::swap(Py, Py1);
			std::swap(Py, Py2);
			std::swap(Pz, Pz1);
			std::swap(Pz, Pz2);
			
			// printf("%x\n", Px.data());
			// printf("%x\n", Px1.data());
			// printf("%x\n====\n", Px2.data());
			
			// if(nstep == 2){
				// exit(1);
			// }
			
			time = nstep * dt;
#ifndef CAVITY_TEST
			//Ez[IDX(ic, jc, kc)] = src(nstep, time);
#endif
			Value m = src(nstep, time);
			for(int  j = 0; j < je; j++){
				for(int k = 0; k < ke;k++){
					Hz[j*kb + k] = profile[j*kb + k]*m;
				}
			}
			
			if((nstep % progress) == 0){
				printf("nstep = %d\n", nstep);
			}
			if(movie_step && ((nstep % movie_step) == 0)){
				DumpFieldArray(Hz, "Hz", nstep);
			}
			for(int i = 0; i < probes.size();i++){
				probes_data[i*(nmax*FIELDS_NUM) + nstep*FIELDS_NUM+0] = Ex[IDX(probes[i].x, probes[i].y, probes[i].z)];
				probes_data[i*(nmax*FIELDS_NUM) + nstep*FIELDS_NUM+1] = Ey[IDX(probes[i].x, probes[i].y, probes[i].z)];
				probes_data[i*(nmax*FIELDS_NUM) + nstep*FIELDS_NUM+2] = Ez[IDX(probes[i].x, probes[i].y, probes[i].z)];
				                                       
				probes_data[i*(nmax*FIELDS_NUM) + nstep*FIELDS_NUM+3] = Hx[IDX(probes[i].x, probes[i].y, probes[i].z)];
				probes_data[i*(nmax*FIELDS_NUM) + nstep*FIELDS_NUM+4] = Hy[IDX(probes[i].x, probes[i].y, probes[i].z)];
				probes_data[i*(nmax*FIELDS_NUM) + nstep*FIELDS_NUM+5] = Hz[IDX(probes[i].x, probes[i].y, probes[i].z)];
			}
		}
	}
	
	/************** Dump probes Data *************************/
	int k = 0;
	printf("Dumping Probes Data:------------------------\n");
	for(auto prb: probes){
		printf("%50s \t %3d %3d %3d\n", prb.name, prb.x, prb.y, prb.z);
		FILE* file = fopen(prb.name, "w");
		if(!file){
			printf("Could not open file\n");
			exit(1);
		}
		for(int i = 0; i < nmax; i++){
			int idx = k*(nmax*FIELDS_NUM) + i*FIELDS_NUM;
			fprintf(file, "%e %e %e %e %e %e\n",
				probes_data[idx + 0], probes_data[idx + 1], probes_data[idx + 2], 
				probes_data[idx + 3], probes_data[idx + 4], probes_data[idx + 5]);
		}
		fclose(file);
		k++;
	}
	
}
void Uniform3DGrid::PrintSimData(const char* filename)
{
	FILE* file = NULL;
	if(!filename){
		file = stdout;
	}else{
		file = fopen(filename, "w");
		if(!file){
			printf("Could not open file %s\n", filename);
			exit(1);
		}
	}
	fprintf(file, "-------------------------------\n");
	fprintf(file, "Grid Size: %d x %d x %d\n", ib, jb, kb);
	fprintf(file, "Simulation Steps: %d\n", nmax);
	fprintf(file, "Time Step: %e\n", dt);
	if(filename && file){
		fclose(file);
	}
}

void Uniform3DGrid::DumpParameterArray(std::vector<Value> param, const char* param_name)
{
	char filename[FILENAME_MAX];
	const char* folder = "params";
	sprintf(filename, "%s/%s", folder, param_name);
	FILE* filemv = fopen(filename, "wb");
	if(!filemv){
		printf("Could not open file %s\n", filename);
		exit(1);
	}
	fwrite(param.data(), sizeof(Uniform3DGrid::Value), param.size(), filemv);
	fclose(filemv);
}

void Uniform3DGrid::DumpFieldArray(std::vector<Value> field, const char* fieldname, int no)
{
	char filename[FILENAME_MAX];
	const char* folder = "snapshots";
	sprintf(filename, "%s/%s%d", folder, fieldname, no);
	FILE* filemv = fopen(filename, "wb");
	if(!filemv){
		printf("Could not open file %s\n", filename);
		exit(1);
	}
	fwrite(field.data(), sizeof(Uniform3DGrid::Value), field.size(), filemv);
	fclose(filemv);
}

void Uniform3DGrid::StoreSimParams(const char* filename)
{
	FILE* file = NULL;	
	if(!filename){
		exit(1);
	}else{
		file = fopen(filename, "w");
		if(!file){
			printf("Could not open file %s\n", filename);
			exit(1);
		}
	}
	fprintf(file, "%d\n", ib);
	fprintf(file, "%d\n", jb);
	fprintf(file, "%d\n", kb);
	fprintf(file, "%d\n", nmax);
	fprintf(file, "%e\n", dt);
	fprintf(file, "%e\n", ddx);
	fprintf(file, "%e\n", ddy);
	fprintf(file, "%e\n", ddz);
	fprintf(file, "%d\n", movie_step);
	fprintf(file, "%d\n", 1); //Row Major true for C++ false for FORTRAN
	fprintf(file, "%d\n", sizeof(Uniform3DGrid::Value)); 

	if(filename && file){
		fclose(file);
	}
}
void DumpSigmaArray(std::vector<Uniform3DGrid::Value>& arr, const char* filename)
{
	FILE* file = fopen(filename, "w");
	if(!file){
		printf("Could not open file %s\n", filename);
		exit(1);
	}
	for(int i = 0; i < arr.size();i++){
		fprintf(file, "%f\n", arr[i]);
	}
	fclose(file);
}
void Uniform3DGrid::FillRange(Material mat, Pred pred, CoordRef ref)
{
	//int is, ie, js, je, ks, ke; Hahahaha
	printf("Filling Material:");
	int xi, yj, zk;
	switch(ref){
		case LeftBottomRef:
			xi = pml_back;
			yj = pml_left;
			zk = pml_bottom;
			break;
		case CenterRef:
			xi = ic;
			yj = jc;
			zk = kc;
			break;
	}
	// is = slimits[0] / ddx;
	// ie = slimits[1] / ddx;
	// js = slimits[2] / ddy;
	// je = slimits[3] / ddy;
	// ks = slimits[4] / ddz;
	// ke = slimits[5] / ddz;
	
	printf("Ref = %d, %d, %d\n", xi, yj, zk);
	
	for(int i = 0; i < ie; ++i){
		for(int  j = 0; j < je; j++){
			for(int k = 0; k < ke;k++){
				Value x, y, z;
				x = (i - xi)*ddx;
				y = (j - yj)*ddy;
				z = (k - zk)*ddz;
				// printf("ref = %e, %e, %e\n", x, y, z);
				//x,y,z now represent the coordinates of the cell corner
				//According to our Yee Cell Ex staggered in X, Ey in Y and so on;
				if(pred(x + 0.5*ddx, y, z)){
					epsx[IDX(i,j,k)] = mat.eps*eps0;
					for(int p = 0; p < npoles; p++){
						C1x[IDX(i,j,k)*npoles + p] = mat.C1[p];
						C2x[IDX(i,j,k)*npoles + p] = mat.C2[p];
						C3x[IDX(i,j,k)*npoles + p] = mat.C3[p];
					}
					// printf("True\n");
					material_map[IDX(i,j,k)] = mat.matNo;
				}
				if(pred(x, y + 0.5*ddy, z)){
					epsy[IDX(i,j,k)] = mat.eps*eps0;
					for(int p = 0; p < npoles; p++){
						C1y[IDX(i,j,k)*npoles + p] = mat.C1[p];
						C2y[IDX(i,j,k)*npoles + p] = mat.C2[p];
						C3y[IDX(i,j,k)*npoles + p] = mat.C3[p];
					}
					//material_map[IDX(i,j,k)] = mat.matNo;
				}
				if(pred(x, y, z + 0.5*ddz)){
					epsz[IDX(i,j,k)] = mat.eps*eps0;
					for(int p = 0; p < npoles; p++){
						C1z[IDX(i,j,k)*npoles + p] = mat.C1[p];
						C2z[IDX(i,j,k)*npoles + p] = mat.C2[p];
						C3z[IDX(i,j,k)*npoles + p] = mat.C3[p];
					}
					//material_map[IDX(i,j,k)] = mat.matNo;
				}
			}
		}
	}
	printf("Done Filling Material");
}
void Uniform3DGrid::DumpMaterialArrays()
{
	DumpParameterArray(epsx, "epsx");
	DumpParameterArray(epsy, "epsy");
	DumpParameterArray(epsz, "epsz");
	
	DumpParameterArray(C1x, "C1x");
	DumpParameterArray(C1y, "C1y");
	DumpParameterArray(C1z, "C1z");
	
	DumpParameterArray(C2x, "C2x");
	DumpParameterArray(C2y, "C2y");
	DumpParameterArray(C2z, "C2z");
	    
	DumpParameterArray(C3x, "C3x");
	DumpParameterArray(C3y, "C3y");
	DumpParameterArray(C3z, "C3z");
	DumpParameterArray(material_map, "material_map");
}