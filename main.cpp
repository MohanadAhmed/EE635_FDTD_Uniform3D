#include <stdio.h>
#include <vector>
#include <fstream>
#include <omp.h>
#include <windows.h>
#include <cfenv>
#include "Uniform3DGrid.h"
#include "timer.h"
/******************************************************************************
		Name: General Polzarization Algorithm
		Author: Mohanad Ahmed
		Date: 01/12/2013
		Note: 	This program implements and test the General Polarization 
		Algorithm.
******************************************************************************/

static const int nTime = 10000;
static const double cLength = 1.2e-7;
static const double delta = 2e-9;
static const float lambda = 750e-9;

static const float thickness = 20e-9;
static const float gap = 20e-9;
static const float w_tri = 75e-9;
static const float width = 2*w_tri + gap;
static const float h_tri1 = 10e-9;
static const float h_tri2 = 40e-9;
static const float slope = (h_tri2 - h_tri1)/(w_tri);
		
static const int pml_cells = 10;
static const double prbDist = cLength/2 - 2*delta;
static const int t0 = 1000;
static const double tp = t0/5;

int main(int argc, char* argv[])
{
	printf("Starting 3D Electromagnetic Solver ....... \n");
	//fesetexceptflag (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
	Timer simTime;
	Uniform3DGrid zgrid(nTime);
	{
		zgrid.SetDomainSize(120e-9, 200e-9, 100e-9);
		zgrid.SetSpatialStep(delta, delta, delta);
		zgrid.SetPMLSize(pml_cells,pml_cells,0,pml_cells,pml_cells,pml_cells);
		zgrid.SetTimeStepFactor(0.9);//0.57735026918962576450914878050196
		zgrid.SetProgressStep(50);
		zgrid.SetSnapshotStep(0);
		zgrid.Allocate();
		
		zgrid.AddProbePoint("probes/source", 0.0e-2, 0.0e-2, 0.0e-2, Uniform3DGrid::CenterRef);
		zgrid.AddProbePoint("probes/probeX", prbDist, 0.0, 0.0, Uniform3DGrid::CenterRef);
		zgrid.AddProbePoint("probes/probeY", 0.0, prbDist, 0.0, Uniform3DGrid::CenterRef);
		zgrid.AddProbePoint("probes/probeZ", 0.0, 0.0, prbDist, Uniform3DGrid::CenterRef);
		
		zgrid.AddProbePoint("probes/bt_center",0,0,0,Uniform3DGrid::CenterRef);
		
		zgrid.AddProbePoint("probes/bt_tip_center",0,-gap/2,0,Uniform3DGrid::CenterRef);
		zgrid.AddProbePoint("probes/bt_tip_edge",0,-gap/2,h_tri1,Uniform3DGrid::CenterRef);
		
		zgrid.AddProbePoint("probes/bt_tail_center",0,-width/2,0,Uniform3DGrid::CenterRef);
		zgrid.AddProbePoint("probes/bt_tail_edge",0,-width/2,h_tri2,Uniform3DGrid::CenterRef);
		
		zgrid.AddProbePoint("probes/bt_center_center",0,-(gap + w_tri)/2,0,Uniform3DGrid::CenterRef);
		zgrid.AddProbePoint("probes/bt_center_edge",0,-(gap + w_tri)/2,(h_tri1 + h_tri2)/2,Uniform3DGrid::CenterRef);
		
		
		auto circle = [](float x, float y, float z){
			if( (x*x + y*y + z*z) < pow(cLength/4, 2)){
				return true;
			}else{ return false;}
		};
		
		auto tester = [](float x, float y, float z){
			//printf("INV\n");
			return true;
		};
		
		auto bowtie_simple = [](float x, float y, float z)
		{
			if( (x < (-thickness/2)) || (x > (thickness/2))) 
				return false;
			else
			{
				if( (y < (-width/2)) || (y > (width/2))) 
					return false;
				else if ( (y > (-gap/2)) && (y < (gap/2))) return false;
				else
				{
					if ((z >= 0) && (y < 0)){
						y = -y;
					}else if ((z < 0) && (y >= 0)){
						z = -z;
					}else if ((z < 0) && (y < 0)){
						z = -z; y = -y;
					}
					if ((z >= 0) && (y >= 0)){
						if( (z < h_tri1))	return true;
						if( (z > h_tri2)) 	return false;
						z = z - h_tri1;
						y = y - gap/2;
						if(z <= slope*y) return true;
						else return false;
					}
					
					return false;
				}
			}
		};
		
		auto bowtie = [](float x, float y, float z)
		{
			if( (x < (-thickness/2)) || (x > (thickness/2))) 
				return false;
			else
			{
				if( (y < (-width/2)) || (y > (width/2))) 
					return false;
				else if ( (y > (-gap/2)) && (y < (gap/2))) return false;
				else
				{
					if (y < 0) y = y + gap/2;
					else y = y - gap/2;
					
					if ((z > 0) && (y > 0)){
						if(z < y) return true;
					}else if ((z > 0) && (y < 0)){
						if(z < -y) return true;
					}else if ((z < 0) && (y > 0)){
						if(-z < y) return true;
					}else if ((z < 0) && (y < 0)){
						if(-z < -y) return true;
					}
					return false;
				}
			}
		};

		Uniform3DGrid::Material mat = zgrid.GetKnownMaterial(Uniform3DGrid::Silver);
		mat.matNo = 3;
		zgrid.FillRange(mat, bowtie_simple , Uniform3DGrid::CenterRef);
		zgrid.DumpMaterialArrays();
		//exit(1);
		zgrid.SetModeProfile([](float x, float y)
		{
			return (float)( expf(-((x*x + y*y)/(25*delta*delta))));
		}, Uniform3DGrid::CenterRef);
		
		
		
		printf("Simulating Small Grid\n");
		simTime.Start();
		static const float freq = (c0/lambda);
		zgrid.Simulate([=](int n, Uniform3DGrid::Value time){
			float z = pow((n - t0)/(1.0*tp), 2);
			return (Uniform3DGrid::Value)((1000*exp(-z))*cos((2*pi*freq*time)));
		});
		simTime.Stop();
		printf("Simulation Time = %Gs\n", simTime.Elapsed()/1000.0);
	}
}