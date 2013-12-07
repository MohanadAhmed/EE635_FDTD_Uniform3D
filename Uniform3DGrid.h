#include <vector>
#include <cmath>
#include <assert.h>

#define pi (3.1415926535897932384626433832795F)
#define eps0 (8.854e-12F)
#define mu0 (4 * pi * 1.0e-7F)
#define c0 (299795637.69321621645649782624619F)
#define h (6.626e-34F)
#define eV (1.6e-19F)
#define h_bar (h/(2*pi))

inline int round(float x){
	return floor(x + 0.5);
}

class Uniform3DGrid
{
public:
	typedef float Value;
private:
	/*********** Fields ***********/
	std::vector<Value>		Ex, Ey, Ez;
	std::vector<Value>		Dx, Dy, Dz;
	std::vector<Value>		Hx, Hy, Hz;
	std::vector<Value>		Bx, By, Bz;
	
	/*********** Split Fields ***********/
	std::vector<Value>		Dxy, Dxz, Dyx, Dyz, Dzx, Dzy;
	std::vector<Value>		Bxy, Bxz, Byx, Byz, Bzx, Bzy;
	
	/*********** Coefficients ***********/
	std::vector<Value>		caDxy, caDxz, caDyx, caDyz, caDzx, caDzy;
	std::vector<Value>		cbDxy, cbDxz, cbDyx, cbDyz, cbDzx, cbDzy;
	std::vector<Value>		daBxy, daBxz, daByx, daByz, daBzx, daBzy;
	std::vector<Value>		dbBxy, dbBxz, dbByx, dbByz, dbBzx, dbBzy;
	
	/*********** Polarization GA Fields ***********/
	std::vector<Value>		Px, Px1, Px2;
	std::vector<Value>		Py, Py1, Py2;
	std::vector<Value>		Pz, Pz1, Pz2;
	
	/*********** Polarization GA Coefficients ***********/
	std::vector<Value>		epsx, epsy, epsz;
	std::vector<Value>		C1x, C1y, C1z;
	std::vector<Value>		C2x, C2y, C2z;
	std::vector<Value>		C3x, C3y, C3z;
	
	/*********** Material Map for Display***********/
	std::vector<Value>		material_map;
	
	/*********** Sources ***************************/
	std::vector<Value>		profile;
	
	int nmax;
	Value xx, yy, zz;
	Value ddx, ddy, ddz;
	Value dt;
	Value dt_factor;
	
	int ib, jb, kb, ic, jc, kc;
	int ie, je, ke;
	int nstep;
	
	int pml_front, pml_right, pml_back;
	int pml_left, pml_top, pml_bottom;
	
	int progress;
	int movie_step;
	
	/******Point Probes*******/
	struct Probe {
		const char* name;
		int x; int y; int z;
	public:
		Probe(const char* name, int x, int y, int z): name(name), x(x), y(y), z(z){};
	};
	std::vector<Probe> probes;
	std::vector<Value> probes_data;
public:
	Uniform3DGrid(int nmax): nmax(nmax){};
	
	void SetDomainSize(Value xx, Value yy, Value zz){
		this->xx = xx; this->yy = yy; this->zz = zz;
	}
	void SetSpatialStep(Value ddx, Value ddy, Value ddz){
		this->ddx = ddx; this->ddy = ddy; this-> ddz = ddz;
	}
	void SetPMLSize(int front, int right, int back, int left, int top, int bottom){
		pml_front = front;
		pml_right = right;
		pml_back = back;
		pml_left = left;
		pml_top = top;
		pml_bottom = bottom;
	}
	void SetTimeStepFactor(Value factor){
		this->dt_factor = factor;
	}
	enum CoordRef{
		LeftBottomRef = 0,
		CenterRef
	};
	void AddProbePoint(const char* name, Value x, Value y, Value z, CoordRef ref = LeftBottomRef){
		int prb_i, prb_j, prb_k;
		switch(ref){
			case LeftBottomRef:
				prb_i = pml_back + round(x/ddx);
				prb_j = pml_left + round(y/ddy);
				prb_k = pml_bottom + round(z/ddz);
				break;
			case CenterRef:
				prb_i = pml_back + round(xx/(2*ddx)) + round(x/(ddx));
				prb_j = pml_left + round(yy/(2*ddy)) + round(y/(ddy));
				prb_k = pml_bottom + round(zz/(2*ddz)) + round(z/(ddz));
				break;
		}
		probes.push_back(Probe(name, prb_i, prb_j, prb_k));
	}
	void SetProgressStep(int progress){
		this->progress = progress;
	}
	void SetSnapshotStep(int step){
		this->movie_step = step;
	}
	typedef Value (*SrcProfile2D)(Value x, Value z);
	void SetModeProfile(SrcProfile2D srcprf, CoordRef ref){
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
		for(int  j = 0; j < je; j++){
			for(int k = 0; k < ke;k++){
				Value x, y, z;
				y = (j - yj)*ddy;
				z = (k - zk)*ddz;
				profile[j*kb + k] = srcprf(y, z);
			}
		}
	}
public:
	/**Materials Handling Section**/
	#define npoles 6
	struct Material
	{
		Value eps;
		Value C1[npoles];
		Value C2[npoles];
		Value C3[npoles];
		Value matNo;
	};
	
	typedef bool (*Pred)(Value x, Value y, Value z);
	void FillRange(Material mat, Pred pred, CoordRef ref = LeftBottomRef);
	enum KnownMaterials
	{
		FreeSpace = 0,
		Silver,
		Gold,
		Copper
	};
	
	Material GetKnownMaterial(KnownMaterials km)
	{
		Material mat;
		
		Material FS = {1, 
			{0.0F,0.0F,0.0F,0.0F,0.0F,0.0F}, 
			{0.0F,0.0F,0.0F,0.0F,0.0F,0.0F}, 
			{0.0F,0.0F,0.0F,0.0F,0.0F,0.0F}, 0};
		
		Value Cu_Wp = 10.83F;
		Value Cu_f[] = {0.575F, 0.061F, 0.104F, 0.723F, 0.638F, 0.000F};
		Value Cu_G[] = {0.030F, 0.378F, 1.056F, 3.213F, 4.305F, 0.000F};
		Value Cu_W[] = {0.000F, 0.291F, 2.957F, 5.300F, 11.18F, 0.000F};
		
		Value Ag_Wp = 9.01F;
		Value Ag_f[] = {0.845F, 0.065F, 0.124F, 0.011F, 0.840F, 5.646F};
		Value Ag_G[] = {0.048F, 3.886F, 0.452F, 0.065F, 0.916F, 2.419F};
		Value Ag_W[] = {0.000F, 0.816F, 4.481F, 8.185F, 9.083F, 20.29F};
		
		Cu_Wp =(eV/(h_bar))*Cu_Wp;
		Ag_Wp =(eV/h_bar)*Ag_Wp;
		// GM_Wp =(eV/h_bar)*GM_Wp;
		for(int  k = 0; k < npoles;k++){
			Cu_G[k] = (eV/(h_bar))*Cu_G[k];
			Cu_W[k] = (eV/(h_bar))*Cu_W[k];
			Ag_G[k] = (eV/(h_bar))*Ag_G[k];
			Ag_W[k] = (eV/(h_bar))*Ag_W[k];
			// GM_G[k] = (eV/(h_bar))*GM_G[k];
			// GM_W[k] = (eV/(h_bar))*GM_W[k];
		}
		
		switch(km)
		{
			case FreeSpace:
				return FS;
			break;
			case Silver:
				mat = this->TranslateLD2Mat(1, Ag_Wp, Ag_f, Ag_G, Ag_W);
			break;
			case Gold:
				
			break;
			case Copper:
				mat = this->TranslateLD2Mat(1, Cu_Wp, Cu_f, Cu_G, Cu_W);
			break;
		}
		return mat;
	}
	
	Material TranslateLD2Mat(Value epsinf, Value Wp, Value fs[], Value Gs[], Value Ws[])
	{
		double a,b,c,d;
		Material mat;

		for(int k = 0; k < npoles; k++){
			a = fs[k] * eps0 * Wp * Wp;
			b = Ws[k] * Ws[k];
			c = Gs[k];
			d = 1;

			mat.C1[k] = (4*d - 2*b*(dt*dt))/(2*d + c * dt);
			mat.C2[k] = (-2*d + c*dt)/(2*d + c*dt);
			mat.C3[k] = (2*a*(dt*dt))/(2*d + c*dt);

			printf("%e\t %e\t %e\n", mat.C1[k], mat.C2[k], mat.C3[k]);
		}
		mat.eps = epsinf;
		return mat;
	}
	
	void Allocate();
	void DumpMaterialArrays();
	typedef Value (*Source)(int n, Value time);
	void Simulate(Source src);
	
	//void AdvanceOneTimeStep();
	//void CurrentTimeStep();
private:
	void PrintSimData(const char* file);
	void DumpParameterArray(std::vector<Value> param, const char* param_name);
	void DumpFieldArray(std::vector<Value> field, const char* fieldname, int no);
	void StoreSimParams(const char* filename);
};