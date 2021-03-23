#ifndef CLASS_SPH
#define CLASS_SPH

class SPH {
    
private:
    
    const double k = 2000.0;           // Gas constant
    const double rho_0 = 1000.0;       // Resting density
    const double mu = 1.0;             // Viscosity
    const double g = 9.81;             // Acceleration due to gravity
    const double e = 0.5;              // Coefficient of restitution
    const double m_init = 1.0;         // Initialised particle mass
    
    unsigned int N;                    // Number of particles
    double h;                          // Radius of influence
    double dt;                         // Time-step delta_t = 10^-4
    double T;                          // Total integration time
    double t;                          // Current time-step
    
    double m;
    double* x;
    double* xLocal;
    double* rho;
    double* rhoInit;
    double* rhoLocal;
    double* p;
    double* pLocal;
    double* v;
    double* vLocal;
    double* a;
    double* aLocal;
    double* F_p;
    double* F_pLocal;
    double* F_v;
    double* F_vLocal;
    double* F_g;
    double* F_gLocal;
    double* r_ij;
    double* v_ij;
    
    // forMPI
    int localN;
    int startN;
    int endN;
    int currRank;
    
public:
    
    SPH(const unsigned int& numOfParticles, const double& timeStep, const double& finalT, const double& radOfInfl);
    
    std::ofstream xOut;
    std::ofstream energyOut;
    
    void getExecCase(unsigned int caseID);
    void calcRij(unsigned int i, unsigned int j);
    void calcVij(unsigned int i, unsigned int j);
    void calcDensityInit();
    void calcAcceleration();
    void generateVInit();
    void getNextParticleVel();
    void getNextParticlePos();
    void applyBC(unsigned int i);
    void setCurrT();
    void writeParticlePosition(std::ofstream& xCoor);
    void writeEnergy(std::ofstream& energyTxt);
    void writeData(std::ofstream& energyTxt);
    void closePPOutputFile(std::ofstream& fileName);
    void calcDensityWithMPI(int localN, unsigned int nProc, int currRank);
    void calcPressureWithMPI(int localN, unsigned int nProc, int currRank);
    void calcPressureForceWithMPI(int localN, unsigned int nProc, int currRank);
    void calcViscousForceWithMPI(int localN, unsigned int nProc, int currRank);
    void calcGravityForceWithMPI(int localN, unsigned int nProc, int currRank);
    void iterate(std::ofstream& xCoor, std::ofstream& energyTxt, std::ofstream& dataTxt, int localN, unsigned int nProc, int currRank);
    
    int getN();
    double scaleMass();
    double getTotalIntTime();
    double getCurrT();
    double getDeltaT();
    double calcKineticEnergy();
    double calcPotentialEnergy();
    double calcTotalEnergy();
    
    ~SPH();
    
};

#endif