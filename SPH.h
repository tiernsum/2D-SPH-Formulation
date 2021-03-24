#ifndef CLASS_SPH
#define CLASS_SPH

class SPH {
    
private:
    
    const double    k      = 2000.0;   // Gas constant
    const double    rho_0  = 1000.0;   // Resting density
    const double    mu     = 1.0;      // Viscosity
    const double    g      = 9.81;     // Acceleration due to gravity
    const double    e      = 0.5;      // Coefficient of restitution
    const double    m_init = 1.0;      // Initialised particle mass
    
    unsigned int    N;                 // Number of particles
    double          h;                 // Radius of influence
    double          dt;                // Time-step
    double          T;                 // Total integration time
    double          t;                 // Current time-step
    
    // Global set of arrays
    double          m;
    double*         x;
    double*         rho;
    double*         rhoInit;
    double*         p;
    double*         v;
    double*         a;
    double*         F_p;
    double*         F_v;
    double*         F_g;
    double*         r_ij;
    double*         v_ij;
    
    // Arrays that are local to each processor
    double*         xLocal;
    double*         rhoLocal;
    double*         pLocal;
    double*         vLocal;
    double*         aLocal;
    double*         F_pLocal;
    double*         F_vLocal;
    double*         F_gLocal;
    
    // for MPI
    unsigned int    localN;
    unsigned int    startN;
    unsigned int    endN;
    unsigned int    currRank;
    
    // for particle set-up
    unsigned int    counter;
    double          xCoor;
    double          yCoor;
    
public:
    
    // Upon initialisation of SPH object
    SPH(const unsigned int& numOfParticles, const double& timeStep, const double& finalT, const double& radOfInfl);
    
    std::ofstream xOut;
    std::ofstream energyOut;
    
    void getExecCase(unsigned int caseID);
    void calcDensityInit();
    double scaleMass();
    
    // Iterative Algorithm
    void iterate(std::ofstream& xCoor, std::ofstream& energyTxt, unsigned int localN, unsigned int nProc, unsigned int currRank);
    
    void calcDensityWithMPI(unsigned int localN, unsigned int nProc, unsigned int currRank);
    void calcPressureWithMPI(unsigned int localN, unsigned int nProc, unsigned int currRank);
    void calcPressureForceWithMPI(unsigned int localN, unsigned int nProc, unsigned int currRank);
    void calcViscousForceWithMPI(unsigned int localN, unsigned int nProc, unsigned int currRank);
    void calcGravityForceWithMPI(unsigned int localN, unsigned int nProc, unsigned int currRank);
    void calcAcceleration();
    void getNextParticleVel();
    void getNextParticlePos();
    
    void calcRij(unsigned int i, unsigned int j);
    void calcVij(unsigned int i, unsigned int j);
    void applyBC(unsigned int i);
    
    // Calculating energy of the system
    double calcKineticEnergy();
    double calcPotentialEnergy();
    double calcTotalEnergy();
    
    // Writing data to files
    void writeParticlePosition(std::ofstream& xCoor);
    void writeEnergy(std::ofstream& energyTxt);
    void closePPOutputFile(std::ofstream& fileName);
    
    // Getters and Setter
    int getN();
    double getDeltaT();
    double getCurrT();
    void setCurrT();
    double getTotalIntTime();
    
    // Deleting arrays on heap to prevent memory leakage
    ~SPH();
    
};

#endif