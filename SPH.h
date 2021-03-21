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
    double* rho;
    double* rhoInit;
    double* p;
    double* v;
    double* a;
    double* F_p;
    double* F_v;
    double* F_g;
    
public:
    
    SPH(const unsigned int& numOfParticles, const double& timeStep, const double& finalT, const double& radOfInfl);
    
    std::ofstream outputPP;
    
    void getExecCase(unsigned int caseID);
    void calcDensity();
    void calcDensityInit();
    void calcPressure();
    void calcPressureForce();
    void calcViscousForce();
    void calcGravityForce();
    void calcAcceleration();
    void generateVInit();
    void getNextParticleVel();
    void getNextParticlePos();
    void applyBC();
    void setCurrT();
    // void createPPOutputFile(std::ofstream& fileName);
    void writeToPPOutputFile(std::ofstream& fileName);
    void closePPOutputFile(std::ofstream& fileName);
    
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