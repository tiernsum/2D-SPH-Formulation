#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "cblas.h"
#include "mpi.h"
#include "SPH.h"
using namespace std;

SPH::SPH(const unsigned int& numOfParticles, const double& timeStep, const double& finalT, const double& radOfInfl) {
    
    cout << "Starting computation..." << endl;
    
    h = radOfInfl;             // Coefficient of restitution
    dt = timeStep;             // Time-step delta_t = 10^-4
    T = finalT;
    N = numOfParticles;
    
    t = 0;
    
    x = new double[2 * N]();
    v = new double[2 * N]();
    a = new double[2 * N]();
    r_ij = new double[2 * N]();
    v_ij = new double[2 * N]();
    rho = new double[N]();
    rhoInit = new double[N]();
    p = new double[N]();
    F_p = new double[2 * N]();
    F_v = new double[2 * N]();
    F_g = new double[2 * N]();
    
    
    getExecCase(N);
    
    calcDensityInit();
    m = scaleMass();
    
}

SPH::~SPH() {
    
    delete[] x;
    delete[] v;
    delete[] r_ij;
    delete[] v_ij;
    delete[] rho;
    delete[] rhoInit;
    delete[] p;
    delete[] F_p;
    delete[] F_v;
    delete[] F_g;
    cout << "Program Terminated..." << endl;
    
}

int SPH::getN() {
    
    return N;
    
}

double SPH::getCurrT() {
    
    return t;
    
}

void SPH::setCurrT() {
    
    t += dt;
    
}

double SPH::getDeltaT() {
    
    return dt;
    
}

double SPH::getTotalIntTime() {
    
    return T;
    
}

void SPH::getExecCase(unsigned int caseID) {
    
    switch(caseID) {
        
        case 1:
            x[0] = 0.5;
            x[1] = 0.5;
            break;
            
        case 2:
            x[0] = 0.5;
            x[1] = 0.5;
            x[2] = 0.5;
            x[3] = h;
            break;
            
        case 3:
            x[0] = 0.5;
            x[1] = 0.5;
            x[2] = 0.495;
            x[3] = h;
            x[4] = 0.505;
            x[5] = h;
            break;
            
        case 4:
            x[0] = 0.505;
            x[1] = 0.5;
            x[2] = 0.515;
            x[3] = 0.5;
            x[4] = 0.51;
            x[5] = 0.45;
            x[6] = 0.5;
            x[7] = 0.45;
            break;
        
        case 8:
            x[0] = 0.075;
            x[1] = 0.05;
            x[2] = 0.1;
            x[3] = 0.05;
            x[4] = 0.12;
            x[5] = 0.05;
            x[6] = 0.15;
            x[7] = 0.05;
            x[8] = 0.075;
            x[9] = 0.15;
            x[10] = 0.1;
            x[11] = 0.15;
            x[12] = 0.125;
            x[13] = 0.15;
            x[14] = 0.15;
            x[15] = 0.15;
            break;
            
        case 900:
            srand(time(0));
            for (unsigned int i = 0; i < 900; ++i) {
                
                x[2*i] = (double)rand()/RAND_MAX * 0.2;
        
                x[2*i+1] = (double)rand()/RAND_MAX * 0.2;
                
            }
            break;
        
        case 1200:
            srand(time(0));
            for (unsigned int i = 0; i < 1200; ++i) {
                
                x[2*i] = (double)rand()/RAND_MAX * 0.2 + 0.1;
        
                x[2*i+1] = (double)rand()/RAND_MAX * 0.3 + 0.3;
                
            }
            break;
            
        case 360:
            srand(time(0));
            for (unsigned int i = 0; i < 360; ++i) {
                
                x[2*i] = (double)rand()/RAND_MAX * 0.1 * cos(i * M_PI / 180.0) + 0.5;
                
                x[2*i+1] = (double)rand()/RAND_MAX * 0.1 * sin(i * M_PI / 180.0) + 0.7;
                
            }
            
    }
    
}

void SPH::calcRij(unsigned int i, unsigned int j) {
    
    r_ij[2*j] = xLocal[2*i] - x[2*j];
    
    r_ij[2*j+1] = xLocal[2*i+1] - x[2*j+1];
    
}

void SPH::calcVij(unsigned int i, unsigned int j) {
    
    v_ij[2*j] = vLocal[2*i] - v[2*j];
    
    v_ij[2*j+1] = vLocal[2*i+1] - v[2*j+1];
    
}

void SPH::calcDensityInit() {
    
    double q;
    
    for (unsigned int i = 0; i < N; ++i) {
        
        for (unsigned int j = 0; j < N; ++j) {
            
            r_ij[2*j] = x[2*i] - x[2*j];
    
            r_ij[2*j+1] = x[2*i+1] - x[2*j+1];
            
            q = sqrt(r_ij[2*j]*r_ij[2*j] + r_ij[2*j+1]*r_ij[2*j+1]) / h;
            
            if (q < 1) {
                
                rhoInit[i] += ((4 * m_init) / (M_PI * h * h)) * ((1 - (q * q)) * (1 - (q * q)) * (1 - (q * q)));
                
            } else {
                
                rhoInit[i] += 0.0;
                
            }
        
        }
        
    }
    
}

double SPH::scaleMass() {
    
    double sumRho = 0.0;
    double scaledMass = 0.0;
    
    for (unsigned int i = 0; i < N; ++i) {
        
        sumRho += rhoInit[i];
        
    }
    
    scaledMass = (N * rho_0) / sumRho;
    
    return scaledMass;
    
}

void SPH::iterate(ofstream& xCoor, ofstream& energyTxt, ofstream& dataTxt, int localN, unsigned int nProc, int currRank) {
    
    while (getCurrT() <= getTotalIntTime()) {
        
        /*for (unsigned int i = 0; i < N; ++i) {
            
            rho[i] = 0.0;
            F_p[2*i] = 0.0;
            F_p[2*i+1] = 0.0;
            F_v[2*i] = 0.0;
            F_v[2*i+1] = 0.0;
            
        }*/
        
        calcDensityWithMPI(localN, nProc, currRank);
        calcPressureWithMPI(localN, nProc, currRank);
        calcPressureForceWithMPI(localN, nProc, currRank);
        calcViscousForceWithMPI(localN, nProc, currRank);
        calcGravityForceWithMPI(localN, nProc, currRank);
        
        if (currRank == 0) {
            
            calcAcceleration();
            
            if (getCurrT() == 0.0) {
                generateVInit();
            } else {
                getNextParticleVel();
            }
            
            // cout << "Particle Location (@t = " << obj1.getCurrT() << "): " << endl;
            
            getNextParticlePos();
            
            writeEnergy(energyTxt);
            
            writeData(dataTxt);
            
        }
        
        setCurrT();
        
    }
    
    writeParticlePosition(xCoor);
    
    /*cout << "Position at t = " << t - dt << endl;
    for (unsigned int i = 0; i < N; ++i) {
        cout << x[2*i] << " " << x[2*i+1] << endl;
    }*/
    
}

void SPH::calcDensityWithMPI(int localN, unsigned int nProc, int currRank) {
    
    MPI_Bcast(x, 2*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
    if (currRank == 0){
            
        for (unsigned int i = 1; i < nProc; ++i) {
                
            startN = i * localN;
            endN = (i+1) * localN;
                
            // cout << "Start N: " << startN << endl;
            // cout << "End N: " << endN << endl;
                
            MPI_Send(&startN, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(&endN, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
                
            // cout << "Rank " << currRank << " sends " << startN << endl;
            // cout << "Rank " << currRank << " sends " << endN << endl;
             
        }
            
        /*cout << "Particle Location (Rank " << currRank << "): " << endl;
        for (unsigned int i = 0; i < N; ++i) {
            cout << x[2*i] << " " << x[2*i+1] << endl;
        }*/
            
        xLocal = new double[2 * localN]();
        rhoLocal = new double[localN]();
            
        startN = currRank * localN;
        endN = (currRank + 1) * localN;
            
        // cout << "Start N (Rank 0): " << startN << endl;
        // cout << "End N (Rank 0): " << endN << endl;
            
        for (unsigned int i = 0; i < localN; ++i) {
                
            xLocal[2*i] = x[2*(i + startN)];
            xLocal[2*i+1] = x[2*(i + startN) + 1];
                
        }
            
        // cout << "Particle Location (Rank " << currRank << "): " << endl;
        /*for (unsigned int i = 0; i < localN; ++i) {
            cout << "(Rank " << currRank << ") " << "xLocal[" << i << "] (t = " << t << ") " << xLocal[2*i] << " " << xLocal[2*i+1] << endl;
        }*/
            
        double q;
    
        for (unsigned int i = 0; i < localN; ++i) {
            
            for (unsigned int j = 0; j < N; ++j) {
                
                calcRij(i, j);
                
                q = sqrt(r_ij[2*j]*r_ij[2*j] + r_ij[2*j+1]*r_ij[2*j+1]) / h;
                
                if (q < 1) {
                    
                    rhoLocal[i] += ((4 * m) / (M_PI * h * h)) * ((1 - (q * q)) * (1 - (q * q)) * (1 - (q * q)));
                    
                } else {
                    
                    rhoLocal[i] += 0.0;
                    
                }
            
            }
            
        }
        
        // for (unsigned int i = 0; i < localN; ++i) {
        //     cout << "(Rank " << currRank << ") " << "rhoLocal[" << i << "] (t = " << t << ") " << rhoLocal[i] << endl;
        // }
        
    } else {
        
        MPI_Recv(&startN, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&endN, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
        // cout << "Rank " << currRank << " receives " << startN << endl;
        // cout << "Rank " << currRank << " receives " << endN << endl;
            
        /*cout << "Particle Location (Rank " << currRank << "): " << endl;
        for (unsigned int i = 0; i < N; ++i) {
            cout << x[2*i] << " " << x[2*i+1] << endl;
        }*/
            
        xLocal = new double[2 * localN]();
        rhoLocal = new double[localN]();
            
        startN = currRank * localN;
        endN = (currRank + 1) * localN;
            
        // cout << "Start N (Rank 1): " << startN << endl;
        // cout << "End N (Rank 1): " << endN << endl;
            
        for (unsigned int i = 0; i < localN; ++i) {
                
            xLocal[2*i] = x[2*(i + startN)];
            xLocal[2*i+1] = x[2*(i + startN) + 1];
                
        }
            
        // cout << "Particle Location (Rank " << currRank << "): " << endl;
        /*for (unsigned int i = 0; i < localN; ++i) {
            cout << "(Rank " << currRank << ") " << "xLocal[" << i << "] (t = " << t << ") " << xLocal[2*i] << " " << xLocal[2*i+1] << endl;
        }*/
        
            
        double q;
    
        for (unsigned int i = 0; i < localN; ++i) {
            
            for (unsigned int j = 0; j < N; ++j) {
                
                calcRij(i, j);
                
                q = sqrt(r_ij[2*j]*r_ij[2*j] + r_ij[2*j+1]*r_ij[2*j+1]) / h;
                
                if (q < 1) {
                    
                    rhoLocal[i] += ((4 * m) / (M_PI * h * h)) * ((1 - (q * q)) * (1 - (q * q)) * (1 - (q * q)));
                    
                } else {
                    
                    rhoLocal[i] += 0.0;
                    
                }
            
            }
            
        }
            
    }
    
    MPI_Gather(rhoLocal, localN, MPI_DOUBLE, rho, localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
    /*if (currRank == 0) {
        cout << "Particle Density at t = " << t << endl;
        for (unsigned int i = 0; i < N; ++i) {
            cout << rho[i] << endl;
        }
    }*/
}

void SPH::calcPressureWithMPI(int localN, unsigned int nProc, int currRank) {
    
    MPI_Bcast(rho, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
    if (currRank == 0){
            
        for (unsigned int i = 1; i < nProc; ++i) {
                
            startN = i * localN;
            endN = (i+1) * localN;
                
            // cout << "Start N: " << startN << endl;
            // cout << "End N: " << endN << endl;
                
            MPI_Send(&startN, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(&endN, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
                
            // cout << "Rank " << currRank << " sends " << startN << endl;
            // cout << "Rank " << currRank << " sends " << endN << endl;
             
        }
            
        pLocal = new double[localN]();
            
        startN = currRank * localN;
        endN = (currRank + 1) * localN;
            
        // cout << "Start N (Rank 0): " << startN << endl;
        // cout << "End N (Rank 0): " << endN << endl;
            
        for (unsigned int i = 0; i < localN; ++i) {
      
            pLocal[i] = k * (rho[i+startN] - rho_0);
            
        } 
        
        // for (unsigned int i = 0; i < localN; ++i) {
        //     cout << "(Rank " << currRank << ") " << "pLocal[" << i << "] (t = " << t << ") " << pLocal[i] << endl;
        // }
        
    } else {
        
        MPI_Recv(&startN, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&endN, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
        // cout << "Rank " << currRank << " receives " << startN << endl;
        // cout << "Rank " << currRank << " receives " << endN << endl;
            
        pLocal = new double[localN]();
            
        startN = currRank * localN;
        endN = (currRank + 1) * localN;
            
        // cout << "Start N (Rank 1): " << startN << endl;
        // cout << "End N (Rank 1): " << endN << endl;
            
        for (unsigned int i = 0; i < localN; ++i) {
      
            pLocal[i] = k * (rho[i+startN] - rho_0);
            
        } 
            
    }
    
    MPI_Gather(pLocal, localN, MPI_DOUBLE, p, localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
    /*if (currRank == 0) {
        cout << "Particle Pressure at t = " << t << endl;
        for (unsigned int i = 0; i < N; ++i) {
            cout << p[i] << endl;
        }
    }*/
}

void SPH::calcPressureForceWithMPI(int localN, unsigned int nProc, int currRank) {
    
    MPI_Bcast(x, 2*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(rho, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(p, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
    if (currRank == 0){
            
        for (unsigned int i = 1; i < nProc; ++i) {
                
            startN = i * localN;
            endN = (i+1) * localN;
                
            // cout << "Start N: " << startN << endl;
            // cout << "End N: " << endN << endl;
                
            MPI_Send(&startN, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(&endN, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
                
            // cout << "Rank " << currRank << " sends " << startN << endl;
            // cout << "Rank " << currRank << " sends " << endN << endl;
             
        }
            
        /*cout << "Particle Location (Rank " << currRank << "): " << endl;
        for (unsigned int i = 0; i < N; ++i) {
            cout << x[2*i] << " " << x[2*i+1] << endl;
        }*/
            
        xLocal = new double[2 * localN]();
        F_pLocal = new double[2 * localN]();
            
        startN = currRank * localN;
        endN = (currRank + 1) * localN;
            
        // cout << "Start N (Rank 0): " << startN << endl;
        // cout << "End N (Rank 0): " << endN << endl;
            
        for (unsigned int i = 0; i < localN; ++i) {
                
            xLocal[2*i] = x[2*(i + startN)];
            xLocal[2*i+1] = x[2*(i + startN) + 1];
                
        }
            
        // cout << "Particle Location (Rank " << currRank << "): " << endl;
        /*for (unsigned int i = 0; i < localN; ++i) {
            cout << "(Rank " << currRank << ") " << "xLocal[" << i << "] (t = " << t << ") " << xLocal[2*i] << " " << xLocal[2*i+1] << endl;
        }*/
            
        double q;
    
        for (unsigned int i = 0; i < localN; ++i) {
        
            for (unsigned int j = 0; j < N; ++j) {
            
                calcRij(i, j);
                
                // cout << "(Rank " << currRank << ") R_ij at t = " << t << ": " << r_ij[2*j] << " " << r_ij[2*j+1] << endl;
                
                q = sqrt(r_ij[2*j]*r_ij[2*j] + r_ij[2*j+1]*r_ij[2*j+1]) / h;
                
                // cout << "(Rank " << currRank << ") q at t = " << t << ": " << q << endl;
                
                if (q < 1 && q != 0) {
                    
                    // cout << "(Rank: " << currRank << ") Rho[" << j << "] at t = " << t << ": " << rho[j] << endl;
                    
                    // cout << "(Rank: " << currRank << ") p[" << j << "] at t = " << t << ": " << p[j] << endl;
                    
                    F_pLocal[2*i] += ((m / rho[j]) * ((p[i+startN] + p[j]) / 2.0)) * ((-30 / (M_PI * h * h * h)) * (r_ij[2*j]) * (((1 - q) * (1 - q)) / q));
                
                    F_pLocal[2*i+1] += ((m / rho[j]) * ((p[i+startN] + p[j]) / 2.0)) * ((-30 / (M_PI * h * h * h)) * (r_ij[2*j+1]) * (((1 - q) * (1 - q)) / q));
                    
                    /*for (unsigned int i = 0; i < localN; ++i) {
                        cout << "(Rank " << currRank << ") " << "F_p[" << i << "] (t = " << t << ") " << F_p[2*i] << " " << F_p[2*i+1] << endl;
                    }*/
                    
                } else {
                    
                    F_pLocal[2*i] += 0.0;
                
                    F_pLocal[2*i+1] += 0.0;
                    
                }
            
            }
            
            F_pLocal[2*i] *= -1.0;
            F_pLocal[2*i+1] *= -1.0;
            
        }
        
        // for (unsigned int i = 0; i < localN; ++i) {
        //     cout << "(Rank " << currRank << ") " << "F_pLocal[" << i << "] (t = " << t << ") " << F_pLocal[2*i] << " " << F_pLocal[2*i+1] << endl;
        // }
        
    } else {
        
        MPI_Recv(&startN, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&endN, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
        // cout << "Rank " << currRank << " receives " << startN << endl;
        // cout << "Rank " << currRank << " receives " << endN << endl;
            
        /*cout << "Particle Location (Rank " << currRank << "): " << endl;
        for (unsigned int i = 0; i < N; ++i) {
            cout << x[2*i] << " " << x[2*i+1] << endl;
        }*/
            
        xLocal = new double[2 * localN]();
        F_pLocal = new double[2 * localN]();
            
        startN = currRank * localN;
        endN = (currRank + 1) * localN;
            
        // cout << "Start N (Rank 1): " << startN << endl;
        // cout << "End N (Rank 1): " << endN << endl;
            
        for (unsigned int i = 0; i < localN; ++i) {
                
            xLocal[2*i] = x[2*(i + startN)];
            xLocal[2*i+1] = x[2*(i + startN) + 1];
                
        }
            
        // cout << "Particle Location (Rank " << currRank << "): " << endl;
        /*for (unsigned int i = 0; i < localN; ++i) {
            cout << "(Rank " << currRank << ") " << "xLocal[" << i << "] (t = " << t << ") " << xLocal[2*i] << " " << xLocal[2*i+1] << endl;
        }*/
        
            
        double q;
    
        for (unsigned int i = 0; i < localN; ++i) {
        
            for (unsigned int j = 0; j < N; ++j) {
            
                calcRij(i, j);
                
                // cout << "(Rank " << currRank << ") R_ij at t = " << t << ": " << r_ij[2*j] << " " << r_ij[2*j+1] << endl;
                
                q = sqrt(r_ij[2*j]*r_ij[2*j] + r_ij[2*j+1]*r_ij[2*j+1]) / h;
                
                // cout << "(Rank " << currRank << ") q at t = " << t << ": " << q << endl;
                
                if (q < 1 && q != 0) {
                    
                    // cout << "(Rank: " << currRank << ") Rho[" << j << "] at t = " << t << ": " << rho[j] << endl;
                    
                    // cout << "(Rank: " << currRank << ") p[" << j << "] at t = " << t << ": " << p[j] << endl;
                    
                    F_pLocal[2*i] += ((m / rho[j]) * ((p[i+startN] + p[j]) / 2.0)) * ((-30 / (M_PI * h * h * h)) * (r_ij[2*j]) * (((1 - q) * (1 - q)) / q));
                
                    F_pLocal[2*i+1] += ((m / rho[j]) * ((p[i+startN] + p[j]) / 2.0)) * ((-30 / (M_PI * h * h * h)) * (r_ij[2*j+1]) * (((1 - q) * (1 - q)) / q));
                    
                    /*for (unsigned int i = 0; i < localN; ++i) {
                        cout << "(Rank " << currRank << ") " << "F_p[" << i << "] (t = " << t << ") " << F_p[2*i] << " " << F_p[2*i+1] << endl;
                    }*/
                    
                } else {
                    
                    F_pLocal[2*i] += 0.0;
                
                    F_pLocal[2*i+1] += 0.0;
                    
                }
            
            }
            
            F_pLocal[2*i] *= -1.0;
            F_pLocal[2*i+1] *= -1.0;
            
        }
            
    }
    
    MPI_Gather(F_pLocal, 2*localN, MPI_DOUBLE, F_p, 2*localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
    /*if (currRank == 0) {
        cout << "Pressure Force at t = " << t << endl;
        for (unsigned int i = 0; i < N; ++i) {
            cout << F_p[2*i] << " " << F_p[2*i+1] << endl;
        }
    }*/
}

void SPH::calcViscousForceWithMPI(int localN, unsigned int nProc, int currRank) {
    
    MPI_Bcast(x, 2*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(rho, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(v, 2*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
    if (currRank == 0){
            
        for (unsigned int i = 1; i < nProc; ++i) {
                
            startN = i * localN;
            endN = (i+1) * localN;
                
            // cout << "Start N: " << startN << endl;
            // cout << "End N: " << endN << endl;
                
            MPI_Send(&startN, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(&endN, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
                
            // cout << "Rank " << currRank << " sends " << startN << endl;
            // cout << "Rank " << currRank << " sends " << endN << endl;
             
        }
            
        /*cout << "Particle Location (Rank " << currRank << "): " << endl;
        for (unsigned int i = 0; i < N; ++i) {
            cout << x[2*i] << " " << x[2*i+1] << endl;
        }*/
            
        xLocal = new double[2 * localN]();
        vLocal = new double[2 * localN]();
        F_vLocal = new double[2 * localN]();
            
        startN = currRank * localN;
        endN = (currRank + 1) * localN;
            
        // cout << "Start N (Rank 0): " << startN << endl;
        // cout << "End N (Rank 0): " << endN << endl;
            
        for (unsigned int i = 0; i < localN; ++i) {
                
            xLocal[2*i] = x[2*(i + startN)];
            xLocal[2*i+1] = x[2*(i + startN) + 1];
            
            vLocal[2*i] = v[2*(i + startN)];
            vLocal[2*i+1] = v[2*(i + startN) + 1];
                
        }
            
        // cout << "Particle Location (Rank " << currRank << "): " << endl;
        /*for (unsigned int i = 0; i < localN; ++i) {
            cout << "(Rank " << currRank << ") " << "xLocal[" << i << "] (t = " << t << ") " << xLocal[2*i] << " " << xLocal[2*i+1] << endl;
            cout << "(Rank " << currRank << ") " << "vLocal[" << i << "] (t = " << t << ") " << vLocal[2*i] << " " << vLocal[2*i+1] << endl;
        }*/
            
        double q;
    
        for (unsigned int i = 0; i < localN; ++i) {
        
            for (unsigned int j = 0; j < N; ++j) {
            
                calcVij(i, j);
                
                calcRij(i, j);
                
                q = sqrt(r_ij[2*j]*r_ij[2*j] + r_ij[2*j+1]*r_ij[2*j+1]) / h;
                
                if (q < 1 && q != 0) {
                    
                    F_vLocal[2*i] += (m / rho[j]) * (v_ij[2*j]) * ((40 / (M_PI * h * h * h * h)) * (1 - q));
                    
                    F_vLocal[2*i+1] += (m / rho[j]) * (v_ij[2*j+1]) * ((40 / (M_PI * h * h * h * h)) * (1 - q));
                    
                } else {
                    
                    F_vLocal[2*i] += 0.0;
                    
                    F_vLocal[2*i+1] += 0.0;
                    
                }
            
            }
            
            F_vLocal[2*i] *= -mu;
                    
            F_vLocal[2*i+1] *= -mu;
            
        }
        
        // for (unsigned int i = 0; i < localN; ++i) {
        //     cout << "(Rank " << currRank << ") " << "F_vLocal[" << i << "] (t = " << t << ") " << F_vLocal[2*i] << " " << F_vLocal[2*i+1] << endl;
        // }
        
    } else {
        
        MPI_Recv(&startN, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&endN, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
        // cout << "Rank " << currRank << " receives " << startN << endl;
        // cout << "Rank " << currRank << " receives " << endN << endl;
            
        /*cout << "Particle Location (Rank " << currRank << "): " << endl;
        for (unsigned int i = 0; i < N; ++i) {
            cout << x[2*i] << " " << x[2*i+1] << endl;
        }*/
            
        xLocal = new double[2 * localN]();
        vLocal = new double[2 * localN]();
        F_vLocal = new double[2 * localN]();
            
        startN = currRank * localN;
        endN = (currRank + 1) * localN;
            
        // cout << "Start N (Rank 1): " << startN << endl;
        // cout << "End N (Rank 1): " << endN << endl;
            
        for (unsigned int i = 0; i < localN; ++i) {
                
            xLocal[2*i] = x[2*(i + startN)];
            xLocal[2*i+1] = x[2*(i + startN) + 1];
            
            vLocal[2*i] = v[2*(i + startN)];
            vLocal[2*i+1] = v[2*(i + startN) + 1];
                
        }
            
        // cout << "Particle Location (Rank " << currRank << "): " << endl;
        /*for (unsigned int i = 0; i < localN; ++i) {
            cout << "(Rank " << currRank << ") " << "xLocal[" << i << "] (t = " << t << ") " << xLocal[2*i] << " " << xLocal[2*i+1] << endl;
        }*/
        
            
        double q;
    
        for (unsigned int i = 0; i < localN; ++i) {
        
            for (unsigned int j = 0; j < N; ++j) {
            
                calcVij(i, j);
                
                calcRij(i, j);
                
                q = sqrt(r_ij[2*j]*r_ij[2*j] + r_ij[2*j+1]*r_ij[2*j+1]) / h;
                
                if (q < 1 && q != 0) {
                    
                    F_vLocal[2*i] += (m / rho[j]) * (v_ij[2*j]) * ((40 / (M_PI * h * h * h * h)) * (1 - q));
                    
                    F_vLocal[2*i+1] += (m / rho[j]) * (v_ij[2*j+1]) * ((40 / (M_PI * h * h * h * h)) * (1 - q));
                    
                } else {
                    
                    F_vLocal[2*i] += 0.0;
                    
                    F_vLocal[2*i+1] += 0.0;
                    
                }
            
            }
            
            F_vLocal[2*i] *= -mu;
                    
            F_vLocal[2*i+1] *= -mu;
            
        }
            
    }
    
    MPI_Gather(F_vLocal, 2*localN, MPI_DOUBLE, F_v, 2*localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
    /*if (currRank == 0) {
        cout << "Viscous Force at t = " << t << endl;
        for (unsigned int i = 0; i < N; ++i) {
            cout << F_v[2*i] << " " << F_v[2*i+1] << endl;
        }
    }*/
}

void SPH::calcGravityForceWithMPI(int localN, unsigned int nProc, int currRank) {
    
    MPI_Bcast(rho, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
    if (currRank == 0){
            
        for (unsigned int i = 1; i < nProc; ++i) {
                
            startN = i * localN;
            endN = (i+1) * localN;
                
            // cout << "Start N: " << startN << endl;
            // cout << "End N: " << endN << endl;
                
            MPI_Send(&startN, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(&endN, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
                
            // cout << "Rank " << currRank << " sends " << startN << endl;
            // cout << "Rank " << currRank << " sends " << endN << endl;
             
        }
            
        /*cout << "Particle Location (Rank " << currRank << "): " << endl;
        for (unsigned int i = 0; i < N; ++i) {
            cout << x[2*i] << " " << x[2*i+1] << endl;
        }*/
            
        F_gLocal = new double[2 * localN]();
            
        startN = currRank * localN;
        endN = (currRank + 1) * localN;
            
        // cout << "Start N (Rank 0): " << startN << endl;
        // cout << "End N (Rank 0): " << endN << endl;
            
        for (unsigned int i = 0; i < localN; ++i) {
        
            F_gLocal[2*i] = 0.0;
            
            F_gLocal[2*i+1] = -rho[i+startN] * g;
            
        }
        
        // for (unsigned int i = 0; i < localN; ++i) {
        //     cout << "(Rank " << currRank << ") " << "F_gLocal[" << i << "] (t = " << t << ") " << F_gLocal[2*i] << " " << F_gLocal[2*i+1] << endl;
        // }
        
    } else {
        
        MPI_Recv(&startN, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&endN, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
        // cout << "Rank " << currRank << " receives " << startN << endl;
        // cout << "Rank " << currRank << " receives " << endN << endl;
            
        /*cout << "Particle Location (Rank " << currRank << "): " << endl;
        for (unsigned int i = 0; i < N; ++i) {
            cout << x[2*i] << " " << x[2*i+1] << endl;
        }*/
            
        F_gLocal = new double[2 * localN]();
            
        startN = currRank * localN;
        endN = (currRank + 1) * localN;
            
        // cout << "Start N (Rank 1): " << startN << endl;
        // cout << "End N (Rank 1): " << endN << endl;
            
        for (unsigned int i = 0; i < localN; ++i) {
        
            F_gLocal[2*i] = 0.0;
            
            F_gLocal[2*i+1] = -rho[i+startN] * g;
            
        }
            
    }
    
    MPI_Gather(F_gLocal, 2*localN, MPI_DOUBLE, F_g, 2*localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
    /*if (currRank == 0) {
        cout << "Gravity Force at t = " << t << endl;
        for (unsigned int i = 0; i < N; ++i) {
            cout << F_g[2*i] << " " << F_g[2*i+1] << endl;
        }
    }*/
}

void SPH::calcAcceleration() {
    
    for (unsigned int i = 0; i < N; ++i) {
        
        a[2*i] = (F_p[2*i] + F_v[2*i] + F_g[2*i]) / rho[i];
        
        a[2*i+1] = (F_p[2*i+1] + F_v[2*i+1] + F_g[2*i+1]) / rho[i];
        
    }
    
    
    /*cout << "Acceleration at t = " << t << endl;
    for (unsigned int i = 0; i < N; ++i) {
        cout << a[2*i] << " " << a[2*i+1] << endl;
    }*/
    
    
}

void SPH::generateVInit() {
    
    for (unsigned int i = 0; i < N; ++i) {
        
        v[2*i] = 0 + (a[2*i] * (dt / 2));
        
        v[2*i+1] = 0 + (a[2*i+1] * (dt / 2));
        
    }
    
}

void SPH::getNextParticleVel() {
    
    for (unsigned int i = 0; i < N; ++i) {
        
        v[2*i] += a[2*i] * dt;
        
        v[2*i+1] += a[2*i+1] * dt;
        
    }
    
    
    /*cout << "Velocity at t = " << t << endl;
    for (unsigned int i = 0; i < N; ++i) {
        cout << v[2*i] << " " << v[2*i+1] << endl;
    }*/
    
    
}

void SPH::getNextParticlePos() {
    
    for (unsigned int i = 0; i < N; ++i) {
        
        x[2*i] += v[2*i] * dt;
        
        x[2*i+1] += v[2*i+1] * dt;
        
        applyBC(i);
        
    }
    
    /*cout << "Position at t = " << t << endl;
    for (unsigned int i = 0; i < N; ++i) {
        cout << x[2*i] << " " << x[2*i+1] << endl;
    }*/
    
}

void SPH::applyBC(unsigned int i) {
    
    if (x[2*i] < h) {
        
        x[2*i] = h;
        
        v[2*i] *= -e;
        
    } else if (x[2*i] > (1 - h)) {
        
        x[2*i] = 1 - h;
        
        v[2*i] *= -e;
            
    }
        
    if (x[2*i+1] < h) {
            
        x[2*i+1] = h;
            
        v[2*i+1] *= -e; 
            
    } else if (x[2*i+1] > (1 - h)) {
            
        x[2*i+1] = 1 - h;
            
        v[2*i+1] *= -e;
            
    }
    
}

void SPH::writeParticlePosition(ofstream& xCoor) {
    
    if (xCoor.is_open()) {
        
        for (unsigned int i = 0; i < N; ++i) {
            
            xCoor << x[2*i] << " " << x[2*i+1] << endl;
            
        }
        
    }
    
}

void SPH::writeEnergy(ofstream& energyTxt) {
    
    if (energyTxt.is_open()) {
        
        energyTxt << getCurrT() << " " << calcKineticEnergy() << " " << calcPotentialEnergy() << " " << calcTotalEnergy() << endl;
        
    }

}

void SPH::writeData(ofstream& dataTxt) {
    
    if (dataTxt.is_open()) {
        
        for (unsigned int i = 0; i < N; ++i) {
            
            dataTxt << x[2*i] << " " << x[2*i+1] << " ";
            
        }
        
        dataTxt << endl;
        
    }
    
}
  
void SPH::closePPOutputFile(ofstream& fileName) {
    
    fileName.close();
    
}

double SPH::calcKineticEnergy() {
    
    double E_k = 0.0;
    
    for (unsigned int i = 0; i < N; ++i) {
        
        E_k += (v[2*i] * v[2*i]) + v[2*i+1] * v[2*i+1];
        
    }
    
    return 0.5 * m * E_k;
    
}

double SPH::calcPotentialEnergy() {
    
    double E_p = 0.0;
    
    for (unsigned int i = 0; i < N; ++i) {
        
        E_p += x[2*i+1];
        
    }
    
    return m * g * E_p;
    
}

double SPH::calcTotalEnergy() {
    
    return calcKineticEnergy() + calcPotentialEnergy();
    
}