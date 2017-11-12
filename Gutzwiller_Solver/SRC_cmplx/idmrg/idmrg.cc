#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include "itensor/all.h"

using namespace itensor;

int main(int argc, char* argv[])
{
    // Parse the input file
    if(argc != 2) 
    {
        printfln("Usage: %s impurity_index",argv[0]);
        return  0;
    }
    std::string simp = argv[1];
    auto input = InputGroup("GDMRG_"+simp+".CTL","input");
    auto N = input.getInt("N");

    // number of particles, default is N   (half filling)
    auto Npart = input.getInt("Npart",N);
    auto Nsweeps = input.getInt("Nsweeps",50);
    auto maxM = input.getInt("maxM",50);
    auto conserveSz = input.getYesNo("ConserveSz",false);
    auto quiet = input.getYesNo("quiet",false);

    auto sites = Hubbard(N,{"ConserveSz", conserveSz});
    auto ampo = AutoMPO(sites);
   
    // input stream
    std::ifstream infile;

    // output file
    std::ofstream oufile("GDMRG_"+simp+".OUT"); 
    
    // Read two-body part
    infile.open("V2E_"+simp+".INP");
    std::string line; 
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        int m1, m2, m3, m4;
        double coulr, couli;
        if (!(iss >> m1 >> m2 >> m3 >> m4 >> coulr >> couli)) 
        {
            break;
        }
        std::string op1 = (m1%2==0) ? "Cdagdn": "Cdagup";
        std::string op2 = (m2%2==0) ? "Cdagdn": "Cdagup";
        std::string op3 = (m3%2==0) ? "Cdn": "Cup";
        std::string op4 = (m4%2==0) ? "Cdn": "Cup";
        ampo += Cplx(coulr, couli), 
             op1, (m1+1)/2, op2, (m2+1)/2, 
             op3, (m3+1)/2, op4, (m4+1)/2;
    }
    infile.close();
  
    // Read one-body part.
    infile.open("H1E_"+simp+".INP");
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        int m1, m2;
        double evalr, evali;
        if (!(iss >> m1 >> m2 >> evalr >> evali)) 
        {
            break;
        }
        std::string op1 = (m1%2==0) ? "Cdagdn": "Cdagup";
        std::string op2 = (m2%2==0) ? "Cdn": "Cup";
        ampo += Cplx(evalr, evali), op1, (m1+1)/2, op2, (m2+1)/2;
    }
    infile.close();
  
    auto H = IQMPO(ampo);
  
    // set initial state
    auto state = InitState(sites);
    for (int i=1; i<=N; ++i) 
    {
        if(i%2==1)
        {
            state.set(i,"Up");
        }
        else
        {
            state.set(i,"Dn");
        }
    }
  
    auto psi = IQMPS(state);
    Print(totalQN(psi));
  
    // setup dmrg calculations
    auto sweeps = Sweeps(Nsweeps);
    sweeps.minm() = 10, 20, 50;
    sweeps.maxm() = 50, maxM;

    // quite delicate, keep it.
    sweeps.cutoff() = 1E-12;
    sweeps.noise() = 1E-5, 1E-6, 1E-10, 1E-12;
  
    // Begin the DMRG calculation
    auto energy = dmrg(psi,H,sweeps,{"Quiet", quiet});
  
    // record final energy by DMRG
    oufile << "#e_embed: " << std::setprecision(12) << energy << std::endl;

    // calculate density matrix
    std::vector<std::vector<Cplx>> dm(2*N, std::vector<Cplx>(2*N, 0.));

    // case 1: same site, same spin 
    for(int j = 1; j <= N; ++j)
    {
        psi.position(j);
        auto up = (dag(prime(psi.A(j),Site))*sites.op("Nup",j)*psi.A(j))
                .cplx();
        auto dn = (dag(prime(psi.A(j),Site))*sites.op("Ndn",j)*psi.A(j))
                .cplx();
        dm[2*j-2][2*j-2] = up;
        dm[2*j-1][2*j-1] = dn;
     }

    // case 2: same site, different spin
    for (int j=1; j <= N; ++j)
    {
        auto s = sites(j);
        auto sp = prime(s);
        IQTensor Op(dag(s), sp);
        Op.set(s(3), sp(2), 1);
        psi.position(j);
        auto res = (dag(prime(psi.A(j), Site))*Op*psi.A(j)).cplx();
        dm[2*j-2][2*j-1] = res;
        dm[2*j-1][2*j-2] = res;
    }

    // case 3: different sites
    for (int i=1; i < N; ++i)
    {
        for (int ii=0; ii <= 1; ++ii)
        {
            for (int j=i+1; j <= N; ++j)
            {
                for (int jj=0; jj <= 1; ++jj)
                {    
                    IQTensor A_i, A_j;
                    int orb_i, orb_j;
                    if (ii == 0) 
                    {
                        A_i = sites.op("Adagup",i);
                        orb_i = 2 * i - 1;
                    } 
                    else 
                    {
                        A_i = sites.op("Adagdn",i);
                        orb_i = 2 * i;
                    } 
                    if (jj == 0) 
                    { 
                        A_j = sites.op("Aup",j);
                        orb_j = 2 * j - 1 ;
                    } 
                    else 
                    {
                        A_j = sites.op("Adn",j);
                        orb_j = 2 * j;
                    }
                    psi.position(i); 
                    auto ir = commonIndex(psi.A(i),psi.A(i+1),Link);
                    auto first_F = sites.op("F",i);
                    auto Corr = dag(prime(prime(psi.A(i),ir, Site), Site))
                            * prime(A_i);
                    if (ii == 0) 
                    {
                        Corr *= first_F;
                        Corr *= psi.A(i);
                    } 
                    else 
                    {
                        Corr *= prime(psi.A(i),Site);
                    }
                    
                    for(int k = i+1; k < j; ++k)
                    {
                        Corr *= psi.A(k);
                        Corr *= sites.op("F",k);
                        Corr *= dag(prime(psi.A(k)));
                    }
                    Corr *= psi.A(j);
                    Corr *= A_j;
                    auto jl = commonIndex(psi.A(j),psi.A(j-1),Link);
  
                    auto last_F = sites.op("F",j); 
                    if (jj == 1) 
                    {
                        Corr *= prime(last_F);
                        Corr *= dag(prime(prime(psi.A(j),jl,Site), Site));
                    } 
                    else 
                    {
                        Corr *= dag(prime(psi.A(j),jl, Site));
                    }
                    auto res = Corr.cplx();
                    dm[orb_i-1][orb_j-1] = res;
                    dm[orb_j-1][orb_i-1] = conj(res);
                }
            }
        }
    }

    for (int i=0; i < 2*N; ++i)
    {
        for (int j=0; j < 2*N; ++j)
        {
            oufile << std::setw(20) << std::setprecision(12) 
                    << dm[i][j].real() 
                    << std::setw(20) << std::setprecision(12) 
                    << dm[i][j].imag(); 
        }
        oufile << std::endl;
    }
    oufile.close();

    return 0;
}
