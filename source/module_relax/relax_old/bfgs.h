#ifndef BFGS_H
#define BFGS_H

#include <vector>
#include <tuple> 
#include<algorithm>
#include<cmath>
#include"module_base/lapack_connector.h"

#include "module_base/matrix.h"
#include "module_base/matrix3.h"
#include "module_cell/unitcell.h"



class BFGS
{
public:
    
    double alpha;//initialize H,diagonal element is alpha
    double maxstep;//every movement smaller than maxstep
    int size;//number of etoms
    
    
    std::vector<double> steplength;
    std::vector<std::vector<double>> H;
    std::vector<double> force0;
    std::vector<std::vector<double>> force;
    std::vector<double> pos0;
    std::vector<std::vector<double>> pos;
    std::vector<std::vector<double>> dpos;

    void allocate(const int _size);//initialize parameters
    void relax_step(ModuleBase::matrix _force,UnitCell& ucell);//
    void PrepareStep(std::vector<std::vector<double>>& force,std::vector<std::vector<double>>& pos,std::vector<std::vector<double>>& H,std::vector<double>& pos0,std::vector<double>& force0,std::vector<double>& steplength,std::vector<std::vector<double>>& dpos);
    void IsRestrain(std::vector<std::vector<double>>& dpos);

private:
    bool sign;
    
    
    void GetPos(UnitCell& ucell,std::vector<std::vector<double>>& pos);
    void Update(std::vector<double> pos, std::vector<double> force,std::vector<std::vector<double>>& H);
    void DetermineStep(std::vector<double> steplength,std::vector<std::vector<double>>& dpos,double maxstep);
    void UpdatePos(UnitCell& ucell);
    

    // matrix method
    std::vector<double> ReshapeMToV(std::vector<std::vector<double>> matrix);
    std::vector<std::vector<double>> MAddM(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b);
    std::vector<double> VSubV(std::vector<double> a, std::vector<double> b);
    std::vector<std::vector<double>> ReshapeVToM(std::vector<double> matrix);
    std::vector<double> DotInMAndV1(std::vector<std::vector<double>> matrix, std::vector<double> vec);
    std::vector<double> DotInMAndV2(std::vector<std::vector<double>> matrix, std::vector<double> vec);
    double DotInVAndV(std::vector<double> vec1, std::vector<double> vec2);
    std::vector<std::vector<double>> OuterVAndV(std::vector<double> a, std::vector<double> b);
    std::vector<std::vector<double>> MPlus(std::vector<std::vector<double>> a, double b);
    std::vector<std::vector<double>> MSubM(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b);
    
    std::tuple<std::vector<double>, std::vector<std::vector<double>>> GetEigenvalueAndEigenVector(std::vector<std::vector<double>> matrix);
};

#endif // BFGS_H
