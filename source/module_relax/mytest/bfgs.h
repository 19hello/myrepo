#ifndef BFGS_H
#define BFGS_H

#include <vector>
#include <tuple> // 确保包含 <tuple> 头文件
#include<algorithm>
#include<cmath>
#include"module_base/lapack_connector.h"

class bfgs
{
public:
    // 公有成员变量
    double alpha;
    int maxstep;
    int size;
    bool sign;
    std::vector<double> steplength;
    std::vector<std::vector<double>> H;
    std::vector<double> force0;
    std::vector<std::vector<double>> force;
    std::vector<double> pos0;
    std::vector<std::vector<double>> pos;
    std::vector<std::vector<double>> dpos;


    // 成员函数声明
    void initialize(int _size);
    bool Step(std::vector<std::vector<double>> _force);
    void PrepareStep();
    void Update(std::vector<double> pos, std::vector<double> force);
    void DetermineStep();

    // 矩阵方法声明
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
