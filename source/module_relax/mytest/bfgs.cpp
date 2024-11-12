#include "bfgs.h"
#include "module_base/matrix3.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"

void bfgs::initialize(int _size) // 初始化H0、H、pos0、force0、force
{
    alpha=70;//relax_scale_force
    maxstep=100; //relax_nmax
    size=_size;
    sign =true;
    H = std::vector<std::vector<double>>(3*size, std::vector<double>(3*size, 0.0));
    for (int i = 0; i < 3*size; ++i) {
        H[i][i] = alpha;  
    }

    pos = std::vector<std::vector<double>> (size, std::vector<double>(3, 0.0));


    pos[0][0]=0;
    pos[0][1]=0;
    pos[0][2]=0;
    pos[1][0]=1.6208698;
    pos[1][1]=1.18863785;
    pos[1][2]=1.45878282;


    pos0 = std::vector<double>(3*size, 0.0);
    dpos = std::vector<std::vector<double>>(size, std::vector<double>(3, 0.0));
    force0 = std::vector<double>(3*size, 0.0);
    force = std::vector<std::vector<double>>(size, std::vector<double>(3, 0.0));
    steplength = std::vector<double>(size, 0.0);
}

bool bfgs::Step(std::vector<std::vector<double>> _force) // 正式开始每一步迭代
{
    std::cout<<"enter Step"<<std::endl;
    GlobalC::ucell.ionic_position_updated = true;
    force = _force;
    //force={{1.83632107,-0.79570368,1.5420276},{-1.71118338,0.67709389,-1.41739503}};

    PrepareStep();
    /*std::cout<<"printdpos"<<std::endl;
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<3;j++)
        {
            std::cout<<dpos[i][j]<<' ';
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;*/
    DetermineStep();
    /*std::cout<<"printdpos"<<std::endl;
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<3;j++)
        {
            std::cout<<dpos[i][j]<<' ';
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;*/
    double move_ion[3*size];
    ModuleBase::zeros(move_ion, size*3);

    for(int iat=0; iat<size; iat++)
    {
        //Cartesian coordinate
        //convert from Angstrom to unit of latvec (Bohr)

        //单位转换
        ModuleBase::Vector3<double> move_ion_cart;
        move_ion_cart.x = dpos[iat][0] / ModuleBase::BOHR_TO_A / GlobalC::ucell.lat0;
        move_ion_cart.y = dpos[iat][1] / ModuleBase::BOHR_TO_A / GlobalC::ucell.lat0;
        move_ion_cart.z = dpos[iat][2] / ModuleBase::BOHR_TO_A / GlobalC::ucell.lat0;
        /*GlobalC::ucell.atoms[0].tau[iat].x+=move_ion_cart.x;
        GlobalC::ucell.atoms[0].tau[iat].y+=move_ion_cart.y;
        GlobalC::ucell.atoms[0].tau[iat].z+=move_ion_cart.z;*/




        /*move_ion_cart.x = dpos[iat][0] ;
        move_ion_cart.y = dpos[iat][1] ;
        move_ion_cart.z = dpos[iat][2] ;*/
        /*std::cout<<"move_ion_cart"<<std::endl;
        std::cout<<move_ion_cart.x<<std::endl;
        std::cout<<move_ion_cart.y<<std::endl;
        std::cout<<move_ion_cart.z<<std::endl;*/



        //convert to Direct coordinate
        //note here the old GT is used

        //坐标转换
        ModuleBase::Vector3<double> move_ion_dr = move_ion_cart * GlobalC::ucell.GT;
        /*std::cout<<"move_ion_dr"<<std::endl;
        std::cout<<move_ion_dr.x<<std::endl;
        std::cout<<move_ion_dr.y<<std::endl;
        std::cout<<move_ion_dr.z<<std::endl;*/

        int it = GlobalC::ucell.iat2it[iat];
        int ia = GlobalC::ucell.iat2ia[iat];
        Atom* atom = &GlobalC::ucell.atoms[it];

        if(atom->mbl[ia].x == 1)
        {
            move_ion[iat * 3] = move_ion_dr.x;
        }
        if(atom->mbl[ia].y == 1)
        {
            move_ion[iat * 3 + 1] = move_ion_dr.y ;
        }
        if(atom->mbl[ia].z == 1)
        {
            move_ion[iat * 3 + 2] = move_ion_dr.z ;
        }
    }
    /*std::cout<<"move_ion"<<std::endl;
    for(int i=0;i<3*size;i++)
    {
        std::cout<<move_ion[i]<<' ';
    }
    std::cout<<std::endl;


    std::cout<<"printtau"<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].tau[0].x<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].tau[0].y<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].tau[0].z<<std::endl;
    std::cout<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].taud[0].x<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].taud[0].y<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].taud[0].z<<std::endl;
    std::cout<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].tau[1].x<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].tau[1].y<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].tau[1].z<<std::endl;
    std::cout<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].taud[1].x<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].taud[1].y<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].taud[1].z<<std::endl;
    std::cout<<std::endl;*/


	GlobalC::ucell.update_pos_taud(move_ion);
    std::cout<<__LINE__<<":"<<GlobalC::ucell.ionic_position_updated<<std::endl;

	// Print the struct
    /*for(int i=0;i<size;i++)
    {
        for(int j=0;j<3;j++)
        {
            move_ion[3*i+j]=(dpos[i][j]);
        }
    }*/



    /*std::cout<<"enterhere"<<std::endl;

    std::cout<<"printtau"<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].tau[0].x<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].tau[0].y<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].tau[0].z<<std::endl;
    std::cout<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].taud[0].x<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].taud[0].y<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].taud[0].z<<std::endl;
    std::cout<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].tau[1].x<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].tau[1].y<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].tau[1].z<<std::endl;
    std::cout<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].taud[1].x<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].taud[1].y<<std::endl;
    std::cout<<GlobalC::ucell.atoms[0].taud[1].z<<std::endl;*/

    pos = MAddM(pos, dpos);
    /*std::cout<<"printpos"<<std::endl;
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<3;j++)
        {
            std::cout<<pos[i][j]<<' ';
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;*/
    double a=0;
    for(int i=0;i<size;i++)
    {
        for(int j=0;j<3;j++)
        {
            double w;
            if(dpos[i][j]>0)
            {
                w=dpos[i][j];
            }
            else
            {
                w=-dpos[i][j];
            }
            if(w>a)
            {
                a=w;
            }
        }
    }
    std::cout<<a<<std::endl;
    bool k=a<0.00001;
    std::cout<<k<<std::endl;
    return k;
}

void bfgs::PrepareStep()
{
    std::vector<double> changedforce = ReshapeMToV(force);
    std::vector<double> changedpos = ReshapeMToV(pos);

    ///输出changedforce和changedpos
    std::cout<<"PrintChangedforce"<<std::endl;
    for(int i=0;i<changedforce.size();i++)
    {
        std::cout<<changedforce[i]<<' ';
    }
    std::cout<<std::endl;
    std::cout<<"PrintChangedpos"<<std::endl;
    for(int i=0;i<changedpos.size();i++)
    {
        std::cout<<changedpos[i]<<' ';
    }
    std::cout<<std::endl;

    //输出H
    std::cout<<"PrintH"<<std::endl;
    for(int i=0;i<3*size;i++)
    {
        for(int j=0;j<3*size;j++)
        {
            std::cout<<H[i][j]<<' ';
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;


    //更新H
    Update(changedpos, changedforce);

    //输出更新后的H
    std::cout<<"UpdateH"<<std::endl;
    for(int i=0;i<3*size;i++)
    {
        for(int j=0;j<3*size;j++)
        {
            std::cout<<H[i][j]<<' ';
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;


    //调用dysev
    std::vector<double> omega(3*size);
    std::vector<double> work(3*size*3*size);
    int lwork=3*size*3*size;
    int info;
    std::vector<double> H_flat;
    for(const auto& row : H)
    {
        H_flat.insert(H_flat.end(), row.begin(), row.end());
    }   
    int value=3*size;
    int* ptr=&value;
    dsyev_("V","U",ptr,H_flat.data(),ptr,omega.data(),work.data(),&lwork,&info);


    //输出V，omega

    std::vector<std::vector<double>> V(3*size, std::vector<double>(3*size, 0.0));
    for(int i = 0; i < 3*size; i++)
    {
        for(int j = 0; j < 3*size; j++)
        {
            V[j][i] = H_flat[3*size*i + j];
        }
    }

    std::cout<<"PrintV"<<std::endl;
    for(int i=0;i<3*size;i++)
    {
        for(int j=0;j<3*size;j++)
        {
            std::cout<<V[i][j]<<' ';
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    std::cout<<"PrintOmega"<<std::endl;
    for(int i=0;i<3*size;i++)
    {
        std::cout<<omega[i]<<' ';
    }
    std::cout<<std::endl;



    std::vector<double> a=DotInMAndV2(V, changedforce);



    std::cout<<"Printa"<<std::endl;
    for(int i=0;i<a.size();i++)
    {
        std::cout<<a[i]<<' ';
    }
    std::cout<<std::endl;

    for(int i = 0; i < a.size(); i++)
    {
        a[i] /= abs(omega[i]);
    }

    std::cout<<"Printa"<<std::endl;
    for(int i=0;i<a.size();i++)
    {
        std::cout<<a[i]<<' ';
    }
    std::cout<<std::endl;

    std::vector<double> tmpdpos = DotInMAndV1(V, a);


    //输出tmpdpos
    std::cout<<"PrintTmpdpos"<<std::endl;
    for(int i=0;i<tmpdpos.size();i++)
    {
        std::cout<<tmpdpos[i]<<' ';
    }
    std::cout<<std::endl;


    dpos = ReshapeVToM(tmpdpos);
    for(int i = 0; i < size; i++)
    {
        double k = 0;
        for(int j = 0; j < 3; j++)
        {
            k += dpos[i][j] * dpos[i][j];
        }
        steplength[i] = sqrt(k);
    }
    pos0 = ReshapeMToV(pos);
    force0 = ReshapeMToV(force);
}

void bfgs::Update(std::vector<double> pos, std::vector<double> force)
{
    if(sign)
    {
        sign=false;
        return;
    }
    //pos0={0,0,0,1.6208698,1.18863785,1.45878282};
    std::vector<double> dpos = VSubV(pos, pos0);
    if(*max_element(dpos.begin(), dpos.end()) < 1e-7)
    {
        return;
    }
    //force0={3.94885101,-2.46796196,3.5609503,-3.82079544,2.33869184,-3.43330592};
    std::vector<double> dforce = VSubV(force, force0);
    double a = DotInVAndV(dpos, dforce);
    std::vector<double> dg = DotInMAndV1(H, dpos);
    double b = DotInVAndV(dpos, dg);
    H = MSubM(H, MPlus(OuterVAndV(dforce, dforce), a));
    H = MSubM(H, MPlus(OuterVAndV(dg, dg), b));
}

void bfgs::DetermineStep()
{
    auto maxsteplength = max_element(steplength.begin(), steplength.end());
    double a = *maxsteplength;
    if(a >= maxstep)
    {
        double scale = maxstep / a;
        for(int i = 0; i < size; i++)
        {
            for(int j=0;j<3;j++)
            {
                dpos[i][j]*=scale;
            }
        }
    }
}

// 矩阵方法

std::vector<double> bfgs::ReshapeMToV(std::vector<std::vector<double>> matrix) // 将n*3矩阵变成3n的数组
{
    int size = matrix.size();
    std::vector<double> result;
    result.reserve(3*size);
    for (const auto& row : matrix) {
        result.insert(result.end(), row.begin(), row.end());
    }
    return result;
}

std::vector<std::vector<double>> bfgs::MAddM(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b) // 矩阵相加
{
    std::vector<std::vector<double>> result = std::vector<std::vector<double>>(a.size(), std::vector<double>(a[0].size(), 0.0));
    for(int i = 0; i < a.size(); i++)
    {
        for(int j = 0; j < a[0].size(); j++)
        {
            result[i][j] = a[i][j] + b[i][j];
        }
    }
    return result;
}

std::vector<double> bfgs::VSubV(std::vector<double> a, std::vector<double> b) // 向量相减
{
    std::vector<double> result = std::vector<double>(a.size(), 0.0);
    for(int i = 0; i < a.size(); i++)
    {
        result[i] = a[i] - b[i];
    }
    return result;
}

std::vector<std::vector<double>> bfgs::ReshapeVToM(std::vector<double> matrix) // 将3n的数组变成n*3矩阵
{
    std::vector<std::vector<double>> result = std::vector<std::vector<double>>(matrix.size() / 3, std::vector<double>(3));
    for(int i = 0; i < result.size(); i++)
    {
        for(int j = 0; j < 3; j++)
        {
            result[i][j] = matrix[i*3 + j];
        }
    }
    return result;
}

std::vector<double> bfgs::DotInMAndV1(std::vector<std::vector<double>> matrix, std::vector<double> vec) // 矩阵与向量相乘
{
    std::vector<double> result(matrix.size(), 0.0);
    for(int i = 0; i < result.size(); i++)
    {
        for(int j = 0; j < vec.size(); j++)
        {
            result[i] += matrix[i][j] * vec[j];
        }
    }
    return result;
}
std::vector<double> bfgs::DotInMAndV2(std::vector<std::vector<double>> matrix, std::vector<double> vec) // 矩阵与向量相乘
{
    std::vector<double> result(matrix.size(), 0.0);
    for(int i = 0; i < result.size(); i++)
    {
        for(int j = 0; j < vec.size(); j++)
        {
            result[i] += matrix[j][i] * vec[j];
        }
    }
    return result;
}

double bfgs::DotInVAndV(std::vector<double> vec1, std::vector<double> vec2) // 向量与向量点乘
{
    double result = 0.0;
    for(int i = 0; i < vec1.size(); i++)
    {
        result += vec1[i] * vec2[i];
    }
    return result;
}

std::vector<std::vector<double>> bfgs::OuterVAndV(std::vector<double> a, std::vector<double> b) // 两个向量外积
{
    std::vector<std::vector<double>> result = std::vector<std::vector<double>>(a.size(), std::vector<double>(b.size(), 0.0));
    for(int i = 0; i < a.size(); i++)
    {
        for(int j = 0; j < b.size(); j++)
        {
            result[i][j] = a[i] * b[j];
        }
    }
    return result;
}

std::vector<std::vector<double>> bfgs::MPlus(std::vector<std::vector<double>> a, double b) // 矩阵除以常数
{
    std::vector<std::vector<double>> result = std::vector<std::vector<double>>(a.size(), std::vector<double>(a[0].size(), 0.0));
    for(int i = 0; i < a.size(); i++)
    {
        for(int j = 0; j < a[0].size(); j++)
        {
            result[i][j] = a[i][j] / b;
        }
    }
    return result;
}

std::vector<std::vector<double>> bfgs::MSubM(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b) // 两个矩阵相减
{
    std::vector<std::vector<double>> result = std::vector<std::vector<double>>(a.size(), std::vector<double>(a[0].size(), 0.0));
    for(int i = 0; i < a.size(); i++)
    {
        for(int j = 0; j < a[0].size(); j++)
        {
            result[i][j] = a[i][j] - b[i][j];
        }
    }
    return result;
}

std::tuple<std::vector<double>, std::vector<std::vector<double>>> bfgs::GetEigenvalueAndEigenVector(std::vector<std::vector<double>> matrix) // 求特征值和特征向量
{
    std::vector<double> omega;
    std::vector<std::vector<double>> V;

    return std::make_tuple(omega, V);
}
