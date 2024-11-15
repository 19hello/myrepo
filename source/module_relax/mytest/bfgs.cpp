#include "bfgs.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/matrix3.h"


void bfgs::initialize(int _size) // initialize H0、H、pos0、force0、force
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
    pos0 = std::vector<double>(3*size, 0.0);
    dpos = std::vector<std::vector<double>>(size, std::vector<double>(3, 0.0));
    force0 = std::vector<double>(3*size, 0.0);
    force = std::vector<std::vector<double>>(size, std::vector<double>(3, 0.0));
    steplength = std::vector<double>(size, 0.0);
    GetPos();
    
}

bool bfgs::Step(std::vector<std::vector<double>> _force,UnitCell& ucell) 
{
    std::cout<<"enter Step"<<std::endl;
    GlobalC::ucell.ionic_position_updated = true;
    force = _force;
    
    //GetPos(ucell);
    PrepareStep();
    DetermineStep();
    /*for(int i=0;i<size;i++)
    {
        for(int j=0;j<3;j++)
        {
            std::cout<<dpos[i][j]<<' ';
        }
        std::cout<<std::endl;
    }*/
    UpdatePos();
    return IsRestrain();
}

void bfgs::GetPos()
{
    int k=0;
    for(int i=0;i<GlobalC::ucell.ntype;i++)
    {
        for(int j=0;j<GlobalC::ucell.atoms[i].na;j++)
        {
            pos[k+j][0]=GlobalC::ucell.atoms[i].tau[j].x*ModuleBase::BOHR_TO_A*GlobalC::ucell.lat0;
            pos[k+j][1]=GlobalC::ucell.atoms[i].tau[j].y*ModuleBase::BOHR_TO_A*GlobalC::ucell.lat0;
            pos[k+j][2]=GlobalC::ucell.atoms[i].tau[j].z*ModuleBase::BOHR_TO_A*GlobalC::ucell.lat0; 
        }
        k+=GlobalC::ucell.atoms[i].na;
    }
}


void bfgs::PrepareStep()
{
    std::vector<double> changedforce = ReshapeMToV(force);
    std::vector<double> changedpos = ReshapeMToV(pos);
    Update(changedpos, changedforce);
    /*for(int i = 0; i < 3*size; i++)
    {
        for(int j = 0; j < 3*size; j++)
        {
            std::cout<<H[i][j]<<' ';
        }
        std::cout<<std::endl;
    }*/
    //call dysev
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
    std::vector<std::vector<double>> V(3*size, std::vector<double>(3*size, 0.0));
    for(int i = 0; i < 3*size; i++)
    {
        for(int j = 0; j < 3*size; j++)
        {
            V[j][i] = H_flat[3*size*i + j];
        }
    }

    /*for(int i=0;i<3*size;i++)
    {
        std::cout<<omega[i]<<' ';
    }
    std::cout<<std::endl;*/
    std::vector<double> a=DotInMAndV2(V, changedforce);
    for(int i = 0; i < a.size(); i++)
    {
        a[i] /= abs(omega[i]);
    }

    std::vector<double> tmpdpos = DotInMAndV1(V, a);


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
    /*std::cout<<"PrintPos"<<std::endl;
    for(int i=0;i<3*size;i++)
    {
        std::cout<<pos[i]<<' ';
    }
    std::cout<<std::endl;
    std::cout<<"PrintPos0"<<std::endl;
    for(int i=0;i<3*size;i++)
    {
        std::cout<<pos0[i]<<' ';
    }
    std::cout<<std::endl;*/
    std::vector<double> dpos = VSubV(pos, pos0);
    if(*max_element(dpos.begin(), dpos.end()) < 1e-7)
    {
        return;
    }
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

void bfgs::UpdatePos()
{
    double move_ion[3*size];
    ModuleBase::zeros(move_ion, size*3);

    for(int iat=0; iat<size; iat++)
    {
        //Cartesian coordinate
        //convert from Angstrom to unit of latvec (Bohr)

        //convert unit
        ModuleBase::Vector3<double> move_ion_cart;
        move_ion_cart.x = dpos[iat][0] / ModuleBase::BOHR_TO_A / GlobalC::ucell.lat0;
        move_ion_cart.y = dpos[iat][1] / ModuleBase::BOHR_TO_A / GlobalC::ucell.lat0;
        move_ion_cart.z = dpos[iat][2] / ModuleBase::BOHR_TO_A / GlobalC::ucell.lat0;

        //convert to Direct coordinate
        //note here the old GT is used

        //convert pos
        ModuleBase::Vector3<double> move_ion_dr = move_ion_cart * GlobalC::ucell.GT;


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
	GlobalC::ucell.update_pos_taud(move_ion);
    pos = MAddM(pos, dpos);
}

bool bfgs::IsRestrain()
{
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
    return a<0.00001;
}
// matrix methods

std::vector<double> bfgs::ReshapeMToV(std::vector<std::vector<double>> matrix) 
{
    int size = matrix.size();
    std::vector<double> result;
    result.reserve(3*size);
    for (const auto& row : matrix) {
        result.insert(result.end(), row.begin(), row.end());
    }
    return result;
}

std::vector<std::vector<double>> bfgs::MAddM(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b) 
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

std::vector<double> bfgs::VSubV(std::vector<double> a, std::vector<double> b) 
{
    std::vector<double> result = std::vector<double>(a.size(), 0.0);
    for(int i = 0; i < a.size(); i++)
    {
        result[i] = a[i] - b[i];
    }
    return result;
}

std::vector<std::vector<double>> bfgs::ReshapeVToM(std::vector<double> matrix) 
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

std::vector<double> bfgs::DotInMAndV1(std::vector<std::vector<double>> matrix, std::vector<double> vec) 
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
std::vector<double> bfgs::DotInMAndV2(std::vector<std::vector<double>> matrix, std::vector<double> vec) 
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

double bfgs::DotInVAndV(std::vector<double> vec1, std::vector<double> vec2) 
{
    double result = 0.0;
    for(int i = 0; i < vec1.size(); i++)
    {
        result += vec1[i] * vec2[i];
    }
    return result;
}

std::vector<std::vector<double>> bfgs::OuterVAndV(std::vector<double> a, std::vector<double> b) 
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

std::vector<std::vector<double>> bfgs::MPlus(std::vector<std::vector<double>> a, double b)
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

std::vector<std::vector<double>> bfgs::MSubM(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b)
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

std::tuple<std::vector<double>, std::vector<std::vector<double>>> bfgs::GetEigenvalueAndEigenVector(std::vector<std::vector<double>> matrix) 
{
    std::vector<double> omega;
    std::vector<std::vector<double>> V;

    return std::make_tuple(omega, V);
}
