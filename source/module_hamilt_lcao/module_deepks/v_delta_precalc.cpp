/// cal_v_delta_precalc : v_delta_precalc is used for training with v_delta label,
//                         which equals gvdm * v_delta_pdm_shell,
//                         v_delta_pdm_shell = overlap * overlap
/// check_v_delta_precalc : check v_delta_precalc

#ifdef __DEEPKS

#include "LCAO_deepks.h"
#include "LCAO_deepks_io.h" // mohan add 2024-07-22
#include "module_base/blas_connector.h"
#include "module_base/constants.h"
#include "module_base/libm/libm.h"
#include "module_base/parallel_reduce.h"
#include "module_hamilt_lcao/module_hcontainer/atom_pair.h"
#include "module_parameter/parameter.h"


// calculates v_delta_precalc[nks,nlocal,nlocal,NAt,NDscrpt] = gvdm * v_delta_pdm_shell;
// v_delta_pdm_shell[nks,nlocal,nlocal,Inl,nm*nm] = overlap * overlap;
// for deepks_v_delta = 1
void LCAO_Deepks::cal_v_delta_precalc(const int nlocal,
    const int nat,
    const UnitCell &ucell,
    const LCAO_Orbitals &orb,
    Grid_Driver &GridD)
{
    ModuleBase::TITLE("LCAO_Deepks", "calc_v_delta_precalc");
    // timeval t_start;
    // gettimeofday(&t_start,NULL);

    this->cal_gvdm(nat);
    const double Rcut_Alpha = orb.Alpha[0].getRcut();
    this->init_v_delta_pdm_shell(1,nlocal);
   
    for (int T0 = 0; T0 < ucell.ntype; T0++)
    {
		Atom* atom0 = &ucell.atoms[T0]; 
        
        for (int I0 =0; I0< atom0->na; I0++)
        {
            const int iat = ucell.itia2iat(T0,I0);
            const ModuleBase::Vector3<double> tau0 = atom0->tau[I0];
            GridD.Find_atom(ucell, atom0->tau[I0] ,T0, I0);

            for (int ad1=0; ad1<GridD.getAdjacentNum()+1 ; ++ad1)
            {
                const int T1 = GridD.getType(ad1);
                const int I1 = GridD.getNatom(ad1);
                const int start1 = ucell.itiaiw2iwt(T1, I1, 0);
                const ModuleBase::Vector3<double> tau1 = GridD.getAdjacentTau(ad1);
				const Atom* atom1 = &ucell.atoms[T1];
				const int nw1_tot = atom1->nw*PARAM.globalv.npol;
				const double Rcut_AO1 = orb.Phi[T1].getRcut(); 

                const double dist1 = (tau1-tau0).norm() * ucell.lat0;
                if (dist1 >= Rcut_Alpha + Rcut_AO1)
                {
                    continue;
                }

				for (int ad2=0; ad2 < GridD.getAdjacentNum()+1 ; ad2++)
				{
					const int T2 = GridD.getType(ad2);
					const int I2 = GridD.getNatom(ad2);
					const int start2 = ucell.itiaiw2iwt(T2, I2, 0);
					const ModuleBase::Vector3<double> tau2 = GridD.getAdjacentTau(ad2);
					const Atom* atom2 = &ucell.atoms[T2];
					const int nw2_tot = atom2->nw*PARAM.globalv.npol;
					
					const double Rcut_AO2 = orb.Phi[T2].getRcut();
                	const double dist2 = (tau2-tau0).norm() * ucell.lat0;

					if (dist2 >= Rcut_Alpha + Rcut_AO2)
					{
						continue;
					}

					for (int iw1=0; iw1<nw1_tot; ++iw1)
					{
						const int iw1_all = start1 + iw1; // this is \mu
						const int iw1_local = pv->global2local_row(iw1_all);
						if(iw1_local < 0) {continue;
}
						const int iw1_0 = iw1/PARAM.globalv.npol;
						for (int iw2=0; iw2<nw2_tot; ++iw2)
						{
							const int iw2_all = start2 + iw2; // this is \nu
							const int iw2_local = pv->global2local_col(iw2_all);
							if(iw2_local < 0) {continue;
}
							const int iw2_0 = iw2/PARAM.globalv.npol;

                            std::vector<double> nlm1 = this->nlm_save[iat][ad1][iw1][0];
                            std::vector<double> nlm2 = this->nlm_save[iat][ad2][iw2][0];

                            assert(nlm1.size()==nlm2.size());
                            int ib=0;
                            for (int L0 = 0; L0 <= orb.Alpha[0].getLmax();++L0)
                            {
                                for (int N0 = 0;N0 < orb.Alpha[0].getNchi(L0);++N0)
                                {
                                    const int inl = this->inl_index[T0](I0, L0, N0);
                                    const int nm = 2*L0+1;
                                    
                                    for (int m1=0; m1<nm; ++m1) // nm = 1 for s, 3 for p, 5 for d
                                    {
                                        for (int m2=0; m2<nm; ++m2) // nm = 1 for s, 3 for p, 5 for d
                                        {
                                            v_delta_pdm_shell[0][iw1_all][iw2_all][inl][m1*nm+m2] += nlm1[ib+m1]*nlm2[ib+m2];
                                        }
                                    }
                                    ib+=nm;
                                }
                            }                            
						}//iw2
					}//iw1
				}//ad2
			}//ad1   
            
        }
    }
#ifdef __MPI
    const int mn_size=(2 * this->lmaxd + 1) * (2 * this->lmaxd + 1);
    for(int inl = 0; inl < this->inlmax; inl++)
    {
        for(int mu = 0; mu < nlocal ; mu++)
        {
            for(int nu=0; nu< nlocal ; nu++)
            {
                Parallel_Reduce::reduce_all(this->v_delta_pdm_shell[0][mu][nu][inl],mn_size);
            }
        }
    }
#endif    
    
    // transfer v_delta_pdm_shell to v_delta_pdm_shell_vector
    
    int nlmax = this->inlmax/nat;
   
    std::vector<torch::Tensor> v_delta_pdm_shell_vector;
    for(int nl = 0; nl < nlmax; ++nl)
    {
        std::vector<torch::Tensor> kuuammv;
        for(int iks = 0; iks < 1; ++iks)
        {
            std::vector<torch::Tensor> uuammv;
            for(int mu = 0; mu < nlocal; ++mu)
            {
                std::vector<torch::Tensor> uammv;
                for(int nu =0 ; nu < nlocal; ++nu)
                {
                    std::vector<torch::Tensor> ammv;
                    for (int iat=0; iat<nat; ++iat)
                    {
                        int inl = iat*nlmax+nl;
                        int nm = 2*this->inl_l[inl]+1;
                        std::vector<double> mmv;
                    
                        for (int m1=0; m1<nm; ++m1) // m1 = 1 for s, 3 for p, 5 for d
                        {
                            for (int m2=0; m2<nm; ++m2) // m1 = 1 for s, 3 for p, 5 for d
                            {
                                mmv.push_back(this->v_delta_pdm_shell[iks][mu][nu][inl][m1*nm+m2]);
                            }
                        }
                        torch::Tensor mm = torch::tensor(mmv, torch::TensorOptions().dtype(torch::kFloat64) ).reshape({nm, nm});    //nm*nm
                        ammv.push_back(mm);
                    }
                    torch::Tensor amm = torch::stack(ammv, 0); 
                    uammv.push_back(amm);                    
                }
                torch::Tensor uamm = torch::stack(uammv, 0); 
                uuammv.push_back(uamm);
            }
            torch::Tensor uuamm = torch::stack(uuammv, 0); 
            kuuammv.push_back(uuamm);
        }
        torch::Tensor kuuamm = torch::stack(kuuammv, 0);  
        v_delta_pdm_shell_vector.push_back(kuuamm);
    }
       
    assert(v_delta_pdm_shell_vector.size() == nlmax);
    
    //einsum for each nl: 
    std::vector<torch::Tensor> v_delta_precalc_vector;
    for (int nl = 0; nl<nlmax; ++nl)
    {
        v_delta_precalc_vector.push_back(at::einsum("kxyamn, avmn->kxyav", {v_delta_pdm_shell_vector[nl], this->gevdm_vector[nl]}));
    }

    this->v_delta_precalc_tensor = torch::cat(v_delta_precalc_vector, -1);
    this->del_v_delta_pdm_shell(1,nlocal);

    //check_v_delta_precalc(nlocal,nat);
    // timeval t_end;
    // gettimeofday(&t_end,NULL);
    // std::cout<<"calculate v_delta_precalc time:\t"<<(double)(t_end.tv_sec-t_start.tv_sec) + (double)(t_end.tv_usec-t_start.tv_usec)/1000000.0<<std::endl;
    return;
}

void LCAO_Deepks::check_v_delta_precalc(const int nat, const int nks,const int nlocal)
{
    std::ofstream ofs("v_delta_precalc.dat");
    ofs << std::setprecision(10);
    for (int iks = 0; iks < nks; ++iks)
    {
        for (int mu = 0; mu < nlocal; ++mu)
        {
            for (int nu = 0; nu < nlocal; ++nu)
            {
                for (int iat = 0;iat < nat;++iat)
                {
                    for(int p=0; p<this->des_per_atom; ++p)
                    {
                        ofs<<this->v_delta_precalc_tensor.index({iks, mu, nu, iat, p }).item().toDouble()<<" ";
                    }
                }
                ofs << std::endl;                
            }
        }
    }
    ofs.close();
}

#endif
