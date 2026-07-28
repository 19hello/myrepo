#include "verlet.h"

#include "md_func.h"
#include "source_base/timer.h"

#ifdef __MPI
#include <mpi.h>
#endif

#include <map>

namespace
{
#ifdef __MPI
struct AndersonAtomKey
{
    int type;
    int type_index;
    int mbl[3];
};

struct AndersonRandomValue
{
    int type;
    int type_index;
    int collision;
    double velocity[3];
};
#endif
}

Verlet::Verlet(const Parameter& param_in, MdCell& mdcell_in) : MD_base(param_in, mdcell_in)
{
}

Verlet::~Verlet()
{
}


void Verlet::setup(ModuleESolver::ESolver* p_esolver, const std::string& global_readin_dir)
{
    ModuleBase::TITLE("Verlet", "setup");
    ModuleBase::timer::start("Verlet", "setup");

    MD_base::setup(p_esolver, global_readin_dir);

    ModuleBase::timer::end("Verlet", "setup");
}


void Verlet::first_half(std::ofstream& ofs)
{
    ModuleBase::TITLE("Verlet", "first_half");
    ModuleBase::timer::start("Verlet", "first_half");

    MD_base::update_vel();
    MD_base::update_pos();

    ModuleBase::timer::end("Verlet", "first_half");
}


void Verlet::second_half()
{
    ModuleBase::TITLE("Verlet", "second_half");
    ModuleBase::timer::start("Verlet", "second_half");

    MD_base::update_vel();
    apply_thermostat();

    ModuleBase::timer::end("Verlet", "second_half");
}


void Verlet::apply_thermostat(void)
{
    double t_target = 0.0;
    t_current = MD_func::current_temp(kinetic, mdcell, frozen_freedom_);

    if (mdp.md_type == "nve")
    {
    }
    else if (mdp.md_thermostat == "rescaling")
    {
        t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_nstep, md_tfirst, md_tlast);
        if (std::abs(t_target - t_current) * ModuleBase::Hartree_to_K > mdp.md_tolerance)
        {
            thermalize(0, t_current, t_target);
        }
    }
    else if (mdp.md_thermostat == "rescale_v")
    {
        if ((step_ + step_rst_) % mdp.md_nraise == 0)
        {
            t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_nstep, md_tfirst, md_tlast);
            thermalize(0, t_current, t_target);
        }
    }
    else if (mdp.md_thermostat == "anderson")
    {
#ifdef __MPI
        if (mdcell.mpi_size() > 1)
        {
            const MPI_Comm comm = mdcell.communicator();
            int rank = 0;
            int size = 1;
            MPI_Comm_rank(comm, &rank);
            MPI_Comm_size(comm, &size);

            std::vector<AndersonAtomKey> local_keys(static_cast<std::size_t>(mdcell.nlocal()));
            std::map<std::pair<int, int>, int> local_indices;
            for (int i = 0; i < mdcell.nlocal(); ++i)
            {
                const LocalAtom& atom = mdcell.owned_atoms()[static_cast<std::size_t>(i)];
                AndersonAtomKey& key = local_keys[static_cast<std::size_t>(i)];
                key.type = atom.type;
                key.type_index = atom.type_index;
                for (int k = 0; k < 3; ++k) key.mbl[k] = atom.mbl[k];
                local_indices[std::make_pair(atom.type, atom.type_index)] = i;
            }

            const MdStruMetadata& metadata = mdcell.stru_metadata();
            std::vector<int> type_offsets(metadata.species.size() + 1, 0);
            for (std::size_t it = 0; it < metadata.species.size(); ++it)
            {
                type_offsets[it + 1] = type_offsets[it] + metadata.species[it].atom_count;
            }
            if (type_offsets.back() != mdcell.nat())
            {
                ModuleBase::WARNING_QUIT("Verlet::apply_thermostat", "MdCell STRU metadata does not match the global atom count.");
            }

            if (rank == 0)
            {
                std::vector<int> owners(static_cast<std::size_t>(mdcell.nat()), -1);
                std::vector<ModuleBase::Vector3<int> > move_flags(static_cast<std::size_t>(mdcell.nat()));
                for (const AndersonAtomKey& key : local_keys)
                {
                    const int index = type_offsets[static_cast<std::size_t>(key.type)] + key.type_index;
                    owners[static_cast<std::size_t>(index)] = 0;
                    move_flags[static_cast<std::size_t>(index)].set(key.mbl[0], key.mbl[1], key.mbl[2]);
                }
                for (int source = 1; source < size; ++source)
                {
                    int count = 0;
                    MPI_Recv(&count, 1, MPI_INT, source, 9600, comm, MPI_STATUS_IGNORE);
                    std::vector<AndersonAtomKey> keys(static_cast<std::size_t>(count));
                    if (count > 0)
                    {
                        MPI_Recv(keys.data(), count * static_cast<int>(sizeof(AndersonAtomKey)), MPI_BYTE, source, 9601, comm, MPI_STATUS_IGNORE);
                    }
                    for (const AndersonAtomKey& key : keys)
                    {
                        const int index = type_offsets[static_cast<std::size_t>(key.type)] + key.type_index;
                        owners[static_cast<std::size_t>(index)] = source;
                        move_flags[static_cast<std::size_t>(index)].set(key.mbl[0], key.mbl[1], key.mbl[2]);
                    }
                }

                for (std::size_t it = 0; it < metadata.species.size(); ++it)
                {
                    const double deviation = sqrt(md_tlast / (metadata.species[it].mass / ModuleBase::AU_to_MASS));
                    for (int ia = 0; ia < metadata.species[it].atom_count; ++ia)
                    {
                        AndersonRandomValue random_value;
                        random_value.type = static_cast<int>(it);
                        random_value.type_index = ia;
                        random_value.collision = static_cast<double>(std::rand()) / RAND_MAX <= 1.0 / mdp.md_nraise;
                        const int index = type_offsets[it] + ia;
                        for (int k = 0; k < 3; ++k)
                        {
                            random_value.velocity[k] = 0.0;
                            if (random_value.collision && move_flags[static_cast<std::size_t>(index)][k])
                            {
                                random_value.velocity[k] = deviation * MD_func::gaussrand();
                            }
                        }
                        const int owner = owners[static_cast<std::size_t>(index)];
                        if (owner == 0)
                        {
                            const int local_index = local_indices[std::make_pair(random_value.type, random_value.type_index)];
                            if (random_value.collision)
                            {
                                LocalAtom& atom = mdcell.mutable_owned_atoms()[static_cast<std::size_t>(local_index)];
                                for (int k = 0; k < 3; ++k) if (atom.mbl[k]) atom.vel[k] = random_value.velocity[k];
                            }
                        }
                        else
                        {
                            MPI_Send(&random_value, sizeof(AndersonRandomValue), MPI_BYTE, owner, 9602, comm);
                        }
                    }
                }
            }
            else
            {
                const int count = static_cast<int>(local_keys.size());
                MPI_Send(&count, 1, MPI_INT, 0, 9600, comm);
                if (count > 0)
                {
                    MPI_Send(local_keys.data(), count * static_cast<int>(sizeof(AndersonAtomKey)), MPI_BYTE, 0, 9601, comm);
                }
                for (int i = 0; i < count; ++i)
                {
                    AndersonRandomValue random_value;
                    MPI_Recv(&random_value, sizeof(AndersonRandomValue), MPI_BYTE, 0, 9602, comm, MPI_STATUS_IGNORE);
                    if (random_value.collision)
                    {
                        const int local_index = local_indices[std::make_pair(random_value.type, random_value.type_index)];
                        LocalAtom& atom = mdcell.mutable_owned_atoms()[static_cast<std::size_t>(local_index)];
                        for (int k = 0; k < 3; ++k) if (atom.mbl[k]) atom.vel[k] = random_value.velocity[k];
                    }
                }
            }
            return;
        }
#endif
        if (my_rank == 0)
        {
            for (LocalAtom& atom : mdcell.mutable_owned_atoms())
            {
                if (static_cast<double>(std::rand()) / RAND_MAX <= 1.0 / mdp.md_nraise)
                {
                    const double deviation = sqrt(md_tlast / atom.mass);
                    for (int k = 0; k < 3; ++k)
                    {
                        if (atom.mbl[k])
                        {
                            atom.vel[k] = deviation * MD_func::gaussrand();
                        }
                    }
                }
            }
        }
    }
    else if (mdp.md_thermostat == "berendsen")
    {
        t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_nstep, md_tfirst, md_tlast);
        thermalize(mdp.md_nraise, t_current, t_target);
    }
    else if (mdp.md_thermostat == "csvr")
    {
        t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_nstep, md_tfirst, md_tlast);
        apply_csvr(t_current, t_target);
    }
    else
    {
        ModuleBase::WARNING_QUIT("Verlet", "No such thermostat!");
    }
}


void Verlet::thermalize(const int& nraise, const double& current_temp, const double& target_temp)
{
    double fac = 0.0;
    if (nraise > 0 && current_temp > 0 && target_temp > 0)
    {
        fac = sqrt(1 + (target_temp / current_temp - 1) / nraise);
    }
    else if (nraise == 0 && current_temp > 0 && target_temp > 0)
    {
        fac = sqrt(target_temp / current_temp);
    }

    for (LocalAtom& atom : mdcell.mutable_owned_atoms()) atom.vel *= fac;
}


void Verlet::apply_csvr(const double& current_temp, const double& target_temp)
{
    // CSVR thermostat: Canonical Sampling through Velocity Rescaling
    // Reference: G. Bussi, D. Donadio, M. Parrinello, J. Chem. Phys. 126, 014101 (2007)

    if (current_temp <= 0.0 || target_temp <= 0.0)
    {
        return;
    }

    // Get degrees of freedom (3N - frozen)
    int ndeg = MD_func::global_dof(mdcell, frozen_freedom_);

    // Calculate kinetic energies
    double kin_energy = current_temp * ndeg * 0.5;  // in Hartree
    double kin_target = target_temp * ndeg * 0.5;   // in Hartree

    // Calculate tau parameter (characteristic time scale / dt)
    double taut = mdp.md_csvr_tau / mdp.md_dt;

    // Calculate decay factor
    double factor = 0.0;
    if (taut > 0.1)
    {
        factor = exp(-1.0 / taut);
    }

    // Generate Gaussian random numbers using MD_func
    double rr = MD_func::gaussrand();

    // Calculate sum of squared Gaussian random numbers (ndeg - 1)
    double sumnoises = 0.0;
    for (int i = 0; i < ndeg - 1; ++i)
    {
        double r = MD_func::gaussrand();
        sumnoises += r * r;
    }

    // CSVR core formula (simplified)
    double factor2 = (1.0 - factor) * kin_target / kin_energy / ndeg;
    double resample = factor + factor2 * (rr * rr + sumnoises) + 2.0 * rr * sqrt(factor * factor2);

    // Ensure non-negative
    resample = std::max(0.0, resample);

    // Calculate scaling factor
    double scale = sqrt(resample);

    // Apply velocity scaling
    for (LocalAtom& atom : mdcell.mutable_owned_atoms()) atom.vel *= scale;
}


void Verlet::print_md(std::ofstream& ofs, const bool& cal_stress)
{
    MD_base::print_md(ofs, cal_stress);
    return;
}


void Verlet::write_restart(const std::string& global_out_dir)
{
    MD_base::write_restart(global_out_dir);
    return;
}


void Verlet::restart(const std::string& global_readin_dir)
{
    MD_base::restart(global_readin_dir);
    return;
}
