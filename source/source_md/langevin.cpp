#include "langevin.h"

#include "md_func.h"
#include "source_base/parallel_common.h"
#include "source_base/timer.h"

#include <map>

namespace
{
#ifdef __MPI
struct LangevinAtomKey
{
    int type;
    int type_index;
};

struct LangevinRandomValue
{
    int type;
    int type_index;
    double value[3];
};
#endif
}

Langevin::Langevin(const Parameter& param_in, MdCell& mdcell_in) : MD_base(param_in, mdcell_in)
{
    /// convert to a.u. unit
    assert(ModuleBase::AU_to_FS!=0.0);

    md_damp = mdp.md_damp / ModuleBase::AU_to_FS;

    assert(mdcell.nlocal() > 0);

    total_force.resize(static_cast<std::size_t>(mdcell.nlocal()));
}


void Langevin::setup(ModuleESolver::ESolver* p_esolver, const std::string& global_readin_dir)
{
    ModuleBase::TITLE("Langevin", "setup");
    ModuleBase::timer::start("Langevin", "setup");

    MD_base::setup(p_esolver, global_readin_dir);

    post_force();

    ModuleBase::timer::end("Langevin", "setup");
    return;
}


void Langevin::first_half(std::ofstream& ofs)
{
    ModuleBase::TITLE("Langevin", "first_half");
    ModuleBase::timer::start("Langevin", "first_half");

    for (int i = 0; i < mdcell.nlocal(); ++i)
    {
        LocalAtom& atom = mdcell.mutable_owned_atoms()[static_cast<std::size_t>(i)];
        for (int k = 0; k < 3; ++k)
        {
            if (atom.mbl[k]) atom.vel[k] += 0.5 * total_force[i][k] * md_dt / atom.mass;
        }
    }
    MD_base::update_pos();

    ModuleBase::timer::end("Langevin", "first_half");
    return;
}


void Langevin::second_half()
{
    ModuleBase::TITLE("Langevin", "second_half");
    ModuleBase::timer::start("Langevin", "second_half");

    post_force();
    for (int i = 0; i < mdcell.nlocal(); ++i)
    {
        LocalAtom& atom = mdcell.mutable_owned_atoms()[static_cast<std::size_t>(i)];
        for (int k = 0; k < 3; ++k)
        {
            if (atom.mbl[k]) atom.vel[k] += 0.5 * total_force[i][k] * md_dt / atom.mass;
        }
    }

    ModuleBase::timer::end("Langevin", "second_half");
    return;
}


void Langevin::print_md(std::ofstream& ofs, const bool& cal_stress)
{
    MD_base::print_md(ofs, cal_stress);
    return;
}


void Langevin::write_restart(const std::string& global_out_dir)
{
    MD_base::write_restart(global_out_dir);
    return;
}


void Langevin::restart(const std::string& global_readin_dir)
{
    MD_base::restart(global_readin_dir);
    return;
}


void Langevin::post_force()
{
    double t_target = MD_func::target_temp(step_ + step_rst_, mdp.md_nstep, md_tfirst, md_tlast);
    total_force.resize(static_cast<std::size_t>(mdcell.nlocal()));

    std::vector<ModuleBase::Vector3<double> > random_values(static_cast<std::size_t>(mdcell.nlocal()));
#ifdef __MPI
    if (mdcell.mpi_size() > 1)
    {
        const MPI_Comm comm = mdcell.communicator();
        int rank = 0;
        int size = 1;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &size);

        std::vector<LangevinAtomKey> local_keys(static_cast<std::size_t>(mdcell.nlocal()));
        std::map<std::pair<int, int>, int> local_indices;
        for (int i = 0; i < mdcell.nlocal(); ++i)
        {
            const LocalAtom& atom = mdcell.owned_atoms()[static_cast<std::size_t>(i)];
            local_keys[static_cast<std::size_t>(i)].type = atom.type;
            local_keys[static_cast<std::size_t>(i)].type_index = atom.type_index;
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
            ModuleBase::WARNING_QUIT("Langevin::post_force", "MdCell STRU metadata does not match the global atom count.");
        }

        if (rank == 0)
        {
            std::vector<int> owners(static_cast<std::size_t>(mdcell.nat()), -1);
            for (const LangevinAtomKey& key : local_keys)
            {
                owners[static_cast<std::size_t>(type_offsets[static_cast<std::size_t>(key.type)] + key.type_index)] = 0;
            }
            for (int source = 1; source < size; ++source)
            {
                int count = 0;
                MPI_Recv(&count, 1, MPI_INT, source, 9500, comm, MPI_STATUS_IGNORE);
                std::vector<LangevinAtomKey> keys(static_cast<std::size_t>(count));
                if (count > 0)
                {
                    MPI_Recv(keys.data(), count * static_cast<int>(sizeof(LangevinAtomKey)), MPI_BYTE, source, 9501, comm, MPI_STATUS_IGNORE);
                }
                for (const LangevinAtomKey& key : keys)
                {
                    owners[static_cast<std::size_t>(type_offsets[static_cast<std::size_t>(key.type)] + key.type_index)] = source;
                }
            }

            for (std::size_t it = 0; it < metadata.species.size(); ++it)
            {
                for (int ia = 0; ia < metadata.species[it].atom_count; ++ia)
                {
                    LangevinRandomValue random_value;
                    random_value.type = static_cast<int>(it);
                    random_value.type_index = ia;
                    for (int k = 0; k < 3; ++k)
                    {
                        random_value.value[k] = static_cast<double>(std::rand()) / RAND_MAX - 0.5;
                    }
                    const int owner = owners[static_cast<std::size_t>(type_offsets[it] + ia)];
                    if (owner == 0)
                    {
                        const int local_index = local_indices[std::make_pair(random_value.type, random_value.type_index)];
                        for (int k = 0; k < 3; ++k) random_values[static_cast<std::size_t>(local_index)][k] = random_value.value[k];
                    }
                    else
                    {
                        MPI_Send(&random_value, sizeof(LangevinRandomValue), MPI_BYTE, owner, 9502, comm);
                    }
                }
            }
        }
        else
        {
            const int count = static_cast<int>(local_keys.size());
            MPI_Send(&count, 1, MPI_INT, 0, 9500, comm);
            if (count > 0)
            {
                MPI_Send(local_keys.data(), count * static_cast<int>(sizeof(LangevinAtomKey)), MPI_BYTE, 0, 9501, comm);
            }
            for (int i = 0; i < count; ++i)
            {
                LangevinRandomValue random_value;
                MPI_Recv(&random_value, sizeof(LangevinRandomValue), MPI_BYTE, 0, 9502, comm, MPI_STATUS_IGNORE);
                const int local_index = local_indices[std::make_pair(random_value.type, random_value.type_index)];
                for (int k = 0; k < 3; ++k) random_values[static_cast<std::size_t>(local_index)][k] = random_value.value[k];
            }
        }
    }
    else
#endif
    {
        for (int i = 0; i < mdcell.nlocal(); ++i)
        {
            for (int k = 0; k < 3; ++k)
            {
                random_values[static_cast<std::size_t>(i)][k] = static_cast<double>(std::rand()) / RAND_MAX - 0.5;
            }
        }
    }

    for (int i = 0; i < mdcell.nlocal(); ++i)
    {
        const LocalAtom& atom = mdcell.owned_atoms()[static_cast<std::size_t>(i)];
        ModuleBase::Vector3<double> fictitious_force = -atom.mass * atom.vel / md_damp;
        for (int j = 0; j < 3; ++j)
        {
            fictitious_force[j] += sqrt(24.0 * t_target * atom.mass / md_damp / md_dt)
                                   * random_values[static_cast<std::size_t>(i)][j];
        }
        total_force[i] = atom.force + fictitious_force;
    }
}
