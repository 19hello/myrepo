#include <regex>
#include <cassert>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <sstream>
#include <stdexcept>
#include <unistd.h>

#include "print_cell.h"
#include "source_cell/md_cell.h"
#include "source_base/formatter.h"
#include "source_base/tool_title.h"
#include "source_base/global_variable.h"

#ifdef __MPI
#include <mpi.h>
#endif

namespace unitcell
{
    void print_tau(Atom* atoms,
                   const std::string& Coordinate,
                   const int ntype,
                   const double lat0,
                   std::ofstream &ofs)
    {
        ModuleBase::TITLE("UnitCell", "print_tau");
        // assert (direct || Coordinate == "Cartesian" || Coordinate == "Cartesian_angstrom"); // this line causes abort in unittest ReadAtomPositionsCACXY.
        // previously there are two if-statements, the first is `if(Coordinate == "Direct")` and the second is `if(Coordinate == "Cartesian" || Coordiante == "Cartesian_angstrom")`
        // however the Coordinate can also be value among Cartesian_angstrom_center_xy, Cartesian_angstrom_center_xz, Cartesian_angstrom_center_yz and Cartesian_angstrom_center_xyz

        // if Coordinate has value one of them, this print_tau will not print anything.
        std::regex pattern("Direct|Cartesian(_angstrom)?(_center_(xy|xz|yz|xyz))?");
        assert(std::regex_search(Coordinate, pattern));
        bool direct = (Coordinate == "Direct");

        //----------------------
        // print atom positions
        //----------------------
        std::string table;
        table += direct? " DIRECT COORDINATES\n": FmtCore::format(" CARTESIAN COORDINATES ( UNIT = %15.8f Bohr )\n", lat0);
        table += FmtCore::format("%5s%19s%19s%19s%8s\n", "atom", "x", "y", "z", "mag");
        for(int it = 0; it < ntype; it++)
        {
            for (int ia = 0; ia < atoms[it].na; ia++)
            {
                const double& x = direct? atoms[it].taud[ia].x: atoms[it].tau[ia].x;
                const double& y = direct? atoms[it].taud[ia].y: atoms[it].tau[ia].y;
                const double& z = direct? atoms[it].taud[ia].z: atoms[it].tau[ia].z;
                table += FmtCore::format("%5s%19.12f%19.12f%19.12f%8.4f\n", 
                                        atoms[it].label, 
                                        x, 
                                        y, 
                                        z, 
                                        atoms[it].mag[ia]); 
            }
        }
        table += "\n";
        ofs << table; 


        // print velocities
        ofs << " ATOMIC VELOCITIES" << std::endl;
        ofs << std::setprecision(12);
        ofs << std::setw(5) << "atom" 
            << std::setw(19) << "vx" 
            << std::setw(19) << "vy" 
            << std::setw(19) << "vz"
            << std::endl;
 
        for(int it = 0; it < ntype; it++)
        {
            for (int ia = 0; ia < atoms[it].na; ia++)
            {
                ofs << std::setw(5) << atoms[it].label;
                ofs << " " << std::setw(18) << atoms[it].vel[ia].x;
                ofs << " " << std::setw(18) << atoms[it].vel[ia].y;
                ofs << " " << std::setw(18) << atoms[it].vel[ia].z;
                ofs << std::endl;
            }
        }
        ofs << std::endl;
        ofs << std::setprecision(6); // return to 6, as original


        return;
    }

    void print_stru_file(const UnitCell& ucell,
                         const Atom*     atoms,
                         const ModuleBase::Matrix3& latvec,
                         const std::string& fn, 
                         const int& nspin,
                         const bool& direct,
                         const bool& vel,
                         const bool& magmom,
                         const bool& orb,
                         const bool& dpks_desc,
                         const int& iproc)
    {
        ModuleBase::TITLE("UnitCell","print_stru_file");
        if (iproc != 0) 
        {
            return; // old: if(GlobalV::MY_RANK != 0) return;
        }
        // ATOMIC_SPECIES
        std::string str = "ATOMIC_SPECIES\n";
        for(int it=0; it<ucell.ntype; it++)
        { 
            str += FmtCore::format("%s %8.4f %s %s\n", 
                                    ucell.atom_label[it], 
                                    ucell.atom_mass[it], 
                                    ucell.pseudo_fn[it], 
                                    ucell.pseudo_type[it]); 
        }
        // NUMERICAL_ORBITAL
        if(orb)
        {
            str += "\nNUMERICAL_ORBITAL\n";
            for(int it = 0; it < ucell.ntype; it++) 
            { 
                str += ucell.orbital_fn[it] + "\n"; 
            }
        }
        // NUMERICAL_DESCRIPTOR
        if(dpks_desc) 
        { 
            str += "\nNUMERICAL_DESCRIPTOR\n" + ucell.descriptor_file + "\n"; 
        }
        // LATTICE_CONSTANT
        str += "\nLATTICE_CONSTANT\n" + FmtCore::format("%-.10f\n", ucell.lat0);
        // LATTICE_VECTORS
        str += "\nLATTICE_VECTORS\n";
        str += FmtCore::format("%20.10f%20.10f%20.10f\n", latvec.e11, latvec.e12, latvec.e13);
        str += FmtCore::format("%20.10f%20.10f%20.10f\n", latvec.e21, latvec.e22, latvec.e23);
        str += FmtCore::format("%20.10f%20.10f%20.10f\n", latvec.e31, latvec.e32, latvec.e33);
        // ATOMIC_POSITIONS
        str += "\nATOMIC_POSITIONS\n";
        const std::string scale = direct? "Direct": "Cartesian";
        int nat_ = 0; // counter iat, for printing out Mulliken magmom who is indexed by iat
        str += scale + "\n";
        for(int it = 0; it < ucell.ntype; it++)
        {
            str += "\n" + ucell.atoms[it].label + " #label\n";
            str += FmtCore::format("%-8.4f #magnetism\n", ucell.magnet.start_mag[it]);
            str += FmtCore::format("%d #number of atoms\n", atoms[it].na);
            for(int ia = 0; ia < atoms[it].na; ia++)
            {
                // output position
                const double& x = direct? atoms[it].taud[ia].x: atoms[it].tau[ia].x;
                const double& y = direct? atoms[it].taud[ia].y: atoms[it].tau[ia].y;
                const double& z = direct? atoms[it].taud[ia].z: atoms[it].tau[ia].z;
                str += FmtCore::format("%20.10f%20.10f%20.10f", x, y, z);
                str += FmtCore::format(" m%2d%2d%2d", atoms[it].mbl[ia].x, atoms[it].mbl[ia].y, atoms[it].mbl[ia].z);
                if (vel) // output velocity
                {
                    str += FmtCore::format(" v%20.10f%20.10f%20.10f", atoms[it].vel[ia].x, atoms[it].vel[ia].y, atoms[it].vel[ia].z);
                }
                if (nspin == 2 && magmom) // output magnetic information
                {
                    str += FmtCore::format(" mag%8.4f", ucell.atom_mulliken[nat_][1]);
                }
                else if (nspin == 4 && magmom) // output magnetic information
                {
                    str += FmtCore::format(" mag%8.4f%8.4f%8.4f", 
                                            ucell.atom_mulliken[nat_][1], 
                                            ucell.atom_mulliken[nat_][2], 
                                            ucell.atom_mulliken[nat_][3]);
                }
                str += "\n";
                nat_++;
            }
        }
        std::ofstream ofs(fn.c_str());
        ofs << str;
        ofs.close();
        return;
    }
}

namespace
{
std::string mdcell_stru_header(const MdCell& cell)
{
    const MdStruMetadata& metadata = cell.stru_metadata();
    std::ostringstream output;
    output << std::fixed << std::setprecision(10);
    output << "ATOMIC_SPECIES\n";
    for (std::size_t it = 0; it < metadata.species.size(); ++it)
    {
        const MdStruSpecies& species = metadata.species[it];
        output << species.label << " " << std::setprecision(4) << species.mass << std::setprecision(10);
        if (!species.pseudo_file.empty()) output << " " << species.pseudo_file;
        if (!species.pseudo_type.empty()) output << " " << species.pseudo_type;
        output << "\n";
    }
    bool has_orbitals = false;
    for (std::size_t it = 0; it < metadata.species.size(); ++it)
        has_orbitals = has_orbitals || !metadata.species[it].orbital_file.empty();
    if (has_orbitals)
    {
        output << "\nNUMERICAL_ORBITAL\n";
        for (std::size_t it = 0; it < metadata.species.size(); ++it)
            output << metadata.species[it].orbital_file << "\n";
    }
    if (!metadata.descriptor_file.empty()) output << "\nNUMERICAL_DESCRIPTOR\n" << metadata.descriptor_file << "\n";
    output << "\nLATTICE_CONSTANT\n" << cell.lat0() << "\n\nLATTICE_VECTORS\n";
    const ModuleBase::Matrix3& lattice = cell.latvec();
    output << lattice.e11 << " " << lattice.e12 << " " << lattice.e13 << "\n";
    output << lattice.e21 << " " << lattice.e22 << " " << lattice.e23 << "\n";
    output << lattice.e31 << " " << lattice.e32 << " " << lattice.e33 << "\n";
    output << "\nATOMIC_POSITIONS\nCartesian\n";
    return output.str();
}

std::string mdcell_type_header(const MdStruSpecies& species)
{
    std::ostringstream output;
    output << "\n" << species.label << " #label\n";
    output << std::fixed << std::setprecision(4) << species.start_mag << " #magnetism\n";
    output << species.atom_count << " #number of atoms\n";
    return output.str();
}

std::string mdcell_atom_line(const LocalAtom& atom)
{
    std::ostringstream output;
    output << std::fixed << std::setprecision(10)
           << atom.cart.x << " " << atom.cart.y << " " << atom.cart.z
           << " m " << atom.mbl.x << " " << atom.mbl.y << " " << atom.mbl.z
           << " v " << atom.vel.x << " " << atom.vel.y << " " << atom.vel.z << "\n";
    return output.str();
}

#ifdef __MPI
bool write_at(const int file, const std::string& data, MPI_Offset offset)
{
    std::size_t written = 0;
    while (written < data.size())
    {
        const ssize_t count = pwrite(file,
                                     data.data() + written,
                                     data.size() - written,
                                     static_cast<off_t>(offset + written));
        if (count <= 0)
        {
            return false;
        }
        written += static_cast<std::size_t>(count);
    }
    return true;
}
#endif
}

namespace mdcell
{
void print_stru_file(const MdCell& cell, const std::string& fn)
{
    const MdStruMetadata& metadata = cell.stru_metadata();
    if (metadata.species.empty()) throw std::runtime_error("MdCell STRU metadata is unavailable.");
    std::vector<std::string> local_atoms(metadata.species.size());
    for (std::size_t iat = 0; iat < cell.owned_atoms().size(); ++iat)
    {
        const LocalAtom& atom = cell.owned_atoms()[iat];
        local_atoms[static_cast<std::size_t>(atom.type)] += mdcell_atom_line(atom);
    }
    const std::string header = mdcell_stru_header(cell);
#ifdef __MPI
    int rank = 0;
    const MPI_Comm comm = cell.communicator();
    MPI_Comm_rank(comm, &rank);
    int header_ok = 1;
    if (rank == 0)
    {
        const int header_file = open(fn.c_str(), O_CREAT | O_TRUNC | O_WRONLY, 0666);
        if (header_file < 0)
        {
            header_ok = 0;
        }
        else
        {
            const bool header_written = write_at(header_file, header, 0);
            const bool header_closed = close(header_file) == 0;
            header_ok = header_written && header_closed;
        }
    }
    MPI_Bcast(&header_ok, 1, MPI_INT, 0, comm);
    if (header_ok == 0)
    {
        throw std::runtime_error("Unable to create MdCell restart STRU file: " + fn + ": " + std::strerror(errno));
    }
    MPI_Barrier(comm);

    const int file = open(fn.c_str(), O_WRONLY);
    int file_ok = file >= 0 ? 1 : 0;
    int all_files_ok = 0;
    MPI_Allreduce(&file_ok, &all_files_ok, 1, MPI_INT, MPI_MIN, comm);
    if (all_files_ok == 0)
    {
        if (file >= 0) close(file);
        throw std::runtime_error("Unable to open MdCell restart STRU file: " + fn + ": " + std::strerror(errno));
    }

    MPI_Offset offset = static_cast<MPI_Offset>(header.size());
    for (std::size_t it = 0; it < metadata.species.size(); ++it)
    {
        const std::string type_header = mdcell_type_header(metadata.species[it]);
        int type_header_ok = 1;
        if (rank == 0) type_header_ok = write_at(file, type_header, offset) ? 1 : 0;
        MPI_Bcast(&type_header_ok, 1, MPI_INT, 0, comm);
        if (type_header_ok == 0)
        {
            close(file);
            throw std::runtime_error("Unable to write MdCell restart STRU type header: " + fn + ": " + std::strerror(errno));
        }
        offset += static_cast<MPI_Offset>(type_header.size());
        const int local_size = static_cast<int>(local_atoms[it].size());
        int type_size = 0;
        int rank_offset = 0;
        MPI_Allreduce(&local_size, &type_size, 1, MPI_INT, MPI_SUM, comm);
        MPI_Exscan(&local_size, &rank_offset, 1, MPI_INT, MPI_SUM, comm);
        if (rank == 0) rank_offset = 0;
        int atom_data_ok = 1;
        if (local_size > 0)
        {
            atom_data_ok = write_at(file, local_atoms[it], offset + rank_offset) ? 1 : 0;
        }
        int all_atom_data_ok = 0;
        MPI_Allreduce(&atom_data_ok, &all_atom_data_ok, 1, MPI_INT, MPI_MIN, comm);
        if (all_atom_data_ok == 0)
        {
            close(file);
            throw std::runtime_error("Unable to write MdCell restart STRU atom data: " + fn + ": " + std::strerror(errno));
        }
        MPI_Barrier(comm);
        offset += type_size;
    }
    if (close(file) != 0)
    {
        throw std::runtime_error("Unable to close MdCell restart STRU file: " + fn + ": " + std::strerror(errno));
    }
#else
    std::ofstream output(fn.c_str());
    output << header;
    for (std::size_t it = 0; it < metadata.species.size(); ++it)
        output << mdcell_type_header(metadata.species[it]) << local_atoms[it];
#endif
}
}
