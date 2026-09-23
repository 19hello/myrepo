#ifndef ESOLVER_NEP_GPU_NEIGHBOR_LIST_H
#define ESOLVER_NEP_GPU_NEIGHBOR_LIST_H

#include <string>
#include <vector>

namespace ModuleESolver
{

class NEP_GPU_NeighborList
{
  public:
    NEP_GPU_NeighborList() = default;
    ~NEP_GPU_NeighborList();

    NEP_GPU_NeighborList(const NEP_GPU_NeighborList&) = delete;
    NEP_GPU_NeighborList& operator=(const NEP_GPU_NeighborList&) = delete;

    bool build(int nall,
               double cutoff,
               const std::vector<double>& position,
               std::vector<int>& neighbor_count,
               std::vector<int>& neighbor_indices,
               int& max_neighbors,
               std::string& error) const;

    bool filter(int nall,
                double cutoff,
                const std::vector<double>& position,
                const std::vector<int>& candidate_count,
                const std::vector<int>& candidate_indices,
                int candidate_max_neighbors,
                std::vector<int>& neighbor_count,
                std::vector<int>& neighbor_indices,
                int& max_neighbors,
                std::string& error) const;

    bool build_device(int nall,
                      double cutoff,
                      const std::vector<double>& position,
                      int& max_neighbors,
                      std::string& error);

    bool filter_device(int nall,
                       double cutoff,
                       const std::vector<double>& position,
                       int candidate_max_neighbors,
                       int& max_neighbors,
                       std::string& error);

    void use_candidate_device();

    bool upload(int nall,
                int max_neighbors,
                const std::vector<int>& neighbor_count,
                const std::vector<int>& neighbor_indices,
                std::string& error);

    const void* device_neighbor_count() const;
    const void* device_neighbor_indices() const;
    int device_nall() const;
    int device_max_neighbors() const;

  private:
    void* candidate_device_count_ = nullptr;
    void* candidate_device_indices_ = nullptr;
    int candidate_device_nall_ = 0;
    int candidate_device_max_neighbors_ = 0;
    void* filtered_device_count_ = nullptr;
    void* filtered_device_indices_ = nullptr;
    int filtered_device_nall_ = 0;
    int filtered_device_max_neighbors_ = 0;
    void* device_neighbor_count_ = nullptr;
    void* device_neighbor_indices_ = nullptr;
    int device_nall_ = 0;
    int device_max_neighbors_ = 0;
};

} // namespace ModuleESolver

#endif
