
#include <Rcpp.h>
#include <random>
#include <numeric>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix perturb_cells_cpp(
    IntegerMatrix count_matrix,
    IntegerVector damaged_cell_selections,
    NumericVector damage_levels,
    IntegerVector non_mito_idx,
    IntegerVector ribo_idx,
    double ribosome_penalty,
    int seed,
    std::string depletion_type
) {
  // Validate depletion type
  if (depletion_type != "abun_random" &&
      depletion_type != "unf_random") {
    stop(
      "depletion_type must be 'abun_random' or 'unf_random'."
    );
  }

  std::mt19937 rng(seed);

  std::unordered_set<int> ribo_set(
    ribo_idx.begin(), ribo_idx.end()
  );

  for (int i = 0; i < damaged_cell_selections.size(); ++i) {
    int cell = damaged_cell_selections[i];

    // Extract non-zero non-mitochondrial gene counts
    std::vector<int> gene_indices;
    std::vector<int> gene_counts;
    int total_count = 0;

    for (int j = 0; j < non_mito_idx.size(); ++j) {
      int gene_idx = non_mito_idx[j];
      int count = count_matrix(gene_idx, cell);

      if (count > 0) {
        gene_indices.push_back(gene_idx);
        gene_counts.push_back(count);
        total_count += count;
      }
    }

    if (total_count == 0) continue;

    // Determine number of transcripts to remove
    double cell_damage_level = damage_levels[cell];

    int total_loss = std::round(
      cell_damage_level * total_count
    );

    total_loss = std::max(
      0, std::min(total_loss, total_count)
    );

    if (total_loss == 0) continue;

    // Expand genes into individual transcript instances
    std::vector<int> transcripts;
    transcripts.reserve(total_count);

    for (int j = 0; j < gene_indices.size(); ++j) {
      for (int k = 0; k < gene_counts[j]; ++k) {
        transcripts.push_back(gene_indices[j]);
      }
    }

    // Store the number of remaining transcripts per gene
    std::unordered_map<int, int> new_counts;

    // If all transcripts are to be removed, avoid sampling
    if (total_loss == total_count) {
      for (int gene_idx : gene_indices) {
        new_counts[gene_idx] = 0;
      }
    }

    // --------------------------------------------------
    // OPTION 1: Abundance-weighted transcript depletion
    // depletion_type = "abun_random"
    // --------------------------------------------------
    else if (depletion_type == "abun_random") {

      // Calculate gene-level weights proportional to abundance
      std::vector<double> gene_weights;

      for (int j = 0; j < gene_indices.size(); ++j) {
        double weight =
          static_cast<double>(gene_counts[j]) / total_count;

        // Apply ribosomal penalty
        if (ribo_set.find(gene_indices[j]) != ribo_set.end()) {
          weight *= ribosome_penalty;
        }

        gene_weights.push_back(weight);
      }

      // Normalise gene weights
      double weight_sum = std::accumulate(
        gene_weights.begin(), gene_weights.end(), 0.0
      );

      if (weight_sum <= 0.0) {
        stop(
          "All transcript weights are zero in a damaged cell."
        );
      }

      for (double& w : gene_weights) {
        w /= weight_sum;
      }

      // Each transcript inherits its gene's weight
      std::vector<double> transcript_weights;
      transcript_weights.reserve(total_count);

      for (int j = 0; j < gene_indices.size(); ++j) {
        for (int k = 0; k < gene_counts[j]; ++k) {
          transcript_weights.push_back(gene_weights[j]);
        }
      }

      // Weighted sampling without replacement
      std::discrete_distribution<> dist(
        transcript_weights.begin(),
        transcript_weights.end()
      );

      std::unordered_set<int> to_remove;

      while (to_remove.size() <
             static_cast<size_t>(total_loss)) {
        int sampled_index = dist(rng);
        to_remove.insert(sampled_index);
      }

      // Count transcripts that remain
      for (int j = 0; j < transcripts.size(); ++j) {
        if (to_remove.find(j) == to_remove.end()) {
          new_counts[transcripts[j]] += 1;
        }
      }
    }

    // --------------------------------------------------
    // OPTION 2: Uniform random transcript depletion
    // depletion_type = "unf_random"
    // --------------------------------------------------
    else if (depletion_type == "unf_random") {

      // Every transcript has equal probability of removal.
      // Ribosomal penalty is intentionally not applied.
      std::shuffle(
        transcripts.begin(), transcripts.end(), rng
      );

      // Count transcripts that remain after random removal
      for (int j = total_loss; j < total_count; ++j) {
        new_counts[transcripts[j]] += 1;
      }
    }

    // Write updated counts back to the matrix
    for (int j = 0; j < non_mito_idx.size(); ++j) {
      int gene_idx = non_mito_idx[j];
      count_matrix(gene_idx, cell) = new_counts[gene_idx];
    }
  }

  return count_matrix;
}