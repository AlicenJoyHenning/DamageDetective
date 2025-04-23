// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include <random>
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
    int seed
) {
  std::mt19937 rng(seed);
  std::unordered_set<int> ribo_set(ribo_idx.begin(), ribo_idx.end());

  for (int i = 0; i < damaged_cell_selections.size(); ++i) {
    int cell = damaged_cell_selections[i];

    // Extract counts for the cell
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

    // Compute number of transcripts to remove
    int barcode_idx = cell;
    double cell_damage_level = damage_levels[barcode_idx];
    int total_loss = std::round(cell_damage_level * total_count);

    // Compute probabilities
    std::vector<double> probabilities;
    for (int j = 0; j < gene_indices.size(); ++j) {
      probabilities.push_back((double)gene_counts[j] / total_count);
    }

    // Penalize ribosomal genes
    for (int j = 0; j < gene_indices.size(); ++j) {
      if (ribo_set.find(gene_indices[j]) != ribo_set.end()) {
        probabilities[j] *= ribosome_penalty;
      }
    }

    // Normalize probabilities
    double prob_sum = std::accumulate(probabilities.begin(), probabilities.end(), 0.0);
    for (double& p : probabilities) p /= prob_sum;

    // Expand transcripts
    std::vector<int> transcripts;
    std::vector<double> expanded_probs;
    for (int j = 0; j < gene_indices.size(); ++j) {
      for (int k = 0; k < gene_counts[j]; ++k) {
        transcripts.push_back(gene_indices[j]);
        expanded_probs.push_back(probabilities[j]);
      }
    }

    // Sample transcripts to remove
    std::discrete_distribution<> dist(expanded_probs.begin(), expanded_probs.end());
    std::unordered_set<int> to_remove;
    while (to_remove.size() < (size_t)total_loss && to_remove.size() < transcripts.size()) {
      int sampled_index = dist(rng);
      to_remove.insert(sampled_index);
    }

    // Reset counts
    std::unordered_map<int, int> new_counts;
    for (int j = 0; j < transcripts.size(); ++j) {
      if (to_remove.find(j) == to_remove.end()) {
        new_counts[transcripts[j]] += 1;
      }
    }

    // Write back to matrix
    for (int j = 0; j < non_mito_idx.size(); ++j) {
      int gene_idx = non_mito_idx[j];
      count_matrix(gene_idx, cell) = new_counts[gene_idx];
    }
  }

  return count_matrix;
}
