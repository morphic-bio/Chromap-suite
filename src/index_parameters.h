#ifndef INDEX_PARAMETERS_H_
#define INDEX_PARAMETERS_H_

namespace chromap {

struct IndexParameters {
  int kmer_size = 17;
  int window_size = 7;
  int num_threads = 1;
  std::string reference_file_path;
  std::string index_output_file_path;
  // Optional binary reference sidecar written together with the index.
  std::string reference_sidecar_path;
};

}  // namespace chromap

#endif  // INDEX_PARAMETERS_H_
