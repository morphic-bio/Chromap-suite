#include <cerrno>
#include <climits>
#include <cstdlib>
#include <iostream>
#include <string>

#include "sequence_batch.h"

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "usage: " << argv[0] << " REFERENCE_SIDECAR THREADS\n";
    return EXIT_FAILURE;
  }

  errno = 0;
  char *end = nullptr;
  const long parsed_threads = std::strtol(argv[2], &end, 10);
  if (errno != 0 || end == argv[2] || *end != '\0' || parsed_threads < 1 ||
      parsed_threads > INT_MAX) {
    std::cerr << "invalid thread count: " << argv[2] << "\n";
    return EXIT_FAILURE;
  }

  chromap::SequenceBatch reference;
  chromap::MaterializedReferenceInfo info;
  std::string error;
  if (!reference.LoadMaterializedReference(
          argv[1], static_cast<int>(parsed_threads), &info, &error)) {
    std::cerr << "materialized reference load failed: " << error << "\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
