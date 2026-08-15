#include <cstdlib>
#include <iostream>

#include "index.h"

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "usage: " << argv[0] << " INDEX\n";
    return EXIT_FAILURE;
  }

  chromap::Index index(argv[1]);
  index.Load();
  return EXIT_SUCCESS;
}
