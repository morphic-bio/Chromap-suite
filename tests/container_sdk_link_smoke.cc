#include "libchromap.h"

int main() {
  chromap::ChromapRunResult (*entrypoint)(
      const chromap::MappingParameters &) = &chromap::RunAtacMapping;
  return entrypoint == nullptr ? 1 : 0;
}
