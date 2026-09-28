#!/usr/bin/env bash
set -euo pipefail
prefix="$(realpath "$1")"
out="$(realpath -m "$2")"
mkdir -p "$out"
cat > "$out/consumer.cc" <<'CPP'
#include "libchromap.h"
int main() {
  chromap::ChromapRunResult (*entry)(const chromap::MappingParameters &) = &chromap::RunAtacMapping;
  return entry == nullptr;
}
CPP
g++ -std=c++11 -fopenmp -I"$prefix/include/chromap-suite" -I"$prefix/include" \
  "$out/consumer.cc" "$prefix/lib/libchromap.a" "$prefix/lib/librapidmacs.a" \
  -lhts -lm -lz -lpthread -ldl -lcurl -lcrypto -lbz2 -llzma -ldeflate -o "$out/consumer"
"$out/consumer"
echo 'PASS: packaged SDK consumer'
