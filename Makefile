CXX=g++

# htslib submodule
HTSLIB_DIR=third_party/htslib
ifeq (,$(wildcard $(HTSLIB_DIR)/htslib/sam.h))
$(error "htslib submodule not initialized. Run: git submodule update --init --recursive")
endif

CXXFLAGS=-std=c++11 -Wall -O3 -fopenmp -msse4.1 -I$(HTSLIB_DIR)
LDFLAGS=-L$(HTSLIB_DIR) -lhts -lm -lz -lpthread -lcurl -lcrypto -lbz2 -llzma -ldeflate

core_cpp_source=sequence_batch.cc index.cc minimizer_generator.cc candidate_processor.cc alignment.cc feature_barcode_matrix.cc ksw.cc draft_mapping_generator.cc mapping_generator.cc mapping_writer.cc overflow_writer.cc overflow_reader.cc bam_sorter.cc chromap.cc
driver_cpp_source=chromap_driver.cc
libchromap_cpp_source=libchromap.cc
runner_cpp_source=chromap_lib_runner.cc
src_dir=src
objs_dir=objs
core_objs=$(patsubst %.cc,$(objs_dir)/%.o,$(core_cpp_source))
driver_objs=$(patsubst %.cc,$(objs_dir)/%.o,$(driver_cpp_source))
libchromap_objs=$(patsubst %.cc,$(objs_dir)/%.o,$(libchromap_cpp_source))
runner_objs=$(patsubst %.cc,$(objs_dir)/%.o,$(runner_cpp_source))

exec=chromap
libchromap=libchromap.a
runner=chromap_lib_runner
peak_caller=chromap_callpeaks
peak_caller_lib_cpp_source=peak_caller/fragment_input.cc peak_caller/frag_compact_store.cc \
	peak_caller/binned_signal.cc \
	peak_caller/call_peaks.cc peak_caller/peak_io.cc peak_caller/bdgpeakcall.cc \
	peak_caller/macs3_frag_peak_pipeline.cc peak_caller/macs3_frag_workspace.cc \
	peak_caller/stage_profile.cc
peak_caller_lib_objs=$(patsubst %.cc,$(objs_dir)/%.o,$(peak_caller_lib_cpp_source))
peak_caller_standalone_cpp_source=chromap_callpeaks.cc
peak_caller_standalone_obj=$(patsubst %.cc,$(objs_dir)/%.o,$(peak_caller_standalone_cpp_source))

ifneq ($(asan),)
	CXXFLAGS+=-fsanitize=address -g
	LDFLAGS+=-fsanitize=address -ldl -g
endif

ifneq ($(LEGACY_OVERFLOW),)
	CXXFLAGS+=-DLEGACY_OVERFLOW
endif

all: dir $(exec) $(peak_caller)
	
dir:
	mkdir -p $(objs_dir) $(objs_dir)/peak_caller

$(exec): $(core_objs) $(driver_objs) $(peak_caller_lib_objs)
	$(CXX) $(CXXFLAGS) $(core_objs) $(driver_objs) $(peak_caller_lib_objs) -o $(exec) $(LDFLAGS)

$(libchromap): $(core_objs) $(libchromap_objs) $(objs_dir)/peak_caller/frag_compact_store.o
	ar rcs $(libchromap) $(core_objs) $(libchromap_objs) $(objs_dir)/peak_caller/frag_compact_store.o

$(runner): $(libchromap) $(runner_objs)
	$(CXX) $(CXXFLAGS) $(runner_objs) $(libchromap) -o $(runner) $(LDFLAGS)

$(peak_caller): $(peak_caller_lib_objs) $(peak_caller_standalone_obj) | dir
	$(CXX) $(CXXFLAGS) $(peak_caller_lib_objs) $(peak_caller_standalone_obj) -o $(peak_caller) $(LDFLAGS)
	
# mkdir ensures objs/peak_caller/ (and any future subdir) exists when only
# a subset of targets (e.g. make chromap_callpeaks) is built.
$(objs_dir)/%.o: $(src_dir)/%.cc
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -I$(src_dir) -c $< -o $@

.PHONY: clean test-unit test-frag-compact-store test-peak-memory-source-100k \
	benchmark-peak-memory-fullset test-peak-100k test-peak-calibration-100k \
	test-peak-input-repr-100k test-peak-pileup-100k test-peak-frag-pileup-100k \
	test-peak-lambda-100k test-peak-score-100k test-peak-bdgpeakcall-100k \
	test-peak-narrowpeak-100k test-peak-integration-100k
clean:
	-rm -rf $(exec) $(libchromap) $(runner) $(peak_caller) $(objs_dir)

# 100K fragment peak-caller benchmark. Pair inputs: CHROMAP_PEAK_RUN_ROOT, or
# FRAGMENTS_TSV_GZ+ATAC_BAM, or same-run auto-pair under CHROMAP_100K_BENCH; RUN_MACS3=0 for internal only.
test-peak-100k: chromap_callpeaks
	RUN_MACS3=0 ./tests/run_peak_caller_100k.sh

# Parameter grid (CALIB_MODE=reduced by default for faster iteration). Override OUTDIR, CALIB_MODE, RUN_MACS3 as needed.
test-peak-calibration-100k: chromap_callpeaks
	CALIB_MODE=reduced RUN_MACS3=0 ./tests/run_peak_caller_calibration_100k.sh

# Phase 1: MACS3 tagAlign vs FRAG vs FRAG --max-count 1 (requires macs3, samtools, bedtools, benchmark pair).
test-peak-input-repr-100k: chromap_callpeaks
	./tests/run_macs3_input_repr_100k.sh

# Phase 2: MACS3 pileup -f FRAG vs chromap_callpeaks --pileup-bdg (requires macs3, benchmark fragments).
test-peak-pileup-100k: chromap_callpeaks
	./tests/run_pileup_parity_100k.sh

# Phase 2B: MACS3 pileup -f FRAG vs chromap_callpeaks --frag-span-pileup-bdg (requires macs3, benchmark fragments).
test-peak-frag-pileup-100k: chromap_callpeaks
	./tests/run_frag_span_pileup_parity_100k.sh

# Phase 3: MACS3 callpeak -f FRAG -B treat + control_lambda vs C++ (requires macs3; RUN_MACS3=0 = fragments-only ok).
test-peak-lambda-100k: chromap_callpeaks
	RUN_MACS3=0 ./tests/run_lambda_parity_100k.sh

# Phase 4: MACS3 bdgcmp ppois + FE vs C++ diagnostic scores (requires macs3 + benchmark fragments).
test-peak-score-100k: chromap_callpeaks
	RUN_MACS3=0 ./tests/run_score_parity_100k.sh

# Phase 5: MACS3 bdgpeakcall vs C++ diagnostic region caller (requires macs3 + benchmark fragments).
test-peak-bdgpeakcall-100k: chromap_callpeaks
	RUN_MACS3=0 ./tests/run_bdgpeakcall_parity_100k.sh

# Phase 6: MACS3 callpeak narrowPeak/summits vs C++ FRAG pipeline (macs3 + benchmark).
test-peak-narrowpeak-100k: chromap_callpeaks
	RUN_MACS3=0 ./tests/run_narrowpeak_parity_100k.sh

# Integrated chromap opt-in peak caller vs chromap_callpeaks + MACS3 BED3 (100K fixture).
test-peak-integration-100k: chromap chromap_callpeaks
	RUN_MACS3=0 ./tests/run_chromap_peak_integration_100k.sh

# Unit test for Y-filtering
test-unit: dir
	@mkdir -p tests
	$(CXX) $(CXXFLAGS) -I$(src_dir) tests/test_y_filter.cc $(src_dir)/sequence_batch.cc -o tests/test_y_filter $(LDFLAGS)
	./tests/test_y_filter

# Unit tests for in-memory fragment packing / accumulator
test-frag-compact-store: dir $(objs_dir)/peak_caller/frag_compact_store.o
	@mkdir -p tests
	$(CXX) $(CXXFLAGS) -I$(src_dir) tests/test_frag_compact_store.cc \
		$(objs_dir)/peak_caller/frag_compact_store.o -o tests/test_frag_compact_store $(LDFLAGS)
	./tests/test_frag_compact_store

# 100K: memory vs file fragment source for integrated MACS3 FRAG peaks
test-peak-memory-source-100k: chromap chromap_callpeaks
	RUN_MACS3=0 ./tests/run_chromap_peak_memory_source_100k.sh

# Full (unsampled) 4-lane ATAC: file vs memory peak source, GNU time -f wall/user/sys/max RSS.
# Tighten BAM sorter with SORT_BAM_RAM=512M (default in script) so compact-store cost is more visible.
# Long-running. Override e.g. FIXTURE_ATAC=.../fixture/atac for a faster smoke, or RUN_SET=file.
benchmark-peak-memory-fullset: chromap chromap_callpeaks
	bash $(CURDIR)/tests/benchmark_chromap_peak_memory_fullset.sh
