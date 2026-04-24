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

ifneq ($(asan),)
	CXXFLAGS+=-fsanitize=address -g
	LDFLAGS+=-fsanitize=address -ldl -g
endif

ifneq ($(LEGACY_OVERFLOW),)
	CXXFLAGS+=-DLEGACY_OVERFLOW
endif

all: dir $(exec) 
	
dir:
	mkdir -p $(objs_dir)

$(exec): $(core_objs) $(driver_objs)
	$(CXX) $(CXXFLAGS) $(core_objs) $(driver_objs) -o $(exec) $(LDFLAGS)

$(libchromap): $(core_objs) $(libchromap_objs)
	ar rcs $(libchromap) $(core_objs) $(libchromap_objs)

$(runner): $(libchromap) $(runner_objs)
	$(CXX) $(CXXFLAGS) $(runner_objs) $(libchromap) -o $(runner) $(LDFLAGS)
	
$(objs_dir)/%.o: $(src_dir)/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean test-unit
clean:
	-rm -rf $(exec) $(libchromap) $(runner) $(objs_dir)

# Unit test for Y-filtering
test-unit: dir
	@mkdir -p tests
	$(CXX) $(CXXFLAGS) -I$(src_dir) tests/test_y_filter.cc $(src_dir)/sequence_batch.cc -o tests/test_y_filter $(LDFLAGS)
	./tests/test_y_filter
