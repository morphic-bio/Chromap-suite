#include "mapping_generator.h"

namespace chromap {

namespace {
thread_local const std::unordered_set<uint32_t> *tls_y_contig_rids = nullptr;
thread_local std::vector<uint32_t> *tls_y_hit_read_ids = nullptr;
}  // namespace

void SetThreadYHitTracking(const std::unordered_set<uint32_t> *y_contig_rids,
                           std::vector<uint32_t> *thread_y_hit_read_ids) {
  tls_y_contig_rids = y_contig_rids;
  tls_y_hit_read_ids = thread_y_hit_read_ids;
}

const std::unordered_set<uint32_t> *GetThreadYContigRids() {
  return tls_y_contig_rids;
}

std::vector<uint32_t> *GetThreadYHitReadIds() {
  return tls_y_hit_read_ids;
}

// For strand, kPositive is 1, kNegative is 0;
template <>
void MappingGenerator<MappingWithoutBarcode>::EmplaceBackSingleEndMappingRecord(
    MappingInMemory &mapping_in_memory,
    std::vector<std::vector<MappingWithoutBarcode>>
        &mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs[mapping_in_memory.rid].emplace_back(
      mapping_in_memory.read_id, mapping_in_memory.GetFragmentStartPosition(),
      mapping_in_memory.GetFragmentLength(), mapping_in_memory.mapq,
      mapping_in_memory.GetStrand(), mapping_in_memory.is_unique,
      /*num_dups=*/1);
}

template <>
void MappingGenerator<MappingWithBarcode>::EmplaceBackSingleEndMappingRecord(
    MappingInMemory &mapping_in_memory,
    std::vector<std::vector<MappingWithBarcode>> &mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs[mapping_in_memory.rid].emplace_back(
      mapping_in_memory.read_id, mapping_in_memory.barcode_key,
      mapping_in_memory.GetFragmentStartPosition(),
      mapping_in_memory.GetFragmentLength(), mapping_in_memory.mapq,
      mapping_in_memory.GetStrand(), mapping_in_memory.is_unique,
      /*num_dups=*/1);
}

template <>
void MappingGenerator<PAFMapping>::EmplaceBackSingleEndMappingRecord(
    MappingInMemory &mapping_in_memory,
    std::vector<std::vector<PAFMapping>> &mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs[mapping_in_memory.rid].emplace_back(
      mapping_in_memory.read_id, std::string(mapping_in_memory.read_name),
      mapping_in_memory.read_length,
      mapping_in_memory.GetFragmentStartPosition(),
      mapping_in_memory.GetFragmentLength(), mapping_in_memory.mapq,
      mapping_in_memory.GetStrand(), mapping_in_memory.is_unique,
      /*num_dups=*/1);
}

template <>
void MappingGenerator<SAMMapping>::EmplaceBackSingleEndMappingRecord(
    MappingInMemory &mapping_in_memory,
    std::vector<std::vector<SAMMapping>> &mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs[mapping_in_memory.rid].emplace_back(
      mapping_in_memory.read_id, std::string(mapping_in_memory.read_name),
      mapping_in_memory.barcode_key, /*num_dups=*/1,
      mapping_in_memory.GetFragmentStartPosition(), mapping_in_memory.rid,
      /*mpos=*/0, /*mrid=*/-1, /*tlen=*/0, 
      mapping_in_memory.SAM_flag, mapping_in_memory.GetStrand(),
      /*is_alt=*/0, mapping_in_memory.is_unique, mapping_in_memory.mapq,
      mapping_in_memory.NM, mapping_in_memory.n_cigar, mapping_in_memory.cigar,
      mapping_in_memory.MD_tag, std::string(mapping_in_memory.read_sequence),
      mapping_in_memory.qual_sequence ? std::string(mapping_in_memory.qual_sequence) : std::string());
}

template <>
void MappingGenerator<PairedEndMappingWithBarcode>::
    EmplaceBackSingleEndMappingRecord(
        MappingInMemory &mapping_in_memory,
        std::vector<std::vector<PairedEndMappingWithBarcode>>
            &mappings_on_diff_ref_seqs) = delete;

template <>
void MappingGenerator<PairedEndAtacDualMapping>::EmplaceBackSingleEndMappingRecord(
    MappingInMemory &mapping_in_memory,
    std::vector<std::vector<PairedEndAtacDualMapping>>
        &mappings_on_diff_ref_seqs) = delete;

template <>
void MappingGenerator<PairedEndMappingWithoutBarcode>::
    EmplaceBackSingleEndMappingRecord(
        MappingInMemory &mapping_in_memory,
        std::vector<std::vector<PairedEndMappingWithoutBarcode>>
            &mappings_on_diff_ref_seqs) = delete;

template <>
void MappingGenerator<PairedPAFMapping>::EmplaceBackSingleEndMappingRecord(
    MappingInMemory &mapping_in_memory,
    std::vector<std::vector<PairedPAFMapping>> &mappings_on_diff_ref_seqs) =
    delete;

template <>
void MappingGenerator<PairsMapping>::EmplaceBackSingleEndMappingRecord(
    MappingInMemory &mapping_in_memory,
    std::vector<std::vector<PairsMapping>> &mappings_on_diff_ref_seqs) = delete;

template <>
void MappingGenerator<SAMMapping>::EmplaceBackPairedEndMappingRecord(
    PairedEndMappingInMemory &paired_end_mapping_in_memory,
    std::vector<std::vector<SAMMapping>> &mappings_on_diff_ref_seqs) {
  int tlen = (int)paired_end_mapping_in_memory.GetFragmentLength();
  for (int i = 0; i < 2; ++i) {
    MappingInMemory &mapping_in_memory = (i == 0 ? paired_end_mapping_in_memory.mapping_in_memory1 :
        paired_end_mapping_in_memory.mapping_in_memory2);
    MappingInMemory &mate_mapping_in_memory = (i == 0 ? paired_end_mapping_in_memory.mapping_in_memory2 :
        paired_end_mapping_in_memory.mapping_in_memory1);
  
    mappings_on_diff_ref_seqs[mapping_in_memory.rid].emplace_back(
      mapping_in_memory.read_id, std::string(mapping_in_memory.read_name),
      mapping_in_memory.barcode_key, /*num_dups=*/1,
      mapping_in_memory.GetFragmentStartPosition(), mapping_in_memory.rid,
      /*mpos=*/mate_mapping_in_memory.GetFragmentStartPosition(), 
      /*mrid=*/mate_mapping_in_memory.rid, 
      /*tlen=*/mapping_in_memory.GetStrand() ? tlen : -tlen, 
      mapping_in_memory.SAM_flag, mapping_in_memory.GetStrand(),
      /*is_alt=*/0, mapping_in_memory.is_unique, mapping_in_memory.mapq,
      mapping_in_memory.NM, mapping_in_memory.n_cigar, mapping_in_memory.cigar,
      mapping_in_memory.MD_tag, std::string(mapping_in_memory.read_sequence),
      mapping_in_memory.qual_sequence ? std::string(mapping_in_memory.qual_sequence) : std::string());
    mapping_in_memory.cigar = nullptr;
    mapping_in_memory.n_cigar = 0;
  }
}

template <>
void MappingGenerator<PairedEndMappingWithoutBarcode>::
    EmplaceBackPairedEndMappingRecord(
        PairedEndMappingInMemory &paired_end_mapping_in_memory,
        std::vector<std::vector<PairedEndMappingWithoutBarcode>>
            &mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs[paired_end_mapping_in_memory.mapping_in_memory1.rid]
      .emplace_back(paired_end_mapping_in_memory.GetReadId(),
                    paired_end_mapping_in_memory.GetFragmentStartPosition(),
                    paired_end_mapping_in_memory.GetFragmentLength(),
                    paired_end_mapping_in_memory.mapq,
                    paired_end_mapping_in_memory.GetStrand(),
                    paired_end_mapping_in_memory.is_unique, /*num_dups=*/1,
                    paired_end_mapping_in_memory.GetPositiveAlignmentLength(),
                    paired_end_mapping_in_memory.GetNegativeAlignmentLength());
}

template <>
void MappingGenerator<PairedEndMappingWithBarcode>::
    EmplaceBackPairedEndMappingRecord(
        PairedEndMappingInMemory &paired_end_mapping_in_memory,
        std::vector<std::vector<PairedEndMappingWithBarcode>>
            &mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs[paired_end_mapping_in_memory.mapping_in_memory1.rid]
      .emplace_back(paired_end_mapping_in_memory.GetReadId(),
                    paired_end_mapping_in_memory.GetBarcode(),
                    paired_end_mapping_in_memory.GetFragmentStartPosition(),
                    paired_end_mapping_in_memory.GetFragmentLength(),
                    paired_end_mapping_in_memory.mapq,
                    paired_end_mapping_in_memory.GetStrand(),
                    paired_end_mapping_in_memory.is_unique, /*num_dups=*/1,
                    paired_end_mapping_in_memory.GetPositiveAlignmentLength(),
                    paired_end_mapping_in_memory.GetNegativeAlignmentLength());
  for (auto *m : {&paired_end_mapping_in_memory.mapping_in_memory1,
                  &paired_end_mapping_in_memory.mapping_in_memory2}) {
    if (m->cigar) {
      free(m->cigar);
      m->cigar = nullptr;
    }
    m->n_cigar = 0;
  }
}

template <>
void MappingGenerator<PairedEndAtacDualMapping>::EmplaceBackPairedEndMappingRecord(
    PairedEndMappingInMemory &paired_end_mapping_in_memory,
    std::vector<std::vector<PairedEndAtacDualMapping>>
        &mappings_on_diff_ref_seqs) {
  PairedEndMappingWithBarcode bed(
      paired_end_mapping_in_memory.GetReadId(),
      paired_end_mapping_in_memory.GetBarcode(),
      paired_end_mapping_in_memory.GetFragmentStartPosition(),
      paired_end_mapping_in_memory.GetFragmentLength(),
      paired_end_mapping_in_memory.mapq,
      paired_end_mapping_in_memory.GetStrand(),
      paired_end_mapping_in_memory.is_unique, /*num_dups=*/1,
      paired_end_mapping_in_memory.GetPositiveAlignmentLength(),
      paired_end_mapping_in_memory.GetNegativeAlignmentLength());
  const int tlen =
      static_cast<int>(paired_end_mapping_in_memory.GetFragmentLength());
  SAMMapping sam_a;
  SAMMapping sam_b;
  for (int i = 0; i < 2; ++i) {
    MappingInMemory &mapping_in_memory =
        (i == 0 ? paired_end_mapping_in_memory.mapping_in_memory1
                : paired_end_mapping_in_memory.mapping_in_memory2);
    MappingInMemory &mate_mapping_in_memory =
        (i == 0 ? paired_end_mapping_in_memory.mapping_in_memory2
                : paired_end_mapping_in_memory.mapping_in_memory1);
    SAMMapping sam(
        mapping_in_memory.read_id, std::string(mapping_in_memory.read_name),
        mapping_in_memory.barcode_key, /*num_dups=*/1,
        mapping_in_memory.GetFragmentStartPosition(), mapping_in_memory.rid,
        /*mpos=*/mate_mapping_in_memory.GetFragmentStartPosition(),
        /*mrid=*/mate_mapping_in_memory.rid,
        /*tlen=*/mapping_in_memory.GetStrand() ? tlen : -tlen,
        mapping_in_memory.SAM_flag, mapping_in_memory.GetStrand(),
        /*is_alt=*/0, mapping_in_memory.is_unique, mapping_in_memory.mapq,
        mapping_in_memory.NM, mapping_in_memory.n_cigar, mapping_in_memory.cigar,
        mapping_in_memory.MD_tag, std::string(mapping_in_memory.read_sequence),
        mapping_in_memory.qual_sequence
            ? std::string(mapping_in_memory.qual_sequence)
            : std::string());
    if (i == 0) {
      sam_a = std::move(sam);
    } else {
      sam_b = std::move(sam);
    }
  }
  mappings_on_diff_ref_seqs[paired_end_mapping_in_memory.mapping_in_memory1.rid]
      .emplace_back(bed, std::move(sam_a), std::move(sam_b));
  paired_end_mapping_in_memory.mapping_in_memory1.cigar = nullptr;
  paired_end_mapping_in_memory.mapping_in_memory1.n_cigar = 0;
  paired_end_mapping_in_memory.mapping_in_memory2.cigar = nullptr;
  paired_end_mapping_in_memory.mapping_in_memory2.n_cigar = 0;
}

template <>
void MappingGenerator<PairedPAFMapping>::EmplaceBackPairedEndMappingRecord(
    PairedEndMappingInMemory &paired_end_mapping_in_memory,
    std::vector<std::vector<PairedPAFMapping>> &mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs[paired_end_mapping_in_memory.mapping_in_memory1.rid]
      .emplace_back(
          paired_end_mapping_in_memory.GetReadId(),
          std::string(
              paired_end_mapping_in_memory.mapping_in_memory1.read_name),
          std::string(
              paired_end_mapping_in_memory.mapping_in_memory2.read_name),
          paired_end_mapping_in_memory.mapping_in_memory1.read_length,
          paired_end_mapping_in_memory.mapping_in_memory2.read_length,
          paired_end_mapping_in_memory.GetFragmentStartPosition(),
          paired_end_mapping_in_memory.GetNegativeAlignmentLength(),
          paired_end_mapping_in_memory.GetFragmentLength(),
          paired_end_mapping_in_memory.GetPositiveAlignmentLength(),
          paired_end_mapping_in_memory.mapq,
          paired_end_mapping_in_memory.mapping_in_memory1.mapq,
          paired_end_mapping_in_memory.mapping_in_memory2.mapq,
          paired_end_mapping_in_memory.GetStrand(),
          paired_end_mapping_in_memory.is_unique, /*num_dups=*/1);
}

template <>
void MappingGenerator<PairsMapping>::EmplaceBackPairedEndMappingRecord(
    PairedEndMappingInMemory &paired_end_mapping_in_memory,
    std::vector<std::vector<PairsMapping>> &mappings_on_diff_ref_seqs) {
  uint8_t strand1 = paired_end_mapping_in_memory.mapping_in_memory1.GetStrand();
  uint8_t strand2 = paired_end_mapping_in_memory.mapping_in_memory2.GetStrand();

  int position1 =
      paired_end_mapping_in_memory.mapping_in_memory1.ref_start_position;
  int position2 =
      paired_end_mapping_in_memory.mapping_in_memory2.ref_start_position;

  if (paired_end_mapping_in_memory.mapping_in_memory1.strand == kNegative) {
    position1 =
        paired_end_mapping_in_memory.mapping_in_memory1.ref_end_position;
  }

  if (paired_end_mapping_in_memory.mapping_in_memory2.strand == kNegative) {
    position2 =
        paired_end_mapping_in_memory.mapping_in_memory2.ref_end_position;
  }

  int rid1 = paired_end_mapping_in_memory.mapping_in_memory1.rid;
  int rid2 = paired_end_mapping_in_memory.mapping_in_memory2.rid;
  const int rid1_rank = pairs_custom_rid_rank_[rid1];
  const int rid2_rank = pairs_custom_rid_rank_[rid2];

  const bool is_rid1_rank_smaller =
      rid1_rank < rid2_rank || (rid1 == rid2 && position1 < position2);
  if (!is_rid1_rank_smaller) {
    std::swap(rid1, rid2);
    std::swap(position1, position2);
    std::swap(strand1, strand2);
  }

  mappings_on_diff_ref_seqs[rid1].emplace_back(
      paired_end_mapping_in_memory.GetReadId(),
      std::string(paired_end_mapping_in_memory.mapping_in_memory1.read_name),
      paired_end_mapping_in_memory.GetBarcode(), rid1, rid2, position1,
      position2, strand1, strand2, paired_end_mapping_in_memory.mapq,
      paired_end_mapping_in_memory.is_unique, /*num_dups=*/1);
}

template <>
void MappingGenerator<MappingWithBarcode>::EmplaceBackPairedEndMappingRecord(
    PairedEndMappingInMemory &paired_end_mapping_in_memory,
    std::vector<std::vector<MappingWithBarcode>> &mappings_on_diff_ref_seqs) =
    delete;

template <>
void MappingGenerator<MappingWithoutBarcode>::EmplaceBackPairedEndMappingRecord(
    PairedEndMappingInMemory &paired_end_mapping_in_memory,
    std::vector<std::vector<MappingWithoutBarcode>>
        &mappings_on_diff_ref_seqs) = delete;

template <>
void MappingGenerator<PAFMapping>::EmplaceBackPairedEndMappingRecord(
    PairedEndMappingInMemory &paired_end_mapping_in_memory,
    std::vector<std::vector<PAFMapping>> &mappings_on_diff_ref_seqs) = delete;

}  // namespace chromap
