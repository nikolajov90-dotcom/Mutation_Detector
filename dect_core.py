from Bio import pairwise2
from Bio.Seq import Seq

# ======================
# PRETTY ALIGNER
# ======================
def pretty_alignment(ref_aln, query_aln):

    match_line = ""

    for r, q in zip(ref_aln, query_aln):
        if r == q:
            match_line += "|"
        elif r == '-' or q == '-':
            match_line += " "
        else:
            match_line += "*"

    return f"{ref_aln}\n{match_line}\n{query_aln}"


# ======================
# CLEANER
# ======================
def seq_cleaner(seq_fasta):
    seq = ""
    for line in seq_fasta.splitlines():
        line = line.strip()
        if line.startswith(">"):
            continue
        line = line.replace(" ", "")
        seq += line
    return seq.upper()


# ======================
# DETEKCIJA MUTACIJA (DNA)
# ======================
def detect_mutations(ref_seq, query_seq):

    alignments = pairwise2.align.globalms(ref_seq, query_seq, 2, -1, -2, -0.5)
    aln = alignments[0]

    ref_aln = aln.seqA
    query_aln = aln.seqB

    mutations = {
        "SNPs": [],
        "insertions": [],
        "deletions": [],
        "frameshift": False
    }

    ref_pos = 0
    query_pos = 0
    current_indel = 0

    for r, q in zip(ref_aln, query_aln):

        # MATCH / SNP
        if r != '-' and q != '-':
            ref_pos += 1
            query_pos += 1

            if r != q:
                mutations["SNPs"].append((ref_pos, r, q))

            # proveri indel region
            if current_indel != 0:
                if current_indel % 3 != 0:
                    mutations["frameshift"] = True
                current_indel = 0

        # DELETION
        elif r != '-' and q == '-':
            ref_pos += 1
            current_indel += 1
            mutations["deletions"].append(ref_pos)

        # INSERTION
        elif r == '-' and q != '-':
            query_pos += 1
            current_indel += 1
            mutations["insertions"].append(query_pos)

    # poslednji indel check
    if current_indel != 0 and current_indel % 3 != 0:
        mutations["frameshift"] = True

    return mutations


# ======================
# PROTEIN EFEKTI
# ======================
def classify_mutations(ref_seq, query_seq):

    effects = []

    ref_seq = Seq(ref_seq)
    query_seq = Seq(query_seq)

    # trim na 3
    min_len = min(len(ref_seq), len(query_seq))
    trim_len = (min_len // 3) * 3

    ref_seq = ref_seq[:trim_len]
    query_seq = query_seq[:trim_len]

    for i in range(0, trim_len, 3):

        ref_codon = ref_seq[i:i+3]
        query_codon = query_seq[i:i+3]

        # skip ako ima gap
        if "-" in ref_codon or "-" in query_codon:
            continue

        ref_aa = ref_codon.translate()
        query_aa = query_codon.translate()

        pos = i // 3 + 1

        if ref_codon != query_codon:

            if ref_aa == query_aa:
                effects.append((pos, "silent", str(ref_codon), str(query_codon)))

            elif query_aa == "*":
                effects.append((pos, "nonsense", str(ref_codon), str(query_codon)))

            else:
                effects.append((pos, "missense", str(ref_codon), str(query_codon)))

    return effects
