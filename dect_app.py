import streamlit as st
import base64
from Bio.Align import PairwiseAligner
from dect_core import detect_mutations, classify_mutations, pretty_alignment

# ======================
# Background
# ======================
def set_bg(img_file):
    with open(img_file, "rb") as f:
        encoded = base64.b64encode(f.read()).decode()

    st.markdown(
        f"""
        <style>
        .stApp {{
            background-image: url("data:image/jpg;base64,{encoded}");
            background-size: cover;
            background-attachment: fixed;
        }}
        </style>
        """,
        unsafe_allow_html=True
    )


set_bg("background.jpg")


# ======================
# UI styling
# ======================
st.markdown("""
<style>
[data-testid="stMainBlockContainer"] {
    background-color: rgba(255,255,255,0.85);
    padding: 2rem;
    border-radius: 20px;
    box-shadow: 0px 0px 30px rgba(0,0,0,0.3);
}
[data-testid="stMainBlockContainer"] * {
    color: black !important;
}
@media (prefers-color-scheme: dark) {
  [data-testid="stMainBlockContainer"] {
    background-color: rgba(30,30,30,0.85);
  }
  [data-testid="stMainBlockContainer"] * {
    color: white !important;
  }
}
</style>
""", unsafe_allow_html=True)


# ======================
# Helpers
# ======================
def seq_cleaner(seq):
    return "".join([line.strip() for line in seq.splitlines() if not line.startswith(">")])


def colored_alignment(alignment):
    ref = alignment.target
    query = alignment.query

    ref_line = ""
    mid_line = ""
    query_line = ""

    for r, q in zip(ref, query):

        if r == q:
            ref_line += f"<span style='color:green'>{r}</span>"
            query_line += f"<span style='color:green'>{q}</span>"
            mid_line += "|"

        elif r == '-' or q == '-':
            ref_line += f"<span style='color:orange'>{r}</span>"
            query_line += f"<span style='color:orange'>{q}</span>"
            mid_line += " "

        else:
            ref_line += f"<span style='color:red'>{r}</span>"
            query_line += f"<span style='color:red'>{q}</span>"
            mid_line += "*"

    return ref_line, mid_line, query_line


# ======================
# UI
# ======================
st.title("Mutation Detector")
st.markdown("***Alat za pronalaženje mutacija u DNK/CDS sekvenci***")
st.markdown("Laboratorija za molekularnu biologiju")
st.markdown("Departman za Biologiju i Ekologiju")
st.markdown("Prirodno-matematički Fakultet")
st.markdown("Univerzitet u Nišu")

seq_input1 = st.text_area("Uneti referentnu sekvencu (fasta/raw):", height=120)

seq_input2 = st.text_area("Uneti query sekvencu (fasta/raw):", height=120)

# ======================
# ALIGNER
# ======================
if seq_input1.strip() and seq_input2.strip():

    ref_seq = seq_cleaner(seq_input1)
    query_seq = seq_cleaner(seq_input2)

    from Bio import pairwise2
    alignments = pairwise2.align.globalms(ref_seq, query_seq, 2, -1, -2, -0.5)
    aln = alignments[0]

    st.markdown("***Prikaz pairwise globalnog poravnanja (Needleman-Wunsch algoritam):***")
    st.code(pretty_alignment(aln.seqA, aln.seqB))



# ======================
# ANALIZA
# ======================
if st.button("Analiziraj"):

    if not seq_input1 or not seq_input2:
        st.warning("Uneti obe sekvence DNK.")
    else:
        ref_seq = seq_cleaner(seq_input1).upper()
        query_seq = seq_cleaner(seq_input2).upper()

        st.markdown("***Prečišćene sekvence:***")
        st.code(f"REF:\n{ref_seq}\n\nQUERY:\n{query_seq}")

        # ======================
        # MUTACIJE
        # ======================
        mutations = detect_mutations(ref_seq, query_seq)

        st.markdown("***DNA mutacije:***")
        st.json(mutations)

        effects = classify_mutations(ref_seq, query_seq)

        st.markdown("***Protein mutacije:***")
        st.write(effects)

