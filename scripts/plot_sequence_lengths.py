import plotly.express as px
from Bio import SeqIO
import numpy as np

lengths = np.array([len(str(record.seq)) for record in SeqIO.parse("seqs/encapsulin_hits_filtered.fasta", "fasta")])

fig = px.histogram(x=lengths, color_discrete_sequence=px.colors.qualitative.Prism)

fig.update_layout(
    width=1200,
    height=600,
    font=dict(size=16),
    template='plotly_white',
    xaxis_title="Length",
    yaxis_title="Number of Sequences"
)

print(sum([length > 900 for length in lengths]))

fig.show()