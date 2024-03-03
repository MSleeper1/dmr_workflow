# need to edit script to accept arguments, plot dmrs and avg methylation for designated region 

# import pandas as pd
# import matplotlib.pyplot as plt

# df = pd.read_csv("LFRYL_region_methylation.txt", delimiter="\t")

# fig, ax = plt.subplots()

# ax1 = df.plot(kind='scatter', x='start', y='colon-merged', color='r', ax=ax)
# ax2 = df.plot(kind='scatter', x='start', y='pancreas-merged', color='g', ax=ax)
# ax3 = df.plot(kind='scatter', x='start', y='lung-merged', color='b', ax=ax)

# ax.set_title("Methylation Percent for 3 Samples")
# ax.set_xlabel("Genome Location")
# ax.set_ylabel("Methylation Percent")

# print(ax1 == ax2 == ax3) #True

#--------------------------------------------------------

# fig, ax = plt.subplots()

# ax1 = df.plot(kind='line', x='start', y='colon-merged', color='r', ax=ax)
# ax2 = df.plot(kind='line', x='start', y='pancreas-merged', color='g', ax=ax)
# ax3 = df.plot(kind='line', x='start', y='lung-merged', color='b', ax=ax)

# ax.set_title("Methylation of tissue samples along FRYL gene on chromosome 4")
# ax.set_xlabel("Location on Chromosome 4")
# ax.set_ylabel("Mean Methylation Fraction")

# print(ax1 == ax2 == ax3) #True

#--------------------------------------------------------
# run in shell to get the genome downloaded
# pyensembl install --release 75 --species homo_sapiens

# import matplotlib.pyplot as plt
# import pyensembl

# # Load the hg19 reference genome
# ensembl = pyensembl.EnsemblRelease(release=75, species='homo_sapiens')

# # Get the FRYL gene
# gene = ensembl.genes_by_name('FRYL')[0]


# # Get the coordinates of the exons and introns
# exon_starts = gene.exon_starts
# exon_ends = gene.exon_ends
# intron_starts = gene.intron_starts
# intron_ends = gene.intron_ends

# # Create a new figure
# fig, ax = plt.subplots()

# # Plot the exons
# for i in range(len(exon_starts)):
#     ax.axvspan(exon_starts[i], exon_ends[i], color='blue', alpha=0.2)

# # Plot the introns
# for i in range(len(intron_starts)):
#     ax.axvspan(intron_starts[i], intron_ends[i], color='red', alpha=0.2)

# # Set the title and axis labels
# ax.set_title('FRYL gene exons and introns')
# ax.set_xlabel('Genomic position (bp)')
# ax.set_ylabel('Region')

# # Show the plot
# plt.show()

#-------------------------

import matplotlib.pyplot as plt
import pyensembl

# Load the hg19 reference genome
ensembl = pyensembl.EnsemblRelease(release=75, species='homo_sapiens')

# Get the gene object for FRYL
gene = ensembl.genes_by_name('FRYL')[0]

# Get the transcripts for the gene
transcripts = gene.transcripts

# Create lists to store the exon and intron coordinates
exon_starts = []
exon_ends = []
intron_starts = []
intron_ends = []

# Loop through each transcript and extract its exons and introns
for transcript in transcripts:
    # Exons
    for i in range(len(transcript.exon_intervals)):
        exon_start, exon_end = transcript.exon_intervals[i]
        exon_starts.append(exon_start)
        exon_ends.append(exon_end)
        # Introns (except for the last exon)
        if i < len(transcript.exon_intervals) - 1:
            intron_start = transcript.exon_intervals[i][1] + 1
            intron_end = transcript.exon_intervals[i+1][0] - 1
            intron_starts.append(intron_start)
            intron_ends.append(intron_end)

# Create a plot with the exons and introns
fig, ax = plt.subplots()
ax.plot(exon_starts, [1]*len(exon_starts), '|', color='black', markersize=20)
ax.plot(intron_starts, [1]*len(intron_starts), '<-', color='black', markersize=5)
ax.set_xlim([gene.start-5000, gene.end+5000])
ax.set_ylim([0,2])
ax.set_title('Exons and introns of FRYL')
ax.set_xlabel('Genomic position (hg19)')
plt.show()%
