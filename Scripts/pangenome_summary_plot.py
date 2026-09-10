import matplotlib.pyplot as plt

categories = ['Core (≥99%)', 'Soft core (95–99%)', 'Shell (15–95%)', 'Cloud (<15%)']
counts = [2769, 336, 3345, 29771]
colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D']

fig, ax = plt.subplots(figsize=(8, 5))
bars = ax.bar(categories, counts, color=colors, edgecolor='black', linewidth=0.5)

# Add count labels on bars
for bar, count in zip(bars, counts):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 500,
            f'{count:,}', ha='center', va='bottom', fontsize=10, fontweight='bold')

ax.set_ylabel('Number of gene clusters')
ax.set_title('Pangenome Composition of 444 E. coli Genomes', fontsize=14)
ax.set_yscale('log')  # log scale to see all categories clearly

plt.tight_layout()
plt.savefig('pangenome_summary_bar.png', dpi=300, bbox_inches='tight')
plt.savefig('pangenome_summary_bar.pdf', bbox_inches='tight')
print("Bar chart saved as pangenome_summary_bar.png and .pdf")
