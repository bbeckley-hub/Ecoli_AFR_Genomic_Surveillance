import pandas as pd
import matplotlib.pyplot as plt

mds_df = pd.read_csv('pangenome_mds_coordinates.csv')
colors = {'ST10': 'blue', 'ST131': 'red', 'ST167': 'green'}

plt.figure(figsize=(10, 8))
for st, color in colors.items():
    subset = mds_df[mds_df['ST'] == st]
    plt.scatter(subset['MDS1'], subset['MDS2'], c=color, label=st, alpha=0.7, s=30)
plt.xlabel('MDS1')
plt.ylabel('MDS2')
plt.title('Pangenome Accessory Gene Jaccard Distances')
plt.legend()
plt.savefig('pangenome_mds_plot.png', dpi=300, bbox_inches='tight')
print("Plot saved to pangenome_mds_plot.png")
