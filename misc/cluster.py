import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import sys

df = pd.read_csv('summary.txt', sep='\t')
df = df.drop(columns=['gene', 'loc', 'str'])
df = df.drop(columns=['introns', 'icost', 'ilen', 'elen'])
#df = df.drop(columns=['acc', 'don', 'emm', 'imm', 'ilen', 'elen'])
df = df.set_index('name')
sns.clustermap(df, metric ="euclidean", standard_scale=1, method="ward")
plt.show()
