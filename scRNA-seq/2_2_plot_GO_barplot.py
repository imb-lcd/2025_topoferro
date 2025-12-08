#!/usr/bin/env python
# coding: utf-8

# In[11]:


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# In[40]:


# Load file
clust = "4"

go_file = "D:/Spatiotemporal_analysis/scRNA-seq/jh_figure/msigdb_df5_cl5_GOBPMF_sel2/cluster_" + clust + "_enrichment_trimmed.txt"
df = pd.read_csv(go_file, sep = "\t")

out_prefix = "D:/Spatiotemporal_analysis/scRNA-seq/jh_figure/msigdb_df5_cl5_GOBPMF_sel2/cluster_" + clust + "_enrich_barplot"


# In[41]:


# add log p adju
df['log.p.adj'] = -np.log(df['p.adjust'])


# In[42]:


# plot bar plot
plt.figure(figsize=(5, 4))
sns.barplot(
    data = df,
    x = "log.p.adj",
    y = "Description",
    color = "blue"
)
plt.xlabel("-log of adjusted p-value")
plt.ylabel("GO Term")
# plt.show()

plt.savefig(out_prefix + ".jpg")
plt.savefig(out_prefix + ".pdf")

