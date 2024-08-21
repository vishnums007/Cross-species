import scanpy as sc
sc.settings.set_figure_params(scanpy=True, dpi=300, dpi_save=300, frameon=True, \
                              vector_friendly=True, fontsize=14, figsize=None, color_map=None, format='pdf', \
                              facecolor=None, transparent=False, ipython_format='png2x')

import pickle
import pandas as pd
import numpy as np
import os
from glob import glob

#set your working directory where the multiple_seed_results folder is present
os.chdir('/DATAFILES/Cross_species/Fulldatasets/Allcombined/')
print(os.getcwd())


## Load the final anndata output file from the best seed that you identified in the first step.    
selectedfile = glob("multiple_seed_results/saturn_results/*_seed_14.h5ad")
adata = sc.read(selectedfile[0])
adata
meta=adata.obs
meta.to_csv("saturnmetadata.csv")
adata.obs_names_make_unique()


## Now, trying incorporate the timepoint information to the saturn output file. For this you have to load the intial anndata of both the species
zebra= sc.read("Zebrafish_SC.h5ad")
zebra
mouse=sc.read("Mouse_SC.h5ad")
mouse
x=mouse.obs["Timepoint"]
x=np.array(x)
y=zebra.obs["Timepoint"]
y=np.array(y)
z=np.concatenate((x,y)) #check the saturn metadata to identify the order in which two species are integrated in the saturn metadata.
adata.obs["Timepoint"]=z #here I am adding mouse first and then zebrafish under the column Timepoint 
adata 

x=mouse.obs["Dataset"]
x=np.array(x)
y=zebra.obs["Dataset"]
y=np.array(y)
z=np.concatenate((x,y)) #check the saturn metadata to identify the order in which two species are integrated in the saturn metadata.
adata.obs["Dataset"]=z #here I am adding mouse first and then zebrafish under the column Timepoint 
adata 

adata.obs["Dataset"].unique()
adata.obs["Timepoint"].unique()
adata.obs["labels"].unique()

#Now we can visualize the integration
#sc.pp.normalize_total(adata)
#sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pl.pca(adata, color="Dataset", title="Dataset", save="pcaplot.pdf")
sc.pl.pca(adata, color="labels", title="Cell Type",save="pcaplot_labelledCelltype.pdf")
sc.pp.neighbors(adata,)
sc.tl.umap(adata)

sc.pl.umap(adata, color="species", title="Species", save="umap_plot_species_labelled.pdf")
sc.pl.umap(adata, color="Dataset", title="Species", save="umap_plot_Dataset_labelled.pdf")
sc.pl.umap(adata, color="labels2", title="Cell Type",save="umapplot_labelledCelltype.pdf",legend_loc='on data',legend_fontsize=3)
sc.pl.umap(adata, color="labels2", title="Cell Type",save="umapplot_labelledCelltype_v2.pdf")
sc.pl.umap(adata, color="labels", title="Cell Type",save="umapplot_labelledCelltypeSpecies.pdf",legend_loc='on data',legend_fontsize=3)
sc.pl.umap(adata, color="labels", title="Cell Type",save="umapplot_labelledCelltypeSpecies_v2.pdf")


##performing subset analysis



#Now, we can upload the corresponding macrogene weight to the anndata


selectedfile = glob("multiple_seed_results/saturn_results/*_seed_14_genes_to_macrogenes.pkl")

with open(selectedfile[0], "rb") as f:
    macrogene_weights = pickle.load(f)
    

adata.uns['weights'] = macrogene_weights


#adata = sc.read("OpticNerveCrush_zebrafishvsmouse_SaturnIntegrated_seed5.h5ad")

##Finding top DE markers for each cluster
sc.tl.rank_genes_groups(adata, groupby="labels", method="wilcoxon")

##Following is the function that will pull out individual genes in a given macrogene number along with corresponding gene weights
def get_scores(macrogene):
    '''
    Given the index of a macrogene number, return the scores by gene for that centroid
    '''
    scores = {}
    for (gene), score in macrogene_weights.items():
        scores[gene] = score[int(macrogene)]
    return scores


# The following loop will save csv files of top DE expressed genes in each cluster along with their gene weights
# Iterate over clusters

for cluster in adata.obs["labels"].unique():
    print(f"Cluster {cluster}:")
    markers = adata.uns["rank_genes_groups"]["names"][cluster]
    for macrogene in markers:
        macrogene_number = f"Macrogene {macrogene}"
        top_markers = pd.DataFrame(get_scores(macrogene).items(), columns=["gene", "weight"])\
                .sort_values("weight", ascending=False)\
                .head(50) #it is selecting top 50
        filename = f"cluster_{cluster}_topmarkers.csv"
        filename=filename.replace("/", "_")
        append_mode = "a" if macrogene != markers[0] else "w"
        
        #Now, writing the markers to the csv file
        with open(filename, mode=append_mode) as file:
            if append_mode == "w":
                file.write(macrogene_number + '\n')
                top_markers.to_csv(file, index=True, header=True)
                print(f"Created new file for Cluster {cluster}: {filename}")

            else:
                file.write(macrogene_number + '\n')
                top_markers.to_csv(file, index=True, header=False, mode='a')
                #print(f"Appended data for Macrogene {macrogene} in Cluster {cluster} to {filename}")

## Making dotplots for topmarkers in each cluster
for cluster in adata.obs["labels"].unique():
    print(f"Making dotplot for {cluster}")
    top_n_markers = 5  # Number of top markers to display
    clusters = [cluster]
    filename = "Dotplot_top5_"+cluster+"_species"+".pdf"
    filename=filename.replace("/", "_")
    sc.pl.rank_genes_groups_dotplot(adata, groupby='labels', groups=clusters, n_genes=top_n_markers, save=filename)

#adding levels to the anndata, this will useful for converting to the seurat object
new_levels=np.unique(sorted(adata.obs['labels2']))
adata.obs['labels2'] = pd.Categorical(adata.obs['labels2'], categories=new_levels)


adata.write("SpinalCords_zebrafishvsmouse_SaturnIntegrated_seed14.h5ad")
