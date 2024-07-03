# Python script for training an scANVI model for seed labeling
# Updated from scvi_train-scArches.py. Customized for BigRef project
# Last updated: 14May2024
# jrose

import scanpy as sc
import os
import scvi
import torch

# REPLACE THE FOLLOWING
###################
data_dir = "/home/Projects/Scharer_sc/scAtlas_ref/analysis/ref_objs"
data_file = "BcellRefAtlas_SLE.h5ad"
scvi_model_file = "BigRef_SLE_scvi_model_seedlbl"
scANVI_model_file = "BigRef_SLE_scANVI_model_seedlbl"
scvi_aData_file="BigRef_SLE_Bcell_AnnDataSCVI_seedlbl.h5ad"
scANVI_aData_file="BigRef_SLE_Bcell_AnnData_scANVI_seedlbl.h5ad"
###################

#Training scvi and scANVI, along with optimized preprocoessing settings. You many need to change preprocessing per individual project
def train_scvi_model(data_dir, data_file, scvi_model_file, scANVI_model_file, scvi_aData_file, scANVI_aData_file):
    scvi.settings.seed = 1990
    print("Last run with scvi-tools version:", scvi.__version__)

    adata_path = os.path.join(data_dir, data_file)
    # load data
    print("Loading data")
    adata = sc.read(adata_path)
    print("Data Loaded")
     # Preprocess the data

    print("Preprocessing")
    sc.pp.filter_cells(adata, min_counts=1)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata,
                            min_mean=0.0125,
                            max_mean=50,
                            min_disp=0.25,
                            subset=True,
                            batch_key="batch"
                           )
    print("Preprocessing Complete")
    # Setup data for scvi
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        #categorical_covariate_keys=["run"],
	    batch_key="batch",
        labels_key="labels"
        #continuous_covariate_keys=["percent_mt"],
    )

    # Define model
    model = scvi.model.SCVI(adata,
        n_latent=30,
        n_layers=2
    )
    model

    # Train model
    print("Starting training")
    model.train()
    print("Training complete")
    #Save trained model & adata
    model.save(scvi_model_file, overwrite=True)
    adata.write(scvi_aData_file)

    print("Starting scANVI")
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
    	model,
    	adata=adata,
    	unlabeled_category="Unknown",
    )
    print("Training scANVI model")
    scanvi_model.train(max_epochs=40)
    print("Training complete")

    print("Saving scANVI model")
    scanvi_model.save(scANVI_model_file)
    adata.write(scANVI_aData_file)


#if __name__ == "__main__":
    # Set up the argument parser
#    parser = argparse.ArgumentParser(description='Run SCVI analysis with given directories and file.')
#    parser.add_argument('data_dir', type=str, nargs='?', default=os.getcwd(), help='Directory containing the data')
#    parser.add_argument('save_dir', type=str, nargs='?', default=os.getcwd(), help='Directory to save the results')
#    parser.add_argument('data_file', type=str, help='Filename of the dataset')
#    parser.add_argument('save_file', type=str, help='Filename to save model as')

    # Parse the command-line arguments
#    args = parser.parse_args()
#^Never got this to work right

train_scvi_model(data_dir, data_file, scvi_model_file, scANVI_model_file, scvi_aData_file, scANVI_aData_file)
