# B Cell scRNAseq Reference Atlas From SLE Patients

This repo contains jupyter notebooks as well as standalone python scripts documenting the creation of a large single-cell RNAseq reference atlas of B cells taken from pubically available data. These cells have been derived from studies of human B cells from patients with systemic lupus erythematosus (SLE) in various ancestry populations, as well as B cells from matched healthy control individuals. 

**Note:** I'm still actively working to update this with all the data and models! I expect I'll get things added bit by bit over time

### Reference Specifications

All reference atlas contain only cells annotated broadly as B cells (i.e. CD19+) or plasmablasts from previously published studies. These references were designed for transfer learning purposes using scanpy/scVI/scArches tools. 

A variety of reference "flavors" have been created which can be used for various purposes depending on your research question of interest:
* SLE-derived B cells (for questions around differences in diseased individuals)
* Matched HC-derived B cells (for comparing SLE vs HC)
* Full Atlas (contains both HC and SLE patient cells)

### Studies Contributing to Atlas

#### [Perez et. al. Science 2022](https://pubmed.ncbi.nlm.nih.gov/35389781/) 

	Large SLE vs HC cohort
	1.2 mill PBMC, 162 SLE, 99 HC 
	Asian and White ancestry population
	151570 B cells


#### [Slight-Webb et. al. JCI insight 2023](https://pubmed.ncbi.nlm.nih.gov/37606045/)

	29 individuals each of White or Black ancestry
	scRNA + CITE-seq
	~20K B cells


### Steps to Reference Creation 

1. Verify and concatinate input datasets

- See the input_datasets and ref_prep directories

2. Seed label using reference studies' cell annotation to expand to entire reference

- I used the SlightWebb annotations as the seeds for the large references due to the use of CITEseq in their annotations and the granularity of B cell subtypes

