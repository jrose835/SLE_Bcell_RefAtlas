# B Cell scRNAseq Reference Atlas From SLE Patients

This repo contains jupyter notebooks as well as standalone python scripts documenting the creation of a large single-cell RNAseq reference atlas of B cells taken from pubically available data. These cells have been derived from studies of human B cells from patients with systemic lupus erythematosus (SLE) in various ancestry populations, as well as B cells from matched healthy control individuals. 

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



