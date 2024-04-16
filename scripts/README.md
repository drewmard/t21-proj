# scRNA-seq scripts

This GitHub repo contains scripts necessary for running scRNA-seq analyses related to the following manuscript:

Marderstein, A.R. et al. (2023). Single-cell multi-omics map of human foetal blood in Down's Syndrome. biorxiv.

### Scripts

#### anndata --> seurat
anndata_to_seurat.R

#### seurat --> sce
seurat_to_sce.R

#### Cell barcode supp table
create_barcode_cluster_supp.R

#### QC of scRNA-seq ##########

#### Cell type clustering ##########

#### Cell type marker genes
cluster_markers_v2.py

#### Cell type marker genes - compile into supp table
create_cluster_markers_supp.R

#### UMAPs, dotplots
plots_v4.py

#### violin, ridge plots (cycling HSCs)
integration_combine.R

#### Cell type freq comparison - data create
create_cellComp_file.py

#### Cell type freq comparison - analysis
cellComp_analysis_v3.R

#### Sample age vs composition
pcw_vs_composition.R

#### Integration - metadata
integration_combine.R

#### Integration/ref query mapping (comparison of cell freq)
ref_vs_separate.R

#### pseudobulks for cycling HSCs
create_cycling_pseudobulks.R

#### pseudobulks from integrated labels
new_pseudobulks.R

#### Visium reference compare
240410_visium_ref_compare.R

#### DE robustness: CellBender - new object create
cellbender_out.py

#### DE robustness: CellBender - pseudobulks DE
230108_cellbender_pseudobulks.R

#### DE robustness: CellBender - analysis of DE results
230108_analysis_cellbender_pseudobulks.R

#### DE robustness: adding/dropping sex/age from DE analysis
sex_de1.R

#### DE robustness: performing DE analysis w/ different analysis choices
testing_de_mods_v2.R

#### DE robustness: analyzing analysis decisions
240226_pb_fold_changes.R

#### DE robustness: impact of chr 21 on log-fold change correlations
spiking_chr21.R

#### T21 v D21 DE analysis results - supp table
create_de_supp.R

#### Evaluating cross-dependent dysregulation of chr21 genes versus other genes (supp)
chr21_gene_analysis.R

#### correlate DE with Muskens et al (methylation)
240130_t21_methylation.R

#### liver v femur DE analysis
pseudobulk_de_v3_partB.062422update.R

#### downsampling T21 to D21 for liver v femur DE analysis
downsample_de_v3.t21_v_healthy.R

#### liver v femur DE: analysis related to it
analy_diff_exp_v4b_makingplots.062422update.R

#### GO enrichment of liver v femur DE
go_enrich_plots.R

#### GxE plots of liver v femur
gxe_plots.R
gxe_plots.makingdata.R

#### cell cycle analyze
cellCycle_analy.R
pull_cell_cycle_stats.R

#### Cycling HSC vs less cycling HSC - DE analysis
de_hsc_vs_cyclingHSC.R

#### Cycling HSC vs less cycling HSC - volcano plot
cycling_HSC_de_volcano_v3.R

#### Cycling HSC vs less cycling HSC - GO enrichments
cycling_HSC_go_enrich.R

#### CellPhoneDB - overall script
run_cpdb.sh

#### CellPhoneDB - input
create_cpdb_input.R

#### CellPhoneDB - python run
input_cellphone_db.py
input_cellphone_db.subset.sh

#### CellPhoneDB - output
cellphonedb_output_analy.R
reformat_cpdb_input.R

#### mitochondria-related data tables for barplots:
mitosox.R
mitotracker.R

#### mitochondria-related barplots:
mito.R

#### somatic enrichment:
291123_somatic_v3.R

#### Chr 21 DE versus genes on other chr:
logfc_chr21_vs_other.R

#### DE on other cell types
240227_diff_exp_other_cell_types.R

#### Various scrna-seq plots mostly related to DE analysis
scrna_seq_plot.R