import scanpy as sc
import scarf

f="/oak/stanford/groups/smontgom/amarder/data/t21/Cellranger/cellranger310_count_HCCVYDSX2_and_HJJ3FDRXY_F15781B_GRCh38-3_0_0/seurat_obj." 
f_in = f + "h5ad"
f_out = f + "zarr"
adata=sc.read_h5ad(f_in)

reader = scarf.H5adReader(
    f_in, 
    cell_ids_key = 'index',               # Where Cell/barcode ids are saved under 'obs' slot
    feature_ids_key = 'name',            # Where gene ids are saved under 'var' slot
    feature_name_key = 'name'  # Where gene names are saved under 'var' slot
)  

# change value of `zarr_fn` to your choice of filename and path
writer = scarf.H5adToZarr(
    reader,
    zarr_fn=f_out
)
writer.dump()
