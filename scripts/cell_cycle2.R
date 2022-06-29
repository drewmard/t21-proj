# module load R/4.1.2

library(Seurat)

cell_type_filename="HSCs_MPPs"

iter=0; df <- list()
for (sampletype in c("Liver","Femur")) {
  # disease_status="Healthy"
  for (disease_status in c("Healthy","DownSyndrome")) {
    iter=iter+1
    print(iter)
    f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cell_subset/",disease_status,'_',sampletype,'_',cell_type_filename,".rds")
    df[[iter]] <- readRDS(file = f.out)
  }
}

dfcombined = merge(df[[1]],c(df[[2]],df[[3]],df[[4]]))

# iter=iter+1
iter=5
cell_type="Cycling HSCs/MPPs";
cell_type_filename = gsub("/","_",cell_type)
disease_status="DownSyndrome"
sampletype="Liver"
f.out <- paste0("/oak/stanford/groups/smontgom/amarder/t21-proj/out/full/data/cell_subset/",disease_status,'_',sampletype,'_',cell_type_filename,".rds")
df[[iter]] <- readRDS(file = f.out)

# dfcombined = merge(dfcombined,df[[5]])
# dfcombined = merge(df[[1]],c(df[[2]],df[[3]],df[[4]],df[[5]]))

dftmp <- do.call(rbind,lapply(1:5,function(i) {df[[i]]@meta.data[,c("leiden_names","environment","organ","PHASE","phase")]}))
# dftmp <- rbind(dfcombined@meta.data[,c("leiden_names","environment","organ","PHASE")],df[[5]]@meta.data[,c("leiden_names","environment","organ","PHASE")])
# dftmp <- dfcombined@meta.data[,c("leiden_names","environment","organ","PHASE")]
dftmp$environment <- factor(dftmp$environment,levels=c("Healthy","Down Syndrome"))
dftmp$leiden_names <- factor(dftmp$leiden_names,levels=c("HSCs/MPPs","Cycling HSCs/MPPs"))
summary(lm(as.numeric(PHASE=="G2M+S")~environment+organ+leiden_names,data=dftmp))
summary(glm(as.numeric(PHASE=="G2M+S")~environment+organ+leiden_names,data=dftmp,family = binomial(link="logit")))

aggregate(as.numeric(PHASE=="G2M+S")~environment+organ+leiden_names,data=dftmp,mean)
aggregate(as.numeric(phase=="G2M")~environment+organ+leiden_names,data=dftmp,mean)
aggregate(as.numeric(phase=="S")~environment+organ+leiden_names,data=dftmp,mean)
aggregate(as.numeric(phase=="G1")~environment+organ+leiden_names,data=dftmp,mean)

# summary(lm(as.numeric(Phase!="G1")~environment+organ,data=dfcombined@meta.data))

tmp = dfcombined[,1:5000]
tmp <- CellCycleScoring(tmp, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
table(tmp@meta.data$Phase,dfcombined@meta.data$Phase[1:5000])

table(dftmp$PHASE,dftmp$leiden_names)

