rm(list=ls())
library(jsonlite)


carpeta_raiz <- "/Users/monicareyes/Documents/school/LCG/semestre5_lcg/bioinformatics/biologia-estructural-isomerasas/results/prediccion_estructural_alphafold/proteinas_alphafold_fijadas"
salida_csv   <- "~/Desktop/valores_alphafold_fijados.csv"

# Aminoácidos estándar -para excluir ligandos, iones, agua
aa <- c("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
        "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL")

# pLDDT - promedio de B_iso_or_equiv
calcular_plddt <- function(cif) {
  lin <- readLines(cif)
  cols <- sub("^_atom_site\\.", "", grep("^_atom_site\\.", lin, value = TRUE))
  atomos <- lin[grepl("^(ATOM|HETATM)\\s", lin)]
  campos <- strsplit(atomos, "\\s+")
  idx_grp <- match("group_PDB", cols)
  idx_res <- match("label_comp_id", cols)
  idx_b   <- match("B_iso_or_equiv", cols)
  grp <- vapply(campos, `[`, "", idx_grp)
  res <- vapply(campos, `[`, "", idx_res)
  b   <- as.numeric(vapply(campos, `[`, "", idx_b))
  es_prot <- grp == "ATOM" & res %in% aa
  mean(b[es_prot])
}

# Una carpeta por secuencia - esto por que asi las guarde
carpetas <- list.dirs(carpeta_raiz, recursive = FALSE, full.names = TRUE)

resultados <- lapply(carpetas, function(d) {
  id <- basename(d)
  
  cif  <- list.files(d, pattern = "_model_0\\.cif$", full.names = TRUE)
  summ <- list.files(d, pattern = "summary_confidences_0\\.json$", full.names = TRUE)
  conf <- list.files(d, pattern = "(full_data|confidences)_0\\.json$", full.names = TRUE)
  conf <- conf[!grepl("summary", conf)]
  
  if (length(cif) != 1 || length(summ) != 1 || length(conf) != 1) {
    warning("Archivos faltantes o ambiguos en: ", id)
    return(data.frame(secuencia = id, pLDDT = NA, pTM = NA, PAE = NA))
  }
  
  plddt <- calcular_plddt(cif)
  ptm   <- fromJSON(summ)$ptm
  pae   <- mean(as.matrix(fromJSON(conf)$pae), na.rm = TRUE)
  
  data.frame(secuencia = id, pLDDT = plddt, pTM = ptm, PAE = pae)
})

tabla <- do.call(rbind, resultados)
write.csv(tabla, salida_csv, row.names = FALSE)
print(tabla)

