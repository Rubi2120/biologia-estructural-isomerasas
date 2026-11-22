import csv
import glob
import os
import json
import csv
import re

CARPETA = "folds_2026_11_21_22_41"   # cambia esto si tu carpeta se llama distinto

with open("binder223_af3_libres.csv", encoding="utf-8") as f:
    filas = list(csv.DictReader(f))

def plddt_promedio(ruta_cif):
    columnas = []
    valores = []
    with open(ruta_cif) as f:
        for linea in f:
            linea = linea.strip()
            if linea.startswith("_atom_site."):
                columnas.append(linea.split(".")[1])
            elif linea.startswith("ATOM"):
                partes = linea.split()
                i = columnas.index("B_iso_or_equiv")
                valores.append(float(partes[i]))
    return sum(valores) / len(valores), len(valores)

def leer_ptm(ruta_summary):
    with open(ruta_summary) as f:
        return json.load(f)["ptm"]


def pae_promedio(ruta_full):
    with open(ruta_full) as f:
        pae = json.load(f)["pae"]
    total = sum(sum(fila) for fila in pae)
    n = sum(len(fila) for fila in pae)
    return total / n, len(pae)

for fila in filas:
    nombre = fila["Identificador"].replace(".", "_").lower()
    cif = glob.glob(f"{CARPETA}/**/*{nombre}_model_0.cif", recursive=True)[0]
    summary = glob.glob(f"{CARPETA}/**/*{nombre}_summary_confidences_0.json", recursive=True)[0]
    full = glob.glob(f"{CARPETA}/**/*{nombre}_full_data_0.json", recursive=True)[0]

    plddt, n_atomos = plddt_promedio(cif)
    ptm = leer_ptm(summary)
    pae, tam = pae_promedio(full)

    fila["pLDDT_promedio"] = round(plddt, 2)
    fila["pTM"] = ptm
    fila["PAE_promedio"] = round(pae, 2)

columnas_salida = list(filas[0].keys())
with open("binder223_af3_libres_metricas.csv", "w", newline="", encoding="utf-8") as f:
    escritor = csv.DictWriter(f, fieldnames=columnas_salida)
    escritor.writeheader()
    escritor.writerows(filas)

print("Guardado:", len(filas), "filas")

RUTA_PDB = os.path.abspath("binder_patulin_223.pdb")
RUTA_LOG = os.path.abspath("matchmaker_log.html")

lineas = [f'open "{RUTA_PDB}"']
for fila in filas:
    nombre = fila["Identificador"].replace(".", "_").lower()
    cif = os.path.abspath(glob.glob(f"{CARPETA}/**/*{nombre}_model_0.cif", recursive=True)[0])
    lineas.append(f"log text === {fila['Identificador']} ===")
    lineas.append(f'open "{cif}"')
    lineas.append("matchmaker #2/A to #1/A")
    lineas.append("close #2")
lineas.append(f'log save "{RUTA_LOG}"')

with open("matchmaker_libres.cxc", "w") as f:
    f.write("\n".join(lineas) + "\n")

print("Escrito matchmaker_libres.cxc con", len(filas), "modelos")

### Leer los resultados de ChimeraX
with open("matchmaker_log.html", encoding="utf-8", errors="ignore") as f:
    texto = re.sub(r"<[^>]+>", "", f.read())

partes = re.split(r"===\s*(.+?)\s*===", texto)
rmsd = {}
for ident, bloque in zip(partes[1::2], partes[2::2]):
    m = re.search(r"RMSD between (\d+) pruned atom pairs is ([\d.]+) angstroms; "
                  r"\(across all (\d+) pairs: ([\d.]+)\)", bloque)
    if m:
        rmsd[ident.strip()] = (float(m.group(2)), float(m.group(4)),
                               int(m.group(1)), int(m.group(3)))

print(len(rmsd), "modelos con RMSD")

for fila in filas:
    r = rmsd[fila["Identificador"]]
    fila["RMSD_pruned"] = r[0]
    fila["RMSD_all"] = r[1]
    fila["pares_pruned"] = r[2]
    fila["pares_all"] = r[3]

columnas_final = list(filas[0].keys())
with open("binder223_af3_libres_final.csv", "w", newline="", encoding="utf-8") as f:
    escritor = csv.DictWriter(f, fieldnames=columnas_final)
    escritor.writeheader()
    escritor.writerows(filas)

print("Guardado final:", len(filas), "filas,", len(columnas_final), "columnas")