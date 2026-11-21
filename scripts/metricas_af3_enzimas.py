import csv
import json
from pathlib import Path
from statistics import mean

BASE = Path.home() / "tarea_biocestructural"
AF3 = Path.home() / "Downloads/folds_2026_11_20_07_34"
ENTRADA = BASE / "entregables/isomerasas_10_secuencias.csv"
SALIDA = BASE / "entregables/isomerasas_10_metricas_af3.csv"

with ENTRADA.open(newline="") as f:
    lector = csv.DictReader(f)
    columnas = lector.fieldnames
    filas = list(lector)

if len(filas) != 10:
    raise ValueError("Se esperaban exactamente 10 secuencias")

resultados = []

for fila in filas:
    identificador = fila["Identificador"]
    carpetas = [
        p.parent for p in AF3.rglob(
            f"fold_{identificador}_model_0.cif"
        )
    ]
    if len(carpetas) != 1:
        raise ValueError(
            f"{identificador}: se encontraron {len(carpetas)} estructuras"
        )
    carpeta = carpetas[0]

    summary_path = carpeta / f"fold_{identificador}_summary_confidences_0.json"
    full_path = carpeta / f"fold_{identificador}_full_data_0.json"
    cif_path = carpeta / f"fold_{identificador}_model_0.cif"

    for archivo in [summary_path, full_path, cif_path]:
        if not archivo.is_file():
            raise FileNotFoundError(archivo)

    with summary_path.open() as f:
        summary = json.load(f)

    with full_path.open() as f:
        full = json.load(f)

    # Las predicciones se solicitaron como una sola cadena proteica.
    chain_ids = full["atom_chain_ids"]
    if len(set(chain_ids)) != 1:
        raise ValueError(
            f"{identificador}: se encontraron varias cadenas atómicas"
        )

    plddt = full["atom_plddts"]
    pae = full["pae"]

    if len(plddt) != len(chain_ids):
        raise ValueError(f"{identificador}: longitud atómica inconsistente")

    if len(pae) != len(fila["Secuencia"]):
        raise ValueError(f"{identificador}: dimensiones PAE inconsistentes")

    if any(len(r) != len(pae) for r in pae):
        raise ValueError(f"{identificador}: matriz PAE no cuadrada")

    nueva = dict(fila)
    nueva["pLDDT_promedio"] = round(mean(plddt), 3)
    nueva["pTM"] = round(float(summary["ptm"]), 4)
    nueva["PAE_promedio_A"] = round(
        mean(valor for r in pae for valor in r), 3
    )

    resultados.append(nueva)

columnas_finales = columnas + [
    "pLDDT_promedio", "pTM", "PAE_promedio_A"
]

with SALIDA.open("w", newline="") as f:
    escritor = csv.DictWriter(f, fieldnames=columnas_finales)
    escritor.writeheader()
    escritor.writerows(resultados)

print("\nRESULTADOS ALPHAFOLD3 — ISOMERASAS")
print("-" * 75)

for r in resultados:
    print(
        f'{r["Identificador"]:22s} '
        f'pLDDT={r["pLDDT_promedio"]:7.3f}  '
        f'pTM={r["pTM"]:.4f}  '
        f'PAE={r["PAE_promedio_A"]:.3f}'
    )

print("\nCSV generado:", SALIDA)
