import csv
import json

with open("binder223_af3_libres.csv", encoding="utf-8") as f:
    filas = list(csv.DictReader(f))

print(len(filas))
print(filas[0]["Identificador"])
print(filas[0]["Secuencia"])

trabajos = []
for fila in filas:
    trabajo = {
        "name": fila["Identificador"].replace(".", "_"),
        "modelSeeds": [],
        "sequences": [{"proteinChain": {"sequence": fila["Secuencia"], "count": 1}}],
        "dialect": "alphafoldserver",
        "version": 1,
    }
    trabajos.append(trabajo)

print(len(trabajos))
print(trabajos[0])

with open("af3_input_libres.json", "w") as f:
    json.dump(trabajos, f, indent=2)

print("Listo: af3_input_libres.json")