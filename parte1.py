import csv

registros = []
with open("team_pi_llms.fasta") as f:
    id_, ppl, seq = None, None, []
    for linea in f:
        linea = linea.strip()
        if linea.startswith(">"):
            if id_ is not None:
                registros.append((id_, "".join(seq), ppl))
            id_ = linea[1:].split()[0]
            ppl = float(linea.split("PPL:")[1])
            seq = []
        elif linea:
            seq.append(linea)
    registros.append((id_, "".join(seq), ppl))

print("Secuencias leídas:", len(registros), "| IDs únicos:", len({r[0] for r in registros}))

registros.sort(key=lambda r: r[2])
sel = registros[:5] + registros[-5:]

def guardar(nombre, filas):
    with open(nombre, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["Identificador", "Secuencia", "Perplejidad"])
        w.writerows(filas)

guardar("isomerasas_10_secuencias.csv", sel)
guardar("isomerasas_1000_todas.csv", registros)

for r in sel:
    print(r[0], r[2])
