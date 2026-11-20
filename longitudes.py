import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def longitudes(archivo):
    L, n = [], 0
    for l in open(archivo):
        l = l.strip()
        if l.startswith(">"):
            if n: L.append(n)
            n = 0
        else:
            n += len(l)
    L.append(n)
    return L

llm = longitudes("team_pi_llms.fasta")
nat = longitudes("naturales_1000.fasta")
print("LLM:", len(llm), "naturales:", len(nat))

fig, ax = plt.subplots(figsize=(7, 4.5))
bins = range(0, max(llm + nat) + 20, 10)
ax.hist(nat, bins=bins, alpha=0.6, color="#2a9d8f", label=f"Naturales (UniProt), n = {len(nat)}")
ax.hist(llm, bins=bins, alpha=0.6, color="#e76f51", label=f"Generadas por LLM, n = {len(llm)}")
ax.set_title("Distribución de longitudes: isomerasas generadas por LLM vs. naturales", fontsize=10)
ax.set_xlabel("Longitud (número de residuos)", fontsize=9)
ax.set_ylabel("Número de secuencias", fontsize=9)
ax.tick_params(labelsize=8)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig("longitudes.png", dpi=200)
