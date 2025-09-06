import pandas as pd
import numpy as np
import glob
import os
import re
import streamlit as st
import matplotlib.pyplot as plt

def list_available_files(folder):
    files = glob.glob(os.path.join(folder, "*.assoc.linear"))
    info = []

    for f in files:
        base = os.path.basename(f).replace(".assoc.linear", "")
        parts = base.split("_")
        if parts[0] != "embedding":
            continue
        try:
            component = int(parts[1])
            group = "_".join(parts[2:])
            info.append({"file": f, "component": component, "group": group})
        except:
            continue

    df = pd.DataFrame(info)
    return df

def load_single_gwas(path):
    df = pd.read_csv(path, sep='\s+')
    df = df[df["P"] > 0].dropna(subset=["P", "CHR", "BP"])
    df["-log10(P)"] = -np.log10(df["P"])
    df["chrom_numeric"] = pd.to_numeric(df["CHR"], errors="coerce")
    df = df.dropna(subset=["chrom_numeric", "BP", "-log10(P)"])
    return df

def manhattan_and_qq(df):
    df = df.copy()
    df = df[df["P"] > 0].dropna(subset=["P", "CHR", "BP"])
    df["-log10(P)"] = -np.log10(df["P"])
    df["CHR"] = pd.to_numeric(df["CHR"], errors="coerce").astype(int)
    df = df.sort_values(["CHR", "BP"])

    # Manhattan: calcular posición acumulada
    chr_offsets = {}
    offset = 0
    ticks = []
    labels = []

    for chrom in sorted(df["CHR"].unique()):
        chr_df = df[df["CHR"] == chrom]
        chr_offsets[chrom] = offset
        max_bp = chr_df["BP"].max()
        midpoint = offset + max_bp / 2
        ticks.append(midpoint)
        labels.append(str(chrom))
        offset += max_bp

    df["genome_pos"] = df.apply(lambda row: row["BP"] + chr_offsets[row["CHR"]], axis=1)

    # QQ plot: valores esperados vs observados
    pvals = df["P"].sort_values()
    expected = -np.log10(np.linspace(1 / len(pvals), 1, len(pvals)))
    observed = -np.log10(pvals.values)

    # Layout en columnas
    col1, col2 = st.columns([2.5, 1])

    with col1:
        st.subheader("Manhattan Plot")
        fig1, ax1 = plt.subplots(figsize=(10, 4))
        ax1.scatter(df["genome_pos"], df["-log10(P)"], c=df["CHR"] % 2, cmap="coolwarm", s=3)
        ax1.axhline(-np.log10(5e-8), color='grey', linestyle='--', linewidth=1)
        ax1.set_xlabel("Chromosome")
        ax1.set_ylabel("-log10(P)")
        ax1.set_xticks(ticks)
        ax1.set_xticklabels(labels)
        st.pyplot(fig1)

    with col2:
        st.subheader("QQ Plot")
        fig2, ax2 = plt.subplots(figsize=(4, 4.5))
        ax2.scatter(expected, observed, s=3, alpha=0.6)
        ax2.plot([0, max(expected)], [0, max(expected)], color="red", linestyle="--")
        ax2.set_xlabel("Expected -log10(P)")
        ax2.set_ylabel("Observed -log10(P)")
        st.pyplot(fig2)

# === Streamlit app ===

st.set_page_config(layout="wide")
st.title("📊 GWAS Dashboard")

# Elegir dataset version
dataset_version = st.selectbox("Seleccionar dataset", [20, 30, 40, 50, 60], index=0)
FOLDER = f"gwas/gwas_outputs_{dataset_version}"
SUMMARY_FOLDER = os.path.join(FOLDER, "summaries")

# Cargar índice de archivos disponibles
available = list_available_files(FOLDER)

tab1, tab2 = st.tabs(["📈 Visualización por archivo", "🌍 Panorama general por embedding"])

with tab1:
    if available.empty:
        st.error("No GWAS files found.")
    else:
        group = st.selectbox("Grupo", sorted(available["group"].unique()))
        component = st.slider("Componente de embedding", 0, 239, 0)
    
        selected_row = available[(available["component"] == component) & (available["group"] == group)]
    
        if selected_row.empty:
            st.warning("No matching file.")
        else:
            path = selected_row.iloc[0]["file"]
            st.success(f"Archivo seleccionado: {os.path.basename(path)}")
    
            df = load_single_gwas(path)
    
            manhattan_and_qq(df)
    
            st.subheader("Top SNPs")
            st.dataframe(df.sort_values("P").head(20))

@st.cache_data
def collect_summaries(summary_folder):
    files = glob.glob(os.path.join(summary_folder, "embedding_*_chr*.csv"))
    summary_rows = []

    pattern = re.compile(r"embedding_(\d+)_white_([^_]+)_(chr\d+)_1Mb\.csv")
    from tqdm import tqdm

    for f in tqdm(files):
        match = pattern.search(os.path.basename(f))
        if not match:
            continue

        component = int(match.group(1))
        group = match.group(2)
        chr_label = match.group(3)

        try:
            df = pd.read_csv(f)
            num_significant = (df["P"] < 5e-8).sum()
            summary_rows.append({
                "component": component,
                "group": group,
                "chromosome": chr_label,
                "n_sig_blocks": num_significant
            })
        except Exception as e:
            print(f"❌ Error in {f}: {e}")

    summary_df = pd.DataFrame(summary_rows)

    if summary_df.empty:
        st.warning("No summary files found.")
    else:
        table = summary_df.groupby(["component", "group"]).agg({"n_sig_blocks": "sum"}).reset_index()
        table = table.sort_values("n_sig_blocks", ascending=False)
    
    return table

with tab2:
    st.header("Resumen por componente y grupo")
    
    print(SUMMARY_FOLDER)
    
    table = collect_summaries(SUMMARY_FOLDER)
    st.dataframe(table, use_container_width=True)
