import math

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

import py3Dmol
from collections.abc import Sequence
from io import StringIO

from Bio.PDB import PDBParser, Superimposer, PDBIO

def plot_pae_from_df(
    res_df: pd.DataFrame,
    rows=None,
    vmin: float = 0.0,
    vmax: float = 30.0,
    cmap: str = "coolwarm",   # blue = niskie PAE, red = wysokie
    figsize_per_plot: float = 4.5,
    show_colorbar: bool = True,
):
    """
    Rysuj macierze PAE dla wybranych wierszy res_df w jednym rzędzie (1 x N)
    i zwróć listę macierzy PAE (np.ndarray).

    Oczekiwane kolumny w res_df:
      - 'pae' (np.ndarray),
      - 'chain_lengths',
      - 'cluster_id' (opcjonalnie, do podpisu),
      - 'job_name', 'ptm', 'iptm', 'ranking_confidence'.
    """

    # --- wybór wierszy -------------------------------------
    if rows is None:
        subset = res_df
    elif isinstance(rows, pd.DataFrame):
        subset = rows
    else:
        subset = res_df.loc[rows]

    if subset.empty:
        raise ValueError("Brak danych do rysowania PAE.")

    pae_matrices = []
    n_plots = len(subset)

    fig, axes = plt.subplots(
        1,
        n_plots,
        figsize=(figsize_per_plot * n_plots, figsize_per_plot),
        squeeze=False,
    )
    axes = axes[0]

    common_im = None

    # --- ewentualny suptitle z job_name --------------------
    unique_jobs = subset["job_name"].dropna().unique()
    if len(unique_jobs) == 1:
        fig.suptitle(str(unique_jobs[0]), fontsize=10)

    # --- rysowanie -----------------------------------------
    for i, row in enumerate(subset.itertuples()):
        pae = getattr(row, "pae", None)
        if pae is None:
            continue

        pae_matrices.append(pae)
        ax = axes[i]

        im = ax.imshow(
            pae,
            origin="lower",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
        common_im = im

        # granice łańcuchów
        chain_lengths = getattr(row, "chain_lengths", None)
        if chain_lengths is not None:
            boundaries = np.cumsum([0] + list(chain_lengths))
            L = pae.shape[0]
            for b in boundaries:
                if 0 < b < L:
                    ax.axhline(b - 0.5, color="lime", linewidth=5)
                    ax.axvline(b - 0.5, color="lime", linewidth=5)

        ax.set_xticks([])
        ax.set_yticks([])

        # --- krótki tytuł: numer klastra + podstawowe score ---
        cluster_id = getattr(row, "cluster_id", None)
        if cluster_id is not None and not pd.isna(cluster_id):
            title = f"clu {int(cluster_id)}"
        else:
            # awaryjnie, gdy brak cluster_id
            model_id = getattr(row, "model_id", "")
            title = model_id if model_id else ""

        extras = []
        for attr, lbl in [
            ("ptm", "pTM"),
            ("iptm", "ipTM"),
            ("ranking_confidence", "rank"),
        ]:
            val = getattr(row, attr, None)
            if val is not None and not pd.isna(val):
                extras.append(f"{lbl}={val:.3f}")

        if extras:
            title += "\n" + ", ".join(extras)

        ax.set_title(title, fontsize=9)

    # --- wspólny colorbar z boku, bez nachodzenia -----------
    if show_colorbar and common_im is not None:
        # zostawiamy miejsce na bar po prawej
        fig.subplots_adjust(right=0.9)
        cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])  # [left, bottom, width, height]
        fig.colorbar(common_im, cax=cbar_ax, label="PAE [Å]")

    plt.tight_layout(rect=[0, 0, 0.9, 1])  # nie używać prawego 10%, tam jest bar
    plt.show()

    return pae_matrices


def _plot_score_violin_grid(
    data: pd.DataFrame,
    group_col: str,
    group_order,
    scores,
    palette,
    n_cols: int = 3,
    figsize_per_panel: float = 4.0,
    x_tick_rotation: int = 45,
    x_label: str = "",
):
    """
    Wspólna funkcja rysująca siatkę violin-plotów dla zadanych metryk.

    Parametry
    ---------
    data : pd.DataFrame
        Dane wejściowe (już przefiltrowane).
    group_col : str
        Nazwa kolumny z grupami na osi X (np. 'job_name' lub 'cluster_label').
    group_order : list
        Kolejność grup na osi X.
    scores : list[str]
        Lista nazw kolumn z metrykami do rysowania.
    palette : dict
        Mapa {group_value: color}.
    n_cols : int
        Liczba kolumn paneli.
    figsize_per_panel : float
        Rozmiar jednego panelu w calach.
    x_tick_rotation : int
        Obrót etykiet na osi X.
    x_label : str
        Opis osi X (wspólny dla wszystkich paneli).

    Zwraca
    -------
    (fig, axes)
    """
    if data.empty:
        raise ValueError("Brak danych do rysowania (data jest puste).")

    n_scores = len(scores)
    n_cols = min(n_cols, n_scores)
    n_rows = math.ceil(n_scores / n_cols)

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(figsize_per_panel * n_cols,
                 figsize_per_panel * n_rows),
        squeeze=False
    )

    last_idx = -1

    for idx, score in enumerate(scores):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]

        sns.violinplot(
            data=data,
            x=group_col,
            y=score,
            order=group_order,
            palette=palette,
            cut=0,
            ax=ax,
        )

        ax.set_title(score)
        ax.set_xlabel(x_label)
        ax.set_ylabel(score)
        ax.tick_params(axis="x", rotation=x_tick_rotation)

        last_idx = idx

    # Ukryj nieużyte osie
    for j in range(last_idx + 1, n_rows * n_cols):
        r = j // n_cols
        c = j % n_cols
        axes[r, c].set_visible(False)

    plt.tight_layout()

    return fig, axes


def plot_job_score_violins(
    res_df: pd.DataFrame,
    job_order,
    scores=None,
    n_cols: int = 3,
    figsize_per_panel: float = 4.0,
):
    """
    Rysuj siatkę violin-plotów metryk AF2 dla wybranych jobów.

    Parametry
    ----------
    res_df : pd.DataFrame
        DataFrame z AF2ScoreParser.obtain_data()
    job_order : list[(str, str)]
        Lista tupli (job_name, color), np.:
            [
                ("BN_BN_domain", "tab:blue"),
                ("BC_BC_domain", "tab:blue"),
                ("BN_BN_exp", "tab:orange"),
                ("BC_BC_exp", "tab:orange"),
            ]
    scores : list[str], optional
        Lista kolumn do plotowania. Domyślnie:
            ['iptm','ptm','ranking_confidence','plddt_mean',
             'inter_pae_mean','intra_pae_mean']
    n_cols : int
        Liczba kolumn paneli.
    figsize_per_panel : float
        Rozmiar jednego panelu (cale).

    Zwraca
    -------
    (fig, axes)
    """

    if scores is None:
        scores = [
            "iptm",
            "ptm",
            "ranking_confidence",
            "plddt_mean",
            "inter_pae_mean",
            "intra_pae_mean",
        ]

    job_names = [j for j, _ in job_order]
    palette = {j: c for j, c in job_order}

    plot_df = res_df[res_df["job_name"].isin(job_names)]

    fig, axes = _plot_score_violin_grid(
        data=plot_df,
        group_col="job_name",
        group_order=job_names,
        scores=scores,
        palette=palette,
        n_cols=n_cols,
        figsize_per_panel=figsize_per_panel,
        x_tick_rotation=45,
        x_label="job_name",
    )

    return fig, axes


def plot_cluster_score_violins(
    res_df: pd.DataFrame,
    job_name: str,
    scores=None,
    n_cols: int = 3,
    figsize_per_panel: float = 4.0,
    sort_clusters_by: str = "size",  # 'size' lub 'id'
):
    """
    Dla zadanego job_name rysuje violin-ploty metryk AF2 per klaster Foldseeka.

    Parametry
    ----------
    res_df : pd.DataFrame
        DataFrame z AF2ScoreParser.obtain_data(), zawierający m.in. kolumny:
        - 'job_name'
        - 'cluster_id'
        - 'iptm', 'ptm', 'ranking_confidence', 'plddt_mean',
          'inter_pae_mean', 'intra_pae_mean' (domyślnie)
    job_name : str
        Nazwa joba, np. 'full_full_exp'.
    scores : list[str], optional
        Lista metryk do porównania między klastrami.
        Domyślnie jak wyżej.
    n_cols : int
        Liczba kolumn paneli.
    figsize_per_panel : float
        Rozmiar jednego panelu (cale).
    sort_clusters_by : {'size', 'id'}
        Jak sortować klastry na osi X:
        - 'size' -> malejąco po liczbie modeli w klastrze,
        - 'id'   -> rosnąco po cluster_id.

    Zwraca
    -------
    (fig, axes)
    """

    if scores is None:
        scores = [
            "iptm",
            "ptm",
            "ranking_confidence",
            "plddt_mean",
            "inter_pae_mean",
            "intra_pae_mean",
        ]

    # wybór jednego joba
    sel = res_df[res_df["job_name"] == job_name].copy()

    # ignorujemy modele bez przydzielonego klastra
    sel = sel.dropna(subset=["cluster_id"])
    if sel.empty:
        raise ValueError(f"Brak danych z przypisanym cluster_id dla job_name={job_name!r}")

    # etykieta klastra jako string (ładniejsze na osi X)
    sel["cluster_label"] = sel["cluster_id"].astype(int).astype(str)

    # kolejność klastrów
    if sort_clusters_by == "size":
        counts = sel["cluster_label"].value_counts()
        cluster_order = counts.sort_values(ascending=False).index.tolist()
    elif sort_clusters_by == "id":
        cluster_order = sorted(sel["cluster_label"].unique(), key=lambda x: int(x))
    else:
        raise ValueError("sort_clusters_by must be 'size' or 'id'.")

    # paleta kolorów – automatycznie, inny kolor na każdy klaster
    n_clusters = len(cluster_order)
    color_list = sns.color_palette("tab20", n_colors=n_clusters)
    palette = {lab: color_list[i] for i, lab in enumerate(cluster_order)}

    fig, axes = _plot_score_violin_grid(
        data=sel,
        group_col="cluster_label",
        group_order=cluster_order,
        scores=scores,
        palette=palette,
        n_cols=n_cols,
        figsize_per_panel=figsize_per_panel,
        x_tick_rotation=45,
        x_label=f"cluster_id ({job_name})",
    )

    return fig, axes
    

def _cas(structure):
    """Zwróć listę atomów CA ze struktury Bio.PDB."""
    return [
        res["CA"]
        for m in structure
        for ch in m
        for res in ch
        if "CA" in res
    ]


def _apply_chainbow_bfactor(structure):
    """
    Ustaw B-factor (tempfactor) jako indeks wzdłuż łańcucha (0–100),
    osobno dla każdego chaina. Dzięki temu można kolorować gradientem po B.
    """
    for model in structure:
        for chain in model:
            # bierzemy tylko reszty z CA, żeby indeks był po „szkielecie”
            residues = [res for res in chain if "CA" in res]
            n = len(residues)
            if n == 0:
                continue

            if n == 1:
                values = [50.0]
            else:
                values = [100.0 * i / (n - 1) for i in range(n)]

            for res, val in zip(residues, values):
                for atom in res:
                    atom.bfactor = val


import os
from collections.abc import Sequence
import py3Dmol
from io import StringIO
from Bio.PDB import PDBParser, Superimposer, PDBIO


def show_pdb_grid(
    pdb_paths,
    style="cartoon",
    color_mode="default",  # "default", "plddt", "chainbow", "chain"
    width_per_view=400,
    height=400,
    desc=None,             # <<< NOWE: lista opisów
):
    # --- przygotowanie listy ścieżek ---
    if isinstance(pdb_paths, str):
        pdb_paths = [pdb_paths]
    elif not isinstance(pdb_paths, Sequence):
        raise TypeError("pdb_paths must be a string or a sequence of strings")

    n = len(pdb_paths)
    if n == 0:
        raise ValueError("pdb_paths is empty")

    # jeśli podano desc, sprawdź długość
    if desc is not None and len(desc) != n:
        raise ValueError("desc must have the same length as pdb_paths")

    parser = PDBParser(QUIET=True)

    # === referencja ===
    ref_struct = parser.get_structure("ref", pdb_paths[0])
    ref_ca = _cas(ref_struct)

    pdb_strings = []
    io = PDBIO()

    # pierwszy – bez zmian lub z chainbowem
    if color_mode == "chainbow":
        _apply_chainbow_bfactor(ref_struct)
        io.set_structure(ref_struct)
        buf = StringIO()
        io.save(buf)
        pdb_strings.append(buf.getvalue())
    else:
        with open(pdb_paths[0]) as f:
            pdb_strings.append(f.read())

    # reszta – dopasowanie CA→CA
    sup = Superimposer()

    for path in pdb_paths[1:]:
        mov_struct = parser.get_structure("mov", path)
        mov_ca = _cas(mov_struct)

        # zakładamy, że listy CA są porównywalne
        sup.set_atoms(ref_ca, mov_ca)
        sup.apply(mov_struct.get_atoms())

        if color_mode == "chainbow":
            _apply_chainbow_bfactor(mov_struct)

        io.set_structure(mov_struct)
        buf = StringIO()
        io.save(buf)

        pdb_strings.append(buf.getvalue())

    # === rysowanie grida ===

    style_map = {
        "cartoon": "cartoon",
        "sticks": "stick",
        "spheres": "sphere",
    }

    if style not in style_map:
        raise ValueError(f"Unknown style: {style}")

    style_key = style_map[style]

    view = py3Dmol.view(
        viewergrid=(1, n),
        width=width_per_view * n,
        height=height,
    )

    for i, pdb_str in enumerate(pdb_strings):
        # dodaj model do odpowiedniej komórki grida
        view.addModel(pdb_str, "pdb", viewer=(0, i))

        # ustaw styl + kolorowanie
        if color_mode == "plddt":
            style_dict = {
                style_key: {
                    "colorscheme": {
                        "prop": "b",
                        "gradient": "roygb",
                        "min": 0,
                        "max": 100,
                    }
                }
            }
        elif color_mode == "chainbow":
            style_dict = {
                style_key: {
                    "colorscheme": {
                        "prop": "b",
                        "gradient": "roygb",
                        "min": 0,
                        "max": 100,
                    }
                }
            }
        elif color_mode == "chain":
            style_dict = {
                style_key: {
                    "colorscheme": "chain",
                }
            }
        elif color_mode == "default":
            style_dict = {
                style_key: {
                    "color": "spectrum",
                }
            }
        else:
            raise ValueError(
                f"Unknown color_mode: {color_mode}. "
                "Use 'default', 'plddt', 'chainbow', or 'chain'."
            )

        view.setStyle(style_dict, viewer=(0, i))
        view.zoomTo(viewer=(0, i))

        # ---------- LABEL Z IDENTYFIKATOREM -----------------
        # tekst: z desc albo domyślnie nazwa pliku
        label_text = (
            str(desc[i]) if desc is not None
            else os.path.basename(pdb_paths[i])
        )

        view.addLabel(
            label_text,
            {
                "fontColor": "black",
                "backgroundColor": "white",
                "backgroundOpacity": 0.7,
                "fontSize": 14,
                "inFront": True,
                "showBackground": True,
                # pozycja: w okolicach pierwszego atomu (0,0,0 zwykle też działa),
                # ważne jest 'inFront', wtedy label "wisi" nad sceną
                "position": {"x": 0, "y": 0, "z": 0},
            },
            viewer=(0, i),
        )

    return view.show()


