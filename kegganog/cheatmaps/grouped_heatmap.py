"""Single-sample grouped three-panel heatmap generator module for KEGGaNOG.

This module categorizes metabolic pathways into distinct biological buckets,
injects structural NaN spacer blocks, and draws a partitioned matrix with left-aligned
functional group metadata tags.
"""

import warnings
from collections.abc import Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm

from .heatmaps_common import (
    GROUPED_PART1_GROUPS,
    GROUPED_PART2_GROUPS,
    GROUPED_PART3_GROUPS,
    insert_split_rows_between_groups,
    save_heatmap_png,
)

function_groups: dict[str, list[str]] = {
    "Carbon fixation": [
        "3-Hydroxypropionate Bicycle",
        "4-Hydroxybutyrate/3-hydroxypropionate",
        "CBB Cycle",
        "gluconeogenesis",
        "rTCA Cycle",
        "RuBisCo",
        "Wood-Ljungdahl",
    ],
    "Carbohydrate metabolism": [
        "Entner-Doudoroff Pathway",
        "glycolysis",
        "Sulfolipid biosynthesis",
        "TCA Cycle",
        "Glyoxylate shunt",
        "Mixed acid: Acetate",
        "Mixed acid: Ethanol, Acetate to Acetylaldehyde",
        "Mixed acid: Ethanol, Acetyl-CoA to Acetylaldehyde (reversible)",
        "Mixed acid: Ethanol, Acetylaldehyde to Ethanol",
        "Mixed acid: Formate",
        "Mixed acid: Formate to CO2 & H2",
        "Mixed acid: Lactate",
        "Mixed acid: PEP to Succinate via OAA, malate & fumarate",
        "alpha-amylase",
        "polyhydroxybutyrate synthesis",
        "starch/glycogen degradation",
        "starch/glycogen synthesis",
    ],
    "Carbon degradation": [
        "beta-glucosidase",
        "beta-N-acetylhexosaminidase",
        "chitinase",
        "D-galacturonate epimerase",
        "D-galacturonate isomerase",
        "diacetylchitobiose deacetylase",
        "glucoamylase",
        "pullulanase",
        "DMS dehydrogenase",
        "DMSP demethylation",
        "Naphthalene degradation to salicylate",
        "alcohol oxidase",
        "basic endochitinase B",
        "bifunctional chitinase/lysozyme",
        "cellulase",
        "dimethylamine/trimethylamine dehydrogenase",
        "oligogalacturonide lyase",
        "pectinesterase",
        "soluble methane monooxygenase",
    ],
    "Nitrogen metabolism": [
        "dissim nitrate reduction",
        "DNRA",
        "nitric oxide reduction",
        "nitrite oxidation",
        "nitrite reduction",
        "nitrogen fixation",
        "nitrous-oxide reduction",
        "ammonia oxidation (amo/pmmo)",
        "hydrazine dehydrogenase",
        "hydrazine synthase",
        "hydroxylamine oxidation",
    ],
    "Sulfur metabolism": [
        "alt thiosulfate oxidation tsdA",
        "dissimilatory sulfate < > APS",
        "dissimilatory sulfite < > APS",
        "dissimilatory sulfite < > sulfide",
        "DMSO reductase",
        "sulfide oxidation",
        "sulfite dehydrogenase",
        "sulfite dehydrogenase (quinone)",
        "sulfur dioxygenase",
        "thiosulfate oxidation",
        "thiosulfate/polysulfide reductase",
        "alt thiosulfate oxidation doxAD",
        "sulfhydrogenase",
        "sulfur assimilation",
        "sulfur disproportionation",
        "sulfur reductase sreABC",
    ],
    "Oxidative phosphorylation": [
        "Cytochrome bd complex",
        "Cytochrome c oxidase",
        "Cytochrome c oxidase, cbb3-type",
        "F-type ATPase",
        "Na-NADH-ubiquinone oxidoreductase",
        "NADH-quinone oxidoreductase",
        "Ubiquinol-cytochrome c reductase",
        "V-type ATPase",
        "Cytochrome aa3-600 menaquinol oxidase",
        "Cytochrome b6/f complex",
        "Cytochrome o ubiquinol oxidase",
        "NAD(P)H-quinone oxidoreductase",
    ],
    "Hydrogen redox": [
        "NAD-reducing hydrogenase",
        "NiFe hydrogenase Hyd-1",
        "Coenzyme B/Coenzyme M regeneration",
        "Coenzyme M reduction to methane",
        "NADP-reducing hydrogenase",
        "NiFe hydrogenase",
        "ferredoxin hydrogenase",
        "hydrogen:quinone oxidoreductase",
        "membrane-bound hydrogenase",
    ],
    "Amino acid metabolism": [
        "arginine",
        "asparagine",
        "glutamine",
        "histidine",
        "lysine",
        "serine",
        "threonine",
        "Serine pathway/formaldehyde assimilation",
        "alanine",
        "aspartate",
        "cysteine",
        "glutamate",
        "glycine",
        "isoleucine",
        "leucine",
        "methionine",
        "phenylalanine",
        "proline",
        "tryptophan",
        "tyrosine",
        "valine",
    ],
    "Vitamin biosynthesis": [
        "cobalamin biosynthesis",
        "riboflavin biosynthesis",
        "thiamin biosynthesis",
        "MEP-DOXP pathway",
        "Retinal biosynthesis",
        "Retinal from apo-carotenals",
        "carotenoids backbone biosynthesis",
        "end-product astaxanthin",
        "end-product myxoxanthophylls",
        "end-product nostoxanthin",
        "end-product zeaxanthin diglucoside",
        "mevalonate pathway",
    ],
    "Cell mobility": ["Chemotaxis", "Flagellum", "Adhesion"],
    "Biofilm formation": [
        "Biofilm PGA Synthesis protein",
        "Biofilm regulator BssS",
        "Colanic acid and Biofilm protein A",
        "Colanic acid and Biofilm transcriptional regulator",
        "Curli fimbriae biosynthesis",
    ],
    "Bacterial secretion systems": [
        "Sec-SRP",
        "Twin Arginine Targeting",
        "Type I Secretion",
        "Type II Secretion",
        "Type III Secretion",
        "Type IV Secretion",
        "Type Vabc Secretion",
        "Type VI Secretion",
    ],
    "Transporters": [
        "transporter: phosphate",
        "transporter: phosphonate",
        "transporter: thiamin",
        "transporter: urea",
        "C-P lyase cleavage PhnJ",
        "CP-lyase complex",
        "CP-lyase operon",
        "bidirectional polyphosphate",
        "transporter: vitamin B12",
    ],
    "Metal transporters": [
        "Cobalt transporter CbiMQ",
        "Cobalt transporter CorA",
        "Copper transporter CopA",
        "Fe-Mn transporter MntH",
        "Ferric iron ABC-type substrate-binding AfuA",
        "Ferrous iron transporter FeoB",
        "Cobalt transporter CbtA",
        "Nickel ABC-type substrate-binding NikA",
    ],
    "Arsenic reduction": ["Arsenic reduction"],
    "Methanogenesis": [
        "Methanogenesis via CO2",
        "Methanogenesis via acetate",
        "Methanogenesis via dimethylamine",
        "Methanogenesis via dimethylsulfide, methanethiol, methylpropanoate",
        "Methanogenesis via methanol",
        "Methanogenesis via methylamine",
        "Methanogenesis via trimethylamine",
    ],
    "Photosynthesis": [
        "Photosystem I",
        "Photosystem II",
        "anoxygenic type-I reaction center",
        "anoxygenic type-II reaction center",
    ],
    "Genetic competence": [
        "Competence factors",
        "Competence-related core components",
        "Competence-related related components",
    ],
    "Miscellaneous": [
        "Soluble methane monooxygenase",
        "Naphthalene degradation to salicylate",
        "alcohol oxidase",
        "DMS dehydrogenase",
        "ferredoxin hydrogenase",
    ],
}


def generate_grouped_heatmap(
    kegg_decoder_file: str,
    output_folder: str,
    dpi: int,
    color: str,
    sample_name: str | None,
    figsize: tuple[float, float] | None = None,
    annot: bool = True,
) -> tuple[plt.Figure, Sequence[plt.Axes]]:
    """Generate a functional grouped three-panel heatmap for a single sample profile.

    Parses tabular data streams, maps specific keys to categorical clusters, introduces
    structural separator empty blocks, and prints labels in rounded layout boxes.

    Args:
        kegg_decoder_file: System disk location pointing to raw text annotation matrices.
        output_folder: Target location identifying active processing directory loops.
        dpi: Target resolution scale bounding the output drawing canvas.
        color: Target string colormap descriptor passed to downstream styling engines.
        sample_name: Target sample identifier string mapping columns within data arrays.
        figsize: Geometric allocation limits (width, height) defining canvas borders.
        annot: Flag indicating whether cellular scalar elements are printed as text labels.

    Returns:
        tuple: Active Matplotlib Figure instance and the coordinate subplots Axes sequence.
    """
    # Parse and serialize disk spreadsheet rows into continuous strings
    with open(kegg_decoder_file, "r") as file:
        lines = file.readlines()

    # Validate sample elements and map numeric matrices under tracked progress
    with tqdm(total=4, desc="Preparing heatmap data") as pbar:
        header = lines[0].strip().split("\t")
        values = lines[1].strip().split("\t")

        target_column = sample_name if sample_name is not None else "Sample"

        data = {"Function": header[1:], target_column: [float(v) for v in values[1:]]}
        df = pd.DataFrame(data)
        pbar.update(1)

        # Lowercase mapping for robust group categorization matching
        function_groups_lower: dict[str, set[str]] = {
            group: {func.lower() for func in funcs}
            for group, funcs in function_groups.items()
        }

        df["Group"] = df["Function"].apply(
            lambda x: next(
                (
                    group
                    for group, funcs in function_groups_lower.items()
                    if x.lower() in funcs
                ),
                "Miscellaneous",
            )
        )

        df = df.sort_values(by=["Group", "Function"]).reset_index(drop=True)

        df["Function"] = pd.Categorical(
            df["Function"], categories=df["Function"], ordered=True
        )

        part1_groups = GROUPED_PART1_GROUPS
        part2_groups = GROUPED_PART2_GROUPS
        part3_groups = GROUPED_PART3_GROUPS

        # Reshape layout spreadsheet structure via split operations
        part1 = df[df["Group"].isin(part1_groups)].reset_index(drop=True)
        pbar.update(1)
        part2 = df[df["Group"].isin(part2_groups)].reset_index(drop=True)
        pbar.update(1)
        part3 = df[df["Group"].isin(part3_groups)].reset_index(drop=True)
        pbar.update(1)

    # Insert artificial NaN spacer rows between category groups
    with tqdm(total=6, desc="Adding split between groups") as pbar:
        part1 = insert_split_rows_between_groups(
            df[df["Group"].isin(part1_groups)], part1_groups
        ).reset_index(drop=True)
        pbar.update(1)
        part2 = insert_split_rows_between_groups(
            df[df["Group"].isin(part2_groups)], part2_groups
        ).reset_index(drop=True)
        pbar.update(1)
        part3 = insert_split_rows_between_groups(
            df[df["Group"].isin(part3_groups)], part3_groups
        ).reset_index(drop=True)
        pbar.update(1)

        part1["Function"] = pd.Categorical(
            part1["Function"], categories=part1["Function"], ordered=True
        )
        pbar.update(1)
        part2["Function"] = pd.Categorical(
            part2["Function"], categories=part2["Function"], ordered=True
        )
        pbar.update(1)
        part3["Function"] = pd.Categorical(
            part3["Function"], categories=part3["Function"], ordered=True
        )
        pbar.update(1)

    # Establish default structural canvas limits if missing
    target_figsize = figsize if figsize is not None else (28.0, 20.0)

    # Initialize structural subplots container canvas layers
    fig, axes_array = plt.subplots(1, 3, figsize=target_figsize)
    axes: Sequence[plt.Axes] = (
        axes_array.tolist() if hasattr(axes_array, "tolist") else axes_array
    )

    cbar_ax = fig.add_axes((0.92, 0.4, 0.02, 0.2))

    plt.subplots_adjust(left=0.15, right=0.85, wspace=0.4)

    # Functional label generator context helpers
    def add_group_labels(
        current_ax: plt.Axes, part_df: pd.DataFrame, group_labels: Sequence[str]
    ) -> None:
        for group in group_labels:
            group_indices = np.where(part_df["Group"] == group)[0]
            if len(group_indices) > 0:
                y_position = float(np.mean(group_indices) + 0.5)
                x_position = -0.075
                current_ax.text(
                    x_position,
                    y_position,
                    group,
                    fontsize=12,
                    ha="right",
                    va="center",
                    weight="bold",
                    bbox={
                        "boxstyle": "round,pad=0.3",
                        "edgecolor": "none",
                        "facecolor": "white",
                    },
                )

    def plot_heatmap(
        part_df: pd.DataFrame,
        group_labels: Sequence[str],
        ax: plt.Axes,
        cbar: bool,
        cbar_axis: plt.Axes | None = None,
        cbar_kws: dict[str, str] | None = None,
    ) -> None:
        value_columns = part_df.columns[1:-1]
        part_df[value_columns] = part_df[value_columns].fillna(0)

        pivot_table = part_df.pivot_table(
            values=value_columns,
            index="Function",
            aggfunc="mean",
            fill_value=0,
            observed=False,
        )

        mask = pivot_table.index.str.startswith("split_")

        sns.heatmap(
            pivot_table,
            cmap=color,
            annot=annot,
            linewidths=0.5,
            ax=ax,
            cbar=cbar,
            cbar_ax=cbar_axis,
            cbar_kws=cbar_kws,
            mask=np.tile(mask[:, None], (1, pivot_table.shape[1])),
        )
        # Derive tick positions directly from the index — get_yticklabels()
        # returns empty text before the figure is rendered, making the old
        # set_visible(False) loop a silent no-op.
        non_split_pos = [
            i + 0.5
            for i, name in enumerate(pivot_table.index)
            if not str(name).startswith("split_")
        ]
        non_split_labels = [
            name for name in pivot_table.index if not str(name).startswith("split_")
        ]
        ax.set_yticks(non_split_pos)
        ax.set_yticklabels(non_split_labels, rotation=0)
        ax.tick_params(axis="y", which="both", left=False)
        add_group_labels(ax, part_df, group_labels)

    # Map underlying statistical matrices chunk-by-chunk onto partitioned views
    with tqdm(total=3, desc="Creating heatmap parts") as pbar:
        plot_heatmap(part1, part1_groups, axes[0], cbar=False)
        axes[0].set_title("Part 1")
        pbar.update(1)

        plot_heatmap(part2, part2_groups, axes[1], cbar=False)
        axes[1].set_title("Part 2")
        pbar.update(1)

        plot_heatmap(
            part3,
            part3_groups,
            axes[2],
            cbar=True,
            cbar_axis=cbar_ax,
            cbar_kws={"label": "Pathway completeness"},
        )
        axes[2].set_title("Part 3")
        pbar.update(1)

        axes[0].set_ylabel("")
        axes[1].set_ylabel("")
        axes[2].set_ylabel("")

        for ax in axes:
            ax.yaxis.tick_right()
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0, va="center", ha="left")

        # Capture transient layout warnings and write final PNG assets to system disk
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", category=UserWarning, message=".*tight_layout.*"
            )
            plt.tight_layout(rect=(0.0, 0.0, 0.9, 1.0))

    save_heatmap_png(output_folder, dpi)

    # Return active layout context wrappers for external saving pipelines
    return fig, axes
