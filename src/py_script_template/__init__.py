import os
import pandas as pd
import numpy as np

from py_script_template.cli import get_cli_argument
from py_script_template.utils import mkdir_if_not_exist, read_single_csv


class Extractor:
    def __init__(self, cli_file_path: str):
        self._load_args(cli_file_path)
        self._load_db_datas()
        self._load_input_file()
        self._calculate_max_values()

        # 对外接口功能
        if self.get_genes_correlation:
            print("使用")
            self.export_all_interactions()
            
        self.export_string_interactions()

    def _load_args(self, cli_file_path: str):
        self.args = get_cli_argument(cli_file_path)
        self.input_file = self.args["input_file"]
        self.db_dir = self.args["db_dir"]
        self.save_dir = self.args["save_dir"]
        self.cut_off = self.args["cut_off"]
        self.get_genes_correlation = self.args["get_genes_correlation"]
        mkdir_if_not_exist(self.save_dir)
        print("Args loaded.")

    def _load_db_datas(self):
        mkdir_if_not_exist(os.path.join("__pycache__", "cache"))

        print("Loading data...")
        self.info_df = read_single_csv(
            input_path=os.path.join(self.db_dir, "9606.protein.info.v12.0.txt"),
            cache_save_path=os.path.join(
                "__pycache__",
                "cache",
                "9606_protein_info_v12_0_parquet",
            ),
            sep="\t",
            dtype={
                "string_protein_id": "string",
                "preferred_name": "string",
                "protein_size": "int32",
                "annotation": "string",
            },
        )
        print(f"[{self.db_dir}/9606.protein.info.v12.0.txt] loaded")

        self.links_df = read_single_csv(
            input_path=os.path.join(self.db_dir, "9606.protein.links.full.v12.0.txt"),
            cache_save_path=os.path.join(
                "__pycache__",
                "cache",
                "9606_protein_links_full_v12_0.parquet",
            ),
            sep=r"\s+",
            dtype={
                "protein1": "string",
                "protein2": "string",
                "neighborhood": "int32",
                "neighborhood_transferred": "int32",
                "fusion": "int32",
                "cooccurence": "int32",
                "homology": "int32",
                "coexpression": "int32",
                "coexpression_transferred": "int32",
                "experiments": "int32",
                "experiments_transferred": "int32",
                "database": "int32",
                "database_transferred": "int32",
                "textmining": "int32",
                "textmining_transferred": "int32",
                "combined_score": "int32",
            },
        )
        print(f"[{self.db_dir}/9606.protein.links.full.v12.0.txt] loaded")

    def _load_input_file(self):
        with open(self.input_file, "r", encoding="utf-8") as r:
            self.genes = [line.strip("\n") for line in r.readlines()]
        print(
            f"Read {len(self.genes)} gene IDs from {self.input_file}, has {self.genes} genes"
        )

    def _calculate_max_values(self):
        self.max_values = {
            "neighborhood": self.links_df["neighborhood"].max(),
            "fusion": self.links_df["fusion"].max(),
            "cooccurence": self.links_df["cooccurence"].max(),
            "homology": self.links_df["homology"].max(),
            "coexpression": self.links_df["coexpression"].max(),
            "experiments": self.links_df["experiments"].max(),
            "database": self.links_df["database"].max(),
            "textmining": self.links_df["textmining"].max(),
            "combined_score": self.links_df["combined_score"].max(),
        }
        print(f"Max values calculated: {self.max_values}")

    def normalize(self, value, max_value, decimals=3):
        if max_value == 0:
            return 0
        normalized_value = value / max_value
        if normalized_value == 1.0:
            normalized_value = 0.900
        truncated_value = np.floor(normalized_value * 1000) / 1000
        final_value = f"%.{decimals}f" % truncated_value
        return 0 if final_value == "0.000" else final_value

    def export_all_interactions(self):
        fieldnames = [
            "node1",
            "node2",
            "node1_string_id",
            "node2_string_id",
            "neighborhood_on_chromosome",
            "gene_fusion",
            "phylogenetic_cooccurrence",
            "homology",
            "coexpression",
            "experimentally_determined_interaction",
            "database_annotated",
            "automated_textmining",
            "combined_score",
        ]

        not_found_genes = []

        for gene in self.genes:
            current_info_df = self.info_df[self.info_df.preferred_name == gene]
            if current_info_df.empty:
                not_found_genes.append(gene)
                continue

            string_ids = current_info_df["string_protein_id"].astype(str).values
            final_datas = []

            for string_id in string_ids:
                current_links_df = self.links_df[
                    self.links_df.protein1 == string_id
                ].copy()
                protein2_string_ids = current_links_df["protein2"].astype(str).values

                for p2_string_id in protein2_string_ids:
                    protein2_info = self.info_df[
                        self.info_df.string_protein_id == p2_string_id
                    ].copy()
                    if protein2_info.empty:
                        continue  # Skip if protein2_info is empty

                    node1 = gene
                    node2 = protein2_info["preferred_name"].astype(str).values[0]
                    link_info = current_links_df[
                        current_links_df.protein2 == p2_string_id
                    ].iloc[0]

                    row = {
                        "node1": node1,
                        "node2": node2,
                        "node1_string_id": string_id,
                        "node2_string_id": p2_string_id,
                        "neighborhood_on_chromosome": self.normalize(
                            link_info.neighborhood, self.max_values["neighborhood"]
                        ),
                        "gene_fusion": self.normalize(
                            link_info.fusion, self.max_values["fusion"]
                        ),
                        "phylogenetic_cooccurrence": self.normalize(
                            link_info.cooccurence, self.max_values["cooccurence"]
                        ),
                        "homology": self.normalize(
                            link_info.homology, self.max_values["homology"]
                        ),
                        "coexpression": self.normalize(
                            link_info.coexpression, self.max_values["coexpression"]
                        ),
                        "experimentally_determined_interaction": self.normalize(
                            link_info.experiments, self.max_values["experiments"]
                        ),
                        "database_annotated": self.normalize(
                            link_info.database, self.max_values["database"]
                        ),
                        "automated_textmining": self.normalize(
                            link_info.textmining, self.max_values["textmining"]
                        ),
                        "combined_score": self.normalize(
                            link_info.combined_score, self.max_values["combined_score"]
                        ),
                    }
                    final_datas.append(row)

            if final_datas:
                final_df = pd.DataFrame(final_datas)
                final_df = final_df.sort_values(by="combined_score", ascending=False)
                final_df = final_df[
                    final_df.combined_score.astype(float) >= float(self.cut_off)
                ]
                save_path = os.path.join(self.save_dir, f"{gene}.tsv")
                final_df.to_csv(
                    save_path,
                    index=False,
                    columns=fieldnames,
                    encoding="utf-8",
                    sep="\t",
                )
                print(f"Gene [{gene}] data saved to: [{save_path}]")

        if not_found_genes:
            with open(
                os.path.join(self.save_dir, "notfound_gene.txt"), "w", encoding="utf-8"
            ) as f:
                for gene in not_found_genes:
                    f.write(f"{gene}\n")
        fieldnames = [
            "node1",
            "node2",
            "node1_string_id",
            "node2_string_id",
            "neighborhood_on_chromosome",
            "gene_fusion",
            "phylogenetic_cooccurrence",
            "homology",
            "coexpression",
            "experimentally_determined_interaction",
            "database_annotated",
            "automated_textmining",
            "combined_score",
        ]

        not_found_genes = []

        for gene in self.genes:
            current_info_df = self.info_df[self.info_df.preferred_name == gene]
            if current_info_df.empty:
                not_found_genes.append(gene)
                continue

            string_ids = current_info_df["string_protein_id"].astype(str).values
            final_datas = []

            for string_id in string_ids:
                current_links_df = self.links_df[
                    self.links_df.protein1 == string_id
                ].copy()
                protein2_string_ids = current_links_df["protein2"].astype(str).values

                for p2_string_id in protein2_string_ids:
                    protein2_info = self.info_df[
                        self.info_df.string_protein_id == p2_string_id
                    ].copy()
                    if protein2_info.empty:
                        continue  # Skip if protein2_info is empty

                    node1 = gene
                    node2 = protein2_info["preferred_name"].astype(str).values[0]
                    link_info = current_links_df[
                        current_links_df.protein2 == p2_string_id
                    ].iloc[0]

                    row = {
                        "node1": node1,
                        "node2": node2,
                        "node1_string_id": string_id,
                        "node2_string_id": p2_string_id,
                        "neighborhood_on_chromosome": self.normalize(
                            link_info.neighborhood, self.max_values["neighborhood"]
                        ),
                        "gene_fusion": self.normalize(
                            link_info.fusion, self.max_values["fusion"]
                        ),
                        "phylogenetic_cooccurrence": self.normalize(
                            link_info.cooccurence, self.max_values["cooccurence"]
                        ),
                        "homology": self.normalize(
                            link_info.homology, self.max_values["homology"]
                        ),
                        "coexpression": self.normalize(
                            link_info.coexpression, self.max_values["coexpression"]
                        ),
                        "experimentally_determined_interaction": self.normalize(
                            link_info.experiments, self.max_values["experiments"]
                        ),
                        "database_annotated": self.normalize(
                            link_info.database, self.max_values["database"]
                        ),
                        "automated_textmining": self.normalize(
                            link_info.textmining, self.max_values["textmining"]
                        ),
                        "combined_score": self.normalize(
                            link_info.combined_score, self.max_values["combined_score"]
                        ),
                    }
                    final_datas.append(row)

            if final_datas:
                final_df = pd.DataFrame(final_datas)
                final_df = final_df[
                    final_df.combined_score.astype(float) >= float(self.cut_off)
                ]
                final_df = final_df.sort_values(by="combined_score", ascending=False)
                save_path = os.path.join(self.save_dir, f"{gene}.tsv")
                final_df.to_csv(
                    save_path,
                    index=False,
                    columns=fieldnames,
                    encoding="utf-8",
                    sep="\t",
                )
                print(f"Gene [{gene}] data saved to: [{save_path}]")

        if not_found_genes:
            with open(
                os.path.join(self.save_dir, "notfound_gene.txt"), "w", encoding="utf-8"
            ) as f:
                for gene in not_found_genes:
                    f.write(f"{gene}\n")

    def export_string_interactions(self):
        fieldnames = [
            "node1",
            "node2",
            "node1_string_id",
            "node2_string_id",
            "neighborhood_on_chromosome",
            "gene_fusion",
            "phylogenetic_cooccurrence",
            "homology",
            "coexpression",
            "experimentally_determined_interaction",
            "database_annotated",
            "automated_textmining",
            "combined_score",
        ]

        valid_interactions = []
        matched_genes = set()

        for gene in self.genes:
            current_info_df = self.info_df[self.info_df.preferred_name == gene]
            if current_info_df.empty:
                continue

            string_ids = current_info_df["string_protein_id"].astype(str).values

            for string_id in string_ids:
                current_links_df = self.links_df[
                    self.links_df.protein1 == string_id
                ].copy()
                protein2_string_ids = current_links_df["protein2"].astype(str).values

                for p2_string_id in protein2_string_ids:
                    protein2_info = self.info_df[
                        self.info_df.string_protein_id == p2_string_id
                    ].copy()
                    if protein2_info.empty:
                        continue

                    node1 = gene
                    node2 = protein2_info["preferred_name"].astype(str).values[0]

                    if node2 not in self.genes:
                        continue

                    matched_genes.add(node1)
                    matched_genes.add(node2)

                    link_info = current_links_df[
                        current_links_df.protein2 == p2_string_id
                    ].iloc[0]
                    row = {
                        "node1": node1,
                        "node2": node2,
                        "node1_string_id": string_id,
                        "node2_string_id": p2_string_id,
                        "neighborhood_on_chromosome": self.normalize(
                            link_info.neighborhood, self.max_values["neighborhood"]
                        ),
                        "gene_fusion": self.normalize(
                            link_info.fusion, self.max_values["fusion"]
                        ),
                        "phylogenetic_cooccurrence": self.normalize(
                            link_info.cooccurence, self.max_values["cooccurence"]
                        ),
                        "homology": self.normalize(
                            link_info.homology, self.max_values["homology"]
                        ),
                        "coexpression": self.normalize(
                            link_info.coexpression, self.max_values["coexpression"]
                        ),
                        "experimentally_determined_interaction": self.normalize(
                            link_info.experiments, self.max_values["experiments"]
                        ),
                        "database_annotated": self.normalize(
                            link_info.database, self.max_values["database"]
                        ),
                        "automated_textmining": self.normalize(
                            link_info.textmining, self.max_values["textmining"]
                        ),
                        "combined_score": self.normalize(
                            link_info.combined_score, self.max_values["combined_score"]
                        ),
                    }
                    valid_interactions.append(row)

        final_df = pd.DataFrame(valid_interactions)
        final_df = final_df[
            (final_df.node1.isin(self.genes) & (final_df.node2.isin(self.genes)))
            & (final_df.combined_score.astype(float) >= float(self.cut_off))
        ]
        final_df = final_df.sort_values(by="combined_score", ascending=False)
        final_df.to_csv(
            os.path.join(self.save_dir, "string_interactions.tsv"),
            index=False,
            columns=fieldnames,
            encoding="utf-8",
            sep="\t",
        )
        print("String interactions saved to string_interactions.tsv")

        # Determine the genes that were not matched
        notmatch_genes = set(self.genes) - matched_genes

        if notmatch_genes:
            with open(
                os.path.join(self.save_dir, "notmatch_gene.txt"), "w", encoding="utf-8"
            ) as f:
                for gene in notmatch_genes:
                    f.write(f"{gene}\n")
            print("Not matched genes saved to notmatch_gene.txt")


def main() -> int:
    Extractor(cli_file_path="./cli_config.toml")

    return 0


if __name__ == "__main__":
    main()
