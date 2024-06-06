import os
import numpy as np
import polars as pl

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
            print("Use --genes_correlation command to extract all interactions genes")
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
        not_found_genes = []
        for gene in self.genes:
            current_info_df = self.info_df.filter(pl.col("preferred_name") == gene)
            if current_info_df.is_empty():
                not_found_genes.append(gene)
                continue

            string_ids = current_info_df["string_protein_id"].to_list()
            final_datas = []

            for string_id in string_ids:
                current_links_df = self.links_df.filter(pl.col("protein1") == string_id)
                protein2_string_ids = current_links_df["protein2"].to_list()

                for p2_string_id in protein2_string_ids:
                    protein2_info = self.info_df.filter(
                        pl.col("string_protein_id") == p2_string_id
                    )
                    if protein2_info.is_empty():
                        continue  # Skip if protein2_info is empty

                    node1 = gene
                    node2 = protein2_info["preferred_name"][0]
                    link_info = current_links_df.filter(
                        pl.col("protein2") == p2_string_id
                    ).to_dicts()[0]

                    row = {
                        "node1": node1,
                        "node2": node2,
                        "node1_string_id": string_id,
                        "node2_string_id": p2_string_id,
                        "neighborhood_on_chromosome": self.normalize(
                            link_info["neighborhood"], self.max_values["neighborhood"]
                        ),
                        "gene_fusion": self.normalize(
                            link_info["fusion"], self.max_values["fusion"]
                        ),
                        "phylogenetic_cooccurrence": self.normalize(
                            link_info["cooccurence"], self.max_values["cooccurence"]
                        ),
                        "homology": self.normalize(
                            link_info["homology"], self.max_values["homology"]
                        ),
                        "coexpression": self.normalize(
                            link_info["coexpression"], self.max_values["coexpression"]
                        ),
                        "experimentally_determined_interaction": self.normalize(
                            link_info["experiments"], self.max_values["experiments"]
                        ),
                        "database_annotated": self.normalize(
                            link_info["database"], self.max_values["database"]
                        ),
                        "automated_textmining": self.normalize(
                            link_info["textmining"], self.max_values["textmining"]
                        ),
                        "combined_score": self.normalize(
                            link_info["combined_score"],
                            self.max_values["combined_score"],
                        ),
                    }
                    final_datas.append(row)

            if final_datas:
                final_df = pl.DataFrame(final_datas)
                final_df = final_df.filter(
                    pl.col("combined_score").cast(pl.Float64) >= float(self.cut_off)
                )
                final_df = final_df.sort("combined_score", descending=True)
                save_path = os.path.join(self.save_dir, f"{gene}.tsv")
                final_df.write_csv(file=save_path, separator="\t")
                print(f"Gene [{gene}] data saved to: [{save_path}]")

        if not_found_genes:
            with open(
                os.path.join(self.save_dir, "notfound_gene.txt"), "w", encoding="utf-8"
            ) as f:
                for gene in not_found_genes:
                    f.write(f"{gene}\n")

    def export_string_interactions(self):
        valid_interactions = []
        matched_genes = set()

        for gene in self.genes:
            current_info_df = self.info_df.filter(pl.col("preferred_name") == gene)
            if current_info_df.is_empty():
                continue

            string_ids = current_info_df["string_protein_id"].to_list()

            for string_id in string_ids:
                current_links_df = self.links_df.filter(pl.col("protein1") == string_id)
                protein2_string_ids = current_links_df["protein2"].to_list()

                for p2_string_id in protein2_string_ids:
                    protein2_info = self.info_df.filter(
                        pl.col("string_protein_id") == p2_string_id
                    )

                    if protein2_info.is_empty():
                        continue  # Skip if protein2_info is empty

                    node1 = gene
                    node2 = protein2_info["preferred_name"][0]

                    if node1 in self.genes and node2 in self.genes:
                        link_info = current_links_df.filter(
                            pl.col("protein2") == p2_string_id
                        ).to_dicts()[0]

                        row = {
                            "node1": node1,
                            "node2": node2,
                            "node1_string_id": string_id,
                            "node2_string_id": p2_string_id,
                            "neighborhood_on_chromosome": self.normalize(
                                link_info["neighborhood"],
                                self.max_values["neighborhood"],
                            ),
                            "gene_fusion": self.normalize(
                                link_info["fusion"], self.max_values["fusion"]
                            ),
                            "phylogenetic_cooccurrence": self.normalize(
                                link_info["cooccurence"], self.max_values["cooccurence"]
                            ),
                            "homology": self.normalize(
                                link_info["homology"], self.max_values["homology"]
                            ),
                            "coexpression": self.normalize(
                                link_info["coexpression"],
                                self.max_values["coexpression"],
                            ),
                            "experimentally_determined_interaction": self.normalize(
                                link_info["experiments"], self.max_values["experiments"]
                            ),
                            "database_annotated": self.normalize(
                                link_info["database"], self.max_values["database"]
                            ),
                            "automated_textmining": self.normalize(
                                link_info["textmining"], self.max_values["textmining"]
                            ),
                            "combined_score": self.normalize(
                                link_info["combined_score"],
                                self.max_values["combined_score"],
                            ),
                        }

                        valid_interactions.append(row)
                        matched_genes.update([node1, node2])

        if valid_interactions:
            final_df = pl.DataFrame(valid_interactions)
            final_df = final_df.filter(
                pl.col("combined_score").cast(pl.Float64) >= float(self.cut_off)
            )
            final_df = final_df.sort("combined_score", descending=True)

            # 保存原始文件
            save_path = os.path.join(self.save_dir, "string_interactions.tsv")
            final_df.write_csv(save_path, separator="\t")
            print(f"String interactions data saved to: [{save_path}]")

            # 创建简短版文件
            unique_pairs = set()
            short_interactions = []

            for row in final_df.to_dicts():
                node1 = row["node1"]
                node2 = row["node2"]
                if (node1, node2) not in unique_pairs and (
                    node2,
                    node1,
                ) not in unique_pairs:
                    unique_pairs.add((node1, node2))
                    short_interactions.append(row)

            short_df = pl.DataFrame(short_interactions)
            short_save_path = os.path.join(
                self.save_dir, "string_interactions_short.tsv"
            )
            short_df.write_csv(short_save_path, separator="\t")
            print(f"Short string interactions data saved to: [{short_save_path}]")

        not_matched_genes = set(self.genes) - matched_genes
        if not_matched_genes:
            with open(
                os.path.join(self.save_dir, "not_matched_genes.txt"),
                "w",
                encoding="utf-8",
            ) as f:
                for gene in not_matched_genes:
                    f.write(f"{gene}\n")


def main() -> int:
    # 获取当前脚本的文件路径
    script_path = __file__

    # 获取当前脚本所在的目录
    script_dir = os.path.dirname(script_path)

    Extractor(cli_file_path=os.path.join(script_dir, "..", "..", "cli_config.toml"))

    return 0
