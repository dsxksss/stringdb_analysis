import pandas as pd


def compare_files(file1, file2):
    # 读取文件
    df1 = pd.read_csv(file1, sep="\t")
    df2 = pd.read_csv(file2, sep="\t")

    # 将node1和node2列组合成一个集合，忽略顺序
    df1["nodes_set"] = df1.apply(
        lambda row: frozenset([row["node1"], row["node2"]]), axis=1
    )
    df2["nodes_set"] = df2.apply(
        lambda row: frozenset([row["node1"], row["node2"]]), axis=1
    )

    # 比较数据
    merged_df = pd.merge(df1, df2, on="nodes_set", suffixes=("_file1", "_file2"))

    # 检查其余列的匹配情况
    columns_to_check = [
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

    for column in columns_to_check:
        merged_df[f"match_{column}"] = (
            merged_df[f"{column}_file1"] == merged_df[f"{column}_file2"]
        )
        merged_df[f"value_{column}_file1"] = merged_df[f"{column}_file1"]
        merged_df[f"value_{column}_file2"] = merged_df[f"{column}_file2"]

    # 获取匹配和不匹配的行
    matching_rows = merged_df[
        merged_df[[f"match_{column}" for column in columns_to_check]].all(axis=1)
    ]
    non_matching_rows = merged_df[
        ~merged_df[[f"match_{column}" for column in columns_to_check]].all(axis=1)
    ]

    return matching_rows, non_matching_rows


# 文件路径
file1 = "./result/HTR2B.tsv"
file2 = "./string_interactions.tsv"

matching_rows, non_matching_rows = compare_files(file1, file2)

# 输出匹配和不匹配的结果
matching_rows[
    [
        col
        for col in matching_rows.columns
        if col.startswith("match_") or col.startswith("value_")
    ]
].to_csv("./matching_rows.tsv", sep="\t", encoding="utf-8", index=False)

non_matching_rows[
    [
        col
        for col in non_matching_rows.columns
        if col.startswith("match_") or col.startswith("value_")
    ]
].to_csv("./non_matching_rows.tsv", sep="\t", encoding="utf-8", index=False)
