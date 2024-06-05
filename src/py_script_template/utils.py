import os
from datetime import datetime
import pandas as pd


def mkdir_if_not_exist(path: str) -> None:
    if not os.path.exists(path):
        os.makedirs(path)


def read_single_csv(
    input_path: str, cache_save_path: str, sep: str, dtype: dict
) -> pd.DataFrame:

    # Check if the cache file exists
    if os.path.exists(cache_save_path):
        return pd.read_parquet(cache_save_path)

    # Determine file size and set chunksize dynamically
    file_size = os.path.getsize(input_path)
    chunksize = calculate_chunksize(file_size)

    # Read the CSV in chunks
    df_chunk = pd.read_csv(
        input_path, chunksize=chunksize, sep=sep, encoding="utf-8", dtype=dtype
    )
    res_chunk = []
    for chunk in df_chunk:
        res_chunk.append(chunk)
    res_df = pd.concat(res_chunk)

    # Save to cache for future use
    res_df.to_parquet(cache_save_path)

    return res_df


def calculate_chunksize(file_size: int) -> int:
    # Approximate desired memory usage per chunk in bytes (e.g., 50MB)
    desired_memory_usage = 50 * 1024 * 1024  # 50MB

    # Estimate the number of rows per chunk based on file size
    # This estimation assumes average row size to be 1KB (tune as needed)
    average_row_size = 1024  # 1KB

    # Calculate the number of chunks needed
    num_chunks = file_size // desired_memory_usage

    # Calculate chunksize
    if num_chunks > 0:
        chunksize = file_size // num_chunks // average_row_size
    else:
        chunksize = (
            file_size // average_row_size
        )  # Handle case when file_size < desired_memory_usage

    return chunksize


def parse_timestamp(timestamp, custom_strfmt="%Y-%m-%d %H:%M:%S") -> str:
    # 将时间戳转换为datetime对象
    dt = datetime.fromtimestamp(timestamp)

    # 格式化日期和时间
    formatted = dt.strftime(custom_strfmt)

    return formatted
