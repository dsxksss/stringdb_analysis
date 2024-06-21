# Stringdb Analysis

stringdb_analysis

## 支持参数

| 参数 | 描述 |
| --- | --- |
| --input_file | 要查找的基因名称 |
| --db_dir | 指定数据库目录, 默认/data/PRG/tools/TCM/dbs/stringdb_analysis |
| --save_dir | 指定导出的数据目录路径, 默认./ |
| --output_file | 指定数据库目录, 默认./ |
| --cut_off | 截断值(0~1之间), 默认0.5(只导出combined_score为0.5以上的数据) |
| --get_genes_correlation| 是否获取输入基因的其他关联数据(Yes/No)大小写均可, 默认No不获取   |

### 使用例子

```bash
python src/stringdb_analysis --input './genes.txt' --save_dir './result'
```
