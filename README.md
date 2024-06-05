# Stringdb Analysis

stringdb_analysis

## 支持参数

| 参数 | 描述 |
| --- | --- |
| --input_file | 要查找的基因名称 |
| --save_dir | 指定导出的数据目录路径, 默认./ |
| --output_file | 指定数据库目录, 默认./ |
| --cut_off | 截断值(0~1之间), 默认0.5(只导出combined_score为0.5以上的数据) |

### 使用例子

```bash
python src/david_analysis --input './genes.txt' --save_dir './result'
```
