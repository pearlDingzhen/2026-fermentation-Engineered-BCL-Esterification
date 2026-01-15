# Lipase 突变评分处理工具

## 功能说明

该模块用于处理 LigandMPNN 原始输出，将 `.pt` 格式转换为可读的 CSV 表格。

## 使用方法

```bash
python process_scores.py -i <输入.pt文件> -o <输出.csv文件>
```

## 输入输出

### 输入格式

LigandMPNN 生成的 `.pt` 文件，包含：
- `sequence`: 野生型氨基酸序列
- `mean_of_probs`: 20种氨基酸的概率分布矩阵

### 输出格式

CSV 表格，包含：
- **数据区 (20行 × 284列)**: 每行一种氨基酸，每列一个残基位点
- **附加行**:
  - `Sequence`: 野生型序列
  - `Max_vaules`: 每个位置的最优突变建议
  - `Mutations`: 建议的突变氨基酸

## 计算方法

### 1. 突变分数计算

对每个残基位点 \( i \) 和目标氨基酸 \( aa \):

\[
\text{Score}_{i \to aa} = \log(P_{aa}) - \log(P_{wt,i})
\]

### 2. 最优突变识别

对于每个残基位点，选择分数最高的突变：

\[
\text{BestMutation}_i = \arg\max_{aa} \text{Score}_{i \to aa}
\]

### 3. 显著性筛选

仅输出分数增量 ≥ 0.4 的突变建议：

\[
\text{if } \max_{aa}(\text{Score}_{i \to aa}) \geq 0.4: \text{report}
\]

## 输出示例

```
# CSV 文件结构
       1     2     3   ...   284
A   0.12   0.05  0.08  ...   0.15
C   0.03   0.21  0.01  ...   0.02
...
Y   0.45   0.12  0.33  ...   0.28

Sequence: HPVFVLVHGAWHGAWCYAHVAAALA...
Max_vaules: A1C:0.45 A2F:0.32 ...
Mutations: C F ...
```

## 代码实现

```python
import argparse
import torch
import pandas as pd
import numpy as np

def main(input_file, output_file):
    # 1. 加载.pt文件
    data = torch.load(input_file)
    sequence = data["sequence"]
    
    # 2. 创建DataFrame
    df_sub = pd.DataFrame(data["mean_of_probs"], index=index)
    
    # 3. 计算野生型参考值
    wt_values = [df_sub.loc[sequence[i], df_sub.columns[j]] 
                 for j, i in enumerate(range(len(sequence)))]
    
    # 4. 计算突变分数（对数差异）
    for col in df_sub.columns:
        df_sub[col] = round(np.log(df_sub[col].astype(float)) 
                          - np.log(wt_values[df_sub.columns.get_loc(col)]), 2)
    
    # 5. 识别最优突变
    max_values = round(df_sub.max(axis=0), 2)
    for column_name in df_sub.columns:
        max_index = np.argmax(df_sub[column_name].values)
        column_number = re.search(r'\d+', column_name).group()
        max_label = sequence[i] + column_number + index[max_index] + ":" + str(max_value)
        
        # 6. 筛选显著性突变（阈值 0.4）
        if float(max_value) >= 0.4:
            print(max_label)
    
    # 7. 保存结果
    df_sub.loc['Sequence'] = list(sequence)
    df_sub.loc['Max_vaules'] = list(max_labels)
    df_sub.loc['Mutations'] = list(mutations)
    df_sub.to_csv(output_file)
```

## 注意事项

1. **对数变换**: 使用自然对数确保数值的可加性
2. **参考点**: 以野生型氨基酸的概率为基准
3. **阈值**: 0.4 的筛选阈值可根据具体需求调整
