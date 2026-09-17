# rs5068 全数据集再审计报告（re-audit）

**日期：** 2026-09-13 | **执行：** @scene#2 专家团队 + 主审实测复核
**目的：** 裁决 A 稿（rs5068→CAD/MI meta OR 0.752）与 B 稿（NT-proBNP→CAD/MI/HF 全 null）的矛盾
**数据：** OpenGWAS API 真实关联落盘（`replication_assoc_v2.json` / `mr_*_assoc.json`），无模拟数据
**暴露：** rs5068 G 等位基因 → NT-proBNP（SCALLOP, N=21,758）：β=0.1452, SE=0.0338（F=18.5）

---

## 1. 主裁决：CAD/MI 家族（7 个结局 GWAS 逐个单 SNP Wald）

| 数据集 | 结局 | OR (95% CI) | p |
|---|---|---|---|
| ieu-a-7 | CAD (CARDIoGRAMplusC4D, N=184k) | 0.969 (0.725–1.295) | 0.83 |
| ieu-a-798 | MI (UKB+C4D) | 1.090 (0.792–1.501) | 0.60 |
| GCST005195 | CAD (van der Harst, N=547k) | 0.769 (0.623–0.948) | 0.014 |
| GCST90013868 | CAD (GeneBANK/REGENIE-UKB) | 0.727 (0.555–0.954) | 0.021 |
| GCST011364 | MI (Hartiala, N=17.5k) | 0.715 (0.518–0.988) | 0.042 |
| finn-b-I9_CHD | CHD (FinnGen R9) | 0.781 (0.576–1.060) | 0.112 |
| finn-b-I9_MI | MI (FinnGen R9) | 0.820 (0.568–1.181) | 0.286 |
| **固定效应合并** | | **0.818 (0.735–0.912)** | **2.7×10⁻⁴** |
| **随机效应合并（DL）** | | **0.819 (0.733–0.914)** | **3.8×10⁻⁴** |

**异质性：Q=6.2 (df=6), p=0.40, I²=4%, τ²=0.0008 —— 几乎无异质性。**

**近似独立子集**（vdH + Hartiala + FinnGen CHD，剔除样本重叠）：**OR=0.759 (0.652–0.884), p=4×10⁻⁴, I²=0%** —— 信号在独立样本中更强。

排除：GCST90038610（Sakaue MI）beta/SE 量纲异常（原稿 a priori 排除规则）。

## 2. 参照家族

| 家族 | 数据集 | 结果 | 判定 |
|---|---|---|---|
| **HF** | HERMES OR 1.20 (0.95–1.52)；FinnGen 1.06 (0.75–1.49)；合并 1.15 (0.95–1.40) p=0.15, I²=0% | 真实 null（非功效问题：HERMES N≈97.7 万） |
| **AF** | Nielsen 1.487 (1.19–1.86)；Roselli 1.302 (1.07–1.59)；合并 **1.381 (1.19–1.60) p=2.2×10⁻⁵**, I²=0% | 风险信号**稳健**，非"仅两个数据集待确认" |
| **Stroke** | MEGASTROKE 0.90 ns；FinnGen 广义 0.63 (0.47–0.85)；心源性栓塞 **2.02 (1.12–3.63)**；I²=**84%** | 真实异质性：保护集中在非心源性/广义表型，心源性栓塞有害（与 AF 风险机制一致） |

## 3. 矛盾的裁决结论

1. **A 稿 OR 0.752 不是数据集挑选的假象**——纳入 B 稿的全部 null 数据集后，7 库合并仍 OR≈0.82、p≈3×10⁻⁴、I²=4%。红队"塌缩到 null"的预判**不成立**。
2. **B 稿的 null 是结局选择所致**：CAD 用了老版 C4D（ieu-a-7，N=184k，方向一致但不显著）与 ieu-a-798（null），且把 HF（真实 null 家族）混入"CAD/MI/HF 全 null"的表述。按家族拆分后，CAD/MI 保护、HF null，两者并不矛盾。
3. **真实图景是三家族分裂**（这比任何一稿现有的叙事都丰富）：
   - CAD/MI（动脉粥样硬化性）：**保护，OR≈0.76–0.82**
   - HF：**null**（1.15, ns）
   - AF + 心源性栓塞卒中：**有害**（1.38 / 2.02）
4. 对 B 稿的影响：B 稿"genetically instrumented NT-proBNP null for CAD/MI"的表述**需要修正**为"保护性（OR≈0.8）但幅度有限"；其"barometer not driver"主结论不受影响（保护方向反而支持"NP 升高非毒性"）。
5. 对 A 稿的影响：headline 应从 0.752 改为全库 0.82（或独立子集 0.76），HF null 须补入安全性叙事，AF/心源性栓塞风险须与保护信号同等权重呈现。

## 4. 样本重叠矩阵（定性，基于联盟构成）

| 数据集对 | 重叠评估 |
|---|---|
| vdH CAD ↔ ieu-a-7 (C4D) | 高：vdH 纳入全部 C4D 样本 |
| vdH CAD ↔ ieu-a-798 (UKB MI) | 高：共享 UKB 病例 |
| GeneBANK (REGENIE/UKB) ↔ ieu-a-798 / vdH | 高：同为 UKB 主体 |
| Hartiala MI ↔ C4D 家族 | 低（独立 MI 联盟） |
| FinnGen ↔ 全部其他 | 无（芬兰孤立人群） |
| HERMES HF ↔ UKB 家族 | 中（HERMES 含 UKB 等 47 队列） |

→ 有效独立样本数约 4–5；独立子集 meta（OR 0.76）已确认结论不依赖重叠样本。

## 5. 方法学附注

- Wald ratio SE 为一阶近似（仅用结局 SE），F=18.5 下弱工具偏倚方向朝 null → 估计保守。
- 回文校验：rs5068 (G/A) 非回文 SNP，无链 ambiguity。
- 复现路径：`rs5068_full_dataset_reaudit.py`（Wald+meta）→ `rs5068_meta_plot_overlap.py`（独立子集+图+矩阵）；全部中间数据在 `replication_assoc_v2.json` / `mr_v2_results_full.csv` 可溯。

## 6. 对投稿策略的回灌（决策树激活）

**重跑结果有利 → 完全合并方案激活**（架构师方案成立）：
- 合并版 headline：rs5068→CAD/MI 全 7 库随机效应 OR 0.82（I²=4%）+ 独立子集 0.76
- "barometer vs driver"弧获得新支柱：NP 升高观察性=坏（代偿标志），遗传性升高=动脉粥样硬化保护 + HF null + AF/心源性栓塞有害 —— 三家族分裂正是"通路干预需择轴而入（NPR3 清除轴优于分泌轴）"的最强论据
- A/B 互斥引用问题消解：统一用本审计的 7 库结果
- 预估合并版 CD 特刊录用概率：≈45–55%（前提：B 的 CVD-free null 叙事同步修复）

**产物清单**
- `rs5068_reaudit/rs5068_full_dataset_reaudit.py` — 主分析脚本
- `rs5068_reaudit/rs5068_meta_plot_overlap.py` — 独立子集/绘图/矩阵脚本
- `rs5068_reaudit/rs5068_wald_all_datasets.csv` — 全数据集 Wald 结果
- `rs5068_reaudit/rs5068_overlap_matrix.csv` — 重叠矩阵
- `rs5068_reaudit/rs5068_forest_CADMI.png` — CAD/MI 森林图
