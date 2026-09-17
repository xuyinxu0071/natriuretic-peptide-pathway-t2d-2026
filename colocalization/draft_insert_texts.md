# Coloc 手稿插入文本草稿（暂存区 — 结果数值待 step5 完成后填充）

## Methods 2.6 追加段（在 mediation MR 小节之后）

**Locus-wide colocalization.** To formally distinguish a shared causal variant from distinct causal variants in linkage at the two lead loci, we conducted Bayesian colocalization (coloc.abf [ref: Giambartolomei et al., PLoS Genet 2014]) between each exposure dataset and each outcome GWAS at the rs5068 (NPPA/NPPB) and rs1421811 (NPR3) loci. The analysis window was ±500 kb around each lead variant (GRCh37). Locus-wide common (minor allele frequency ≥ 1%) biallelic SNVs with allele frequencies were obtained from dbSNP build 151 (GRCh37 p13, common-variant release) via tabix-guided regional extraction. Association statistics (beta, SE, effect allele, allele frequency, sample size) for every variant were retrieved from the same OpenGWAS datasets used in the MR analysis [exposures: eQTLGen cis-eQTL for NPPA and NPR3 (whole blood), Evangelou et al. systolic blood pressure; outcomes: CARDIoGRAMplusC4D (ieu-a-7), van der Harst et al. CAD (ebi-a-GCST005195), Nielsen et al. atrial fibrillation (ebi-a-GCST006061), HERMES heart failure (ebi-a-GCST009541), MEGASTROKE any stroke (ebi-a-GCST005838) and cardioembolic stroke (ebi-a-GCST006910), Hartiala et al. MI (ebi-a-GCST011364), and FinnGen R9 MI (finn-b-I9_MI)]. Alleles were harmonized to the dbSNP reference/alternate orientation (strand and effect-allele alignment; palindromic variants and frequency-discordant variants excluded), and analyses were restricted to variants present in both the exposure and outcome dataset. Quantitative exposures were analyzed as type = "quant" (eQTL effects on SD-scale expression, sdY = 1; SBP with sdY estimated from the effect-size distribution), binary outcomes as type = "cc" with case fractions from the source GWAS. We report posterior probabilities for hypotheses H0–H4, with H4 (a single shared causal variant) as the colocalization hypothesis. As sensitivity analyses, all models were re-run with the colocalization prior p12 varied one order of magnitude in each direction (1×10⁻⁴ and 1×10⁻⁶). Analyses used R 4.5.2 with the coloc package; full per-variant inputs and outputs are provided in the Supplement (Table S11).

## Results 新增小节（3.3 之后，作为 3.4 前插入 — 编号需顺延调整）

**Locus-wide colocalization.** [待填充: 每对 PP.H4 数值 + nsnps + 敏感性]

## 4.3 修改（替换 "will require locus-wide conditional and colocalization analysis, which we have defined as the next analytic step"）

[待填充]

## 4.5 Strengths 追加

[待填充]

## 4.6 Limitations 重写段（替换 "Formal colocalization analysis ... do not present them as such."）

[待填充]

## 5. Conclusions 修改

[待填充]

## Abstract 修改（Conclusions 句）

[待填充]
