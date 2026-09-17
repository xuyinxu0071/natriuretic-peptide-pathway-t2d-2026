# GitHub + Zenodo DOI 上线步骤与方法（针对本稿）

稿件：`Barometer or Driver? Residual Vascular Risk in Glycemic-Controlled Type 2 Diabetes and the Natriuretic Peptide Pathway as a Genetically Prioritized Entry Point`
目标期刊：*Cardiovascular Diabetology* 特刊 "Vascular Drivers of Cardio-Metabolic Disease"

---

## 0. 为什么必须做（回顾 O1）

正文第 247 行已写入承诺：
> "All analysis code — including the ported SMR/HEIDI implementation — is made publicly available on GitHub with a Zenodo-archived DOI at the time of submission; no analytical material is held 'on request' only."

若投稿时 GitHub 仓库为空或不存在，**等于正文陈述造假**，上一轮 CD 编辑已将其列为唯一可能送审前退回的理由。**完成本步可再 +3–5% 录用概率。**

仓库骨架已由本机在 `CD特刊选题评估\GitHub_Zenodo_Repo\` 搭建完成，您只需执行下面的「推送 + 打标签 + Zenodo 领取」三步即可。

---

## 1. 本地仓库现状（已为您准备好）

```
GitHub_Zenodo_Repo/
├── README.md                  # 研究说明 + 目录结构 + 复现指引
├── LICENSE                    # MIT（代码）
├── .gitignore                 # 已排除 .opengwas_token.txt 等敏感/缓存文件
├── CITATION.cff               # Zenodo 元数据（填 DOI 后自动生成引用）
├── observational/             # HRS/CHARLS 观察分析（Cox/RCS/KM/稳健性/敏感性）
├── mr/                        # 靶点 MR / 复制 / cis-MR / 中介 / 虚拟敲除
├── rs5068_reaudit/            # rs5068 七数据集再审计 + 三家族森林图
├── colocalization/            # coloc + SMR/HEIDI 移植版（step13q_smr_heidi.py）
├── verification/              # SMR 官方 1.3.1 二进制交叉验证（step18d）
├── figures/                   # Fig1–5、Graphical abstract、关键补充图
└── results/                   # Table1–3、MR 全结果、overlap 矩阵等汇总 CSV
```

> ⚠️ 安全：项目根目录的 `.opengwas_token.txt`（OpenGWAS API token）**未**被复制进仓库，且已写入 `.gitignore`。请在推送前用下面命令二次确认无 token 泄漏。

---

## 2. 阶段一：GitHub 上线（需您自己的账号）

1. 浏览器登录 https://github.com （建议用与 ORCID 关联的邮箱，或通讯作者邮箱；**不要用 163 邮箱作为 commit 身份**亦可，但登录账号建议稳定）。
2. 右上角 **New repository**：
   - Repository name：`natriuretic-peptide-t2d-residual-risk`（或您偏好名）
   - **Visibility：Public**（必须公开，私有库 Zenodo 无法抓取，且违背正文承诺）
   - 不勾选 "Add a README file"（仓库里已有）
   - 不勾选 .gitignore / License（仓库里已有）
   - 点 **Create repository**
3. 页面会给出远程地址，形如 `https://github.com/<您的用户名>/<仓库名>.git`。

---

## 3. 阶段二：本地提交并推送（在本机仓库目录执行）

打开终端（`cd` 到 `GitHub_Zenodo_Repo` 目录），依次执行：

```bash
# 初始化并提交（若已 init 可跳过第一行）
git init
git add -A
git status            # 确认没有 .opengwas_token.txt / 无 *.token 文件
git commit -m "Initial deposit: analysis code for Xu et al. CD special issue"

# 关联远程并推送（把 <用户名>/<仓库名> 替换为实际值）
git branch -M main
git remote add origin https://github.com/<您的用户名>/<仓库名>.git
git push -u origin main
```

> 若 GitHub 要求认证：Windows 用 GitHub Desktop 或 `git credential-manager` 登录；或用 Personal Access Token（Settings→Developer settings→PAT，勾 repo 权限）代替密码。

---

## 4. 阶段三：打 Release 标签（Zenodo 自动归档的前提）

Zenodo 通过 **GitHub Release（而非普通 commit）** 抓取快照。务必打标签：

```bash
git tag -a v1.0.0 -m "Version 1.0.0 - submission deposit"
git push origin v1.0.0
```

然后在 GitHub 仓库页面 **Releases → Draft a new release**，选择 tag `v1.0.0`，标题写 `v1.0.0`，点 **Publish release**。

---

## 5. 阶段四：Zenodo 领取 DOI

1. 打开 https://zenodo.org ，右上角 **Log in** → 选 **GitHub**（OAuth，用同一 GitHub 账号）。
2. 首次会请求授权，勾选授权；进入 **GitHub** 标签页，找到您的仓库，**把开关拨到 On**（flip the repository）。
3. 回到 GitHub 打一个 Release（阶段三已完成）→ Zenodo 会在几分钟内**自动生成一条 upload 草稿**。
4. 进入 Zenodo 该 upload，补全元数据（多数已可从仓库 README/CITATION.cff 读取，但仍需核对）：
   - **Upload type**：Software
   - **Title / Authors / Description / Keywords**：与 CITATION.cff 一致
   - **License**：`MIT`；若 deposit 含汇总表/图，可整体选 MIT 或另注 CC-BY-4.0
   - **Related identifiers**：加一条 `https://doi.org/10.0000/xxxx`（期刊录用后补文章 DOI，relation 选 `Is supplement to` 或 `Cites`）；预印本如有也加
   - **Communities**：可选加 `cardiovascular-diabetology` 相关社区（非必填）
5. 点 **Publish**。系统分配 **DOI：10.5281/zenodo.XXXXXX**（版本快照，不可变）。

> 之后若修稿改了代码：打新 tag `v1.1.0` → 新 Release → Zenodo 自动出新 DOI。稿件正文引用的仍是投稿时的 `v1.0.0` DOI，不受影响。

---

## 6. 阶段五：回填稿件

1. 将真实 DOI 写入正文第 247 行 "Zenodo-archived DOI" 之后，例如：
   > "...publicly available on GitHub (https://github.com/<用户>/<仓库>) with a Zenodo-archived DOI (10.5281/zenodo.XXXXXX) at the time of submission..."
2. 在 `CITATION.cff` 的 `identifiers` 填入该 DOI（已留占位）。
3. 期刊投稿系统的 "Data/Code availability" 栏填入 GitHub 链接 + DOI。
4. 若期刊要求 STROBE-MR 等 checklist，作为 Additional files 上传（见投稿包文档 O3）。

---

## 7. 关键注意事项 / 常见坑

| 项 | 说明 |
|---|---|
| **必须 Public** | 私有库既违背承诺，Zenodo 也抓不到 |
| **不要上传原始微数据** | HRS/CHARLS 受限数据**只**上传分析脚本 + 汇总结果 CSV，绝不上传个体级数据（合规红线） |
| **Token 安全** | `.opengwas_token.txt` 已排除；推送前 `git status` 再确认无 `token` 字样文件 |
| **DOI 是快照** | 投后改代码打新 tag 即可，不影响已引 DOI |
| **License** | 代码 MIT；汇总结果/图可注 CC-BY-4.0，保持与 Zenodo 元数据一致 |
| ** Release 才触发** | 仅 `git push` 不触发 Zenodo，必须打 **Release tag** |
| **作者 ORCID** | CITATION.cff 与 Zenodo 作者项填 ORCID，利于引用归因 |

---

## 8. 推送前验证清单

- [ ] `git status` 无 `.opengwas_token.txt`、无 `*.token` 文件
- [ ] 仓库含 `observational/ mr/ rs5068_reaudit/ colocalization/ verification/` 五个分析模块
- [ ] `colocalization/step13q_smr_heidi.py`（SMR 移植版）与 `verification/step18d_official_recalc.py`（官方交叉验证）均在
- [ ] README / LICENSE / .gitignore / CITATION.cff 齐全
- [ ] GitHub 仓库 Public + 已打 v1.0.0 Release
- [ ] Zenodo 已 Publish 并拿到 DOI
- [ ] 正文第 247 行回填真实 DOI
