# cancer-egfr-landscape

**临床突变频率 × 结构覆盖图谱** —— 一个学习驱动的计算药物设计（CADD）开源项目的第二个子仓库。

> 📘 **没有生物/化学背景？** 看 [`docs/explainer.md`](docs/explainer.md) —— 写给程序员的讲解，每个生物概念都配了编程类比，把三个模块每一步在干什么、踩了什么坑、为什么这么做，从头到尾讲了一遍。

## ⚠️ 项目定位（必读）

这是一个**学习驱动的开源计算工具项目**，产出的是**计算假说 + 可复现流程**，**不是药物发现成果，不构成任何医疗建议**。

- 本仓库所有"某突变位点是否有结构覆盖"的结论，含义仅仅是"RCSB PDB 上是否存在公开的、解析了抑制剂结合复合物的晶体结构"，**不代表该位点"有药可用"或"无药可用"**——现实中可能存在未公开、专利保护、或尚未结晶的药物设计工作。
- 任何计算输出都需要湿实验（体外/细胞/动物实验）验证才有科学意义。
- 如果你不是本项目作者、偶然看到这个仓库：请不要把这里的任何数字当作临床决策依据。

作者背景：10+ 年经验的程序员，**没有生物学/药物化学背景**，正在通过实操学习 CADD 基础流程。姊妹仓库 `target-egfr-c797s`（单靶点结构分析与虚拟筛选）本地没有找到，`03_docking_validation/` 没能直接复用它的脚本，而是按同样的思路（RDKit 处理配体、AutoDock Vina 做对接）独立实现了一遍，用真实 PDB 结构 + 真实跑起来的 Vina 走完了整条流程。

## 这个仓库做什么

把两类信息拼在一起看：

1. **临床上 EGFR 到底哪些位点突变得多**（来自真实患者测序数据，cBioPortal）
2. **这些位点里，哪些已经有公开的、带抑制剂共晶的晶体结构可以拿来做计算研究**（来自 RCSB PDB）

拼出来的图谱能回答一个很朴素的工程问题："当前计算药物设计社区的公开结构资源，覆盖了多少真实患者身上出现的耐药突变？还有哪些高频突变位点是'结构空白'，理论上更值得优先补做对接验证？"

### 目前跑出来的具体结论

用真实数据把三个模块全部跑了一遍，链条是这样的：

1. 在 TCGA 肺腺癌队列里统计 EGFR 突变，最高频的两个是 **L858R(31.5%)** 和 **19号外显子缺失/19del(27.4%)**。
2. 去 RCSB PDB 查这两个位点有没有"抑制剂共晶结构"：L858R 有（找到多个真实结构），**19del 没有**——公开数据库里唯一能查到的 19del 结构（7TVD）是不含任何药物的裸结构。
3. 于是把 19del 这个"结构空白"当成目标，用 AutoDock Vina 探索奥希替尼能不能对接进这个受体：
   - **先做阳性对照**：把 osimertinib 对接回它真实结合的 6LUD 结构，检验这套流程本身靠不靠谱——对接出的最佳姿势跟晶体结构里的真实配体位置只差 **0.45 埃**，说明流程本身是可信的。
   - **再对 7TVD(19del) 做探索性对接**：最佳姿势 affinity = **-7.357 kcal/mol**（阳性对照那边是 -8.196 kcal/mol，作为参照）。

也就是说，这个项目目前实际做出来的东西是：**用真实公开数据，找到了一个高频耐药突变（19del）在公开结构资源里是空白的，用一套经过阳性对照验证过的对接流程，对它跑出了一个探索性的结合假说**——但这仍然只是计算假说，不是"能结合"的证据，需要湿实验验证才有意义。

## 目录结构

```
cancer-egfr-landscape/
├── README.md                              本文件
├── docs/
│   └── methodology.md                     数据来源、聚合方法、已知局限性
├── 01_clinical_data/
│   ├── fetch_egfr_mutations.py            用 bravado 调用 cBioPortal API v3，拉取真实患者 EGFR 突变记录
│   ├── aggregate_mutation_frequency.py    把原始突变记录聚合成"位点 → 频率"表
│   └── results/                           输出目录（不入库真实患者数据的完整记录，只保留聚合后的频率表）
├── 02_structure_mapping/
│   └── map_structure_coverage.py          用 RCSB Search API + Data API 比对突变位点 × PDB 结构覆盖
├── 03_docking_validation/                 对模块2找出的结构空白位点(19del/7TVD)做探索性对接
│   ├── 1_fetch_structures.py              下载真实 PDB 结构(7TVD apo受体 + 6LUD配体参考)
│   ├── 2_prepare_ligand.py                RDKit生成奥希替尼3D构象 + meeko转PDBQT
│   ├── 3_prepare_receptor.py              CA原子叠合定位对接框(含独立交叉验证) + 生成受体PDBQT
│   ├── 4_dock_vina.sh                     AutoDock Vina 对接(需要本地装好 vina, 见文件内注释)
│   ├── 5_analyze_results.py               解析 Vina 日志, 按结合能排序
│   └── 6_self_dock_control.py             阳性对照: 把 osimertinib 对接回它真实结合的 6LUD, 验证流程本身可靠
└── requirements.txt
```

## 快速开始

```bash
# 需要 Python 3.10+ (meeko 用了 match 语句, 3.9 会直接语法错误)
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# AutoDock Vina 单独装, 不在 requirements.txt 里(pip装不了, 见下方"安装 Vina 踩坑记录")
conda create -n vina-env --override-channels -c conda-forge vina -y
conda activate vina-env   # 后续跑 4_dock_vina.sh / 6_self_dock_control.py 前需要 vina 在 PATH 里
                          # (如果同时还要用 .venv 里的 python/meeko, 把两边 bin 目录都加进 PATH 即可,
                          #  不需要真的同时 activate 两个环境)

# 第一步：拉取真实患者突变数据(需要能访问 cbioportal.org 的网络环境)
python 01_clinical_data/fetch_egfr_mutations.py

# 第二步：聚合成按位点统计的频率表
python 01_clinical_data/aggregate_mutation_frequency.py

# 第三步：突变位点 x PDB结构覆盖比对(默认只处理 count>=2 的位点, 见脚本内注释)
python 02_structure_mapping/map_structure_coverage.py

# 第四步：对结构空白位点做探索性对接
python 03_docking_validation/1_fetch_structures.py
python 03_docking_validation/2_prepare_ligand.py
python 03_docking_validation/3_prepare_receptor.py
bash 03_docking_validation/4_dock_vina.sh
python 03_docking_validation/5_analyze_results.py
python 03_docking_validation/6_self_dock_control.py   # 阳性对照, 建议先跑这个再相信上面的结果
```

### 安装 Vina 踩坑记录

- `pip install vina` 会失败：`Boost library location was not found`，这个 PyPI 包要从源码编译，需要先装 Boost。
- homebrew 没有现成的 `autodock-vina` 配方。
- 真正跑通的装法：`conda create --override-channels -c conda-forge vina`，装的是预编译好的 1.2.7，不用编译。`--override-channels` 是为了跳过 Anaconda `defaults` 渠道需要额外接受服务条款的问题，只用不受影响的 conda-forge 渠道。
- Vina 1.2.7 的命令行**没有 `--log` 参数**（比一些教程写的旧版本用法少这个选项），日志改用 shell 的 `tee` 自己存。

产出：`01_clinical_data/results/egfr_mutation_frequency.csv`（位点 | 出现次数 | 频率%）、`02_structure_mapping/results/structure_coverage.csv`（位点 | 频率 | 是否结构覆盖 | 对应PDB）、`03_docking_validation/results/docking_ranked.csv`（对接姿势按结合能排序）。

如果当前网络环境访问不了 cBioPortal/RCSB，对应脚本会明确报错退出，**不会**用编造的数字顶替——这种情况下 `aggregate_mutation_frequency.py` 和 `5_analyze_results.py` 都支持 `--demo` 参数，用字段 schema 与真实数据完全一致、但明确标注为"演示数据"的样例数据跑通逻辑本身。

## 数据来源

- **cBioPortal REST API v3**：`https://www.cbioportal.org/api/v3/api-docs`，默认查询队列 `luad_tcga_pan_can_atlas_2018`（TCGA 肺腺癌，PanCancer Atlas 版本）
- **RCSB PDB**：`https://www.rcsb.org`（`02_structure_mapping/`、`03_docking_validation/` 使用）

详细方法与已知局限性见 [`docs/methodology.md`](docs/methodology.md)。

## 当前进度

- [x] `01_clinical_data`：临床突变频率统计。用真实 TCGA LUAD PanCancer Atlas 队列数据跑通，L858R(31.5%) > 19del(27.4%) 两大主力突变，排序方向符合已知临床流行病学。
- [x] `02_structure_mapping`：突变位点 × 结构覆盖比对（"结构覆盖"定义为**要求该结构同时解析了抑制剂结合复合物**，而不只是出现过该突变位点）。真实结果：L858R、T790M 有覆盖；**19del 只找到 apo（无配体）结构 7TVD，按定义不算覆盖**——这是当前最值得关注的"公开结构空白"。
- [x] `03_docking_validation`：全流程真实跑通，包括 Vina 本身。真实下载 7TVD/6LUD、RDKit 生成奥希替尼 3D 构象、CA 原子叠合算出对接框（用已知口袋残基做了独立交叉验证，仅差 1.72 埃）、`conda install -c conda-forge vina` 装上 Vina 1.2.7、跑了阳性对照（osimertinib 对接回 6LUD，姿势质心跟真实晶体配体只差 0.45 埃）、对 7TVD(19del) 做了探索性对接（最佳姿势 -7.357 kcal/mol）。
