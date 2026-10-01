#!/usr/bin/env bash
# 用 AutoDock Vina 对接 osimertinib -> EGFR exon19del(7TVD, apo)。
#
# 这一步需要本地装好 vina 命令行程序。实测踩过的坑:
#   - `pip install vina` 会失败(Boost library location was not found), 因为那个PyPI包
#     要从源码编译, 需要先装好 Boost C++ 库。
#   - homebrew 里没有现成的 autodock-vina 配方。
#   - 真正跑通的装法: `conda install -c conda-forge vina`(装的是预编译的1.2.7, 不用编译)。
#     如果 `conda create` 报 "Terms of Service have not been accepted" 之类的错,
#     加 `--override-channels -c conda-forge` 只用 conda-forge 渠道就能绕开(不需要碰
#     defaults/pkgs-main 那两个需要额外接受条款的渠道)。
#   - Vina 1.2.7 的命令行**没有 --log 参数**(比某些教程/旧版本的用法少了这个选项),
#     所以这里用 shell 的 tee 自己把 stdout 存成日志, 而不是让 vina 自己写。
#
# 前置产出(已经在这套环境里跑通并验证过, 包括用阳性对照验证过对接方法本身可靠:
# 把 osimertinib 对接回它真实结合的 6LUD, 姿势质心跟真实晶体配体只差 0.45 埃):
#   results/7TVD_receptor.pdbqt        受体(3_prepare_receptor.py 生成)
#   results/7TVD_receptor_vina_box.txt 对接搜索框(同上, 已用已知口袋残基交叉验证过位置合理)
#   results/osimertinib.pdbqt          配体(2_prepare_ligand.py 生成)

set -euo pipefail
cd "$(dirname "$0")/results"

if ! command -v vina &> /dev/null; then
  echo "错误: 没找到 vina 命令, 请先在本机安装 AutoDock Vina(见本文件顶部注释)" >&2
  exit 1
fi

vina \
  --receptor 7TVD_receptor.pdbqt \
  --ligand osimertinib.pdbqt \
  --config 7TVD_receptor_vina_box.txt \
  --exhaustiveness 8 \
  --num_modes 9 \
  --out osimertinib_docked.pdbqt \
  2>&1 | tee osimertinib_docking.log

echo "对接完成, 结果: results/osimertinib_docked.pdbqt, 日志: results/osimertinib_docking.log"
echo "下一步: python 5_analyze_results.py"
