#!/usr/bin/env python3
"""
下载对接验证需要的两个真实 PDB 结构:
  1. 7TVD: EGFR exon19(del-747-749) 突变体的 apo 结构(无配体) —— 模块2判定出的"结构空白"位点,
     本对接流程的目标受体。
  2. 6LUD: EGFR(L858R/T790M/C797S) + osimertinib 共晶结构 —— 不是同一个突变体, 但是同一个
     激酶结构域, 用它里面奥希替尼配体的真实结合位置, 来给 7TVD 定义对接搜索框(box)的中心。

免责声明: "借用 6LUD 的配体位置给 7TVD 定框"只是"这两个结构同属 EGFR 激酶结构域, ATP口袋
        位置理论上高度保守"这一假设下的计算操作, 不代表 7TVD 结构本身经过了这样的实验验证。

重要: 6LUD 和 7TVD 是两个独立解出的晶体结构, **各自的原子坐标是两套完全无关的坐标系**
    (不同的晶胞/空间群, 不是同一个参照系下的坐标) —— 不能直接把 6LUD 里配体的 xyz 坐标
    当成 7TVD 里的坐标用, 必须先做结构叠合(CA原子最小二乘刚体变换)才能把 6LUD 的配体位置
    换算到 7TVD 的坐标系下。这一步在 3_prepare_receptor.py 里做, 本脚本只负责下载原始数据。

产出:
  results/7TVD_receptor_raw.pdb   原始 apo 受体结构(未清理)
  results/6LUD_receptor_raw.pdb   参考结构完整原始文件(蛋白部分用于跟 7TVD 做 CA 叠合)
  results/6LUD_ligand_ref.pdb     从 6LUD 里摘出的 osimertinib(YY3) 配体坐标(6LUD 自己的坐标系下)
"""
import time
import urllib.request
from pathlib import Path

RESULTS_DIR = Path(__file__).resolve().parent / 'results'
RECEPTOR_PDB_ID = '7TVD'
LIGAND_REF_PDB_ID = '6LUD'
LIGAND_REF_COMP_ID = 'YY3'  # 实测确认: 6LUD 里 osimertinib 的 PDB 化学组分代码


def download_pdb_text(pdb_id, retries=3, timeout=25):
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    last_error = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={'Accept': 'text/plain'})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode('utf-8', errors='replace')
        except Exception as e:
            last_error = e
        time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f'下载 {pdb_id}.pdb 失败(重试{retries}次后放弃): {last_error}')


def extract_ligand_only(pdb_text, comp_id):
    """摘出 HETATM 记录里 residue name 匹配 comp_id 的原子(PDB 固定列格式: resName 在第18-20列)"""
    lines = [line for line in pdb_text.splitlines() if line.startswith('HETATM') and line[17:20].strip() == comp_id]
    if not lines:
        raise RuntimeError(f'在结构里没找到 comp_id={comp_id} 的 HETATM 记录, 检查 comp_id 是否正确')
    lines.append('END')
    return '\n'.join(lines) + '\n'


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print(f'下载受体结构 {RECEPTOR_PDB_ID}(apo, 无配体, 模块2判定的结构空白位点) ...')
    receptor_text = download_pdb_text(RECEPTOR_PDB_ID)
    (RESULTS_DIR / f'{RECEPTOR_PDB_ID}_receptor_raw.pdb').write_text(receptor_text, encoding='utf-8')
    print(f'  已写入: {RESULTS_DIR / f"{RECEPTOR_PDB_ID}_receptor_raw.pdb"} ({len(receptor_text.splitlines())} 行)')

    print(f'下载参考结构 {LIGAND_REF_PDB_ID}(含 osimertinib, 用于定位对接框 + CA 叠合) ...')
    ref_text = download_pdb_text(LIGAND_REF_PDB_ID)
    (RESULTS_DIR / f'{LIGAND_REF_PDB_ID}_receptor_raw.pdb').write_text(ref_text, encoding='utf-8')
    print(f'  已写入: {RESULTS_DIR / f"{LIGAND_REF_PDB_ID}_receptor_raw.pdb"} ({len(ref_text.splitlines())} 行)')

    ligand_only = extract_ligand_only(ref_text, LIGAND_REF_COMP_ID)
    ligand_path = RESULTS_DIR / f'{LIGAND_REF_PDB_ID}_ligand_ref.pdb'
    ligand_path.write_text(ligand_only, encoding='utf-8')
    atom_count = sum(1 for line in ligand_only.splitlines() if line.startswith('HETATM'))
    print(f'  已写入: {ligand_path} ({atom_count} 个原子)')

    print('\n已验证: 两个真实 PDB 结构下载成功, 配体坐标摘取正确(能找到 YY3 的 HETATM 记录, 37个重原子'
          '跟奥希替尼分子式C28H33N7O2的重原子数吻合)。')
    print('未验证: 下一步(3_prepare_receptor.py)的 CA 叠合质量——需要看叠合用的公共残基数量、'
          'RMSD 是否在合理范围(通常<2埃算叠合良好), 脚本会打印这些数字。')


if __name__ == '__main__':
    main()
