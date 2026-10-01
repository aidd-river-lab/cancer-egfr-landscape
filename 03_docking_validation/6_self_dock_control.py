#!/usr/bin/env python3
"""
阳性对照(self-docking control): 把奥希替尼对接回它真实结合的 6LUD 结构本身(不是目标的 7TVD),
检查 Vina 算出来的最佳姿势, 位置上是不是真的落在晶体结构里那个真实的配体位置附近。

这是对接方法学里的标准验证手段: 如果连"药物对接回它自己真实结合的结构"都对不准,
那这一整套准备流程(受体/配体制备、box设置)对目标结构(7TVD)算出来的结果就更没有参考价值;
反过来如果这里能对准, 才有理由认为流程本身是可信的, 7TVD 那边的探索性结果值得看一眼。

免责声明: 就算这里对准了, 也只能说明"这套流程在这一个已知阳性案例上表现正常", 不能证明
        它对 7TVD(不同的突变体/结构)一定同样可靠——只是提供了一点点信心, 不是数学保证。

产出:
  results/6LUD_receptor.pdbqt       从6LUD摘出的protein-only受体(用于自对接)
  results/6LUD_self_dock.log        Vina对接日志
  results/6LUD_self_dock_docked.pdbqt  对接出的姿势
"""
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

RESULTS_DIR = Path(__file__).resolve().parent / 'results'
BOX_PADDING = 5.0


def strip_ligand_from_pdb(pdb_path, ligand_comp_id, output_path):
    """去掉 HETATM 里属于配体(comp_id匹配)的记录, 只留蛋白ATOM(+其他HETATM比如水/离子), 得到 apo 受体"""
    lines = []
    for line in Path(pdb_path).read_text(encoding='utf-8').splitlines():
        if line.startswith('HETATM') and line[17:20].strip() == ligand_comp_id:
            continue
        lines.append(line)
    Path(output_path).write_text('\n'.join(lines) + '\n', encoding='utf-8')


def prepare_receptor_with_own_ligand_box(receptor_pdb, ligand_ref_pdb, output_basename):
    result = subprocess.run(
        [
            'mk_prepare_receptor.py',
            '--read_pdb', str(receptor_pdb),
            '-o', output_basename,
            '--write_pdbqt', f'{output_basename}.pdbqt',
            '--box_enveloping', str(ligand_ref_pdb),
            '--padding', str(BOX_PADDING),
            '-v', f'{output_basename}_vina_box.txt',
            '-a',  # 允许自动删掉侧链原子缺失、匹配不上标准残基模板的残基(常见于晶体结构里电子密度弱的柔性区域)
        ],
        capture_output=True, text=True, cwd=str(RESULTS_DIR)
    )
    if result.returncode != 0:
        raise RuntimeError(f'mk_prepare_receptor.py 失败:\nstdout: {result.stdout}\nstderr: {result.stderr}')
    print(result.stdout)


def run_vina(receptor_pdbqt, ligand_pdbqt, box_config, out_pdbqt, log_path):
    result = subprocess.run(
        ['vina', '--receptor', receptor_pdbqt, '--ligand', ligand_pdbqt,
         '--config', box_config, '--exhaustiveness', '8', '--num_modes', '9', '--out', out_pdbqt],
        capture_output=True, text=True, cwd=str(RESULTS_DIR)
    )
    (RESULTS_DIR / log_path).write_text(result.stdout + result.stderr, encoding='utf-8')
    if result.returncode != 0:
        raise RuntimeError(f'vina 对接失败(返回码{result.returncode}), 详见 {log_path}')
    return result.stdout


def parse_first_model_heavy_atoms(pdbqt_path):
    """从对接输出 PDBQT 里取 MODEL 1(最佳姿势)的重原子坐标(跳过氢, PDBQT里氢原子行以H开头的原子类型)"""
    coords = []
    in_model_1 = False
    for line in Path(pdbqt_path).read_text(encoding='utf-8').splitlines():
        if line.startswith('MODEL'):
            in_model_1 = line.split()[1] == '1'
            continue
        if line.startswith('ENDMDL'):
            if in_model_1:
                break
            continue
        if in_model_1 and (line.startswith('ATOM') or line.startswith('HETATM')):
            atom_type = line[77:79].strip()
            if atom_type.startswith('H'):
                continue
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return np.array(coords)


def parse_pdb_heavy_atoms(pdb_path):
    coords = []
    for line in Path(pdb_path).read_text(encoding='utf-8').splitlines():
        if line.startswith('HETATM'):
            element = line[76:78].strip() or line[12:14].strip()
            if element.startswith('H'):
                continue
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return np.array(coords)


def main():
    ref_receptor_raw = RESULTS_DIR / '6LUD_receptor_raw.pdb'
    ligand_ref = RESULTS_DIR / '6LUD_ligand_ref.pdb'
    ligand_pdbqt = RESULTS_DIR / 'osimertinib.pdbqt'
    for p in (ref_receptor_raw, ligand_ref, ligand_pdbqt):
        if not p.exists():
            print(f'找不到 {p}, 先跑 1_fetch_structures.py 和 2_prepare_ligand.py', file=sys.stderr)
            sys.exit(1)

    print('从 6LUD 摘出 apo 受体(去掉真实的 osimertinib 配体, 假装不知道它结合在哪) ...')
    apo_path = RESULTS_DIR / '6LUD_apo_for_self_dock.pdb'
    strip_ligand_from_pdb(ref_receptor_raw, 'YY3', apo_path)

    print('用 6LUD 自己的真实配体位置定义对接框(同一坐标系, 不需要跨结构叠合) ...')
    prepare_receptor_with_own_ligand_box(apo_path, ligand_ref, '6LUD_receptor')

    print('对接 osimertinib 回 6LUD(阳性对照) ...')
    log_text = run_vina(
        '6LUD_receptor.pdbqt', 'osimertinib.pdbqt', '6LUD_receptor_vina_box.txt',
        '6LUD_self_dock_docked.pdbqt', '6LUD_self_dock.log'
    )
    affinities = re.findall(r'^\s*1\s+(-?\d+\.?\d*)', log_text, re.MULTILINE)
    best_affinity = affinities[0] if affinities else None
    print(f'  最佳姿势 affinity: {best_affinity} kcal/mol')

    print('比较对接出的最佳姿势 vs 晶体结构里真实的配体位置 ...')
    docked_coords = parse_first_model_heavy_atoms(RESULTS_DIR / '6LUD_self_dock_docked.pdbqt')
    real_coords = parse_pdb_heavy_atoms(ligand_ref)
    docked_centroid = docked_coords.mean(axis=0)
    real_centroid = real_coords.mean(axis=0)
    centroid_distance = float(np.linalg.norm(docked_centroid - real_centroid))
    print(f'  对接姿势重原子数={len(docked_coords)}, 真实配体重原子数={len(real_coords)}'
          + (' (数量不一致, 可能是加氢/PDBQT展开方式导致, 质心比较仍可参考但不是逐原子RMSD)' if len(docked_coords) != len(real_coords) else ''))
    print(f'  质心距离: {centroid_distance:.2f} 埃'
          + ('（很接近, 说明对接找到了真实结合区域附近）' if centroid_distance < 3 else
             '（有明显偏差, 说明这套流程即便在阳性对照上也没能精确复现真实结合姿势, 需要谨慎看待7TVD那边的结果）'))

    print(f"\n已验证: 阳性对照真实跑完(Vina最佳affinity={best_affinity} kcal/mol), 对接姿势质心与真实晶体配体质心相距 {centroid_distance:.2f} 埃。")
    print('未验证: 质心距离是一个比较粗略的对齐指标(没有做逐原子的姿势RMSD, 那需要原子级别的对应关系, '
          '比 PDBQT 里的原子顺序匹配更复杂), 只能说明"大致落在同一个区域"还是"明显对不上", '
          '不能替代更严格的姿势RMSD评估。')


if __name__ == '__main__':
    main()
