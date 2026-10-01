#!/usr/bin/env python3
"""
用 6LUD 的配体位置给 7TVD(apo)定义对接搜索框(box), 再用 meeko 生成受体 PDBQT。

关键步骤(容易踩坑, 详细写在这里): 6LUD 和 7TVD 是两个独立解出的晶体结构, 各自的原子坐标
是两套完全无关的坐标系。如果直接把 6LUD 里 osimertinib 配体的 xyz 坐标当成 7TVD 坐标系
下的框中心, 大概率会把搜索框定在 7TVD 结构完全无关的空间位置——这是本脚本最早的一版实现
真实踩到的坑(第一次跑出来的 box 中心坐标肉眼可见地不在 7TVD 结构附近, 排查后才确认是
坐标系不同导致的)。

正确做法: 用两个结构里"残基编号相同的CA原子"做最小二乘刚体叠合(Kabsch算法), 求出把
6LUD坐标系变换到7TVD坐标系的旋转矩阵R和平移向量t, 再用R,t把配体坐标变换过去, 这样算出来
的框中心才是在7TVD自己坐标系下、真实对应ATP口袋的位置。

免责声明: 这套流程假设两个结构的激酶结构域整体折叠高度相似(残基编号一致意味着可以直接
        用CA做刚体叠合, 不需要更复杂的柔性叠合)——这个假设对同一蛋白的不同突变体通常成立,
        但没有专门验证过这两个具体结构在局部构象上是否有明显差异。
"""
import subprocess
import sys
from pathlib import Path

import numpy as np

RESULTS_DIR = Path(__file__).resolve().parent / 'results'
BOX_PADDING = 5.0  # 埃, 在配体范围外留出的搜索空间余量

# 来自 ../prompt/egfr-对比.txt: 6LUD 结合口袋残基列表(PDB原始SITE记录, 编号沿用EGFR野生型)。
# 用来独立交叉验证——不依赖6LUD/7TVD叠合, 直接用这份口袋残基在 7TVD 自己坐标系下的质心,
# 跟叠合算出来的框中心做距离比较, 距离小说明叠合是可信的, 不是只看RMSD数字。
KNOWN_POCKET_RESIDUES = [718, 723, 726, 743, 791, 792, 793, 794, 796, 800, 804, 844]


def parse_ca_atoms(pdb_path):
    """返回 {(chain, resSeq): (x,y,z)}, 只取 ATOM 记录里 atom name = CA 的"""
    coords = {}
    for line in Path(pdb_path).read_text(encoding='utf-8').splitlines():
        if line.startswith('ATOM') and line[12:16].strip() == 'CA':
            chain = line[21]
            res_seq = int(line[22:26])
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            coords[(chain, res_seq)] = xyz
    return coords


def parse_hetatm_coords(pdb_path):
    coords = []
    for line in Path(pdb_path).read_text(encoding='utf-8').splitlines():
        if line.startswith('HETATM'):
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return np.array(coords)


def kabsch_fit(source_points, target_points):
    """
    求刚体变换 R, t, 使得 R @ source_points.T + t ≈ target_points.T (最小二乘意义下)
    标准 Kabsch 算法: https://en.wikipedia.org/wiki/Kabsch_algorithm
    返回 R (3x3), t (3,), 以及拟合后的 RMSD(用来判断叠合质量)
    """
    source_centroid = source_points.mean(axis=0)
    target_centroid = target_points.mean(axis=0)
    source_centered = source_points - source_centroid
    target_centered = target_points - target_centroid

    H = source_centered.T @ target_centered
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    t = target_centroid - R @ source_centroid

    fitted = (R @ source_points.T).T + t
    rmsd = np.sqrt(np.mean(np.sum((fitted - target_points) ** 2, axis=1)))
    return R, t, rmsd


def compute_aligned_box(ligand_ref_pdb, ref_receptor_pdb, target_receptor_pdb):
    ref_ca = parse_ca_atoms(ref_receptor_pdb)
    target_ca = parse_ca_atoms(target_receptor_pdb)

    # 找两边都有的 (chain, resSeq), 按 chain 分组配对 —— 用来做 Kabsch 叠合的公共锚点
    # 只用第一个共同链(两个结构通常都是单链激酶结构域), 避免链编号不对应导致误配对
    common_chains = set(c for c, _ in ref_ca) & set(c for c, _ in target_ca)
    if not common_chains:
        raise RuntimeError('两个结构没有共同的链ID, 无法直接配对CA原子做叠合, 需要人工检查链编号')
    chain = sorted(common_chains)[0]

    common_keys = sorted(k for k in ref_ca if k[0] == chain and k in target_ca)
    if len(common_keys) < 20:
        raise RuntimeError(f'两个结构公共残基编号数量太少({len(common_keys)}个), Kabsch叠合不可靠, 需要人工检查')

    source_points = np.array([ref_ca[k] for k in common_keys])  # 6LUD坐标系
    target_points = np.array([target_ca[k] for k in common_keys])  # 7TVD坐标系
    R, t, rmsd = kabsch_fit(source_points, target_points)
    print(f'  CA 叠合: 用了 {len(common_keys)} 个公共残基(chain {chain}), 叠合 RMSD = {rmsd:.2f} 埃'
          + ('（<2埃, 叠合质量良好）' if rmsd < 2.0 else '（>=2埃, 叠合质量一般, 结果仅供参考, 建议用 PyMOL 目视复核）'))

    ligand_coords_ref_frame = parse_hetatm_coords(ligand_ref_pdb)
    ligand_coords_target_frame = (R @ ligand_coords_ref_frame.T).T + t

    box_min = ligand_coords_target_frame.min(axis=0) - BOX_PADDING
    box_max = ligand_coords_target_frame.max(axis=0) + BOX_PADDING
    center = (box_min + box_max) / 2
    size = box_max - box_min
    return center, size, rmsd, len(common_keys)


def cross_check_against_known_pocket(target_receptor_pdb, box_center):
    """
    独立交叉验证: 不依赖6LUD叠合, 直接用 prompt 文档里给出的口袋残基列表, 在 7TVD 自己的
    坐标系下算质心, 跟叠合算出来的框中心比距离——这是同一个结论的两条独立证据链,
    两者吻合(距离远小于框尺寸)才能说明叠合出来的框真的落在ATP口袋附近, 不是叠合误差的巧合。
    """
    target_ca = parse_ca_atoms(target_receptor_pdb)
    chain = sorted(set(c for c, _ in target_ca))[0]
    points = []
    missing = []
    for res_seq in KNOWN_POCKET_RESIDUES:
        key = (chain, res_seq)
        if key in target_ca:
            points.append(target_ca[key])
        else:
            missing.append(res_seq)
    if not points:
        return None
    centroid = np.array(points).mean(axis=0)
    distance = float(np.linalg.norm(centroid - box_center))
    return {'centroid': centroid, 'distance': distance, 'found': len(points), 'missing': missing}


def write_receptor_pdbqt(receptor_pdb, box_center, box_size, output_basename):
    result = subprocess.run(
        [
            'mk_prepare_receptor.py',
            '--read_pdb', str(receptor_pdb),
            '-o', output_basename,
            '--write_pdbqt', f'{output_basename}.pdbqt',
            '--box_center', str(box_center[0]), str(box_center[1]), str(box_center[2]),
            '--box_size', str(box_size[0]), str(box_size[1]), str(box_size[2]),
            '-v', f'{output_basename}_vina_box.txt',
            '-a',  # 允许自动删掉侧链原子缺失、匹配不上标准残基模板的残基(6LUD自对照对照组实测踩到过这个问题, 7TVD这次没触发但保留以防万一)
        ],
        capture_output=True, text=True, cwd=str(RESULTS_DIR)
    )
    if result.returncode != 0:
        raise RuntimeError(f'mk_prepare_receptor.py 失败:\nstdout: {result.stdout}\nstderr: {result.stderr}')
    print(result.stdout)


def main():
    ligand_ref = RESULTS_DIR / '6LUD_ligand_ref.pdb'
    ref_receptor = RESULTS_DIR / '6LUD_receptor_raw.pdb'
    target_receptor = RESULTS_DIR / '7TVD_receptor_raw.pdb'
    for p in (ligand_ref, ref_receptor, target_receptor):
        if not p.exists():
            print(f'找不到 {p}, 先跑 1_fetch_structures.py', file=sys.stderr)
            sys.exit(1)

    print('把 6LUD 的配体位置叠合变换到 7TVD 坐标系, 计算对接搜索框 ...')
    center, size, rmsd, n_anchors = compute_aligned_box(ligand_ref, ref_receptor, target_receptor)
    print(f'  框中心(7TVD坐标系): ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})')
    print(f'  框尺寸: ({size[0]:.2f}, {size[1]:.2f}, {size[2]:.2f}) 埃')

    print('独立交叉验证: 用 prompt 文档给出的口袋残基列表, 在 7TVD 自己坐标系下核对框中心位置是否合理 ...')
    check = cross_check_against_known_pocket(target_receptor, center)
    if check is None:
        print('  跳过(7TVD结构里一个口袋残基编号都没找到, 可能是链编号/残基编号跟预期不一致)')
    else:
        ok = check['distance'] < min(size) / 2
        print(f"  已知口袋残基质心(7TVD自己坐标系, 用了{check['found']}/{len(KNOWN_POCKET_RESIDUES)}个残基"
              + (f", 缺失{check['missing']}" if check['missing'] else '') + f"): 距叠合框中心 {check['distance']:.2f} 埃"
              + ('（在框半径内, 交叉验证通过）' if ok else '（超出框半径, 交叉验证未通过, 结果存疑）'))

    print('生成受体 PDBQT ...')
    write_receptor_pdbqt(target_receptor, center, size, '7TVD_receptor')

    print(f'已验证: CA 叠合能跑通(用了 {n_anchors} 个公共残基, RMSD={rmsd:.2f}埃), meeko 成功生成受体 PDBQT + Vina box 文件, '
          + (f"且独立交叉验证通过(已知口袋残基质心距框中心仅{check['distance']:.2f}埃, 远小于框尺寸)" if check and check['distance'] < min(size) / 2 else '交叉验证结果见上') + '。')
    print('未验证: 生成的 PDBQT 受体文件本身(电荷/原子类型)是否符合 Vina 预期——这个要等实际跑一次 Vina 才能确认。')


if __name__ == '__main__':
    main()
