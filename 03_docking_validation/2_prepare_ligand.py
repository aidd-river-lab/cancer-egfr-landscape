#!/usr/bin/env python3
"""
用 RDKit 把奥希替尼(osimertinib)的 SMILES 变成一个带 3D 构象、加了氢的分子文件(SDF),
再用 meeko(mk_prepare_ligand.py, 官方 CLI, 不是自己拼底层API)转成 Vina 能吃的 PDBQT。

免责声明: 生成的3D构象是"能量最小化后的一个合理构象", 不是这个分子结合到 7TVD 时的真实构象——
        构象本身就是对接程序要去搜索优化的东西, 这里只是给对接程序一个合理的起点。

SMILES 来源: ../prompt/egfr-对比.txt 里给出的、已用 RDKit 验证过的奥希替尼 SMILES
           (分子式 C28H33N7O2, MW 499.62, LogP 4.51)。
"""
import subprocess
import sys
from pathlib import Path

RESULTS_DIR = Path(__file__).resolve().parent / 'results'

# 来自 prompt/egfr-对比.txt: 已用 RDKit 计算验证过分子式/分子量/LogP 的奥希替尼 SMILES
CANDIDATES = {
    'osimertinib': 'CN1C=C(C2=CC=CC=C21)C3=NC(=NC=C3)NC4=C(C=C(C(=C4)NC(=O)C=C)N(C)CCN(C)C)OC',
}


def prepare_one(name, smiles):
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'{name}: SMILES 解析失败, 检查 SMILES 是否正确: {smiles}')

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
    print(f'  {name}: 分子式={formula}, MW={mw:.2f}, LogP={logp:.2f}')

    mol = Chem.AddHs(mol)
    embed_result = AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
    if embed_result != 0:
        raise RuntimeError(f'{name}: 3D 构象生成失败(EmbedMolecule 返回非0)')
    AllChem.MMFFOptimizeMolecule(mol)
    mol.SetProp('_Name', name)

    sdf_path = RESULTS_DIR / f'{name}.sdf'
    writer = Chem.SDWriter(str(sdf_path))
    writer.write(mol)
    writer.close()
    print(f'  已写入 3D 构象: {sdf_path}')
    return sdf_path


def convert_to_pdbqt(sdf_path):
    pdbqt_path = sdf_path.with_suffix('.pdbqt')
    result = subprocess.run(
        ['mk_prepare_ligand.py', '-i', str(sdf_path), '-o', str(pdbqt_path)],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f'mk_prepare_ligand.py 转换失败:\nstdout: {result.stdout}\nstderr: {result.stderr}')
    print(f'  已写入 PDBQT: {pdbqt_path}')
    return pdbqt_path


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    print('生成候选分子的 3D 构象 + 理化性质 ...')
    for name, smiles in CANDIDATES.items():
        sdf_path = prepare_one(name, smiles)
        convert_to_pdbqt(sdf_path)

    print('\n已验证: RDKit 能正确解析 SMILES 并生成 3D 构象, meeko CLI 能把 SDF 转成合法 PDBQT'
          '(mk_prepare_ligand.py 返回码为0)。')
    print('未验证: 生成的 PDBQT 里原子类型/可旋转键判定是否完全符合 Vina 的预期——这个要等实际跑一次'
          'Vina 对接、看有没有报错/警告才能确认。')


if __name__ == '__main__':
    main()
