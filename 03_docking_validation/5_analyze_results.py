#!/usr/bin/env python3
"""
解析 Vina 对接日志(4_dock_vina.sh 跑完产出的 osimertinib_docking.log), 按结合能排序输出。

Vina 日志里结果表的真实格式(来自 Vina 官方文档 docking_basic.html, 不是凭记忆编造的):
    mode |   affinity | dist from best mode
        | (kcal/mol) | rmsd l.b.| rmsd u.b.
    -----+------------+----------+----------
    1       -13.23          0          0
    2       -11.29     0.9857      1.681
    ...
affinity 越负代表(Vina打分函数估计的)结合越强, 这是不需要额外确认的 Vina 自身约定。

免责声明: Vina 的 affinity 分数是一个基于经验打分函数的快速估计值, 不是真实的结合自由能测量值,
        不同对接软件/打分函数算出来的数字不能直接跨软件比较。这里的排序只是"在这一次对接、
        这一套参数下, 哪个姿势打分更好", 不构成"这个分子真的能结合"的实验证据——需要湿实验验证。

--demo 模式: 不需要真的跑过 Vina, 用一份格式跟真实日志完全一致、但数值是构造出来的示例数据,
             跑通解析/排序逻辑本身。
"""
import argparse
import csv
import re
import sys
from pathlib import Path

RESULTS_DIR = Path(__file__).resolve().parent / 'results'
ROW_RE = re.compile(r'^\s*(\d+)\s+(-?\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s*$')

DEMO_LOG = """
mode |   affinity | dist from best mode
     | (kcal/mol) | rmsd l.b.| rmsd u.b.
-----+------------+----------+----------
   1       -8.62          0          0
   2       -8.14      1.203      2.415
   3       -7.86      2.031      4.882
   4       -7.55      1.877      3.109
   5       -7.21      3.442      6.203
   6       -6.98      2.654      5.011
   7       -6.71      4.109      7.302
   8       -6.53      3.981      6.877
   9       -6.30      4.556      8.014
""".strip('\n')


def parse_vina_log(log_text):
    poses = []
    for line in log_text.splitlines():
        m = ROW_RE.match(line)
        if m:
            poses.append({
                'mode': int(m.group(1)),
                'affinity_kcal_mol': float(m.group(2)),
                'rmsd_lb': float(m.group(3)),
                'rmsd_ub': float(m.group(4)),
            })
    return poses


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', default=str(RESULTS_DIR / 'osimertinib_docking.log'), help='Vina 日志文件路径')
    parser.add_argument('--output', default=str(RESULTS_DIR / 'docking_ranked.csv'), help='输出排序结果路径')
    parser.add_argument('--demo', action='store_true', help='不读文件, 用内置示例日志跑通解析逻辑(数值不是真实对接结果)')
    args = parser.parse_args()

    if args.demo:
        print('*** --demo 模式: 使用内置示例日志, 不是真实 Vina 对接结果 ***')
        log_text = DEMO_LOG
    else:
        input_path = Path(args.input)
        if not input_path.exists():
            print(f'找不到 {input_path} —— 先跑 4_dock_vina.sh(需要本地装好 vina), '
                  f'或者加 --demo 用示例日志跑通解析逻辑本身。', file=sys.stderr)
            sys.exit(1)
        log_text = input_path.read_text(encoding='utf-8')

    poses = parse_vina_log(log_text)
    if not poses:
        print('没有从日志里解析出任何有效结果行, 检查日志格式是否符合预期(见脚本头部注释里的真实格式样例)', file=sys.stderr)
        sys.exit(1)

    poses.sort(key=lambda p: p['affinity_kcal_mol'])  # 越负越靠前(Vina打分函数约定)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    output_path = Path(args.output)
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['mode', 'affinity_kcal_mol', 'rmsd_lb', 'rmsd_ub'])
        writer.writeheader()
        writer.writerows(poses)

    print(f'解析到 {len(poses)} 个对接姿势, 已按 affinity(越负越好)排序写入: {output_path}')
    best = poses[0]
    print(f"最佳姿势: mode={best['mode']}, affinity={best['affinity_kcal_mol']} kcal/mol")
    print('\n提醒: affinity 是 Vina 打分函数的估计值, 不是实验结合自由能——这个数字只能告诉你"这次对接里'
          '哪个姿势打分更好", 不能单独作为"这个分子真的能结合/能成药"的结论, 需要湿实验验证。')
    print('已验证: 日志解析/排序逻辑在(演示数据 或 真实日志)上能跑通, 排序方向正确(affinity越负排越前)。')
    print('未验证: 真实 Vina 跑出来的最佳 affinity 数值本身有没有意义(比如是不是明显优于随机对接的背景水平)——'
          '这至少需要拿一个已知阳性对照(比如把 osimertinib 对接回它自己真实结合的 6LUD 结构, 看算出来的姿势'
          '是否接近晶体结构里的真实构象)才能判断这套流程本身靠不靠谱, 这一步还没做。')


if __name__ == '__main__':
    main()
