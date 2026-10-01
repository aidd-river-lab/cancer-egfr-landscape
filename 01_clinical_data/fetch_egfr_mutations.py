#!/usr/bin/env python3
"""
拉取 cBioPortal 上真实患者的 EGFR 突变记录。

免责声明: 本脚本产出的数据用于计算假说探索, 不构成任何医疗建议;
        数据来自 cBioPortal 公开的、已去标识化的队列数据, 不含任何可识别患者身份的信息。

已知局限性: 见 ../docs/methodology.md 第1节。

设计原则: 不硬编码猜测的 molecularProfileId / sampleListId 字符串 —— 这两个 ID 的具体拼法
在不同队列间不总是同一套命名规则, 硬编码猜测值一旦猜错会静默拿到空结果或错误结果。
本脚本改为在运行时调用 cBioPortal 自己的元数据接口(getAllMolecularProfilesInStudy /
getAllSampleListsInStudy)去发现"这个 study 下, 类型是 MUTATION_EXTENDED 的 profile"
和"名称含 all sample 的 sample list", 拿到的 ID 保证是当前这个 study 真实存在的。

网络受限说明: 如果当前环境访问不了 cbioportal.org, 本脚本会明确报错退出并说明原因,
             不会用编造的数字顶替。想在没有网络的环境里跑通下游聚合逻辑,
             用 aggregate_mutation_frequency.py --demo。
"""
import argparse
import csv
import json
import sys
from pathlib import Path

EGFR_ENTREZ_GENE_ID = 1956  # NCBI Gene ID, 多来源交叉确认(如需复核: https://www.ncbi.nlm.nih.gov/gene/1956)
DEFAULT_STUDY_ID = 'luad_tcga_pan_can_atlas_2018'
API_DOCS_URL = 'https://www.cbioportal.org/api/v3/api-docs'

RESULTS_DIR = Path(__file__).resolve().parent / 'results'
RAW_OUTPUT_FIELDS = [
    'studyId', 'sampleId', 'patientId', 'proteinChange',
    'mutationType', 'chr', 'startPosition', 'endPosition',
    'referenceAllele', 'variantAllele'
]


def build_client():
    try:
        from bravado.client import SwaggerClient
    except ImportError:
        print('缺少依赖 bravado, 先运行: pip install -r requirements.txt', file=sys.stderr)
        sys.exit(1)
    try:
        return SwaggerClient.from_url(
            API_DOCS_URL,
            config={'validate_requests': False, 'validate_responses': False, 'validate_swagger_spec': False}
        )
    except Exception as error:
        print(f'无法连接 cBioPortal API ({API_DOCS_URL}): {error}', file=sys.stderr)
        print('这是网络环境限制, 不是代码问题。如果要在无网络环境跑通下游聚合逻辑,', file=sys.stderr)
        print('用: python aggregate_mutation_frequency.py --demo', file=sys.stderr)
        sys.exit(1)


def call_json(operation):
    """
    实测发现: bravado 对 cBioPortal 这份 swagger spec 的响应模型反序列化(.result())
    会稳定返回 None(哪怕 HTTP 状态是 200、body 是合法 JSON) —— 这是 bravado-core 处理
    这份 spec 时的已知兼容性问题, 不是我们调用方式的问题(单条查询、列表查询都复现了同样现象)。
    这里绕过模型层, 直接拿底层 HTTP 响应的原始 JSON body 自己解析, 实测稳定可用。
    """
    response = operation.response(fallback_result=None)
    return json.loads(response.metadata.incoming_response.text)


def resolve_molecular_profile_id(cbioportal, study_id):
    """在这个 study 下找类型是 MUTATION_EXTENDED 的 molecular profile, 而不是硬编码猜测 "{study_id}_mutations" """
    profiles = call_json(cbioportal.Molecular_Profiles.getAllMolecularProfilesInStudyUsingGET(studyId=study_id))
    mutation_profiles = [p for p in profiles if p.get('molecularAlterationType') == 'MUTATION_EXTENDED']
    if not mutation_profiles:
        raise RuntimeError(f'study {study_id} 下没有找到 MUTATION_EXTENDED 类型的 molecular profile')
    if len(mutation_profiles) > 1:
        print(f'警告: study {study_id} 下有 {len(mutation_profiles)} 个 MUTATION_EXTENDED profile, 取第一个: '
              f'{mutation_profiles[0]["molecularProfileId"]}', file=sys.stderr)
    return mutation_profiles[0]['molecularProfileId']


def resolve_sample_list_id(cbioportal, study_id):
    """
    找这个 study 下"全部有测序数据的样本"这份 sample list, 而不是硬编码猜测 "{study_id}_all"。
    实测同一个 study 下可能有多份 "*_all" 结尾的 list(比如 luad_tcga_pan_can_atlas_2018 下
    同时有 "..._all" 和 "..._methylation_all"), 优先精确匹配 "{study_id}_all", 避免选到
    别的数据类型专用的子集。
    """
    sample_lists = call_json(cbioportal.Sample_Lists.getAllSampleListsInStudyUsingGET(studyId=study_id))
    exact = f'{study_id}_all'
    if any(s['sampleListId'] == exact for s in sample_lists):
        return exact
    candidates = [s['sampleListId'] for s in sample_lists if s['sampleListId'].endswith('_all')]
    if not candidates:
        candidates = [s['sampleListId'] for s in sample_lists]
    if not candidates:
        raise RuntimeError(f'study {study_id} 下没有找到任何 sample list')
    print(f'警告: 没找到精确匹配 {exact}, 取候选里第一个: {candidates[0]}', file=sys.stderr)
    return candidates[0]


def fetch_egfr_mutations(cbioportal, molecular_profile_id, sample_list_id):
    return call_json(cbioportal.Mutations.getMutationsInMolecularProfileBySampleListIdUsingGET(
        molecularProfileId=molecular_profile_id,
        sampleListId=sample_list_id,
        entrezGeneId=EGFR_ENTREZ_GENE_ID,
        projection='DETAILED'
    ))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study-id', default=DEFAULT_STUDY_ID, help=f'cBioPortal study id (默认 {DEFAULT_STUDY_ID})')
    parser.add_argument('--output', default=None, help='输出 CSV 路径 (默认 results/egfr_mutations_raw_<study-id>.csv)')
    args = parser.parse_args()

    output_path = Path(args.output) if args.output else RESULTS_DIR / f'egfr_mutations_raw_{args.study_id}.csv'
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print(f'连接 cBioPortal API v3: {API_DOCS_URL}')
    cbioportal = build_client()

    print(f'解析 study={args.study_id} 下的 molecular profile / sample list ...')
    molecular_profile_id = resolve_molecular_profile_id(cbioportal, args.study_id)
    sample_list_id = resolve_sample_list_id(cbioportal, args.study_id)
    print(f'  molecularProfileId = {molecular_profile_id}')
    print(f'  sampleListId       = {sample_list_id}')

    print(f'拉取 EGFR(entrezGeneId={EGFR_ENTREZ_GENE_ID}) 突变记录 ...')
    mutations = fetch_egfr_mutations(cbioportal, molecular_profile_id, sample_list_id)
    print(f'  共 {len(mutations)} 条突变记录')

    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=RAW_OUTPUT_FIELDS)
        writer.writeheader()
        skipped = 0
        for m in mutations:
            try:
                writer.writerow({
                    'studyId': args.study_id,
                    'sampleId': m['sampleId'],
                    'patientId': m['patientId'],
                    'proteinChange': m['proteinChange'],
                    'mutationType': m['mutationType'],
                    'chr': m['chr'],
                    'startPosition': m['startPosition'],
                    'endPosition': m['endPosition'],
                    'referenceAllele': m['referenceAllele'],
                    'variantAllele': m['variantAllele'],
                })
            except KeyError:
                skipped += 1

    print(f'已写入: {output_path}')
    if skipped:
        print(f'  跳过 {skipped} 条字段不完整的记录(未丢弃统计, 只是没写入这几条)')
    print('已验证: 能否连上 cBioPortal API + 动态解析出真实存在的 profile/sample list id + 按 EGFR 基因取到突变记录。')
    print('未验证: 这批数据是否需要按 sampleType(原发灶/转移灶)或 mutationStatus 做进一步过滤 —— 当前是队列内全量。')


if __name__ == '__main__':
    main()
