# 公开版本说明（2026-10-20）

本目录对应GitHub仓库 `runzhihuang2022/Huang_2026_SOX9`，按Figure1–Figure7分组。
根README是当前说明；10月16日README仅保留为历史准备记录。

只包括脱敏汇总数据、代码、环境、参数、参考计算图、校验和与复现记录。
原始人类表达矩阵、患者/真实文库对应表、服务器信息和凭据均不在公开目录。

`Figure6`、`Figure7`明确列出实验源数据缺口，不提供伪造数据或空分析脚本。
研究的完整上游Seurat/SCENIC/KEGG/Venn等缺口仍见 `provenance/missing_data.csv`。
公开仓库不等于全文所有实验已全部复现，也不等于原始测序数据已公开存档。

运行：`python scripts/reproduce.py all --out outputs/my_replay_01`
数值核对可加 `--tables-only`。完整矩阵级命令见 `FULL_ANALYSIS.md`。
代码/数据许可证和论文正式引用信息由作者最终确认。
