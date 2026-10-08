# 引用指南（中文）

如果您在研究中使用了 PhyloDater，请引用本软件。关联论文正式发表后，也请一并引用。

## 软件引用（中文格式）

> 曾子超（Zeng, Zichao）. PhyloDater：多软件并行系统发育定年平台. GitHub 仓库. https://github.com/ZengZichao/phylodater 。ORCID: 0000-0001-6553-970X

## 软件引用（BibTeX 格式）

```bibtex
@software{zeng2026phylodater,
  author = {Zeng, Zichao},
  title = {PhyloDater: A multi-software parallel platform for phylogenetic molecular dating},
  year = {2026},
  url = {https://github.com/ZengZichao/phylodater},
  version = {0.1.0},
  doi = {10.5281/zenodo.23067037},
  note = {ORCID: 0000-0001-6553-970X}
}
```

## 归档软件引用（Zenodo，版本 DOI）

v0.1.0 发布版已在 Zenodo 归档并注册了可解析的 DOI。需要"锁定到具体版本"的引用时，
请用这一条：

> Zeng, Zichao. (2026). PhyloDater (Version 0.1.0) [Software]. Zenodo.
> https://doi.org/10.5281/zenodo.23067037

| 标识符 | 值 | 指向 |
|---|---|---|
| 版本 DOI（本次发布） | `10.5281/zenodo.23067037` | 已归档的 v0.1.0 |
| 概念 DOI（任意版本） | `10.5281/zenodo.23066788` | PhyloDater 全部版本 |
| 源码标签 | https://github.com/ZengZichao/phylodater/releases/tag/v0.1.0 | GitHub release |

两个 DOI 均从 Zenodo 记录本身回读核验（`https://zenodo.org/api/records/23067037`，
2026-10-08 获取）：许可证 MIT、版本 `v0.1.0`、关键词表与
[CITATION.cff](CITATION.cff) 一致。

## 关联论文（预印本 / 稿件）

**目前没有已发表论文，因此这里不给任何论文条目。** 上面那个是**软件归档** DOI，
不能替代论文引用——它指向归档代码，不指向同行评议的正文。本文件也刻意不放
`10.XXXX/XXXX` 这类论文占位符：GitHub、Zenodo、Zotero 这些引用消费者
分不清占位符与真标识符，一个假 DOI 比"没有 DOI"更糟。

论文正式出现后，请在此处（并且只在此处）补上真实的期刊、年份、卷期页码与 DOI，
同时在 [CITATION.cff](CITATION.cff) 里加上对应的 `preferred-citation` 块。

## 机器可读引用

如果不想手工复制以上内容，可以直接使用项目提供的机器可读引用文件 [CITATION.cff](CITATION.cff)。该文件采用标准 CFF 格式，Zenodo、GitHub 的 CITATION.cff 解析器和 Zotero 都能直接读取。

## 定年引擎的引用

一个定年结果是 PhyloDater 与真正算出它的那个引擎共同产出的，两者都要引用。
各引擎支持（以及在渲染中丢失）哪些校准语义，按方法逐条记在
[docs/TESTING_CN.md](docs/TESTING_CN.md) §8，并在
[docs/TEST_RESULTS_CN.md](docs/TEST_RESULTS_CN.md) 中给出实测。

## 许可证

PhyloDater 基于 MIT 许可证发布，详见 [LICENSE](LICENSE) 文件。

---

> 📚 本文档的英文权威版见 [README.md](README.md) 中的《Citation》一节与 [CITATION.cff](CITATION.cff)。
> The canonical English citation information is in [README.md](README.md) (Citation section) and [CITATION.cff](CITATION.cff).
