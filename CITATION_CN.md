> **Language:** [English (canonical)](CITATION.md) · 中文（本文档）

# 引用指南

如果您在研究中使用了 PhyloDater，请引用本软件。关联论文正式发表后，也请一并引用。

## 软件引用（纯文本）

> Zeng, Zichao, & Wang, Yinzhao. PhyloDater: a unified platform integrating
> multiple methods for phylogenetic molecular dating. GitHub repository.
> https://github.com/ZengZichao/phylodater
> ORCID: 0000-0001-6553-970X (Z. Zeng); 0000-0001-8715-1386 (Y. Wang)

以上以 `CITATION.cff` 登记的英文题名为准。

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

## 软件引用（BibTeX 格式）

```bibtex
@software{zeng2026phylodater,
  author = {Zeng, Zichao and Wang, Yinzhao},
  title = {PhyloDater: a unified platform integrating multiple methods for phylogenetic molecular dating},
  year = {2026},
  url = {https://github.com/ZengZichao/phylodater},
  version = {0.1.0},
  doi = {10.5281/zenodo.23067037},
  note = {ORCID: 0000-0001-6553-970X}
}
```

## 关联论文（预印本 / 稿件）

**目前没有已发表论文，因此这里不给任何论文条目。** 上面那个是**软件归档** DOI，
不能替代论文引用——它指向归档代码，不指向同行评议的正文。本文件也刻意不放
`10.XXXX/XXXX` 这类论文占位符：机器可读与人类可读的引用消费者
（GitHub、Zenodo、Zotero）分不清占位符与真标识符，而在一条可被引用的记录里，
一个假 DOI 比"没有 DOI"更糟。

论文正式出现后，请在此处（并且只在此处）补上真实的期刊、年份、卷、期、页码与 DOI，
同时在 [CITATION.cff](CITATION.cff) 里加上对应的 `preferred-citation` 块。

## 机器可读引用

如果不想手工复制以上内容，可以直接使用项目提供的机器可读引用文件
[CITATION.cff](CITATION.cff)。该文件采用标准 CFF 格式，Zenodo、GitHub 的
CITATION.cff 解析器和 Zotero 都能直接读取。

> **哪份引用是权威的？** [CITATION.cff](CITATION.cff) 是机器可读记录；本文件与
> [CITATION.md](CITATION.md) 是同一信息的英文/中文人类可读版本。当版本号、发布
> 日期或作者列表变化时，请保持这三份文件同步。

## 定年引擎的引用

一个定年结果是 PhyloDater 与真正算出它的那个引擎共同产出的，两者都要引用。
各引擎支持（以及在渲染中丢失）哪些校准语义，按方法逐条记在
[docs/TESTING_CN.md](docs/TESTING_CN.md) §8，并在
[docs/TEST_RESULTS_CN.md](docs/TEST_RESULTS_CN.md) 中给出实测。

## 许可证

PhyloDater 基于 MIT 许可证发布，详见 [LICENSE](LICENSE) 文件。可选依赖可能附带
其他许可条款，详见 [THIRD_PARTY_LICENSES_CN.md](THIRD_PARTY_LICENSES_CN.md)。

---

> 📚 本引用指南的英文权威版是 [CITATION.md](CITATION.md)，它是上述信息及
> [CITATION.cff](CITATION.cff) 的人类可读英文版本。
