# 引用指南（中文）

如果您在研究中使用了 PhyloDater，请引用本软件。关联论文正式发表后，也请一并引用。

## 软件引用（中文格式）

> 曾子超（Zeng, Zichao）. PhyloDater：多软件并行系统发育定年平台. GitHub 仓库. https://github.com/ZengZichao/phylodater 。DOI: 10.5281/zenodo.23067037 (https://doi.org/10.5281/zenodo.23067037)。ORCID: 0000-0001-6553-970X

## 软件引用（BibTeX 格式）

```bibtex
@software{zeng2026phylodater,
  author = {Zeng, Zichao},
  title = {PhyloDater: A multi-software parallel platform for phylogenetic molecular dating},
  year = {2026},
  url = {https://github.com/ZengZichao/phylodater},
  doi = {10.5281/zenodo.23067037},
  version = {0.1.0},
  note = {ORCID: 0000-0001-6553-970X}
}
```

## 关联论文（预印本 / 稿件）

**目前没有已发表论文，因此这里不给任何条目。** PhyloDater 目前唯一注册的 DOI
是上方引用中的 Zenodo 软件 DOI（`10.5281/zenodo.23067037`）——它标识的是软件
记录，不是论文。论文位置也刻意不放 `10.XXXX/XXXX` 这类占位：GitHub、Zenodo、
Zotero 这些引用消费者分不清占位符与真标识符，一个假 DOI 比"没有 DOI"更糟。

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
