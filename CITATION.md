> **Language:** English (canonical) · [中文](CITATION_CN.md)

# Citation guide

If PhyloDater contributes to your research, please cite the software. Once the
companion manuscript is published, please cite that too.

## Software citation (plain text)

> Zeng, Zichao. PhyloDater: A multi-software parallel platform for
> phylogenetic molecular dating. GitHub repository.
> https://github.com/ZengZichao/phylodater
> DOI: 10.5281/zenodo.23067037 (https://doi.org/10.5281/zenodo.23067037)
> ORCID: 0000-0001-6553-970X

## Software citation (BibTeX)

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

## Companion paper (preprint / manuscript)

**There is none yet, so none is quoted here.** The only DOI registered for
PhyloDater is the Zenodo software DOI (`10.5281/zenodo.23067037`) in the
citation above — it identifies the software record, not a paper. The article
slot deliberately carries no `10.XXXX/XXXX` placeholder: machine-readable and
human-readable citation consumers (GitHub, Zenodo, Zotero) cannot tell a
placeholder from a real identifier, and a fake DOI in a citable record is worse
than no DOI.

Once the manuscript appears, add here — and only here — the real journal, year,
volume, issue, pages and DOI, and add the matching `preferred-citation` block to
[CITATION.cff](CITATION.cff).

## Machine-readable citation

Rather than copying the blocks above by hand, use the machine-readable
[CITATION.cff](CITATION.cff): it is standard CFF, readable directly by Zenodo,
GitHub's CITATION.cff parser and Zotero.

> **Which citation is authoritative?** [CITATION.cff](CITATION.cff) is the
> machine-readable record; this file and
> [CITATION_CN.md](CITATION_CN.md) are the human-readable English/Chinese
> versions of the same information. Keep all three in sync when the version,
> release date or author list changes.

## Citing the dating engines

A dated result is produced jointly by PhyloDater and the engine that estimated it.
Cite both: the calibration semantics each engine supports (and loses) are documented
method by method in [docs/TESTING.md](docs/TESTING.md) §8 and measured in
[docs/TEST_RESULTS.md](docs/TEST_RESULTS.md).

## Licence

PhyloDater is released under the MIT License — see [LICENSE](LICENSE). Optional
dependencies may carry other terms; see
[THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md).
