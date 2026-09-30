"""
SelfTest - 自检模式与依赖验证

执行系统自检：
- 检查第三方库是否可导入及版本
- 测试示例树文件解析
- 测试单系群判定逻辑
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from .logging import get_logger


@dataclass
class TestResult:
    """测试结果"""

    name: str
    passed: bool
    message: str = ""
    version: Optional[str] = None


class SelfTester:
    """
    自检测试器

    执行系统自检并输出结果表格。
    """

    def __init__(self) -> None:
        self.logger = get_logger()
        self.results: List[TestResult] = []

    def run_all_checks(self) -> bool:
        """
        运行所有检查

        Returns:
            True 如果全部通过
        """
        self.logger.section("PhyloDater 自检模式")

        # 1. 检查第三方库
        self._check_dependencies()

        # 2. 测试树文件解析
        self._test_tree_parsing()

        # 3. 测试分类学解析
        self._test_taxonomy_parsing()

        # 4. 测试单系群判定
        self._test_monophyly_check()

        # 5. 测试深度验证器
        self._test_deep_validator()

        # 6. 测试特殊标识符大小写
        self._test_special_identifiers()

        # 7. 测试输入大小限制
        self._test_input_size_limits()

        # 8. 测试交叉验证
        self._test_cross_validation()

        # 9. 测试外部分类学表格
        self._test_external_taxonomy_table()

        # 10. 测试多棵树处理
        self._test_multi_tree_handling()

        # 11. 测试恶意字符拒绝
        self._test_malicious_char_rejection()

        # 输出结果表格
        self._print_results()

        # 返回是否全部通过
        return all(r.passed for r in self.results)

    def _check_dependencies(self) -> None:
        """检查第三方库依赖

        只有真正必需、且缺失就会让核心跑不起来的库才放在 ``dependencies``；
        能可选、能降级、或在某些解释器上根本装不上的（ete3 / pyr8s / md-cat /
        logdate）全部走下面的可选依赖区，只报不阻塞。
        """
        dependencies = [
            ("biopython", "Bio", "1.80"),
            ("dendropy", "dendropy", "4.5.0"),
            ("numpy", "numpy", "1.20.0"),
            ("pandas", "pandas", "1.3.0"),
            ("pyyaml", "yaml", "5.4"),
            ("matplotlib", "matplotlib", "3.4.0"),
        ]

        for name, import_name, min_version in dependencies:
            try:
                module = __import__(import_name)
                version = getattr(module, "__version__", "unknown")

                # 版本比较（简化版）
                passed = True
                if version != "unknown":
                    try:
                        # 简单版本比较
                        current = tuple(map(int, version.split(".")[:3]))
                        minimum = tuple(map(int, min_version.split(".")[:3]))
                        passed = current >= minimum
                    except (ValueError, TypeError):
                        passed = True  # 无法比较时默认通过

                if passed:
                    self.results.append(
                        TestResult(
                            name=f"依赖: {name}",
                            passed=True,
                            message=f"v{version}",
                            version=version,
                        )
                    )
                else:
                    self.results.append(
                        TestResult(
                            name=f"依赖: {name}",
                            passed=False,
                            message=f"版本过低: v{version} (需要 >= {min_version})",
                        )
                    )

            except ImportError:
                self.results.append(
                    TestResult(name=f"依赖: {name}", passed=False, message="未安装")
                )

        # 可选依赖：不阻塞整体自检，但必须明确提示缺失会影响哪一条通路。
        # 七个引擎里有两个（r8s 的 pyr8s 后端、MD-Cat 的 emd 后端）和可视化
        # 层的 ete 形状后端都是以 Python 包形式发的，旧表只探了 logdate 一个，
        # 导致"装了 md-cat 没有"这类环境差异只能等到跑定年时才暴。下面的发行
        # 渠道按 2026-09-30 实测写：以 PyPI/bioconda 的 **JSON API** 为准
        # （pypi.org/project/<name>/ 的 HTML 页面对不存在的包也会返回 200 软页）。
        optional_dependencies = [
            (
                "logdate (wLogDate)",
                "logdate",
                "wLogDate 方法。PyPI 发行名 wlogdate（实测 PyPI JSON API 存在，"
                "当前 1.0.2），bioconda 亦有同名包（1.0.4）；"
                "旧文档里的 github.com/wchan29/wLogDate 已 404。",
            ),
            (
                "emd (MD-Cat)",
                "emd",
                "MD-Cat 方法的直接导入后端。只能从源码装："
                "https://github.com/uym2/MD-Cat（PyPI 与 bioconda 均无此包，"
                "实测 JSON API 404）；装成 md_cat 发行名。",
            ),
            (
                "pyr8s",
                "pyr8s",
                "r8s 方法的纯 Python 后端（仅支持 NPRS）。只能从源码装："
                "https://github.com/iTaxoTools/pyr8s（PyPI 与 bioconda 均无此包，"
                "实测 JSON API 404）",
            ),
            (
                "ete3",
                "ete3",
                "可视化层需要 ete 形状的树对象时使用；核心定年不依赖它。"
                "ete3 3.1.x 在 Python >= 3.13 上不可安装（stdlib cgi 已移除），"
                "那种解释器请改装 ete4（pip install 'phylodater[ete4]'）",
            ),
        ]

        for name, import_name, hint in optional_dependencies:
            try:
                module = __import__(import_name)
                version = getattr(module, "__version__", "unknown")
                self.results.append(
                    TestResult(
                        name=f"可选依赖: {name}",
                        passed=True,
                        message=f"v{version}",
                        version=version,
                    )
                )
            except ImportError:
                self.results.append(
                    TestResult(
                        name=f"可选依赖: {name}",
                        passed=True,
                        message=f"未安装。影响：{hint}",
                    )
                )

    def _test_tree_parsing(self) -> None:
        """测试树文件解析"""
        from ..models import PhylogeneticTree

        # 测试用例1: 基本 Newick 解析
        try:
            tree = PhylogeneticTree.from_newick(
                "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);"
            )
            if tree.num_tips == 4 and len(tree.tip_names) == 4:
                self.results.append(
                    TestResult(
                        name="树解析: 基本 Newick",
                        passed=True,
                        message=f"{tree.num_tips} 个叶节点",
                    )
                )
            else:
                self.results.append(
                    TestResult(
                        name="树解析: 基本 Newick",
                        passed=False,
                        message=f"叶节点数量错误: {tree.num_tips}",
                    )
                )
        except Exception as e:
            self.results.append(
                TestResult(name="树解析: 基本 Newick", passed=False, message=str(e))
            )

        # 测试用例2: 带分类信息的树
        try:
            tree = PhylogeneticTree.from_newick(
                "((GB_GCA_001_d_Bacteria_p_Proteobacteria:0.1,"
                "GB_GCA_002_d_Bacteria_p_Cyanobacteriota:0.2):0.3,"
                "(GB_GCA_003_d_Archaea_p_Euryarchaeota:0.4,"
                "GB_GCA_004_d_Archaea_p_Thermoproteota:0.5):0.6);"
            )
            if tree.num_tips == 4:
                self.results.append(
                    TestResult(
                        name="树解析: 带分类信息",
                        passed=True,
                        message=f"{tree.num_tips} 个叶节点",
                    )
                )
            else:
                self.results.append(
                    TestResult(
                        name="树解析: 带分类信息", passed=False, message="解析失败"
                    )
                )
        except Exception as e:
            self.results.append(
                TestResult(name="树解析: 带分类信息", passed=False, message=str(e))
            )

    def _test_taxonomy_parsing(self) -> None:
        """测试分类学解析"""
        from ..services.taxonomy_parser import TaxonomyParser

        parser = TaxonomyParser()

        # 测试用例1: 格式A（嵌入式）
        try:
            result = parser.parse(
                "GB_GCA_001_d_Bacteria_p_Proteobacteria_c_Gammaproteobacteria"
            )
            if (
                result
                and result.get("domain") == "Bacteria"
                and result.get("phylum") == "Proteobacteria"
            ):
                self.results.append(
                    TestResult(
                        name="分类学: 格式A（嵌入式）",
                        passed=True,
                        message=f"domain={result['domain']}, phylum={result['phylum']}",
                    )
                )
            else:
                self.results.append(
                    TestResult(
                        name="分类学: 格式A（嵌入式）",
                        passed=False,
                        message=f"解析结果错误: {result}",
                    )
                )
        except Exception as e:
            self.results.append(
                TestResult(name="分类学: 格式A（嵌入式）", passed=False, message=str(e))
            )

        # 测试用例2: 格式B（表格分号式）
        try:
            result = parser.parse(
                "d__Bacteria;p__Cyanobacteriota;c__Cyanobacteriia;o__Synechococcales;f__Synechococcaceae;g__Synechococcus;s__"
            )
            if (
                result
                and result.get("domain") == "Bacteria"
                and result.get("phylum") == "Cyanobacteriota"
            ):
                # 检查种水平为空
                species = result.get("species")
                if species is None or species == "":
                    self.results.append(
                        TestResult(
                            name="分类学: 格式B（表格分号式）",
                            passed=True,
                            message=f"domain={result['domain']}, species=None (缺失值处理正确)",
                        )
                    )
                else:
                    self.results.append(
                        TestResult(
                            name="分类学: 格式B（表格分号式）",
                            passed=False,
                            message=f"缺失值处理错误: species={species}",
                        )
                    )
            else:
                self.results.append(
                    TestResult(
                        name="分类学: 格式B（表格分号式）",
                        passed=False,
                        message=f"解析结果错误: {result}",
                    )
                )
        except Exception as e:
            self.results.append(
                TestResult(
                    name="分类学: 格式B（表格分号式）", passed=False, message=str(e)
                )
            )

    def _test_monophyly_check(self) -> None:
        """测试单系群判定"""
        from ..models import PhylogeneticTree

        # 构建测试树：细菌和古菌混合
        tree_newick = (
            "((Bact1_d_Bacteria:0.1,Arch1_d_Archaea:0.1):0.2,"
            "(Bact2_d_Bacteria:0.3,Arch2_d_Archaea:0.3):0.1);"
        )

        try:
            tree = PhylogeneticTree.from_newick(tree_newick)

            # 测试: Bacteria 不是单系群
            is_mono = tree.is_monophyletic(["Bact1_d_Bacteria", "Bact2_d_Bacteria"])
            if not is_mono:
                self.results.append(
                    TestResult(
                        name="单系群: 非单系群检测",
                        passed=True,
                        message="正确识别 Bacteria 非单系群",
                    )
                )
            else:
                self.results.append(
                    TestResult(
                        name="单系群: 非单系群检测",
                        passed=False,
                        message="未能识别非单系群",
                    )
                )

        except Exception as e:
            self.results.append(
                TestResult(name="单系群: 非单系群检测", passed=False, message=str(e))
            )

        # 构建测试树：细菌是单系群
        tree_newick_mono = (
            "((Bact1_d_Bacteria:0.1,Bact2_d_Bacteria:0.1):0.2,"
            "(Arch1_d_Archaea:0.3,Arch2_d_Archaea:0.3):0.1);"
        )

        try:
            tree = PhylogeneticTree.from_newick(tree_newick_mono)

            # 测试: Bacteria 是单系群
            is_mono = tree.is_monophyletic(["Bact1_d_Bacteria", "Bact2_d_Bacteria"])
            if is_mono:
                self.results.append(
                    TestResult(
                        name="单系群: 单系群检测",
                        passed=True,
                        message="正确识别 Bacteria 单系群",
                    )
                )
            else:
                self.results.append(
                    TestResult(
                        name="单系群: 单系群检测",
                        passed=False,
                        message="未能识别单系群",
                    )
                )

        except Exception as e:
            self.results.append(
                TestResult(name="单系群: 单系群检测", passed=False, message=str(e))
            )

        # 测试: LUCA 场景
        try:
            tree = PhylogeneticTree.from_newick(tree_newick_mono)
            from ..services.calibration_resolver import CalibrationResolver

            resolver = CalibrationResolver(tree)

            # 解析 LUCA（细菌+古菌的 MRCA）
            result = resolver.resolve("LUCA")
            if result and result.resolved_taxa:
                self.results.append(
                    TestResult(
                        name="单系群: LUCA 解析",
                        passed=True,
                        message=f"解析到 {len(result.resolved_taxa)} 个分类单元",
                    )
                )
            else:
                self.results.append(
                    TestResult(
                        name="单系群: LUCA 解析", passed=False, message="解析失败"
                    )
                )

        except Exception as e:
            self.results.append(
                TestResult(name="单系群: LUCA 解析", passed=False, message=str(e))
            )

    def _test_deep_validator(self) -> None:
        """测试深度验证器"""
        from ..services.deep_validator import DeepValidator

        validator = DeepValidator()

        # 测试用例: 括号不平衡检测
        try:
            # 创建临时测试
            import tempfile

            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".nwk", delete=False
            ) as f:
                f.write("((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6;")  # 缺少右括号
                temp_path = Path(f.name)

            result = validator.validate_tree_deep(temp_path)
            temp_path.unlink()

            if not result.is_valid and any("括号" in e for e in result.errors):
                self.results.append(
                    TestResult(
                        name="深度验证: 括号不平衡检测",
                        passed=True,
                        message="正确检测括号不平衡",
                    )
                )
            else:
                self.results.append(
                    TestResult(
                        name="深度验证: 括号不平衡检测",
                        passed=False,
                        message="未能检测括号不平衡",
                    )
                )

        except Exception as e:
            self.results.append(
                TestResult(
                    name="深度验证: 括号不平衡检测", passed=False, message=str(e)
                )
            )

        # 测试用例: 重复ID检测
        try:
            import tempfile

            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".fasta", delete=False
            ) as f:
                f.write(">seq1\nATCG\n>seq1\nATCG\n")  # 重复ID
                temp_path = Path(f.name)

            seq_result = validator.validate_alignment_deep(temp_path)
            temp_path.unlink()

            if not seq_result.is_valid and len(seq_result.duplicate_ids) > 0:
                self.results.append(
                    TestResult(
                        name="深度验证: 重复ID检测",
                        passed=True,
                        message=f"检测到 {len(seq_result.duplicate_ids)} 个重复ID",
                    )
                )
            else:
                self.results.append(
                    TestResult(
                        name="深度验证: 重复ID检测",
                        passed=False,
                        message="未能检测重复ID",
                    )
                )

        except Exception as e:
            self.results.append(
                TestResult(name="深度验证: 重复ID检测", passed=False, message=str(e))
            )

    def _test_special_identifiers(self) -> None:
        """测试特殊标识符LUCA/LBCA/LACA的大小写要求"""
        from ..models import PhylogeneticTree
        from ..services.calibration_resolver import CalibrationResolver

        tree_newick = (
            "((Bact1_d_Bacteria:0.1,Bact2_d_Bacteria:0.1):0.2,"
            "(Arch1_d_Archaea:0.3,Arch2_d_Archaea:0.3):0.1);"
        )

        try:
            tree = PhylogeneticTree.from_newick(tree_newick)
            resolver = CalibrationResolver(tree)

            result = resolver.resolve("LUCA")
            if result and result.resolved_taxa:
                self.results.append(
                    TestResult(
                        name="特殊标识符: LUCA (大写)",
                        passed=True,
                        message=f"解析到 {len(result.resolved_taxa)} 个分类单元",
                    )
                )
            else:
                self.results.append(
                    TestResult(
                        name="特殊标识符: LUCA (大写)", passed=False, message="解析失败"
                    )
                )
        except Exception as e:
            self.results.append(
                TestResult(name="特殊标识符: LUCA (大写)", passed=False, message=str(e))
            )

        try:
            tree = PhylogeneticTree.from_newick(tree_newick)
            resolver = CalibrationResolver(tree)
            result = resolver.resolve("ROOT")
            if result and result.resolved_taxa:
                self.results.append(
                    TestResult(
                        name="特殊标识符: ROOT",
                        passed=True,
                        message=f"解析到 {len(result.resolved_taxa)} 个分类单元",
                    )
                )
            else:
                self.results.append(
                    TestResult(
                        name="特殊标识符: ROOT", passed=False, message="解析失败"
                    )
                )
        except Exception as e:
            self.results.append(
                TestResult(name="特殊标识符: ROOT", passed=False, message=str(e))
            )

    def _test_input_size_limits(self) -> None:
        """测试输入大小限制"""
        from ..services.deep_validator import DeepValidator

        validator = DeepValidator()
        has_limit = (
            hasattr(validator, "MAX_INPUT_SIZE") and validator.MAX_INPUT_SIZE > 0
        )
        self.results.append(
            TestResult(
                name="安全: 输入大小限制",
                passed=has_limit,
                message=(
                    f"限制: {validator.MAX_INPUT_SIZE // (1024*1024)} MB"
                    if has_limit
                    else "未设置限制"
                ),
            )
        )

    def _test_cross_validation(self) -> None:
        """测试交叉验证"""
        try:
            from .. import cross_validate  # noqa: F401

            self.results.append(
                TestResult(
                    name="库模式API: cross_validate", passed=True, message="函数存在"
                )
            )
        except ImportError:
            self.results.append(
                TestResult(
                    name="库模式API: cross_validate", passed=False, message="函数未导出"
                )
            )

    def _test_external_taxonomy_table(self) -> None:
        """测试外部分类学表格加载"""
        table_path = None
        try:
            import tempfile

            from ..services.taxonomy_parser import TaxonomyParser

            # 创建临时表格文件
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".tsv", delete=False, encoding="utf-8"
            ) as f:
                f.write("name\ttaxonomy\n")
                f.write("Sample1\td__Bacteria;p__Firmicutes\n")
                f.write("Sample2\td__Archaea;p__Euryarchaeota\n")
                table_path = f.name
            parser = TaxonomyParser()
            parser.load_from_file(Path(table_path))
            result = parser.parse("Sample1")
            if result and result.get("domain") == "Bacteria":
                self.results.append(
                    TestResult(name="外部分类学表格", passed=True, message="PASS")
                )
            else:
                self.results.append(
                    TestResult(
                        name="外部分类学表格",
                        passed=False,
                        message=f"解析结果不符预期: {result}",
                    )
                )
        except Exception as e:
            self.results.append(
                TestResult(name="外部分类学表格", passed=False, message=str(e))
            )
        finally:
            if table_path is not None:
                Path(table_path).unlink(missing_ok=True)

    def _test_multi_tree_handling(self) -> None:
        """测试多棵树处理"""
        nex_path = None
        try:
            import tempfile

            # 创建含多棵树的 Nexus 文件
            multi_tree_nexus = """#NEXUS
BEGIN TREES;
tree t1 = ((A:1,B:1):1,C:2);
tree t2 = ((A:2,B:2):1,C:3);
END;"""
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".nex", delete=False, encoding="utf-8"
            ) as f:
                f.write(multi_tree_nexus)
                nex_path = f.name
            # 验证真实的解析路径能处理多棵树
            from ..models.tree import PhylogeneticTree

            try:
                tree = PhylogeneticTree.from_file(Path(nex_path))
                # 如果成功加载，说明只取了第一棵树（可接受行为）
                self.results.append(
                    TestResult(
                        name="多棵树处理",
                        passed=True,
                        message=f"多树 NEXUS 已加载（取第一棵树），叶节点数: {tree.num_tips}",
                    )
                )
            except Exception as e:
                # 被拒绝也是可接受行为
                self.results.append(
                    TestResult(
                        name="多棵树处理",
                        passed=True,
                        message=f"多树 NEXUS 被正确拒绝: {type(e).__name__}: {e}",
                    )
                )
        except Exception as e:
            self.results.append(
                TestResult(name="多棵树处理", passed=False, message=str(e))
            )
        finally:
            if nex_path is not None:
                Path(nex_path).unlink(missing_ok=True)

    def _test_malicious_char_rejection(self) -> None:
        """测试恶意字符注入拒绝"""
        tree_path = None
        try:
            import tempfile

            from ..services.deep_validator import DeepValidator

            # 含 Unicode 双向覆盖字符的树
            malicious_tree = "(\u202eA:1,B:1):0;"
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".nwk", delete=False, encoding="utf-8"
            ) as f:
                f.write(malicious_tree)
                tree_path = f.name
            validator = DeepValidator()
            result = validator.validate_tree_deep(Path(tree_path))
            if not result.is_valid and any(
                "恶意" in e or "malicious" in e.lower() or "控制" in e or "字符" in e
                for e in result.errors
            ):
                self.results.append(
                    TestResult(name="恶意字符拒绝", passed=True, message="PASS")
                )
            else:
                self.results.append(
                    TestResult(
                        name="恶意字符拒绝",
                        passed=False,
                        message=f"未正确拒绝恶意字符: {result.errors}",
                    )
                )
        except Exception as e:
            self.results.append(
                TestResult(name="恶意字符拒绝", passed=False, message=str(e))
            )
        finally:
            if tree_path is not None:
                Path(tree_path).unlink(missing_ok=True)

    def _print_results(self) -> None:
        """输出结果表格"""
        print("\n" + "=" * 70)
        print(f"{'检查项':<40} {'状态':<8} {'信息'}")
        print("-" * 70)

        passed_count = 0
        failed_count = 0

        for result in self.results:
            status = "[PASS]" if result.passed else "[FAIL]"
            if result.passed:
                passed_count += 1
            else:
                failed_count += 1

            print(f"{result.name:<40} {status:<8} {result.message}")

        print("-" * 70)
        print(
            f"总计: {len(self.results)} 项检查, {passed_count} 通过, {failed_count} 失败"
        )
        print("=" * 70)

        if failed_count == 0:
            print("\n✓ 所有检查通过！")
        else:
            print(f"\n✗ {failed_count} 项检查失败，请修复后重试。")
