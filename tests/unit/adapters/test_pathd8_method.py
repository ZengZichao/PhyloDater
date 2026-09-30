"""
PATHd8Method 适配器单元测试

测试 PATHd8 定年方法的各项功能

说明：本文件里以 ``REAL_*`` 命名的夹具全部是随项目提供的 PATHd8 二进制
（``PhyloDater-参考软件/原始代码-PATHd8/PATHd8``，Mach-O arm64）在
本机上真实跑出来的产物，仅为写入源码去掉了行尾空格，未做任何改写。
之所以坚持用真实产物做夹具，是因为审阅报告 B-11 的病灶正是
"测试用的假想输出与上游实际输出不符"：上游报告约束改写的那句 printf
被注释掉了（meta.c:36-43），所以任何基于 "conflict"/"adjusted"
关键词的测试都会自我实现。
"""

import shutil
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from phylodater.adapters.pathd8_method import PATHd8Method
from phylodater.core.exceptions import (
    CalibrationError,
    ExecutionError,
    ResultParsingError,
)
from phylodater.infrastructure.configuration import CommonConfig, PATHd8Config
from phylodater.models import (
    CalibrationPoint,
    FixedAgeConstraint,
    MaximumAgeConstraint,
    PhylogeneticTree,
    SoftLowerBoundConstraint,
    UniformAgeConstraint,
)

NEWICK = (
    "((human:0.03,(chimp:0.01,bonobo:0.01):0.02):0.05,"
    "(gorilla:0.06,(orangutan:0.04,sumatran:0.04):0.02):0.02);"
)
TIPS = ["human", "chimp", "bonobo", "gorilla", "orangutan", "sumatran"]

# ---------------------------------------------------------------------------
# 真实产物夹具
# ---------------------------------------------------------------------------

# 一次完全正常、约束全部生效的运行（返回码 1，"Calculation finished."）
REAL_CLEAN_OUTPUT = """

************************************************************************************************
*  d 8   C A L C U L A T I O N                                                                 *
************************************************************************************************

Number of informative fixnodes:   1
Number of informative minnodes:   1
Number of informative maxnodes:   1

d8 tree    : ((human:6.000000,(chimp:2.000000,bonobo:2.000000):4.000000)HumanChimp:4.000000,(gorilla:7.500000,(orangutan:5.000000,sumatran:5.000000):2.500000):2.500000)GreatApes:0;

Ancestor of          Ancestor of          Name                        Age  #Terminals              MPL         Rate *      minage      maxage
human                sumatran             GreatApes                10.000           6           40.000       4.000000         8.0        10.0
human                bonobo               HumanChimp                6.000           3           15.000       2.500000         6.0         6.0
chimp                bonobo               -                         2.000           2            5.000       2.500000           -           -
gorilla              sumatran             -                         7.500           3           30.000       4.000000           -           -
orangutan            sumatran             -                         5.000           2           20.000       4.000000           -           -


  *) Rate = MPL / Age


************************************************************************************************
*  M P L  C A L C U L A T I O N                                                                *
************************************************************************************************

Clock test confidence: 0.950000
Clock tests          : 5   (one for each node)
Accepted             : 5
Rejected             : 0

MPL tree   : ((human:15.000000,(chimp:5.000000,bonobo:5.000000):10.000000)HumanChimp:25.000000,(gorilla:30.000000,(orangutan:20.000000,sumatran:20.000000):10.000000):10.000000)GreatApes:0;

Ancestor of          Ancestor of          Name                            MPL             #Terminals      Clock test: Acc/Rej
human                sumatran             GreatApes                40.000 +/- 7.226                6      Acc
human                bonobo               HumanChimp               15.000 +/- 5.263                3      Acc
chimp                bonobo               -                         5.000 +/- 3.097                2      Acc
gorilla              sumatran             -                        30.000 +/- 6.847                3      Acc
orangutan            sumatran             -                        20.000 +/- 6.193                2      Acc


************************************************************************************************
*  E N D  C A L C U L A T I O N                                                                *
************************************************************************************************
"""

# 上游把冗余（不相容）的 min 约束静默丢掉后的真实产物：写出去 2 条 minage，
# 只报回 1 个 informative minnode，RootMin 那一行两列都是 "-"。
# 注意全文既没有 "conflict" 也没有 "adjusted"——这正是 B-11 的死分支。
REAL_SILENT_REWRITE_OUTPUT = """

************************************************************************************************
*  d 8   C A L C U L A T I O N                                                                 *
************************************************************************************************

Number of informative fixnodes:   1
Number of informative minnodes:   1
Number of informative maxnodes:   0

d8 tree    : ((human:6.000000,(chimp:2.000000,bonobo:2.000000):4.000000)HumanChimp:10.000000,(gorilla:12.000000,(orangutan:8.000000,sumatran:8.000000):4.000000)GorMin:4.000000)RootMin:0;

Ancestor of          Ancestor of          Name                        Age  #Terminals              MPL         Rate *      minage      maxage
human                sumatran             RootMin                  16.000           6           40.000       2.500000           -           -
human                bonobo               HumanChimp                6.000           3           15.000       2.500000         6.0         6.0
chimp                bonobo               -                         2.000           2            5.000       2.500000           -           -
gorilla              sumatran             GorMin                   12.000           3           30.000       2.500000        10.0           -
orangutan            sumatran             -                         8.000           2           20.000       2.500000           -           -


  *) Rate = MPL / Age


************************************************************************************************
*  M P L  C A L C U L A T I O N                                                                *
************************************************************************************************

Clock test confidence: 0.950000
Clock tests          : 5   (one for each node)
Accepted             : 5
Rejected             : 0

MPL tree   : ((human:15.000000,(chimp:5.000000,bonobo:5.000000):10.000000)HumanChimp:25.000000,(gorilla:30.000000,(orangutan:20.000000,sumatran:20.000000):10.000000)GorMin:10.000000)RootMin:0;

Ancestor of          Ancestor of          Name                            MPL             #Terminals      Clock test: Acc/Rej
human                sumatran             RootMin                  40.000 +/- 7.226                6      Acc
human                bonobo               HumanChimp               15.000 +/- 5.263                3      Acc
chimp                bonobo               -                         5.000 +/- 3.097                2      Acc
gorilla              sumatran             GorMin                   30.000 +/- 6.847                3      Acc
orangutan            sumatran             -                        20.000 +/- 6.193                2      Acc


************************************************************************************************
*  E N D  C A L C U L A T I O N                                                                *
************************************************************************************************
"""

# NUM_FIX==0 的真实产物：上游清掉用户全部 min/max、把根钉在 1.0，
# 但照常产出可解析的 d8 tree（B-12 的灾难形态）。
REAL_NO_FIXNODE_OUTPUT = """

************************************************************************************************
*  d 8   C A L C U L A T I O N                                                                 *
************************************************************************************************



No fixnodes were defined in

\t\t        PATHd8.infile

Any user given constraint is ignored and root age is fixed to 1

Number of informative fixnodes:   1
Number of informative minnodes:   0
Number of informative maxnodes:   0

d8 tree    : ((human:0.375000,(chimp:0.125000,bonobo:0.125000):0.250000):0.625000,(gorilla:0.750000,(orangutan:0.500000,sumatran:0.500000):0.250000):0.250000)GreatApes:0;

Ancestor of          Ancestor of          Name                        Age  #Terminals              MPL         Rate *      minage      maxage
human                sumatran             GreatApes                 1.000           6           40.000      40.000000         1.0         1.0
human                bonobo               -                         0.375           3           15.000      40.000000           -           -
chimp                bonobo               -                         0.125           2            5.000      40.000000           -           -
gorilla              sumatran             -                         0.750           3           30.000      40.000000           -           -
orangutan            sumatran             -                         0.500           2           20.000      40.000000           -           -


  *) Rate = MPL / Age
"""

# 同一节点上两个互相矛盾的 fixage：上游 printf 到标准输出后
# error("fill_fixage_data", ERR_FORMAT) → **exit(0)**，结果文件 0 字节。
REAL_FATAL_STDOUT = """


The node with mrca human and gorilla is defined to be
1. A fixnode with fixage 5.000000
2. A fixnode with fixage 20.000000
This is a contradiction.
"""


def _make_tree():
    tree = Mock(spec=PhylogeneticTree)
    tree.newick = NEWICK
    tree.tip_names = list(TIPS)
    tree.num_tips = len(TIPS)
    return tree


def _make_cal(
    name,
    constraint,
    pair=("human", "chimp"),
    is_root=False,
):
    cal = Mock(spec=CalibrationPoint)
    cal.name = name
    cal.age_constraint = constraint
    cal.mrca_leaf_pair = pair
    cal.is_root_node = is_root
    cal.resolved_taxa = list(pair) if pair else list(TIPS)
    return cal


class TestPATHd8Method:
    """PATHd8Method 测试类"""

    @pytest.fixture
    def temp_dir(self):
        """创建临时目录"""
        temp_path = Path(tempfile.mkdtemp(prefix="pathd8_test_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def pathd8_config(self):
        """创建 PATHd8 配置"""
        return PATHd8Config()

    @pytest.fixture
    def common_config(self):
        """创建通用配置"""
        return CommonConfig(nthreads=4)

    @pytest.fixture
    def mock_tree(self):
        """创建模拟树"""
        tree = Mock(spec=PhylogeneticTree)
        tree.newick = "((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);"
        tree.tip_names = ["human", "chimp", "mouse", "rat"]
        tree.num_tips = 4
        tree.without_internal_labels = Mock(return_value=tree)
        tree.write = Mock()
        return tree

    @pytest.fixture
    def mock_calibrations(self):
        """创建模拟校准点"""
        calib1 = Mock(spec=CalibrationPoint)
        calib1.name = "Primates"
        calib1.resolved_taxa = ["human", "chimp"]
        calib1.mrca_leaf_pair = ("human", "chimp")
        calib1.is_root_node = False
        calib1.age_constraint = FixedAgeConstraint(fixed_age=6.0)

        return [calib1]

    def test_init_creates_work_directory(self, temp_dir, pathd8_config, common_config):
        """测试初始化创建工作目录"""
        method = PATHd8Method(pathd8_config, temp_dir, common_config=common_config)

        assert method.work_dir.exists()
        assert method.method_name == "pathd8"

    def test_method_name_property(self, temp_dir, pathd8_config, common_config):
        """测试方法名称属性"""
        method = PATHd8Method(pathd8_config, temp_dir, common_config=common_config)
        assert method.method_name == "pathd8"

    def test_common_config_usage(self, temp_dir, pathd8_config, common_config):
        """测试 common_config 使用"""
        method = PATHd8Method(pathd8_config, temp_dir, common_config=common_config)
        assert method.common_config.nthreads == 4

    @patch("phylodater.adapters.pathd8_method.ProcessRunner")
    def test_validate_environment_success(
        self, mock_runner_class, temp_dir, pathd8_config, common_config
    ):
        """测试环境验证成功"""
        mock_runner = Mock()
        mock_runner.check_executable.return_value = True
        mock_runner.get_version.return_value = "1.0"
        mock_runner_class.return_value = mock_runner

        method = PATHd8Method(pathd8_config, temp_dir, common_config=common_config)
        result = method.validate_environment()

        assert result is True

    @patch("phylodater.adapters.pathd8_method.ProcessRunner")
    def test_validate_environment_failure(
        self, mock_runner_class, temp_dir, pathd8_config, common_config
    ):
        """测试环境验证失败"""
        mock_runner = Mock()
        mock_runner.check_executable.return_value = False
        mock_runner_class.return_value = mock_runner

        method = PATHd8Method(pathd8_config, temp_dir, common_config=common_config)
        result = method.validate_environment()

        assert result is False

    def test_requires_fixed_age_constraint(
        self, temp_dir, pathd8_config, common_config
    ):
        """测试需要固定年龄约束

        注意：本测试原先只是给自己造的 Mock 断言它自己的属性
        （`assert has_fixed is True`），从未调用适配器，属于自我实现型测试。
        现在真正调用 prepare_inputs 走两条分支。
        """
        method = PATHd8Method(pathd8_config, temp_dir, common_config=common_config)

        # 没有任何固定年龄约束 → 必须拒绝
        without_fixed = [
            _make_cal(
                "Primates",
                UniformAgeConstraint(min_age=6.0, max_age=8.0),
                pair=("human", "chimp"),
            )
        ]
        with pytest.raises(CalibrationError, match="at least one fixed age"):
            method.prepare_inputs(_make_tree(), without_fixed)

        # 有固定年龄约束且能落到 MRCA 上 → 应当写出 fixage
        with_fixed = [
            _make_cal("Primates", FixedAgeConstraint(fixed_age=10.0)),
        ]
        method.prepare_inputs(_make_tree(), with_fixed)
        assert "mrca: human, chimp, fixage=10.0;" in (
            method._input_files["infile"].read_text()
        )

    def test_work_dir_isolation(self, temp_dir, pathd8_config, common_config):
        """测试工作目录隔离"""
        method1 = PATHd8Method(pathd8_config, temp_dir, common_config=common_config)
        method2 = PATHd8Method(
            pathd8_config, temp_dir / "other", common_config=common_config
        )

        assert method1.work_dir != method2.work_dir


class TestPATHd8Config:
    """PATHd8 配置测试"""

    def test_default_config(self):
        """测试默认配置"""
        config = PATHd8Config()
        assert config.pathd8_bin == "PATHd8"

    def test_custom_config(self):
        """测试自定义配置"""
        config = PATHd8Config(pathd8_bin="/usr/local/bin/PATHd8")
        assert config.pathd8_bin == "/usr/local/bin/PATHd8"

    def test_empty_bin_validation(self):
        """测试空路径验证"""
        with pytest.raises(ValueError):
            PATHd8Config(pathd8_bin="")


class TestWrittenConstraintLedger:
    """B-12 前置：约束清单必须从**写出去的文件**里回读，而不是从内存列表猜"""

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="pathd8_ledger_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def method(self, temp_dir):
        return PATHd8Method(
            PATHd8Config(), temp_dir, common_config=CommonConfig(nthreads=4)
        )

    def test_ledger_is_parsed_from_the_real_infile(self, method, temp_dir):
        """写盘后回读：kind/值/节点名都要还原出来"""
        method.prepare_inputs(
            _make_tree(),
            [
                _make_cal("HumanChimp", FixedAgeConstraint(fixed_age=6.0)),
                _make_cal(
                    "GreatApes",
                    UniformAgeConstraint(min_age=8.0, max_age=10.0),
                    pair=("human", "gorilla"),
                ),
            ],
        )
        kinds = [
            (e["kind"], e["value"], e["name"]) for e in method._written_constraints
        ]
        assert kinds == [
            ("fixage", 6.0, "HumanChimp"),
            ("minage", 8.0, "GreatApes"),
            ("maxage", 10.0, "GreatApes"),
        ]
        assert method._constraint_report["directives_written"] == {
            "fixage": 1,
            "minage": 1,
            "maxage": 1,
        }
        assert method._constraint_report["constrained_nodes_written"] == 2

    def test_guard_rejects_a_file_that_lost_its_fixage_line(self, method):
        """守卫看的是文件：把 fixage 行删掉后必须拒绝（模拟 B-6 那类上游丢点）"""
        method.prepare_inputs(
            _make_tree(),
            [
                _make_cal("HumanChimp", FixedAgeConstraint(fixed_age=6.0)),
                _make_cal(
                    "GreatApes",
                    UniformAgeConstraint(min_age=8.0, max_age=10.0),
                    pair=("human", "gorilla"),
                ),
            ],
        )
        infile = method._input_files["infile"]
        infile.write_text(
            "\n".join(
                line for line in infile.read_text().splitlines() if "fixage" not in line
            )
            + "\n"
        )
        assert "fixage" not in infile.read_text()

        with pytest.raises(CalibrationError, match="NO 'fixage=' directive"):
            method._guard_written_constraints(infile, method._calibrations, [])

    def test_expected_bounds_follow_upstream_merge_rules(self, method):
        """同节点多条同向边界：min 取最大、max 取最小；fixage 同时钉住两列"""
        entries = method._parse_written_constraints(
            "mrca: A, B, minage=5.0;\n"
            "mrca: A, B, minage=9.0;\n"
            "mrca: A, B, maxage=20.0;\n"
            "mrca: A, B, maxage=14.0;\n"
            "name of mrca: A, B, name=N1;\n"
            "mrca: C, D, fixage=7.5;\n"
            "name of mrca: C, D, name=N2;\n"
        )
        bounds = method._expected_bounds(entries)
        assert bounds["N1"] == {"min": 9.0, "max": 14.0}
        assert bounds["N2"] == {"min": 7.5, "max": 7.5}


class TestPathd8B12Guards:
    """B-12：'全部约束被丢弃、根被固定为 1' 的路径必须在写盘后被拦住"""

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="pathd8_b12_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def method(self, temp_dir):
        return PATHd8Method(
            PATHd8Config(), temp_dir, common_config=CommonConfig(nthreads=4)
        )

    def test_report_bypass_path_raises(self, method):
        """报告里那条完整通路：唯一的固定年龄校准缺 mrca_leaf_pair → 必须拒绝

        旧实现只检查内存列表（`any(isinstance(..., FixedAgeConstraint))`），
        守卫通过后该点因缺 mrca_leaf_pair 被跳过，文件里一行 fixage 都没有，
        于是 PATHd8 清空全部 min/max 并把根钉成 1.0，适配器照样返回"Ma"。
        """
        calibrations = [
            _make_cal(
                "HumanChimp",
                FixedAgeConstraint(fixed_age=6.0),
                pair=None,
            ),
            _make_cal(
                "GreatApes",
                UniformAgeConstraint(min_age=8.0, max_age=10.0),
                pair=("human", "gorilla"),
            ),
        ]
        with pytest.raises(CalibrationError) as exc:
            method.prepare_inputs(_make_tree(), calibrations)

        message = str(exc.value)
        # 直接引用上游 meta.c:30 那句话，让用户知道单位已经变了
        assert (
            "Any user given constraint is ignored and root age is fixed to 1" in message
        )
        assert "HumanChimp" in message  # 点名是哪个校准没落盘
        assert "fixage" in message

    def test_root_range_calibration_is_written_as_minmax(self, method):
        """根上的区间校准必须真的落盘（PATHd8 接受根 mrca 的 minage/maxage）。

        旧实现只处理 FixedAgeConstraint，其余根校准一律"不写文件、仅告警"，
        根年龄于是完全由 MPL 从树长推导——基准测试里实测到请求 [100,110] 而
        返回 153 Ma。已用 PATHd8 二进制确认根 mrca 的 min/max 会被计入
        "informative minnodes/maxnodes"，因此这里要求它真正写进 infile。
        """
        method.prepare_inputs(
            _make_tree(),
            [
                _make_cal(
                    "RootBound",
                    UniformAgeConstraint(min_age=8.0, max_age=10.0),
                    pair=None,
                    is_root=True,
                ),
                _make_cal("HumanChimp", FixedAgeConstraint(fixed_age=6.0)),
            ],
        )
        text = method._input_files["infile"].read_text()
        assert "fixage=6.0" in text
        assert "minage=8.0" in text
        assert "maxage=10.0" in text
        # 区间被真正写出，就不应再有"未写入"的告警
        recorded = " ".join(str(w) for w in method._warnings)
        assert "was NOT written" not in recorded

    def test_second_root_calibration_is_dropped_with_a_warning(self, method):
        """根只能被约束一次：第二个根校准必须点名告警，不得静默丢弃。"""
        method.prepare_inputs(
            _make_tree(),
            [
                _make_cal(
                    "RootBound",
                    UniformAgeConstraint(min_age=8.0, max_age=10.0),
                    pair=None,
                    is_root=True,
                ),
                _make_cal(
                    "RootBound2",
                    MaximumAgeConstraint(max_age=12.0),
                    pair=None,
                    is_root=True,
                ),
                _make_cal("HumanChimp", FixedAgeConstraint(fixed_age=6.0)),
            ],
        )
        recorded = " ".join(str(w) for w in method._warnings)
        assert "RootBound2" in recorded
        assert "NOT written" in recorded

    def test_real_no_fixnode_output_is_refused(self, method):
        """输出侧兜底：上游自报'约束全部作废'时不得把相对年龄当 Ma 返回"""
        method.prepare_inputs(
            _make_tree(),
            [
                _make_cal("HumanChimp", FixedAgeConstraint(fixed_age=6.0)),
                _make_cal(
                    "GreatApes",
                    UniformAgeConstraint(min_age=8.0, max_age=10.0),
                    pair=("human", "gorilla"),
                ),
            ],
        )
        method._input_files["outfile"].write_text(REAL_NO_FIXNODE_OUTPUT)

        with pytest.raises(CalibrationError, match="root age is fixed to 1"):
            method.parse_results()

    def test_unnormalized_root_age_is_refused_without_the_message(self, method):
        """第二道兜底：根年龄≈1 且用户从没给过这个量级 → 疑似未归一化

        夹具取自真实产物，仅删去 "No fixnodes were defined / Any user given
        constraint is ignored" 两行，用来单独验证根年龄判据本身有效。
        """
        trimmed = "\n".join(
            line
            for line in REAL_NO_FIXNODE_OUTPUT.splitlines()
            if "fixnodes were defined" not in line
            and "user given constraint" not in line
        )
        assert "Any user given constraint" not in trimmed
        method.prepare_inputs(
            _make_tree(),
            [
                _make_cal("HumanChimp", FixedAgeConstraint(fixed_age=600.0)),
                _make_cal(
                    "GreatApes",
                    UniformAgeConstraint(min_age=800.0, max_age=1000.0),
                    pair=("human", "gorilla"),
                ),
            ],
        )
        method._input_files["outfile"].write_text(trimmed)

        with pytest.raises(CalibrationError, match="un-normalized relative timescale"):
            method.parse_results()


class TestPathd8ConstraintReconciliation:
    """B-11：用真实产物核对"上游有没有改写我的约束"

    旧实现扫的是 "conflict"/"adjusted" 两个关键词，而上游对应的那句 printf
    被注释掉了（meta.c:36-43），所以那个分支在真实输出上永远不触发。
    下面两个断言把这件事钉住。
    """

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="pathd8_recon_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def method(self, temp_dir):
        return PATHd8Method(
            PATHd8Config(), temp_dir, common_config=CommonConfig(nthreads=4)
        )

    @pytest.mark.parametrize(
        "fixture_name, fixture",
        [
            ("REAL_CLEAN_OUTPUT", REAL_CLEAN_OUTPUT),
            ("REAL_SILENT_REWRITE_OUTPUT", REAL_SILENT_REWRITE_OUTPUT),
            ("REAL_NO_FIXNODE_OUTPUT", REAL_NO_FIXNODE_OUTPUT),
        ],
    )
    def test_upstream_never_says_conflict_or_adjusted(self, fixture_name, fixture):
        """B-11 的病灶本身：真实输出里根本没有那两个关键词（夹具自证）"""
        lowered = fixture.lower()
        assert "conflict" not in lowered, fixture_name
        assert "adjusted" not in lowered, fixture_name

    def test_dead_keyword_branch_is_gone_from_the_source(self):
        """源码里不得再留 'conflict'/'adjusted' 关键词检测（防止假安全感复生）"""
        from phylodater.adapters import pathd8_method as module

        source = Path(module.__file__).read_text(encoding="utf-8")
        code = "\n".join(
            line for line in source.splitlines() if not line.lstrip().startswith("#")
        )
        assert '"conflict" in' not in code
        assert '"adjusted" in' not in code

    def test_clean_run_produces_no_false_alarm(self, method):
        """约束全部生效时不得无病呻吟"""
        method.prepare_inputs(
            _make_tree(),
            [
                _make_cal("HumanChimp", FixedAgeConstraint(fixed_age=6.0)),
                _make_cal(
                    "GreatApes",
                    UniformAgeConstraint(min_age=8.0, max_age=10.0),
                    pair=("human", "gorilla"),
                ),
            ],
        )
        method._input_files["outfile"].write_text(REAL_CLEAN_OUTPUT)

        result = method.parse_results()

        report = result.metadata["constraint_reconciliation"]
        assert report["verified"] is True
        assert report["voided_nodes"] == []
        assert report["unlocated_nodes"] == []
        assert report["directive_count_deficit"] == 0
        assert report["informative_nodes_reported"] == {
            "fixage": 1,
            "minage": 1,
            "maxage": 1,
        }
        assert result.node_ages["HumanChimp"].mean_age == pytest.approx(6.0)
        assert result.node_ages["GreatApes"].mean_age == pytest.approx(10.0)
        assert not [
            w for w in result.warnings if "did not carry all constraints" in str(w)
        ]

    def test_silently_dropped_constraint_is_named(self, method):
        """真实改写场景：点名 RootMin 的边界被上游丢掉，并留在元数据里"""
        method.prepare_inputs(
            _make_tree(),
            [
                _make_cal("HumanChimp", FixedAgeConstraint(fixed_age=6.0)),
                _make_cal(
                    "RootMin",
                    SoftLowerBoundConstraint(min_age=10.0),
                    pair=("human", "gorilla"),
                ),
                _make_cal(
                    "GorMin",
                    SoftLowerBoundConstraint(min_age=10.0),
                    pair=("gorilla", "orangutan"),
                ),
            ],
        )
        method._input_files["outfile"].write_text(REAL_SILENT_REWRITE_OUTPUT)

        result = method.parse_results()

        report = result.metadata["constraint_reconciliation"]
        assert report["verified"] is True
        assert len(report["voided_nodes"]) == 1
        assert "RootMin" in report["voided_nodes"][0]
        # 写出去 2 条 minage，上游只报回 1 个 informative minnode
        assert report["directives_written"]["minage"] == 2
        assert report["informative_nodes_reported"]["minage"] == 1
        assert report["directive_count_deficit"] == 1

        flagged = [
            w for w in result.warnings if "did not carry all constraints" in str(w)
        ]
        assert len(flagged) == 1
        assert "RootMin" in str(flagged[0])

    def test_two_calibrations_sharing_one_mrca_are_counted_at_input_time(self, method):
        """两个校准指向同一个 MRCA：PATHd8 一个节点只留一个名字，必须当场说"""
        method.prepare_inputs(
            _make_tree(),
            [
                _make_cal("HumanChimp", FixedAgeConstraint(fixed_age=6.0)),
                _make_cal(
                    "Ghost",
                    UniformAgeConstraint(min_age=8.0, max_age=10.0),
                    pair=("human", "gorilla"),
                ),
                _make_cal(
                    "GreatApes",
                    UniformAgeConstraint(min_age=8.0, max_age=10.0),
                    pair=("human", "gorilla"),
                ),
            ],
        )
        infile_text = method._input_files["infile"].read_text()
        # 3 个校准只落出 2 个被命名的节点（PATHd8 后写的名覆盖先写的）
        assert infile_text.count("name of mrca: human, gorilla") == 2
        assert method._constraint_report["calibrations_supplied"] == 3
        assert method._constraint_report["constrained_nodes_written"] == 2
        assert any(
            "only 2 node(s) received a directive" in str(w) for w in method._warnings
        )

    def test_unlocatable_node_name_is_reported_not_fatal(self, method):
        """写出去的节点名在上游表里找不到：要说话，但不下死刑"""
        method.prepare_inputs(
            _make_tree(),
            [
                _make_cal("HumanChimp", FixedAgeConstraint(fixed_age=6.0)),
                _make_cal(
                    "GreatApes",
                    UniformAgeConstraint(min_age=8.0, max_age=10.0),
                    pair=("human", "gorilla"),
                ),
            ],
        )
        # REAL_SILENT_REWRITE_OUTPUT 里那一批节点名是 RootMin/GorMin/HumanChimp，
        # 没有 GreatApes → 该节点无法定位，但已定位的 HumanChimp 完全吻合
        method._input_files["outfile"].write_text(REAL_SILENT_REWRITE_OUTPUT)

        result = method.parse_results()

        report = result.metadata["constraint_reconciliation"]
        assert report["voided_nodes"] == []
        assert any("GreatApes" in item for item in report["unlocated_nodes"])
        assert not [w for w in result.warnings if "EVERY constraint" in str(w)]
        assert any("not locatable" in str(w) for w in result.warnings)

    def test_missing_table_is_disclosed_as_unverified(self, method):
        """核验不了就明说'未核验'，不得让人以为已经检查过"""
        method.prepare_inputs(
            _make_tree(),
            [_make_cal("HumanChimp", FixedAgeConstraint(fixed_age=6.0))],
        )
        method._input_files["outfile"].write_text(
            "\nd8 tree    : (human:6.000000,chimp:6.000000)HumanChimp:0;\n"
        )

        result = method.parse_results()

        report = result.metadata["constraint_reconciliation"]
        assert report["verified"] is False
        assert "table was not found" in report["note"]
        assert any("could NOT be performed" in str(w) for w in result.warnings)


class TestPathd8FatalOutput:
    """C-10：致命信息必须中止，而不是 warning 之后继续解析"""

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="pathd8_fatal_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def method(self, temp_dir):
        m = PATHd8Method(
            PATHd8Config(), temp_dir, common_config=CommonConfig(nthreads=4)
        )
        m.prepare_inputs(
            _make_tree(),
            [_make_cal("HumanChimp", FixedAgeConstraint(fixed_age=6.0))],
        )
        return m

    def test_fatal_stdout_aborts(self, method):
        method._process_stdout = REAL_FATAL_STDOUT
        with pytest.raises(ExecutionError, match="This is a contradiction"):
            method.parse_results()

    def test_exit_code_zero_is_not_trusted(self, method):
        """上游用 ERR_FORMAT==0 当退出码：返回码不得被当成功判据"""
        method._process_stdout = REAL_FATAL_STDOUT
        with pytest.raises(ExecutionError, match="exits with code 0"):
            method._check_fatal_upstream_output()

    def test_unexpected_returncode_still_aborts(self, method):
        with patch("phylodater.adapters.pathd8_method.ProcessRunner") as runner_cls:
            runner = Mock()
            runner.run.return_value = Mock(returncode=7, stdout="", stderr="boom")
            runner_cls.return_value = runner
            with pytest.raises(ExecutionError, match="returncode=7"):
                method.execute()

    def test_returncode_one_is_not_claimed_as_auto_adjustment(self, method):
        """实测：完全正常的一次 PATHd8 运行就是返回 1，不得据此宣称'年龄被调整过'"""
        with patch("phylodater.adapters.pathd8_method.ProcessRunner") as runner_cls:
            runner = Mock()
            runner.run.return_value = Mock(
                returncode=1, stdout="Calculation finished.", stderr=""
            )
            runner_cls.return_value = runner
            assert method.execute() is True

        recorded = " ".join(str(w) for w in method._warnings)
        assert "may have been automatically adjusted" not in recorded
        assert "ordinary successful run" in recorded

    def test_empty_result_file_is_an_error(self, method):
        method._input_files["outfile"].write_text("")
        with pytest.raises(ResultParsingError, match="empty"):
            method.parse_results()


class TestPathd8SequenceLength:
    """C-12：缺比对时不再伪造 `Sequence length = 1000`"""

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="pathd8_seqlen_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def method(self, temp_dir):
        return PATHd8Method(
            PATHd8Config(), temp_dir, common_config=CommonConfig(nthreads=4)
        )

    def test_no_alignment_omits_the_directive_and_says_so(self, method):
        method.prepare_inputs(
            _make_tree(),
            [_make_cal("HumanChimp", FixedAgeConstraint(fixed_age=6.0))],
            alignment_path=None,
        )
        text = method._input_files["infile"].read_text()
        assert "Sequence length" not in text
        assert any(
            "NOT rescaled by a fabricated default of 1000" in str(w)
            for w in method._warnings
        )

        # 缺省该指令时上游照常完成（已实测），年龄与标定值一致
        method._input_files["outfile"].write_text(REAL_CLEAN_OUTPUT)
        result = method.parse_results()
        assert result.metadata["sequence_length"] is None
        assert result.metadata["sequence_length_source"] == (
            "not_provided_branch_lengths_used_as_given"
        )
        assert result.node_ages["HumanChimp"].mean_age == pytest.approx(6.0)

    def test_alignment_value_is_written_verbatim(self, method, temp_dir):
        alignment = temp_dir / "aln.fasta"
        alignment.write_text(">human\nACGTACGTAC\n>chimp\nACGTACGTAC\n")
        method.prepare_inputs(
            _make_tree(),
            [_make_cal("HumanChimp", FixedAgeConstraint(fixed_age=6.0))],
            alignment_path=alignment,
        )
        text = method._input_files["infile"].read_text()
        assert "Sequence length = 10;" in text
        assert method._sequence_length_meta["sequence_length_source"] == "alignment"
