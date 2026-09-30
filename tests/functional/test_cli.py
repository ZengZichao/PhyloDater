"""
模块一：CLI 与子命令架构功能测试

验证子命令路由、通用参数、dating 参数解析及退出码。
"""

import pytest

from tests.functional.conftest import run_phylodater_cli


class TestCLISubcommands:
    """CLI 子命令路由测试"""

    def test_cli_no_subcommand_shows_help(self):
        """CLI-001: 无子命令时显示帮助并退出码 2"""
        result = run_phylodater_cli([])
        assert result.returncode == 2
        assert "usage" in (result.stdout + result.stderr).lower()

    def test_cli_dating_help(self):
        """CLI-002: dating 子命令帮助"""
        result = run_phylodater_cli(["dating", "--help"])
        assert result.returncode == 0
        assert "usage" in result.stdout.lower()
        assert "--tree" in result.stdout

    def test_cli_check_runs_self_test(self):
        """CLI-003: check 子命令执行自检"""
        result = run_phylodater_cli(["check"])
        assert result.returncode in [0, 1]
        assert "[PASS]" in result.stdout or "[FAIL]" in result.stdout

    def test_cli_fold_help(self):
        """CLI-004: fold 子命令帮助

        C-36：`fold` 有完整参数面但处理函数只 raise NotImplementedError。
        帮助本身必须**在用户决定调用之前**就说清"尚未实现"，否则用户要跑到
        最后一行才发现功能不存在（而手稿软件表还把它当作已交付能力，P0-7）。
        """
        result = run_phylodater_cli(["fold", "--help"])
        assert result.returncode == 0
        assert "usage" in result.stdout.lower()
        assert "--output" in result.stdout
        combined = result.stdout + result.stderr
        assert "尚未实现" in combined
        assert "not yet implemented" in combined.lower()

    def test_cli_top_level_help_marks_fold_unimplemented(self):
        """CLI-004b: 主 --help 的子命令清单里 fold 不得伪装成已交付能力（C-36）。"""
        result = run_phylodater_cli(["--help"])
        assert result.returncode == 0
        lines = result.stdout.splitlines()
        fold_lines = [ln for ln in lines if ln.strip().startswith("fold")]
        assert fold_lines, "fold 仍应出现在子命令清单里（隐藏不是诚实做法）"
        assert all(
            "尚未实现" in ln or "not yet implemented" in ln.lower() for ln in fold_lines
        ), fold_lines

    def test_cli_unknown_subcommand(self):
        """CLI-005: 未知子命令"""
        result = run_phylodater_cli(["unknown"])
        assert result.returncode == 2
        assert "unknown" in (result.stdout + result.stderr).lower()

    def test_cli_version(self):
        """CLI-012: 版本信息"""
        result = run_phylodater_cli(["--version"])
        assert result.returncode == 0
        assert "phylodater" in result.stdout.lower()


class TestCLIGeneralArguments:
    """通用参数测试"""

    def test_cli_log_level_debug(self):
        """CLI-010: DEBUG 日志级别"""
        result = run_phylodater_cli(["--log-level", "DEBUG", "check"])
        assert result.returncode in [0, 1]
        assert "DEBUG" in result.stdout or "DEBUG" in result.stderr

    def test_cli_log_file(self, tmp_path):
        """CLI-011: 日志文件"""
        log_file = tmp_path / "test.log"
        result = run_phylodater_cli(["--log-file", str(log_file), "check"])
        assert result.returncode in [0, 1]
        assert log_file.exists()
        content = log_file.read_text(encoding="utf-8")
        assert len(content) > 0


class TestCLIDatingArguments:
    """dating 子命令参数测试"""

    def test_cli_dating_missing_required_args(self):
        """CLI-020: 缺少必填参数"""
        result = run_phylodater_cli(["dating"])
        assert result.returncode == 2
        assert "--tree" in (result.stdout + result.stderr)

    def test_cli_dating_invalid_method(
        self, normal_tree, normal_alignment, normal_calibration, tmp_path
    ):
        """CLI-021: 无效方法名"""
        output_dir = tmp_path / "out"
        result = run_phylodater_cli(
            [
                "dating",
                "-t",
                str(normal_tree),
                "-s",
                str(normal_alignment),
                "-c",
                str(normal_calibration),
                "-o",
                str(output_dir),
                "--method",
                "invalid_method",
            ]
        )
        assert result.returncode == 2
        assert "invalid_method" in (result.stdout + result.stderr)

    def test_cli_dating_method_args_unknown_key_rejected(
        self, normal_tree, normal_alignment, normal_calibration, tmp_path
    ):
        """CLI-022a: ``--method-args`` 的未知参数键必须被拒绝，不得静默丢弃。

        ``Configuration.set_cli_override`` 按各方法配置 dataclass 的字段白名单
        校验参数名，未知键抛 ``ValueError("Unknown parameter ...")``，CLI 把它
        归为输入数据错误（``EXIT_DATA_ERROR`` = 3）。

        这是**正确**的严格性，与审阅报告反复批评的失效形态相对立：C-27 要求
        "任何被丢弃的片段一律 raise"（对照 B-24/B-26 的"入口只验形状、深层
        降级为日志"）。PATHd8 本就没有任何可调参数（见 ``dating --help`` 的
        epilog："pathd8: (无额外参数)"），因此旧断言里的 ``returncode in [0, 1]``
        实际是在奖励"用户以为传了参数、其实一个都没生效"这一行为。
        """
        output_dir = tmp_path / "out"
        result = run_phylodater_cli(
            [
                "dating",
                "-t",
                str(normal_tree),
                "-s",
                str(normal_alignment),
                "-c",
                str(normal_calibration),
                "-o",
                str(output_dir),
                "--method",
                "pathd8",
                "--method-args",
                "key1=value1,key2=value2",
            ]
        )
        combined = result.stdout + result.stderr
        assert result.returncode == 3, combined
        assert "Unknown parameter" in combined
        # 必须指名被拒绝的具体键，而不是笼统报"参数错误"
        assert "key1" in combined
        assert "pathd8" in combined

    def test_cli_dating_method_args_unknown_key_suggests_closest(
        self, normal_tree, normal_alignment, normal_calibration, tmp_path
    ):
        """CLI-022b: 拼错的参数名要给出可执行的纠正建议。

        ``set_cli_override`` 用 ``difflib.get_close_matches`` 提示相近字段；
        拒绝如果只是"不行"而不说"你想要的是不是这个"，用户仍会去猜。
        """
        output_dir = tmp_path / "out"
        result = run_phylodater_cli(
            [
                "dating",
                "-t",
                str(normal_tree),
                "-s",
                str(normal_alignment),
                "-c",
                str(normal_calibration),
                "-o",
                str(output_dir),
                "--method",
                "pathd8",
                "--method-args",
                "pathd8_bin2=PATHd8",
            ]
        )
        combined = result.stdout + result.stderr
        assert result.returncode == 3, combined
        assert "pathd8_bin" in combined

    def test_cli_dating_method_args_known_key_accepted(
        self, normal_tree, normal_alignment, normal_calibration, tmp_path
    ):
        """CLI-022c: 合法的参数键照常接受（严格白名单不误伤真实通路）。

        ``pathd8_bin`` 是 ``PATHd8Config`` 的真实字段，传它应顺利通过配置层，
        之后是否跑成取决于本机是否装了 PATHd8 —— 因此允许 0（成功）或
        1（外部软件不可用/运行失败），但不能是 2（用法错误）或 3（数据错误）。
        """
        output_dir = tmp_path / "out"
        result = run_phylodater_cli(
            [
                "dating",
                "-t",
                str(normal_tree),
                "-s",
                str(normal_alignment),
                "-c",
                str(normal_calibration),
                "-o",
                str(output_dir),
                "--method",
                "pathd8",
                "--method-args",
                "pathd8_bin=PATHd8",
            ]
        )
        combined = result.stdout + result.stderr
        assert "Unknown parameter" not in combined
        assert result.returncode in [0, 1], combined

    def test_cli_dating_method_args_invalid_format(
        self, normal_tree, normal_alignment, normal_calibration, tmp_path
    ):
        """CLI-023: method-args 格式错误"""
        output_dir = tmp_path / "out"
        result = run_phylodater_cli(
            [
                "dating",
                "-t",
                str(normal_tree),
                "-s",
                str(normal_alignment),
                "-c",
                str(normal_calibration),
                "-o",
                str(output_dir),
                "--method",
                "pathd8",
                "--method-args",
                "invalid_format",
            ]
        )
        assert result.returncode == 2
        assert "key=value" in (result.stdout + result.stderr)

    @pytest.mark.parametrize("threads_value", ["0", "abc"])
    def test_cli_dating_threads_validation(
        self, threads_value, normal_tree, normal_alignment, normal_calibration, tmp_path
    ):
        """CLI-024/025: threads 正整数校验"""
        output_dir = tmp_path / "out"
        result = run_phylodater_cli(
            [
                "dating",
                "-t",
                str(normal_tree),
                "-s",
                str(normal_alignment),
                "-c",
                str(normal_calibration),
                "-o",
                str(output_dir),
                "--method",
                "pathd8",
                "--threads",
                threads_value,
            ]
        )
        assert result.returncode == 2
        combined = result.stdout + result.stderr
        # 0 触发“必须是正整数”，abc 触发“必须是整数”
        assert any(k in combined for k in ["正整数", "整数", "integer"])

    @pytest.mark.parametrize("mode", ["split", "first", "last", "random", "error"])
    def test_cli_dating_multi_tree_mode_valid(
        self, mode, normal_tree, normal_alignment, normal_calibration, tmp_path
    ):
        """CLI-026: 多棵树模式"""
        output_dir = tmp_path / "out"
        result = run_phylodater_cli(
            [
                "dating",
                "-t",
                str(normal_tree),
                "-s",
                str(normal_alignment),
                "-c",
                str(normal_calibration),
                "-o",
                str(output_dir),
                "--method",
                "pathd8",
                "--multi-tree-mode",
                mode,
            ]
        )
        assert result.returncode in [0, 1]

    def test_cli_dating_multi_tree_mode_invalid(
        self, normal_tree, normal_alignment, normal_calibration, tmp_path
    ):
        """CLI-027: 无效多棵树模式"""
        output_dir = tmp_path / "out"
        result = run_phylodater_cli(
            [
                "dating",
                "-t",
                str(normal_tree),
                "-s",
                str(normal_alignment),
                "-c",
                str(normal_calibration),
                "-o",
                str(output_dir),
                "--method",
                "pathd8",
                "--multi-tree-mode",
                "invalid",
            ]
        )
        assert result.returncode == 2

    def test_cli_fold_normal_call(self, normal_tree, tmp_path):
        """CLI-028: fold 调用（功能未实现应返回运行时错误码而非假成功）。

        §17: ``_run_fold`` 已显式抛出 ``NotImplementedError``，``main`` 捕获后
        返回 ``EXIT_RUNTIME_ERROR``(=1)，不再返回成功（假阳性）。

        C-36 追加两点：未实现声明必须在回显解析参数**之前**给出（旧实现的日志
        前几行与已实现子命令毫无区别），并且不得留下任何看似产物的文件。
        """
        output_dir = tmp_path / "out"
        result = run_phylodater_cli(
            ["fold", "-t", str(normal_tree), "-o", str(output_dir)]
        )
        combined = result.stdout + result.stderr
        assert result.returncode == 1, combined
        assert "未实现" in combined

        warning_at = combined.find("fold 子命令【实验性")
        echo_at = combined.find("折叠级别")
        assert warning_at != -1, combined
        assert echo_at != -1, combined
        assert warning_at < echo_at, "未实现声明必须先于参数回显出现"

        # 只允许日志文件存在，不得有"折叠树"产物
        produced = sorted(p.name for p in output_dir.iterdir())
        assert all(name.endswith(".log") for name in produced), produced

    def test_cli_fold_unknown_rank_is_still_usage_error(self, normal_tree, tmp_path):
        """CLI-028b: 参数面照常校验（标注未实现不等于放松解析）。"""
        result = run_phylodater_cli(
            [
                "fold",
                "-t",
                str(normal_tree),
                "-o",
                str(tmp_path / "out"),
                "--fold-rank",
                "kingdom",
            ]
        )
        assert result.returncode == 2
        assert "--fold-rank" in (result.stdout + result.stderr)

    def test_cli_fold_missing_output(self, normal_tree):
        """CLI-029: fold 缺少必填参数"""
        result = run_phylodater_cli(["fold", "-t", str(normal_tree)])
        assert result.returncode == 2
        assert "--output" in (result.stdout + result.stderr)

    def test_cli_dating_empty_calibration_exits_code_3(
        self, normal_tree, normal_alignment, tmp_path
    ):
        """CLI-030: 空校准文件退出码 3"""
        empty_cal = tmp_path / "empty.yaml"
        empty_cal.write_text("calibrations: []\n", encoding="utf-8")
        output_dir = tmp_path / "out"
        result = run_phylodater_cli(
            [
                "dating",
                "-t",
                str(normal_tree),
                "-s",
                str(normal_alignment),
                "-c",
                str(empty_cal),
                "-o",
                str(output_dir),
                "--method",
                "pathd8",
            ]
        )
        assert result.returncode == 3
        assert (
            "未加载" in result.stdout
            or "calibrat" in (result.stdout + result.stderr).lower()
        )

    def test_cli_dating_taxonomy_table_without_name_column_exits_3(
        self, normal_tree, normal_alignment, normal_calibration, tmp_path
    ):
        """CLI-031 (B-26): 分类学表整份加载失败时必须终止，不得报"全匹配"。

        这份 TSV 有 ``domain``/``phylum`` 等级别列，却没有可识别的名称列
        （``group_id`` 不在 ``NAME_COLUMN_ALIASES`` 里）。修复前：
        ``load_from_file`` 返回 0 且只打一条 INFO，CLI 丢弃该返回值，
        ``check_label_consistency`` 在"表为空"时把**全部**叶节点算成已匹配，
        于是"分类学表↔树标签一致性检查"呈现为 100% 通过。

        现在这条路径上有两道闸（STEP 1 的结构校验与 STEP 4 的加载/一致性核验），
        任一道都必须以数据错误码退出。
        """
        bad_table = tmp_path / "taxonomy.tsv"
        bad_table.write_text(
            "group_id\tdomain\tphylum\tclass\n"
            "G1\tBacteria\tProteobacteria\tGamma\n"
            "G2\tBacteria\tProteobacteria\tGamma\n",
            encoding="utf-8",
        )
        output_dir = tmp_path / "out"
        result = run_phylodater_cli(
            [
                "dating",
                "-t",
                str(normal_tree),
                "-s",
                str(normal_alignment),
                "-c",
                str(normal_calibration),
                "-o",
                str(output_dir),
                "--method",
                "pathd8",
                "--taxonomy-file",
                str(bad_table),
            ]
        )
        combined = result.stdout + result.stderr
        assert result.returncode == 3, combined
        # 必须点名"这张表不可用"，而不是给出任何"已核对/全部匹配"的表述
        assert "分类学" in combined or "taxonomy" in combined.lower()
        assert "100%" not in combined
        assert "tree tips not found in taxonomy table: []" not in combined

    def test_cli_dating_taxonomy_zero_match_exits_3(
        self, normal_tree, normal_alignment, normal_calibration, tmp_path
    ):
        """CLI-032 (B-26): 表能加载但与树标签零交集，同样属于数据错误。

        表格结构完全合法（含 ``name`` 与级别列，能过 STEP 1 校验），但条目
        是 ``SpA/SpB``，而树是 ``human/chimp/dog/mouse/rat``。零匹配意味着
        "表格与这棵树不是同一批样本"，继续跑等于把一项没有核对成功的检查
        当作已核对。
        """
        mismatched = tmp_path / "taxonomy.tsv"
        mismatched.write_text(
            "name\tdomain\tphylum\n"
            "SpA\tBacteria\tProteobacteria\n"
            "SpB\tBacteria\tFirmicutes\n",
            encoding="utf-8",
        )
        output_dir = tmp_path / "out"
        result = run_phylodater_cli(
            [
                "dating",
                "-t",
                str(normal_tree),
                "-s",
                str(normal_alignment),
                "-c",
                str(normal_calibration),
                "-o",
                str(output_dir),
                "--method",
                "pathd8",
                "--taxonomy-file",
                str(mismatched),
            ]
        )
        combined = result.stdout + result.stderr
        assert result.returncode == 3, combined
        assert "零匹配" in combined
        assert "已加载 2 条分类学条目" in combined

    def test_cli_dating_taxonomy_matched_reports_counts(
        self, normal_tree, normal_alignment, normal_calibration, tmp_path
    ):
        """CLI-033 (B-26): 有交集时一致性检查照常放行，并把匹配数说清楚。

        防止修复"过头"变成"只要有未匹配项就一律终止"：树里 5 个叶，表里
        给出其中 2 个 + 1 个表外条目，应当继续运行并逐条 warning。
        """
        table = tmp_path / "taxonomy.tsv"
        table.write_text(
            "name\tdomain\tphylum\n"
            "human\tEukaryota\tChordata\n"
            "mouse\tEukaryota\tChordata\n"
            "ghost\tEukaryota\tChordata\n",
            encoding="utf-8",
        )
        output_dir = tmp_path / "out"
        result = run_phylodater_cli(
            [
                "dating",
                "-t",
                str(normal_tree),
                "-s",
                str(normal_alignment),
                "-c",
                str(normal_calibration),
                "-o",
                str(output_dir),
                "--method",
                "pathd8",
                "--taxonomy-file",
                str(table),
            ]
        )
        combined = result.stdout + result.stderr
        assert "零匹配" not in combined
        assert "一致性检查" in combined
        assert "3 tree tips not found in taxonomy table" in combined
        assert "ghost" in combined
        # pathd8 是否可用取决于本机环境，但绝不能是数据/用法错误
        assert result.returncode in [0, 1], combined

    @pytest.mark.parametrize("reroot", ["midpoint", "outgroup", "none"])
    def test_cli_dating_reroot_valid(
        self, reroot, normal_tree, normal_alignment, normal_calibration, tmp_path
    ):
        """PIP-048/049/050: 有效 reroot 策略（不依赖外部软件）"""
        output_dir = tmp_path / "out"
        args = [
            "dating",
            "-t",
            str(normal_tree),
            "-s",
            str(normal_alignment),
            "-c",
            str(normal_calibration),
            "-o",
            str(output_dir),
            "--method",
            "pathd8",
            "--reroot",
            reroot,
        ]
        if reroot == "outgroup":
            args.extend(["--outgroup", "human"])
        result = run_phylodater_cli(args)
        assert result.returncode in [0, 1]

    def test_cli_dating_reroot_mad_external_dependency(
        self, normal_tree, normal_alignment, normal_calibration, tmp_path
    ):
        """PIP-050: MAD 定根依赖外部 mad 程序"""
        import shutil

        pytest = __import__("pytest")
        if not shutil.which("mad"):
            pytest.skip("MAD 程序未安装")
        output_dir = tmp_path / "out"
        result = run_phylodater_cli(
            [
                "dating",
                "-t",
                str(normal_tree),
                "-s",
                str(normal_alignment),
                "-c",
                str(normal_calibration),
                "-o",
                str(output_dir),
                "--method",
                "pathd8",
                "--reroot",
                "mad",
            ]
        )
        assert result.returncode in [0, 1]

    def test_cli_dating_reroot_invalid(
        self, normal_tree, normal_alignment, normal_calibration, tmp_path
    ):
        """PIP-047: 未知 reroot 策略"""
        output_dir = tmp_path / "out"
        result = run_phylodater_cli(
            [
                "dating",
                "-t",
                str(normal_tree),
                "-s",
                str(normal_alignment),
                "-c",
                str(normal_calibration),
                "-o",
                str(output_dir),
                "--method",
                "pathd8",
                "--reroot",
                "invalid",
            ]
        )
        assert result.returncode == 2
