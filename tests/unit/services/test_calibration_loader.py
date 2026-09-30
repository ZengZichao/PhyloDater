"""
CalibrationLoader 服务单元测试

测试校准配置加载器的各项功能。

审阅项 B-24 / B-19 / C-37 / C-38 的行为契约变化都在本文件里钉住：

* 内容级校验（必填键、类型放置层级、min/max 放置层级、非空条目）；
* ``.nan`` / ``.inf`` 这类非有限年龄在 loader/validator 边界就被拒，
  不再被 ``except (KeyError, ValueError)`` 压成静默丢点；
* 条目不是映射时不再抛裸 ``AttributeError``；
* 缺 ``type`` 不再猜测为 ``fixed``；
* ``soft_lower`` 的 YAML 键 ``p``/``c`` 映射到改名后的
  ``offset_fraction``/``cauchy_scale``；
* 丢点必须汇报"请求 N / 生效 M"。
"""

import shutil
import tempfile
from pathlib import Path

import pytest
import yaml

from phylodater.core.exceptions import CalibrationError
from phylodater.infrastructure.safe_io import safe_writer
from phylodater.models import (
    FixedAgeConstraint,
    MaximumAgeConstraint,
    SoftLowerBoundConstraint,
    UniformAgeConstraint,
)
from phylodater.services.calibration_loader import (
    CalibrationLoader,
    hard_errors,
    validate_calibration_entry,
    validate_calibration_payload,
)


def _write(temp_dir, name, text):
    path = temp_dir / name
    path.write_text(text, encoding="utf-8")
    return path


class TestCalibrationLoader:
    """CalibrationLoader 测试类"""

    @pytest.fixture
    def temp_dir(self):
        """创建临时目录"""
        temp_path = Path(tempfile.mkdtemp(prefix="calib_loader_test_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def loader(self):
        """创建加载器实例"""
        return CalibrationLoader()

    def test_load_valid_yaml(self, temp_dir, loader):
        """测试加载有效的 YAML 文件"""
        config_data = {
            "calibrations": [
                {
                    "name": "Primates",
                    "mrca_pair": ["human", "chimp"],
                    "constraint": {
                        "type": "uniform",
                        "min": 6.0,
                        "max": 8.0,
                    },
                }
            ]
        }

        config_file = temp_dir / "calibrations.yaml"
        with safe_writer(config_file) as f:
            yaml.dump(config_data, f)

        calibrations = loader.load(config_file)

        assert len(calibrations) == 1
        assert calibrations[0].name == "Primates"

    def test_load_multiple_calibrations(self, temp_dir, loader):
        """测试加载多个校准点"""
        config_data = {
            "calibrations": [
                {
                    "name": "Primates",
                    "mrca_pair": ["human", "chimp"],
                    "constraint": {"type": "uniform", "min": 6.0, "max": 8.0},
                },
                {
                    "name": "Rodents",
                    "mrca_pair": ["mouse", "rat"],
                    "constraint": {"type": "soft_lower", "min": 10.0},
                },
                {
                    "name": "Root",
                    "is_root": True,
                    "constraint": {"type": "maximum", "max": 100.0},
                },
            ]
        }

        config_file = temp_dir / "calibrations.yaml"
        with safe_writer(config_file) as f:
            yaml.dump(config_data, f)

        calibrations = loader.load(config_file)

        assert len(calibrations) == 3

    def test_load_file_not_found(self, temp_dir, loader):
        """测试文件不存在时抛出异常"""
        nonexistent_file = temp_dir / "nonexistent.yaml"

        with pytest.raises(FileNotFoundError):
            loader.load(nonexistent_file)

    def test_load_empty_file(self, temp_dir, loader):
        """测试加载空文件"""
        config_file = temp_dir / "empty.yaml"
        config_file.write_text("")

        calibrations = loader.load(config_file)

        assert calibrations == []

    def test_load_missing_name(self, temp_dir, loader):
        """测试缺少名称的校准点被跳过，且以 error 级留痕（不再是 warning）"""
        config_data = {
            "calibrations": [
                {
                    "mrca_pair": ["human", "chimp"],
                    "constraint": {"type": "uniform", "min": 6.0, "max": 8.0},
                }
            ]
        }

        config_file = temp_dir / "no_name.yaml"
        with safe_writer(config_file) as f:
            yaml.dump(config_data, f)

        seen = _record_logger(loader)
        calibrations = loader.load(config_file)

        assert len(calibrations) == 0
        errors = [message for level, message in seen if level == "error"]
        assert errors, "缺 name 必须 error 级留痕"
        assert any("name" in message for message in errors)
        # 丢点必须被对账（B-24）
        assert any(
            "requested 1" in message for message in errors
        ), "必须汇报 '请求 N / 生效 M'"

    def test_parse_uniform_constraint(self, loader):
        """测试解析均匀分布约束"""
        constraint_data = {
            "type": "uniform",
            "min": 5.0,
            "max": 10.0,
        }

        constraint = loader._parse_constraint(constraint_data)

        assert isinstance(constraint, UniformAgeConstraint)
        assert constraint.min_age == 5.0
        assert constraint.max_age == 10.0

    def test_parse_fixed_constraint(self, loader):
        """测试解析固定年龄约束"""
        constraint_data = {
            "type": "fixed",
            "age": 100.0,
        }

        constraint = loader._parse_constraint(constraint_data)

        assert isinstance(constraint, FixedAgeConstraint)
        assert constraint.fixed_age == 100.0

    def test_parse_soft_lower_bound_constraint(self, loader):
        """测试解析软下界约束"""
        constraint_data = {
            "type": "soft_lower",
            "min": 10.0,
        }

        constraint = loader._parse_constraint(constraint_data)

        assert isinstance(constraint, SoftLowerBoundConstraint)
        assert constraint.min_age == 10.0

    def test_soft_lower_yaml_keys_map_to_renamed_fields(self, loader):
        """YAML 的 ``p``/``c`` 必须映射到改名后的 ``offset_fraction``/``cauchy_scale``。

        审阅项 B-2：constraints 的归属方把 ``tail_prob``/``tail_shape`` 改名为
        ``offset_fraction``/``cauchy_scale``（前者是 L() 槽 2 的相对偏移 P，
        后者是槽 3 的 Cauchy 尺度）。加载器若继续用旧名调用，会触发
        DeprecationWarning 并绕一圈兼容换算。
        """
        with warnings_as_errors():
            constraint = loader._parse_constraint(
                {"type": "soft_lower", "min": 100.0, "p": 0.25, "c": 3.0}
            )

        assert isinstance(constraint, SoftLowerBoundConstraint)
        assert constraint.offset_fraction == 0.25
        assert constraint.cauchy_scale == 3.0
        # 未给 tail_prob 时保持上游默认（槽 4 一直是 0.025）
        assert constraint.tail_prob == pytest.approx(0.025)

    def test_soft_lower_tail_prob_key_fills_mcmctree_slot4(self, loader):
        """新增可选键 ``tail_prob`` 写进 L() 的第 4 槽（旧实现从不填）。"""
        constraint = loader._parse_constraint(
            {"type": "soft_lower", "min": 100.0, "tail_prob": 0.05}
        )

        assert constraint.tail_prob == pytest.approx(0.05)

    def test_parse_maximum_constraint(self, loader):
        """测试解析最大年龄约束"""
        constraint_data = {
            "type": "maximum",
            "max": 50.0,
        }

        constraint = loader._parse_constraint(constraint_data)

        assert isinstance(constraint, MaximumAgeConstraint)
        assert constraint.max_age == 50.0

    def test_parse_unknown_constraint_type(self, loader):
        """测试未知约束类型返回 None（并 error 级留痕）"""
        constraint_data = {
            "type": "unknown_type",
        }

        constraint = loader._parse_constraint(constraint_data)

        assert constraint is None

    def test_parse_constraint_missing_type_is_not_guessed(self, loader):
        """C-38：缺 ``type`` 不再静默按 ``fixed`` 处理。

        旧行为让"忘写/写错层级的 type"被翻译成"缺 age 键"，用户会去补一个
        ``age:`` 而真正的问题（min/max 放错层级）依旧存在。
        """
        assert loader._parse_constraint({"age": 50.0}) is None

        messages = validate_entry_errors({"name": "A", "constraint": {"age": 50.0}})
        assert any("缺少 'type'" in m for m in messages)

    def test_nan_age_is_rejected_at_loader(self, temp_dir, loader):
        """B-19/B-24：``min: .nan``（YAML 规范字面量）必须在 loader 边界被拒。"""
        text = (
            "calibrations:\n"
            "  - name: TestClade\n"
            "    mrca_pair: [human, chimp]\n"
            "    constraint:\n"
            "      type: uniform\n"
            "      min: .nan\n"
            "      max: 100\n"
        )
        config_file = _write(temp_dir, "nan.yaml", text)

        seen = _record_logger(loader)
        calibrations = loader.load(config_file)

        assert calibrations == []
        assert any(
            "NaN" in message for level, message in seen
        ), "NaN 不能再被 except (KeyError, ValueError) 压成静默丢点"

    def test_inf_age_is_rejected_at_loader(self, temp_dir, loader):
        """B-19：``.inf`` 同样在边界被拒（constraints 层也会拒，这里是第二道门）。"""
        text = (
            "calibrations:\n"
            "  - name: TestClade\n"
            "    mrca_pair: [human, chimp]\n"
            "    constraint:\n"
            "      type: fixed\n"
            "      age: .inf\n"
        )
        config_file = _write(temp_dir, "inf.yaml", text)

        seen = _record_logger(loader)
        assert loader.load(config_file) == []
        assert any("无穷" in message or "inf" in message for _, message in seen)

    def test_scalar_entry_does_not_raise_attribute_error(self, temp_dir, loader):
        """C-37：条目写成标量时不得再以裸 AttributeError 崩到用户面前。"""
        config_file = _write(temp_dir, "scalar.yaml", "calibrations:\n  - TestClade\n")

        seen = _record_logger(loader)
        calibrations = loader.load(config_file)  # 旧版本在这里直接 AttributeError

        assert calibrations == []
        # 不只是"不崩"：还得点名**哪一条**、**该写成什么**
        errors = [message for level, message in seen if level == "error"]
        assert any("#1" in message for message in errors), errors
        assert any("必须是一个映射" in message for message in errors), errors
        assert any("constraint" in message for message in errors), errors

    def test_empty_entry_does_not_raise_attribute_error(self, temp_dir, loader):
        """C-37：``- `` 空条目同样必须被翻译成人话而不是 NoneType AttributeError。"""
        config_file = _write(
            temp_dir, "empty_entry.yaml", "calibrations:\n  -\n  - name: B\n"
        )

        seen = _record_logger(loader)
        assert loader.load(config_file) == []

        errors = [message for level, message in seen if level == "error"]
        assert any("#1" in message and "映射" in message for message in errors), errors

    def test_partial_loss_is_reported_as_error(self, temp_dir, loader):
        """B-24：N 个条目里丢了 k 个（k<N）必须 error 级对账，而不是两行 warning。"""
        text = (
            "calibrations:\n"
            "  - name: Good\n"
            "    mrca_pair: [human, chimp]\n"
            "    constraint: {type: fixed, age: 100}\n"
            "  - name: Bad\n"
            "    mrca_pair: [mouse, rat]\n"
            "    constraint: {type: uniform, min: 90}\n"  # 缺 max
        )
        config_file = _write(temp_dir, "partial.yaml", text)

        seen = _record_logger(loader)
        calibrations = loader.load(config_file)

        assert len(calibrations) == 1
        assert any("requested 2" in (message or "") for level, message in seen)
        assert any(level == "error" for level, _ in seen)

    def test_strict_mode_raises(self, temp_dir, loader):
        """``strict=True``：任何丢点立刻失败（供流水线/批量入口选用）。"""
        text = (
            "calibrations:\n"
            "  - name: Bad\n"
            "    constraint: {type: uniform, min: 90}\n"
        )
        config_file = _write(temp_dir, "strict.yaml", text)

        with pytest.raises(CalibrationError):
            loader.load(config_file, strict=True)

    def test_load_with_root_calibration(self, temp_dir, loader):
        """测试加载根节点校准（合法：根校准没有 mrca_pair）"""
        config_data = {
            "calibrations": [
                {
                    "name": "Root",
                    "is_root": True,
                    "constraint": {"type": "uniform", "min": 80.0, "max": 120.0},
                }
            ]
        }

        config_file = temp_dir / "root_calib.yaml"
        with safe_writer(config_file) as f:
            yaml.dump(config_data, f)

        calibrations = loader.load(config_file)

        assert len(calibrations) == 1
        assert calibrations[0].is_root_node is True
        assert calibrations[0].mrca_leaf_pair is None


class TestCalibrationLoaderEdgeCases:
    """CalibrationLoader 边界情况测试"""

    @pytest.fixture
    def loader(self):
        return CalibrationLoader()

    def test_empty_calibrations_list(self, loader):
        """测试空校准列表"""
        constraint = loader._parse_constraint(
            {"type": "uniform", "min": 5.0, "max": 10.0}
        )
        assert constraint is not None

    def test_unknown_key_is_rejected_not_ignored(self, loader):
        """B-24：拼错/放错层级的键以前被静默忽略，"你以为生效的参数实际没生效"。

        例如把 ``min``/``max`` 写在条目层、或把 ``type`` 写在条目层。
        """
        messages = validate_entry_errors(
            {
                "name": "TestClade",
                "constraint": {
                    "type": "uniform",
                    "min": 5.0,
                    "max": 10.0,
                    "extra_field": "ignored",
                },
            }
        )

        assert hard_errors(messages), "未知键必须报错"
        assert any("extra_field" in m for m in messages)
        assert (
            loader._parse_constraint(
                {"type": "uniform", "min": 5.0, "max": 10.0, "extra_field": "x"}
            )
            is None
        )

    def test_min_max_at_entry_level_is_detected(self):
        """B-24 用例 2：min/max 直接写在条目层（没有 constraint 包裹）。"""
        messages = validate_entry_errors(
            {"name": "TestClade", "type": "uniform", "min": 60.0, "max": 90.0}
        )

        assert any("缺少 'constraint' 块" in m for m in messages)

    def test_type_at_entry_level_is_detected(self):
        """B-24 用例 1：``type`` 放错层级（写在条目层）。"""
        messages = validate_entry_errors(
            {"name": "TestClade", "type": "uniform_bound", "constraint": {"age": 100}}
        )

        assert any("缺少 'constraint' 块" in m or "缺少 'type'" in m for m in messages)

    def test_reversed_interval_is_rejected_not_swapped(self):
        """min >= max 直接报错，不静默调转（与 constraints 层的自律一致）。"""
        messages = validate_entry_errors(
            {
                "name": "A",
                "constraint": {"type": "uniform", "min": 100.0, "max": 50.0},
            }
        )

        assert any("'min'" in m and "'max'" in m for m in hard_errors(messages))

    def test_age_above_earth_age_is_rejected(self):
        """年龄必须 <= 4600 Ma。"""
        messages = validate_entry_errors(
            {"name": "A", "constraint": {"type": "fixed", "age": 9000}}
        )

        assert any("4600" in m for m in hard_errors(messages))

    def test_boolean_age_is_rejected(self):
        """``age: yes`` 会被 YAML 解析成 True；bool 不是年龄。"""
        messages = validate_entry_errors(
            {"name": "A", "constraint": {"type": "fixed", "age": True}}
        )

        assert hard_errors(messages)

    def test_node_block_must_be_mapping(self):
        """``node: auto``（标量）会让 NodeDefinition.from_dict 崩 AttributeError。"""
        messages = validate_entry_errors(
            {
                "name": "A",
                "constraint": {"type": "fixed", "age": 10},
                "node": "auto",
            }
        )

        assert any("node" in m for m in hard_errors(messages))

    def test_empty_node_block_is_rejected(self):
        """``node:``（值为 None）同样必须报错而不是崩溃。"""
        messages = validate_entry_errors(
            {
                "name": "A",
                "constraint": {"type": "fixed", "age": 10},
                "node": None,
            }
        )

        assert hard_errors(messages)

    def test_mrca_pair_wrong_length_rejected(self):
        """``mrca_pair`` 必须是两个叶节点名。"""
        messages = validate_entry_errors(
            {
                "name": "A",
                "constraint": {"type": "fixed", "age": 10},
                "mrca_pair": ["only-one"],
            }
        )

        assert any("mrca_pair" in m for m in hard_errors(messages))

    def test_payload_with_no_entries_is_invalid(self):
        """``calibrations: []`` 必须在入口就被判为无效（而不是等到零值守卫）。"""
        messages = validate_calibration_payload({"calibrations": []}, source="x.yaml")

        assert hard_errors(messages)

    def test_all_eight_types_have_schemas(self, loader):
        """八类合法约束都在契约表里，且都能被解析出来。"""
        samples = {
            "fixed": {"type": "fixed", "age": 100.0},
            "uniform": {"type": "uniform", "min": 90.0, "max": 110.0},
            "soft_lower": {"type": "soft_lower", "min": 90.0},
            "maximum": {"type": "maximum", "max": 110.0},
            "soft_bounds": {
                "type": "soft_bounds",
                "min": 90.0,
                "max": 110.0,
                "p": 0.025,
            },
            "gamma": {"type": "gamma", "alpha": 2.0, "beta": 0.5},
            "skew_normal": {
                "type": "skew_normal",
                "location": 100.0,
                "scale": 10.0,
                "shape": 1.0,
            },
            "skew_t": {
                "type": "skew_t",
                "location": 100.0,
                "scale": 10.0,
                "shape": 1.0,
                "df": 5.0,
            },
        }

        assert set(samples) == {
            "fixed",
            "uniform",
            "soft_lower",
            "maximum",
            "soft_bounds",
            "gamma",
            "skew_normal",
            "skew_t",
        }
        for name, block in samples.items():
            constraint = loader._parse_constraint(block)
            assert constraint is not None, name

    def test_skew_t_df_must_exceed_two(self, loader):
        """df <= 2 时 t 分布方差发散，必须在校验层就拒。"""
        messages = validate_entry_errors(
            {
                "name": "A",
                "constraint": {
                    "type": "skew_t",
                    "location": 100.0,
                    "scale": 10.0,
                    "shape": 1.0,
                    "df": 2.0,
                },
            }
        )

        assert any("df" in m for m in hard_errors(messages))


# --------------------------------------------------------------------------- #
# 审阅项 C-37 / C-38：确认性契约（B-24 那趟修复已改掉行为，这里把形状钉死）
# --------------------------------------------------------------------------- #
class TestNonMappingEntryAndMissingTypeDisclosure:
    """C-37 / C-38 的确认性测试。

    * **C-37**：非映射条目（标量 / 空 / 列表）必须以带条目序号的**人话**报错，
      而不是 ``AttributeError: 'str' object has no attribute 'get'``；
      ``strict=True`` 时抛 :class:`CalibrationError`（同样不是 AttributeError）。
    * **C-38**：缺 ``type`` 不再猜成 ``fixed``。报告实测三条不同成因（``type``
      放在条目层、``min``/``max`` 放在条目层、完全没有 ``constraint``）都只报一句
      ``Failed to parse constraint: 'age'``，把用户支到错误的改法上；现在三种成因
      各自点名真实成因，且**没有一种**再报"缺 age"。
    """

    @pytest.fixture
    def loader(self):
        return CalibrationLoader()

    def test_list_entry_is_reported_as_bad_shape(self):
        errors = validate_calibration_entry(["name", "TestClade"], 2)

        assert any("#2" in message for message in errors), errors
        assert any("必须是一个映射" in message for message in errors), errors
        assert any("list" in message for message in errors), errors

    def test_strict_mode_raises_calibration_error_not_attribute_error(
        self, temp_dir_factory, loader
    ):
        """报告实测的两份形状合法 YAML（标量条目 / 空条目）在 strict 下的形态。"""
        for name, text in (
            ("scalar.yaml", "calibrations:\n  - TestClade\n"),
            ("nested.yaml", "calibrations:\n  -\n  - name: B\n"),
        ):
            config_file = _write(temp_dir_factory, name, text)
            _record_logger(loader)

            with pytest.raises(CalibrationError) as excinfo:
                loader.load(config_file, strict=True)

            assert "#1" in str(excinfo.value), (name, str(excinfo.value))

    def test_parse_calibration_survives_non_mapping_input(self, loader):
        """绕开 ``load`` 直调同样不得抛裸 AttributeError。"""
        for bad in ("TestClade", None, 42, ["name", "A"]):
            with pytest.raises(CalibrationError):
                loader._parse_calibration(bad, index=3)

    def test_three_same_symptom_miswrites_are_told_apart(self):
        """C-38：三种成因必须给出三种指向真实成因的文案。"""
        type_at_entry = validate_entry_errors(
            {"name": "A", "type": "uniform", "min": 5.0, "max": 8.0}
        )
        no_constraint = validate_entry_errors({"name": "A", "mrca_pair": ["x", "y"]})
        type_missing_in_block = validate_entry_errors(
            {"name": "A", "constraint": {"min": 5.0, "max": 8.0}}
        )

        # 1) type 放错层级 / 2) 根本没有 constraint → 都点名"缺 constraint 块"
        assert any("缺少 'constraint' 块" in m for m in type_at_entry), type_at_entry
        assert any("缺少 'constraint' 块" in m for m in no_constraint), no_constraint
        # 并把实际收到的键列出来，用户才知道自己写成了什么
        assert any(
            "['max', 'min', 'name', 'type']" in m for m in type_at_entry
        ), type_at_entry
        # 3) constraint 块里缺 type → 点名缺 type + 已知键
        assert any(
            "缺少 'type'" in m for m in type_missing_in_block
        ), type_missing_in_block
        assert any("已知键" in m for m in type_missing_in_block), type_missing_in_block

        # 三条都不许再报"缺 age"（旧文案把用户支到错误改法上）
        for messages in (type_at_entry, no_constraint, type_missing_in_block):
            assert not any("必填键 'age'" in m for m in messages), messages

    def test_missing_type_never_falls_back_to_fixed(self, loader):
        """``type`` 缺失即无效：不能构造出任何约束对象（更不会是 fixed）。"""
        assert loader._parse_constraint({"min": 5.0, "max": 8.0}) is None
        assert loader._parse_constraint({"age": 50.0}) is None

        messages = validate_entry_errors(
            {"name": "A", "constraint": {"min": 5.0, "max": 8.0}}
        )
        assert hard_errors(messages)

    def test_bad_shape_entry_enters_the_reconciliation(self, temp_dir_factory, loader):
        """坏形状条目也必须计入"请求 N / 生效 M"（B-24 与 C-37 的交集）。"""
        config_file = _write(
            temp_dir_factory,
            "mixed.yaml",
            "calibrations:\n"
            "  - name: Good\n"
            "    mrca_pair: [human, chimp]\n"
            "    constraint: {type: fixed, age: 100}\n"
            "  - JustACladeName\n",
        )

        seen = _record_logger(loader)
        calibrations = loader.load(config_file)

        assert len(calibrations) == 1
        errors = [message for level, message in seen if level == "error"]
        assert any(
            "requested 2" in message and "1 were applied" in message
            for message in errors
        ), errors
        assert any("#2" in message for message in errors), errors


# --------------------------------------------------------------------------- #
# 测试辅助
# --------------------------------------------------------------------------- #
@pytest.fixture
def temp_dir_factory():
    """模块级临时目录（供上面的独立测试类使用）。"""
    path = Path(tempfile.mkdtemp(prefix="calib_loader_case_"))
    yield path
    shutil.rmtree(path, ignore_errors=True)


class _LogRecorder:
    """把 logger 的调用收集成 (level, message) 列表。"""

    def __init__(self):
        self.records = []

    def _record(self, level):
        def handler(message=None, *args, **kwargs):
            self.records.append((level, str(message)))

        return handler

    def __getattr__(self, name):
        if name in ("info", "warning", "error", "success", "debug", "step"):
            return self._record(name)
        raise AttributeError(name)


def _record_logger(loader):
    recorder = _LogRecorder()
    loader.logger = recorder
    return recorder.records


class warnings_as_errors:
    """上下文：把 DeprecationWarning 升级为异常（用于确认不再走旧字段名）。"""

    def __enter__(self):
        import warnings

        self._warnings = warnings
        self._cm = warnings.catch_warnings()
        self._cm.__enter__()
        warnings.simplefilter("error", DeprecationWarning)
        return self

    def __exit__(self, *exc_info):
        self._cm.__exit__(*exc_info)
        return False


def validate_entry_errors(entry):
    """内容级校验单条目（与加载器/CLI 校验器共用同一份契约表）。"""
    return validate_calibration_entry(entry, 1)
