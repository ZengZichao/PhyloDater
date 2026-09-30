"""
Pipeline 集成测试

测试 DatingPipeline 的完整工作流
"""

import shutil
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from phylodater.core.pipeline import DatingPipeline, PipelineConfig
from phylodater.models import CalibrationPoint, DatingResult, PhylogeneticTree


class TestDatingPipelineIntegration:
    """DatingPipeline 集成测试类"""

    @pytest.fixture
    def temp_dir(self):
        """创建临时目录"""
        temp_path = Path(tempfile.mkdtemp(prefix="pipeline_test_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def pipeline_config(self, temp_dir):
        """创建流水线配置"""
        return PipelineConfig(
            methods=["mcmctree", "lsd2"],
            output_dir=temp_dir / "output",
            threads=2,
            skip_failed=True,
        )

    @pytest.fixture
    def mock_tree(self):
        """创建模拟树"""
        tree = Mock(spec=PhylogeneticTree)
        tree.newick = "((human:0.1,chimp:0.1):0.2,(mouse:0.3,rat:0.3):0.1);"
        tree.tip_names = ["human", "chimp", "mouse", "rat"]
        tree.is_rooted.return_value = True
        tree.is_binary.return_value = True
        tree.validate_branch_lengths.return_value = (True, [])
        tree.get_rooting_info.return_value = {
            "is_rooted": True,
            "likely_artifact": False,
        }
        return tree

    @pytest.fixture
    def mock_calibrations(self):
        """创建模拟校准点"""
        calib = Mock(spec=CalibrationPoint)
        calib.name = "Primates"
        calib.mrca_pair = ["human", "chimp"]
        calib.is_root = False
        calib.node = "Node1"
        calib.resolved_taxa = ["human", "chimp"]
        calib.age_constraint = None
        return [calib]

    def test_pipeline_initialization(self, pipeline_config):
        """测试流水线初始化"""
        pipeline = DatingPipeline(pipeline_config)

        assert pipeline.config == pipeline_config
        assert pipeline.results == {}

    @patch("phylodater.core.pipeline.TreeValidator")
    @patch("phylodater.core.pipeline.CalibrationResolver")
    def test_pipeline_run_success(
        self,
        mock_resolver_class,
        mock_validator_class,
        pipeline_config,
        mock_tree,
        mock_calibrations,
        temp_dir,
    ):
        """测试流水线成功运行"""
        # 设置 mock
        mock_validator = Mock()
        mock_validator.full_validation.return_value = Mock(
            is_valid=True, errors=[], warnings=[]
        )
        mock_validator_class.return_value = mock_validator

        mock_resolver = Mock()
        mock_resolver.resolve_all.return_value = mock_calibrations
        mock_resolver_class.return_value = mock_resolver

        # 创建比对文件
        alignment_file = temp_dir / "alignment.fasta"
        alignment_file.write_text(
            ">human\nATCG\n>chimp\nATCG\n>mouse\nATCG\n>rat\nATCG\n"
        )

        pipeline = DatingPipeline(pipeline_config)

        # Mock 方法执行
        with patch.object(pipeline, "_execute_methods"):
            with patch.object(pipeline, "_generate_reports"):
                results = pipeline.run(mock_tree, alignment_file, mock_calibrations)

        assert isinstance(results, dict)

    def test_pipeline_config_defaults(self):
        """测试流水线配置默认值"""
        config = PipelineConfig()

        assert "mcmctree" in config.methods or "treepl" in config.methods
        assert config.threads >= 1
        assert config.skip_failed is True

    def test_pipeline_config_custom(self, temp_dir):
        """测试自定义流水线配置"""
        config = PipelineConfig(
            methods=["lsd2", "pathd8"],
            output_dir=temp_dir / "custom_output",
            threads=8,
            skip_failed=False,
        )

        assert config.methods == ["lsd2", "pathd8"]
        assert config.threads == 8
        assert config.skip_failed is False


class TestPipelineInputValidation:
    """流水线输入验证测试"""

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="pipeline_val_test_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def pipeline(self, temp_dir):
        config = PipelineConfig(output_dir=temp_dir / "output")
        return DatingPipeline(config)

    @pytest.fixture
    def mock_tree(self):
        tree = Mock(spec=PhylogeneticTree)
        tree.tip_names = ["human", "chimp", "mouse", "rat"]
        return tree

    @patch("phylodater.core.pipeline.TreeValidator")
    def test_input_validation_pass(
        self, mock_validator_class, pipeline, mock_tree, temp_dir
    ):
        """测试输入验证通过"""
        mock_validator = Mock()
        mock_validator.full_validation.return_value = Mock(
            is_valid=True, errors=[], warnings=[]
        )
        mock_validator_class.return_value = mock_validator

        alignment_file = temp_dir / "test.fasta"
        alignment_file.write_text(">human\nATCG\n")

        # 不应抛出异常
        pipeline._validate_inputs(mock_tree, alignment_file)

    @patch("phylodater.core.pipeline.TreeValidator")
    def test_input_validation_fail(
        self, mock_validator_class, pipeline, mock_tree, temp_dir
    ):
        """测试输入验证失败"""
        from phylodater.core.exceptions import PhyloDaterError

        mock_validator = Mock()
        mock_validator.full_validation.return_value = Mock(
            is_valid=False, errors=["Tree is not rooted"], warnings=[]
        )
        mock_validator_class.return_value = mock_validator

        alignment_file = temp_dir / "test.fasta"
        alignment_file.write_text(">human\nATCG\n")

        with pytest.raises(PhyloDaterError):
            pipeline._validate_inputs(mock_tree, alignment_file)


class TestPipelineCalibrationResolution:
    """流水线校准点解析测试"""

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="pipeline_calib_test_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def pipeline(self, temp_dir):
        config = PipelineConfig(output_dir=temp_dir / "output")
        return DatingPipeline(config)

    @patch("phylodater.core.pipeline.CalibrationResolver")
    def test_calibration_resolution(self, mock_resolver_class, pipeline, temp_dir):
        """测试校准点解析 - 预解析的校准点"""
        mock_tree = Mock(spec=PhylogeneticTree)
        mock_cal1 = Mock(spec=CalibrationPoint)
        mock_cal1.resolved_taxa = ["human", "chimp"]
        mock_cal1.age_constraint = None
        mock_cal2 = Mock(spec=CalibrationPoint)
        mock_cal2.resolved_taxa = ["mouse", "rat"]
        mock_cal2.age_constraint = None
        mock_calibrations = [mock_cal1, mock_cal2]

        mock_resolver = Mock()
        mock_resolver.validate_ancestry.return_value = None
        mock_resolver_class.return_value = mock_resolver

        result = pipeline._resolve_calibrations(mock_tree, mock_calibrations)

        # Pre-resolved calibrations should be returned as-is
        assert result == mock_calibrations
        mock_resolver.validate_ancestry.assert_called_once_with(mock_calibrations)

    def test_empty_calibrations_raises_calibration_error(self, pipeline, temp_dir):
        """测试空校准列表抛出清晰的 CalibrationError"""
        from phylodater.core.exceptions import CalibrationError

        mock_tree = Mock(spec=PhylogeneticTree)

        with pytest.raises(
            CalibrationError, match="No valid calibration points provided"
        ):
            pipeline._resolve_calibrations(mock_tree, [])


class TestPipelineResultsHandling:
    """流水线结果处理测试"""

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="pipeline_results_test_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    @pytest.fixture
    def pipeline(self, temp_dir):
        config = PipelineConfig(output_dir=temp_dir / "output")
        return DatingPipeline(config)

    def test_results_storage(self, pipeline):
        """测试结果存储"""
        mock_result = Mock(spec=DatingResult)
        mock_result.method_name = "mcmctree"

        pipeline.results["mcmctree"] = mock_result

        assert "mcmctree" in pipeline.results
        assert pipeline.results["mcmctree"] == mock_result

    def test_multiple_results(self, pipeline):
        """测试多个结果"""
        result1 = Mock(spec=DatingResult)
        result1.method_name = "mcmctree"

        result2 = Mock(spec=DatingResult)
        result2.method_name = "lsd2"

        pipeline.results["mcmctree"] = result1
        pipeline.results["lsd2"] = result2

        assert len(pipeline.results) == 2
        assert "mcmctree" in pipeline.results
        assert "lsd2" in pipeline.results


class TestPipelineEdgeCases:
    """流水线边界情况测试"""

    @pytest.fixture
    def temp_dir(self):
        temp_path = Path(tempfile.mkdtemp(prefix="pipeline_edge_test_"))
        yield temp_path
        shutil.rmtree(temp_path, ignore_errors=True)

    def test_empty_methods_list(self, temp_dir):
        """测试空方法列表"""
        with pytest.raises(ValueError, match="methods list cannot be empty"):
            PipelineConfig(methods=[], output_dir=temp_dir / "output")

    def test_single_method(self, temp_dir):
        """测试单方法"""
        config = PipelineConfig(methods=["mcmctree"], output_dir=temp_dir / "output")
        pipeline = DatingPipeline(config)

        assert pipeline.config.methods == ["mcmctree"]

    def test_output_dir_creation(self, temp_dir):
        """测试输出目录创建"""
        output_dir = temp_dir / "new_output_dir"
        assert not output_dir.exists()

        config = PipelineConfig(output_dir=output_dir)
        pipeline = DatingPipeline(config)

        # 配置已设置，但目录可能在运行时才创建
        assert pipeline.config.output_dir == output_dir
