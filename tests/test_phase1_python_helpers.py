import ast
import hashlib
import json
from pathlib import Path
import re
import sys
import unittest


REPO = Path(__file__).resolve().parents[1]
MODULES = (
    REPO / "modules/local/barcode_gate_metrics/main.nf",
    REPO / "modules/local/tres_report_html/main.nf",
)
CONTAINER = (
    "docker.io/library/python:3.12.13-bookworm@"
    "sha256:3cd9086bdb30f7c9bc08a3fa621d9842e0d3f6f9291aeb4677e0547817c10b12"
)


class Phase1PythonHelperArchitectureTests(unittest.TestCase):
    def test_synthetic_fixture_manifest_is_complete_and_valid(self):
        fixture_root = REPO / "tests/fixtures/python_helpers"
        manifest = json.loads(
            (fixture_root / "manifest.json").read_text(encoding="utf-8")
        )
        self.assertEqual(manifest["asset_type"], "synthetic-test-only")
        self.assertEqual(manifest["license"], "MIT")
        declared = set(manifest["files"])
        actual = {
            path.relative_to(fixture_root).as_posix()
            for path in fixture_root.rglob("*")
            if path.is_file() and path.name not in {"README.md", "manifest.json"}
        }
        self.assertEqual(declared, actual)
        for relative_path, expected_sha256 in manifest["files"].items():
            observed = hashlib.sha256(
                (fixture_root / relative_path).read_bytes()
            ).hexdigest()
            self.assertEqual(observed, expected_sha256, relative_path)

    def test_runtime_imports_are_standard_library_or_staged_repository_modules(self):
        files = [
            REPO / "bin/write_barcode_gate_metrics.py",
            REPO / "bin/run_split_reads_dna.py",
            REPO / "bin/tresflow_fastq_utils.py",
            REPO / "bin/render_tres_report.py",
            *sorted((REPO / "lib/tresflow_qc").glob("*.py")),
        ]
        repository_modules = {
            "run_split_reads_dna",
            "tresflow_fastq_utils",
            "tresflow_qc",
        }
        standard_library = set(getattr(sys, "stdlib_module_names", ())) or {
            "__future__", "argparse", "base64", "collections", "csv",
            "dataclasses", "datetime", "gzip", "hashlib", "html", "itertools",
            "json", "math", "os", "pathlib", "re", "shutil", "subprocess",
            "sys", "tempfile", "time", "typing", "unicodedata",
        }
        unexpected = set()
        for path in files:
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    imported = {alias.name.split(".", 1)[0] for alias in node.names}
                elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                    imported = {node.module.split(".", 1)[0]}
                else:
                    continue
                unexpected.update(imported - standard_library - repository_modules)
        self.assertEqual(unexpected, set())

    def test_conda_environment_contains_only_exact_python(self):
        environment = (
            REPO / "modules/local/python_helpers/environment.yml"
        ).read_text(encoding="utf-8")
        self.assertIn("- python=3.12.13", environment)
        dependency_lines = [
            line.strip()
            for line in environment.splitlines()
            if re.match(r"^\s{2}-\s", line) and "conda-forge" not in line and "nodefaults" not in line
        ]
        self.assertEqual(dependency_lines, ["- python=3.12.13"])

    def test_both_processes_share_the_digest_pinned_multiarch_image(self):
        for module in MODULES:
            text = module.read_text(encoding="utf-8")
            self.assertIn(CONTAINER, text)
            self.assertIn(f"docker://{CONTAINER}", text)
            self.assertIn('conda "${moduleDir}/../python_helpers/environment.yml"', text)
            self.assertNotRegex(text, r"python:[^'\"]+@(?!sha256:)")

    def test_task_commands_are_decoupled_from_host_runtime_and_project_paths(self):
        for module in MODULES:
            text = module.read_text(encoding="utf-8")
            self.assertIn("python3", text)
            self.assertNotIn("PYTHON3_BIN", text)
            self.assertNotIn("runtimeShellExports", text)
            self.assertNotIn("runtime_env_prefix", text)
            self.assertNotIn("/home/annan", text)
            self.assertNotIn("${projectDir}", text)
        self.assertIn('path helperScripts, stageAs: "helper/bin/*"', MODULES[0].read_text())
        report_module = MODULES[1].read_text()
        self.assertIn('path reportScript, stageAs: "helper/bin/render_tres_report.py"', report_module)
        self.assertIn('path reportPackage, stageAs: "helper/lib/tresflow_qc"', report_module)

    def test_workflows_supply_every_repository_helper_as_a_path_input(self):
        rna = (REPO / "subworkflows/local/rna_core/main.nf").read_text()
        dna = (REPO / "subworkflows/local/dna_core/main.nf").read_text()
        workflow = (REPO / "workflows/treseq.nf").read_text()
        for text in (rna, dna):
            self.assertIn("write_barcode_gate_metrics.py", text)
            self.assertIn("run_split_reads_dna.py", text)
            self.assertIn("tresflow_fastq_utils.py", text)
            self.assertIn("BARCODE_GATE_METRICS(ch_", text)
        self.assertIn("bin/render_tres_report.py", workflow)
        self.assertIn("lib/tresflow_qc", workflow)
        self.assertIn("TRES_REPORT_HTML(\n        ch_tres_report_input,", workflow)

    def test_process_test_uses_an_unresolvable_legacy_environment(self):
        harness = (REPO / "tests/phase1_python_helpers.nf").read_text()
        self.assertIn("phase1-must-not-be-used", harness)
        self.assertIn("BARCODE_GATE_METRICS(", harness)
        self.assertIn("TRES_REPORT_HTML(", harness)

    def test_execution_profiles_are_explicitly_test_only(self):
        config = (REPO / "tests/phase1_python_helpers.config").read_text()
        self.assertIn("not production pipeline profiles", config)
        self.assertIn("phase1_helpers_docker", config)
        self.assertIn("docker.enabled = true", config)
        self.assertIn("phase1_helpers_conda", config)
        self.assertIn("conda.useMicromamba = true", config)
        self.assertIn("phase1_helpers_apptainer", config)
        self.assertIn("apptainer.autoMounts = true", config)


if __name__ == "__main__":
    unittest.main()
