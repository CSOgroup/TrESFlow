"""Configured runtime tools must not inherit an unrelated Python interpreter."""
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from test_multi_fastq_inputs import REPO, find_nextflow_jar


@pytest.mark.parametrize("renderer", ["host", "module"])
def test_env_shebang_uses_configured_runtime(tmp_path, renderer):
    java, jar = shutil.which("java"), find_nextflow_jar()
    if not java or jar is None:
        pytest.skip("Java and an offline Nextflow distribution are required")
    runtime = tmp_path / "configured runtime's prefix"
    bindir = runtime / "bin"
    bindir.mkdir(parents=True)
    (bindir / "python").symlink_to(sys.executable)
    probe = bindir / "probe"
    probe.write_text('#!/usr/bin/env python\nimport json, os, sys\n'
                     'print(json.dumps([sys.executable, os.environ["PATH"]]))\n')
    probe.chmod(0o755)
    foreign = tmp_path / "foreign"
    foreign.mkdir()
    (foreign / "python").write_text("#!/bin/sh\nexit 97\n")
    (foreign / "python").chmod(0o755)
    tmpdir = tmp_path / "task temporary files"
    source = REPO / ("lib/RuntimeSupport.groovy" if renderer == "host" else
                     "modules/local/runtime_support/main.nf")
    groovy = '''
def parameters = [runtime_env_prefix: args[2], runtime_tmpdir: args[3]]
if (args[0] == 'host') {
    def support = new GroovyClassLoader().parseClass(new File(args[1]))
    print support.shellExports(parameters)
} else {
    def support = new GroovyShell().parse(new File(args[1]))
    print support.runtimeShellExports(parameters)
}
'''
    rendered = subprocess.run([java, "-cp", str(jar), "groovy.ui.GroovyMain", "-e", groovy, "--",
                               renderer, str(source), str(runtime), str(tmpdir)],
                              check=True, text=True, capture_output=True, timeout=30).stdout
    inherited = f"{foreign}:/usr/bin:/bin"
    result = subprocess.run(["/bin/bash", "-euo", "pipefail", "-c",
                             rendered + '\n"$RUNTIME_BIN_DIR/probe"'],
                            env=dict(os.environ, PATH=inherited), text=True,
                            capture_output=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    interpreter, path = json.loads(result.stdout)
    assert Path(interpreter).parent == bindir
    assert path == f"{bindir}:{inherited}"
    assert tmpdir.is_dir()
