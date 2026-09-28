#!/usr/bin/env python3
"""Regression tests for packaging failure gates and publication dependencies."""
import argparse
import importlib.util
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import yaml

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("release", ROOT / "scripts/release/release.py")
release = importlib.util.module_from_spec(spec)
spec.loader.exec_module(release)


class PackagingGates(unittest.TestCase):
    def args(self, temp):
        return argparse.Namespace(version="v" + release.suite_version(), jobs=2,
                                  out_dir=Path(temp) / "packages", test_dir=Path(temp) / "tests")

    def test_version_mismatch_is_rejected(self):
        with self.assertRaises(ValueError):
            release.version_info("v999.0.0")

    def test_prerelease_debian_ordering(self):
        version = release.suite_version()
        self.assertEqual(release.version_info(f"v{version}-rc1"),
                         (f"v{version}-rc1", version, f"{version}~rc1"))

    def test_stale_output_is_rejected_before_build(self):
        with tempfile.TemporaryDirectory() as temp, patch.object(release, "run") as run:
            args = self.args(temp)
            args.out_dir.mkdir()
            (args.out_dir / "old.deb").write_bytes(b"stale")
            with self.assertRaises(ValueError):
                release.build(args)
            run.assert_not_called()

    def test_build_failure_blocks_packaging(self):
        with tempfile.TemporaryDirectory() as temp, \
                patch.object(release, "run", side_effect=subprocess.CalledProcessError(1, "make")), \
                patch.object(release, "stage_release") as stage:
            with self.assertRaises(subprocess.CalledProcessError):
                release.build(self.args(temp))
            stage.assert_not_called()

    def test_test_failure_blocks_packaging(self):
        with tempfile.TemporaryDirectory() as temp, \
                patch.object(release, "run", side_effect=[None, subprocess.CalledProcessError(1, "tests")]), \
                patch.object(release, "stage_release") as stage:
            with self.assertRaises(subprocess.CalledProcessError):
                release.build(self.args(temp))
            stage.assert_not_called()

    def test_sdk_failure_blocks_artifact_creation(self):
        with tempfile.TemporaryDirectory() as temp, \
                patch.object(release, "run", side_effect=[None, None, None,
                    subprocess.CalledProcessError(1, "sdk")]), \
                patch.object(release, "stage_release", return_value={"source_revision": "a" * 40}), \
                patch.object(release, "build_tarball") as tarball, \
                patch.object(release, "build_deb") as deb:
            with self.assertRaises(subprocess.CalledProcessError):
                release.build(self.args(temp))
            tarball.assert_not_called()
            deb.assert_not_called()


class PublicationGates(unittest.TestCase):
    def setUp(self):
        self.workflow = yaml.safe_load((ROOT / ".github/workflows/release.yml").read_text())
        self.jobs = self.workflow["jobs"]

    def ancestors(self, job):
        needs = self.jobs[job].get("needs", [])
        needs = [needs] if isinstance(needs, str) else needs
        return set(needs).union(*(self.ancestors(parent) for parent in needs))

    def test_all_publishers_require_package_and_runtime_gates(self):
        for job in ("docker-image", "github-release"):
            self.assertTrue({"build-packages", "runtime-checks", "source-package"} <= self.ancestors(job))
            self.assertNotIn("always()", self.jobs[job].get("if", ""))
        self.assertIn("docker-image", self.ancestors("github-release"))

    def test_manual_runs_cannot_publish(self):
        self.assertIn("github.event_name == 'push'", self.jobs["github-release"]["if"])
        pushes = [s for s in self.jobs["docker-image"]["steps"] if "docker push" in s.get("run", "")]
        self.assertEqual(len(pushes), 1)
        self.assertIn("github.event_name == 'push'", pushes[0]["if"])
        self.assertIn("vars.RELEASE_PUSH_IMAGE == 'true'", pushes[0]["if"])

    def test_image_is_tested_before_push_without_rebuilding(self):
        steps = self.jobs["docker-image"]["steps"]
        build = next(i for i, step in enumerate(steps) if step.get("uses", "").startswith("docker/build-push-action"))
        smoke = next(i for i, step in enumerate(steps) if "run_container_smoke.sh" in step.get("run", ""))
        push = next(i for i, step in enumerate(steps) if "docker push" in step.get("run", ""))
        self.assertLess(build, smoke)
        self.assertLess(smoke, push)
        self.assertIs(steps[build]["with"]["push"], False)
        self.assertIs(steps[build]["with"]["load"], True)
        self.assertIn("local/chromap-suite:release-check", steps[smoke]["run"])
        self.assertIn("docker tag local/chromap-suite:release-check", steps[push]["run"])

    def test_runtime_matrix_covers_both_baselines_and_forward_compatibility(self):
        rows = self.jobs["runtime-checks"]["strategy"]["matrix"]["include"]
        self.assertEqual({(r["baseline"], r["runtime"]) for r in rows},
                         {("ubuntu22.04", "22.04"), ("ubuntu22.04", "24.04"), ("ubuntu24.04", "24.04")})


if __name__ == "__main__":
    unittest.main()
