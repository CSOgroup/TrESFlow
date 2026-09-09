import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def load_cleanup_module():
    spec = importlib.util.spec_from_file_location(
        "cleanup_intermediate_fastqs",
        REPO_ROOT / "bin/cleanup_intermediate_fastqs.py",
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


CLEANUP = load_cleanup_module()


class IntermediateFastqCleanupTests(unittest.TestCase):
    def test_producer_output_cleanup_unlinks_only_explicit_regular_fastqs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            work_dir = Path(tmpdir) / "work"
            task_dir = work_dir / "ab" / "task"
            task_dir.mkdir(parents=True)
            r1 = task_dir / "sample_R1.fastq"
            r2 = task_dir / "sample_R2.fq.gz"
            unrelated = task_dir / "keep.fastq"
            for path in (r1, r2, unrelated):
                path.write_bytes(b"reads")

            targets = CLEANUP.validate_targets(
                work_dir, [str(r1), str(r2)], "producer-output"
            )
            results = CLEANUP.unlink_targets(targets)

            self.assertEqual([status for _, status, _ in results], ["removed", "removed"])
            self.assertFalse(r1.exists())
            self.assertFalse(r2.exists())
            self.assertTrue(unrelated.exists())

    def test_staged_symlink_is_unlinked_without_touching_work_target(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            work_dir = Path(tmpdir) / "work"
            producer_dir = work_dir / "aa" / "producer"
            consumer_dir = work_dir / "bb" / "consumer"
            producer_dir.mkdir(parents=True)
            consumer_dir.mkdir(parents=True)
            producer = producer_dir / "sample.fastq"
            producer.write_bytes(b"reads")
            staged = consumer_dir / "sample.fastq"
            staged.symlink_to(producer)

            targets = CLEANUP.validate_targets(work_dir, [str(staged)], "staged-input")
            CLEANUP.unlink_targets(targets)

            self.assertFalse(staged.exists())
            self.assertFalse(staged.is_symlink())
            self.assertEqual(producer.read_bytes(), b"reads")

    def test_staged_copy_and_hardlink_entries_are_safe_to_unlink(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            work_dir = Path(tmpdir) / "work"
            producer_dir = work_dir / "aa" / "producer"
            consumer_dir = work_dir / "bb" / "consumer"
            producer_dir.mkdir(parents=True)
            consumer_dir.mkdir(parents=True)
            producer = producer_dir / "sample.fastq"
            producer.write_bytes(b"reads")
            staged_copy = consumer_dir / "copy.fastq"
            staged_copy.write_bytes(producer.read_bytes())
            staged_hardlink = consumer_dir / "hardlink.fastq"
            staged_hardlink.hardlink_to(producer)

            targets = CLEANUP.validate_targets(
                work_dir,
                [str(staged_copy), str(staged_hardlink)],
                "staged-input",
            )
            CLEANUP.unlink_targets(targets)

            self.assertFalse(staged_copy.exists())
            self.assertFalse(staged_hardlink.exists())
            self.assertEqual(producer.read_bytes(), b"reads")

    def test_outside_path_refusal_happens_before_any_unlink(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            work_dir = root / "work"
            task_dir = work_dir / "aa" / "task"
            task_dir.mkdir(parents=True)
            inside = task_dir / "inside.fastq"
            outside = root / "source.fastq"
            inside.write_bytes(b"inside")
            outside.write_bytes(b"source")

            with self.assertRaisesRegex(ValueError, "outside the controlled work directory"):
                CLEANUP.validate_targets(
                    work_dir,
                    [str(inside), str(outside)],
                    "producer-output",
                )

            self.assertEqual(inside.read_bytes(), b"inside")
            self.assertEqual(outside.read_bytes(), b"source")

    def test_staged_symlink_to_source_fastq_is_refused(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            work_dir = root / "work"
            task_dir = work_dir / "aa" / "task"
            task_dir.mkdir(parents=True)
            source = root / "source.fastq"
            source.write_bytes(b"source")
            staged = task_dir / "source.fastq"
            staged.symlink_to(source)

            with self.assertRaisesRegex(ValueError, "target is outside the work directory"):
                CLEANUP.validate_targets(work_dir, [str(staged)], "staged-input")

            self.assertTrue(staged.is_symlink())
            self.assertEqual(source.read_bytes(), b"source")

    def test_producer_cleanup_is_idempotent_but_rejects_symlinks_and_non_fastqs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            work_dir = Path(tmpdir) / "work"
            task_dir = work_dir / "aa" / "task"
            task_dir.mkdir(parents=True)
            missing = task_dir / "missing.fastq"

            targets = CLEANUP.validate_targets(
                work_dir, [str(missing)], "producer-output"
            )
            self.assertEqual(CLEANUP.unlink_targets(targets)[0][1], "already_missing")

            target = task_dir / "target.fastq"
            target.write_bytes(b"reads")
            symlink = task_dir / "link.fastq"
            symlink.symlink_to(target)
            with self.assertRaisesRegex(ValueError, "producer-output cleanup of a symlink"):
                CLEANUP.validate_targets(work_dir, [str(symlink)], "producer-output")

            non_fastq = task_dir / "metrics.tsv"
            non_fastq.write_text("keep\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "non-FASTQ"):
                CLEANUP.validate_targets(work_dir, [str(non_fastq)], "producer-output")

            self.assertTrue(target.exists())
            self.assertTrue(symlink.is_symlink())
            self.assertTrue(non_fastq.exists())

    def test_target_replacement_between_validation_and_unlink_is_refused(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            work_dir = Path(tmpdir) / "work"
            task_dir = work_dir / "aa" / "task"
            task_dir.mkdir(parents=True)
            target = task_dir / "sample.fastq"
            replacement = task_dir / "replacement.fastq"
            target.write_bytes(b"original")
            replacement.write_bytes(b"replacement")

            targets = CLEANUP.validate_targets(
                work_dir, [str(target)], "producer-output"
            )
            target.unlink()
            replacement.rename(target)

            with self.assertRaisesRegex(RuntimeError, "changed after validation"):
                CLEANUP.unlink_targets(targets)

            self.assertEqual(target.read_bytes(), b"replacement")


if __name__ == "__main__":
    unittest.main()
