"""The download of the model files (:mod:`t1sim.models`).

No network: the files come from a local folder through a ``file://`` URL,
which runs the same code as a download apart from the transport.
"""

import hashlib
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

_SRC = Path(__file__).resolve().parents[1] / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from t1sim import models  # noqa: E402


class TestTheDownload(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        tmp = Path(self._tmp.name)
        self.source = tmp / "source"
        self.target = tmp / "target"
        self.source.mkdir()
        files = {}
        for k in range(3):
            content = f"template {k}".encode() * 100
            (self.source / f"file{k}.bin").write_bytes(content)
            files[f"file{k}.bin"] = hashlib.sha256(content).hexdigest()
        patches = [
            mock.patch.object(models, "MODEL_FILES", files),
            mock.patch.dict(models.os.environ, {"T1SIM_MODEL_URL": self.source.as_uri()}),
        ]
        for patch in patches:
            patch.start()
            self.addCleanup(patch.stop)
        self.addCleanup(self._tmp.cleanup)

    def test_it_fetches_every_file_and_leaves_no_temporary_one(self):
        self.assertEqual(models.missing_models(self.target), list(models.MODEL_FILES))
        models.download_models(self.target, verbose=False)
        self.assertEqual(sorted(p.name for p in self.target.iterdir()),
                         sorted(models.MODEL_FILES))
        self.assertEqual(models.missing_models(self.target, verify=True), [])

    def test_a_file_that_is_there_is_not_fetched_again(self):
        models.download_models(self.target, verbose=False)
        with mock.patch.object(models, "_fetch", side_effect=AssertionError("fetched")):
            self.assertEqual(models.ensure_models(self.target, verbose=False), self.target)

    def test_only_the_missing_file_is_fetched(self):
        models.download_models(self.target, verbose=False)
        (self.target / "file1.bin").unlink()
        with mock.patch.object(models, "_fetch", wraps=models._fetch) as fetch:
            models.ensure_models(self.target, verbose=False)
        self.assertEqual([Path(c.args[0]).name for c in fetch.call_args_list], ["file1.bin"])

    def test_a_wrong_hash_is_an_error_and_keeps_nothing(self):
        (self.source / "file2.bin").write_bytes(b"corrupted")
        with self.assertRaisesRegex(RuntimeError, "file2.bin.*SHA-256"):
            models.download_models(self.target, verbose=False)
        self.assertEqual(sorted(p.name for p in self.target.iterdir()),
                         ["file0.bin", "file1.bin"])

    def test_a_changed_file_is_found_only_by_verifying(self):
        models.download_models(self.target, verbose=False)
        (self.target / "file0.bin").write_bytes(b"changed")
        self.assertEqual(models.missing_models(self.target), [])
        self.assertEqual(models.missing_models(self.target, verify=True), ["file0.bin"])

    def test_a_failed_download_names_the_url(self):
        (self.source / "file0.bin").unlink()
        with self.assertRaisesRegex(RuntimeError, "Failed to download .*file0.bin"):
            models.download_models(self.target, verbose=False)
        self.assertEqual(list(self.target.iterdir()), [])


class TestTheModelFolder(unittest.TestCase):
    def test_the_environment_variable_wins(self):
        with mock.patch.dict(models.os.environ, {"T1SIM_MODEL_DIR": "/some/where"}):
            self.assertEqual(models.model_dir(), Path("/some/where"))

    def test_a_source_checkout_keeps_them_in_the_repository(self):
        with mock.patch.dict(models.os.environ, clear=False) as env:
            env.pop("T1SIM_MODEL_DIR", None)
            self.assertEqual(models.model_dir(), _SRC.parent / "models")

    def test_an_installed_package_uses_the_user_cache(self):
        with mock.patch.dict(models.os.environ, clear=False) as env, \
             mock.patch.object(models, "_repository_root", return_value=None):
            env.pop("T1SIM_MODEL_DIR", None)
            self.assertEqual(models.model_dir(), models._user_cache_dir() / "models")

    def test_the_real_files_are_pinned_to_a_commit(self):
        self.assertEqual(len(models.MODEL_FILES), 5)
        self.assertIn(models.MODEL_SOURCE_COMMIT, models._DEFAULT_MODEL_URL)
        for digest in models.MODEL_FILES.values():
            self.assertRegex(digest, "^[0-9a-f]{64}$")


if __name__ == "__main__":
    unittest.main()
