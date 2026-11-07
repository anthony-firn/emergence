"""WASM reuse must reject mismatched binaries and modified Rust sources."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import build


class WasmReuse(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.package = Path(self.temp.name)
        (self.package / 'emergence_web.wasm').write_bytes(b'known-release')
        self.meta = {'source': 'https://github.com/anthony-firn/emergence.git',
                     'commit': 'a' * 40, 'wasm_sha256': hashlib.sha256(b'known-release').hexdigest()}
        self.write()

    def write(self):
        (self.package / 'BUILD.json').write_text(json.dumps(self.meta))

    def test_matching_release_preserves_original_wasm_commit(self):
        self.meta['wasm_commit'] = 'b' * 40
        self.write()
        with patch('build.subprocess.run') as diff, patch('build.subprocess.check_output', return_value=b''):
            wasm, commit = build.reusable_wasm(self.package)
            self.assertEqual(wasm.read_bytes(), b'known-release')
            self.assertEqual(commit, 'b' * 40)
            self.assertIn('b' * 40, diff.call_args.args[0])

    def test_corrupt_binary_rejected_before_git_check(self):
        (self.package / 'emergence_web.wasm').write_bytes(b'corrupt')
        with self.assertRaisesRegex(ValueError, 'checksum'):
            build.reusable_wasm(self.package)

    def test_dirty_release_rejected(self):
        self.meta['modified'] = True
        self.write()
        with self.assertRaisesRegex(ValueError, 'clean release'):
            build.reusable_wasm(self.package)

    def test_local_rust_edits_rejected(self):
        with patch('build.subprocess.run'), patch('build.subprocess.check_output', return_value=b' M fluid/src.rs'):
            with self.assertRaisesRegex(ValueError, 'local changes'):
                build.reusable_wasm(self.package)

    def test_changed_rust_between_commits_rejected(self):
        with patch('build.subprocess.run', side_effect=build.subprocess.CalledProcessError(1, 'git diff')):
            with self.assertRaises(build.subprocess.CalledProcessError):
                build.reusable_wasm(self.package)


if __name__ == '__main__':
    unittest.main()
