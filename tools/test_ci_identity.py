#!/usr/bin/env python3
"""Challenge version and source/PR/run binding without building application code."""
import json
from pathlib import Path
import tempfile
import unittest
from ci_identity import ci_metadata


class IdentityTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        (self.root / "Cargo.toml").write_text('[package]\nversion = "0.1.0-dev"\n')
        self.sha = "a" * 40
        self.manifest = {"source_version": "0.1.0-dev", "git_head": self.sha, "build_id": "b" * 64}
        self.discovery = {"version": "0.1.0-dev", "build": dict(self.manifest)}
        self.env = {"GITHUB_SHA": self.sha, "GITHUB_RUN_ID": "123", "GITHUB_RUN_ATTEMPT": "2",
                    "GITHUB_EVENT_NAME": "pull_request"}
        self.event = {"pull_request": {"head": {"sha": "c" * 40}, "base": {"sha": "d" * 40}}}

    def value(self):
        return ci_metadata(self.root, self.manifest, self.discovery, self.env, self.event, self.sha)

    def test_event_payload_mismatch_rejected(self):
        self.env["GITHUB_EVENT_NAME"] = "push"
        with self.assertRaises(ValueError): self.value()

    def test_invalid_pr_head_rejected(self):
        self.event["pull_request"]["head"]["sha"] = "not-a-sha"
        with self.assertRaises(ValueError): self.value()

    def test_repository_lock_matches_single_source(self):
        from ci_identity import source_version
        root = Path(__file__).resolve().parents[1]
        block = (root / "Cargo.lock").read_text().split('name = "irreducible"', 1)[1].split("[[package]]", 1)[0]
        self.assertIn(f'version = "{source_version(root)}"', block)

    def test_pr_merge_distinct_from_head(self):
        value = self.value()
        self.assertNotEqual(value["checked_out_sha"], value["pr_head_sha"])
        self.assertEqual(value["checkout_role"], "synthetic_pr_merge")
        self.assertEqual(value["build_id"], self.manifest["build_id"])

    def test_attempt_distinguishes_execution_without_changing_build(self):
        old = self.value()
        self.env["GITHUB_RUN_ATTEMPT"] = "3"
        new = self.value()
        self.assertNotEqual(old["ci_identity"], new["ci_identity"])
        self.assertEqual(old["build_id"], new["build_id"])

    def test_main_has_no_pr_head(self):
        self.event = {}
        self.env["GITHUB_EVENT_NAME"] = "push"
        self.assertIsNone(self.value()["pr_head_sha"])

    def test_version_mismatch_rejected(self):
        self.discovery["version"] = "0.1.0"
        with self.assertRaises(ValueError): self.value()

    def test_checkout_mismatch_rejected(self):
        self.env["GITHUB_SHA"] = "e" * 40
        with self.assertRaises(ValueError): self.value()

    def test_build_mismatch_rejected(self):
        self.discovery["build"]["build_id"] = "f" * 64
        with self.assertRaises(ValueError): self.value()


if __name__ == "__main__":
    unittest.main()
