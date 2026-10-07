#!/usr/bin/env python3
"""Exercise routing, snapshots, structure checks and the model-pin hook in temporary repositories."""

import json
import os
import re
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
import unittest.mock

import vareffect_harness as harness

SOURCE = Path(__file__).resolve().parents[1]
CLEAN_ENV = {k: v for k, v in os.environ.items() if not k.startswith(("CLAUDE_CODE_SUBAGENT_MODEL", "ANTHROPIC_DEFAULT_OPUS_MODEL"))}


class RepoTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.git("init", "-q")
        self.git("config", "user.name", "Fixture")
        self.git("config", "user.email", "fixture@example.invalid")
        self.git("config", "core.excludesFile", os.devnull)
        self.write("README.md", "base\n")
        self.git("add", ".")
        self.git("commit", "-qm", "fixture base")
        self.base = self.git("rev-parse", "HEAD").strip()

    def git(self, *args):
        return subprocess.check_output(["git", "-C", str(self.root), *args], text=True)

    def write(self, path, text):
        p = self.root / path
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
        return p

    def routed(self, verbose=True):
        result = harness.route(self.root, self.base, verbose=verbose)
        return result, {f["path"]: f for f in result["files"]}

    # -- routing ---------------------------------------------------------------------------------

    def test_annotation_paths_route_the_specialist(self):
        for path in ["vareffect/src/consequence/snv.rs", "vareffect/src/hgvs_p.rs", "vareffect/tests/locate_tests.rs",
                     "vareffect-cli/src/builders/transcript_models/mod.rs", "vareffect-cli/src/csq.rs",
                     "vareffect-cli/vareffect_build.toml", "vareffect/VEP_DIVERGENCES.md"]:
            with self.subTest(path=path):
                rules = harness.route_path(path)
                self.assertIn(harness.SPECIALIST, rules["reviewers"])
                self.assertEqual(rules["tier_floor"], "T2")
        self.assertIn("concordance", harness.route_path("vareffect/src/locate/indel.rs")["verification"])
        self.assertIn("concordance", harness.route_path("vareffect/tests/data/expected.tsv")["verification"])

    def test_scaffolding_and_style_references_stay_engineering_only(self):
        for path in ["vareffect-cli/src/init.rs", "vareffect-cli/src/check.rs", "README.md", "vareffect/README.md",
                     "CHANGELOG.md", ".claude/skills/idiomatic-rust/SKILL.md", "Makefile", "deny.toml"]:
            with self.subTest(path=path):
                rules = harness.route_path(path)
                self.assertEqual(rules["reviewers"], [harness.ENGINEERING])
                self.assertEqual(rules["tier_floor"], "T1")

    def test_harness_authority_routes_the_specialist(self):
        for path in [".claude/INVARIANTS.md", ".claude/settings.json", ".claude/hooks/pin-reviewer-models.py",
                     ".claude/agents/vareffect-reviewer.md", ".claude/skills/vareffect-review/SKILL.md",
                     "scripts/vareffect_harness.py", "CLAUDE.md"]:
            with self.subTest(path=path):
                rules = harness.route_path(path)
                self.assertIn(harness.SPECIALIST, rules["reviewers"])
                self.assertIn("harness", rules["verification"])

    def test_verification_categories(self):
        self.assertIn("store-rebuild", harness.route_path("vareffect-cli/src/builders/reference_genome.rs")["verification"])
        self.assertIn("store-rebuild", harness.route_path("vareffect/src/types.rs")["verification"])
        self.assertIn("public-contract", harness.route_path("vareffect/src/lib.rs")["verification"])
        release = harness.route_path(".github/workflows/release.yml")
        self.assertEqual((release["tier_floor"], "release" in release["verification"]), ("T2", True))
        self.assertIn("supply-chain", harness.route_path("Cargo.lock")["verification"])
        self.assertIn("build-config", harness.route_path(".gitignore")["verification"])

    def test_new_ground_truth_under_tests_data_is_routed_not_ignored(self):
        shutil.copy2(SOURCE / ".gitignore", self.root / ".gitignore")
        self.write("vareffect/tests/data/vep_ground_truth.tsv", "chrom\tpos\n")
        self.write("vareffect/tests/data/run_mismatches.log", "dump\n")
        self.write("data/vareffect/GRCh38.bin", "genome\n")
        _, files = self.routed()
        self.assertIn(harness.SPECIALIST, files["vareffect/tests/data/vep_ground_truth.tsv"]["reviewers"])
        self.assertNotIn("vareffect/tests/data/run_mismatches.log", files)
        self.assertNotIn("data/vareffect/GRCh38.bin", files)

    def test_route_covers_staged_unstaged_untracked_and_deleted(self):
        self.write("vareffect/src/codon.rs", "fn a() {}\n")
        self.write("vareffect-cli/src/init.rs", "fn b() {}\n")
        self.git("add", ".")
        self.git("commit", "-qm", "sources")
        self.base = self.git("rev-parse", "HEAD").strip()
        self.write("staged.txt", "staged\n")
        self.git("add", "staged.txt")
        (self.root / "vareffect/src/codon.rs").unlink()
        self.write("vareffect-cli/src/init.rs", "fn b() { }\n")
        self.write("untracked.md", "new\n")
        result, files = self.routed()
        self.assertEqual(set(files), {"staged.txt", "vareffect/src/codon.rs", "vareffect-cli/src/init.rs", "untracked.md"})
        self.assertIn(harness.SPECIALIST, files["vareffect/src/codon.rs"]["reviewers"])
        self.assertEqual(result["tier_floor"], "T2")

    def test_unsafe_code_raises_the_floor_by_content_even_when_removed(self):
        self.write("vareffect-cli/src/init.rs", "fn f() { unsafe { core::hint::unreachable_unchecked() } }\n")
        self.git("add", ".")
        self.git("commit", "-qm", "unsafe")
        self.base = self.git("rev-parse", "HEAD").strip()
        self.write("vareffect-cli/src/init.rs", "fn f() {}\n")
        _, files = self.routed()
        self.assertIn("unsafe", files["vareffect-cli/src/init.rs"]["verification"])
        self.assertEqual(files["vareffect-cli/src/init.rs"]["tier_floor"], "T2")

    def test_every_unsafe_form_is_detected(self):
        for text in ["unsafe { x() }", "unsafe fn f() {}", "unsafe impl Send for X {}", "unsafe trait T {}",
                     'unsafe extern "C" {}', "#[unsafe(no_mangle)]"]:
            self.assertTrue(harness.UNSAFE.search(text), text)
        for text in ["// not unsafe here", "let unsafe_count = 1;", "fn is_unsafe() {}"]:
            self.assertFalse(harness.UNSAFE.search(text), text)

    def test_ci_and_dependencies_raise_the_floor_or_require_rebuilds(self):
        self.assertEqual(harness.route_path(".github/workflows/ci.yml")["tier_floor"], "T2")
        for path in ["Cargo.toml", "Cargo.lock", "vareffect/Cargo.toml", "vareffect-cli/Cargo.toml"]:
            rules = harness.route_path(path)
            self.assertIn("store-rebuild", rules["verification"], path)
            self.assertIn(harness.SPECIALIST, rules["reviewers"], path)
            self.assertEqual(rules["tier_floor"], "T2", path)
        self.assertIn("store-rebuild", harness.route_path("vareffect/src/chrom.rs")["verification"])
        self.assertIn("public-contract", harness.route_path("vareffect/src/consequence/mod.rs")["verification"])

    def test_staged_rename_keeps_old_and_new_paths(self):
        self.write("vareffect/src/codon.rs", "fn a() {}\n")
        self.git("add", ".")
        self.git("commit", "-qm", "codon")
        self.base = self.git("rev-parse", "HEAD").strip()
        self.git("mv", "vareffect/src/codon.rs", "neutral.txt")
        _, files = self.routed()
        self.assertIn("vareffect/src/codon.rs", files)
        self.assertIn("neutral.txt", files)
        self.assertIn(harness.SPECIALIST, files["vareffect/src/codon.rs"]["reviewers"])

    def test_explicit_comparisons_have_distinct_semantics(self):
        self.git("checkout", "-qb", "side")
        self.write("side.txt", "side\n")
        self.git("add", ".")
        self.git("commit", "-qm", "side")
        self.git("checkout", "-q", "-")
        self.write("main.txt", "main\n")
        self.git("add", ".")
        self.git("commit", "-qm", "main")
        main = self.git("rev-parse", "HEAD").strip()
        two = harness.route(self.root, f"{main}..side", verbose=True)
        three = harness.route(self.root, f"{main}...side", verbose=True)
        self.assertEqual(two["comparison"]["start"], main)
        self.assertEqual(three["comparison"]["start"], self.base)
        self.assertIn("main.txt", {f["path"] for f in two["files"]})
        self.assertNotIn("main.txt", {f["path"] for f in three["files"] if any(c["source"] == "comparison" for c in f["changes"])})

    def test_unrelated_histories_fail_merge_base(self):
        default = self.git("rev-parse", "--abbrev-ref", "HEAD").strip()
        self.git("checkout", "-q", "--orphan", "other")
        self.git("commit", "-q", "--allow-empty", "-m", "orphan")
        with self.assertRaises(harness.HarnessError):
            harness.route(self.root, f"other...{default}")

    def test_commits_after_the_comparison_end_still_route(self):
        self.write("vareffect/src/codon.rs", "fn a() {}\n")
        self.git("add", ".")
        self.git("commit", "-qm", "later")
        result = harness.route(self.root, f"{self.base}..{self.base}", verbose=True)
        files = {f["path"]: f for f in result["files"]}
        self.assertEqual(files["vareffect/src/codon.rs"]["changes"][0]["source"], "after-comparison")

    def test_unreadable_rust_routes_the_unsafe_check(self):
        rules = harness.route_path("vareffect-cli/src/init.rs")
        with unittest.mock.patch.object(harness, "content_texts", return_value=([], True)):
            routed = harness.route_content(self.root, "vareffect-cli/src/init.rs", rules)
        self.assertEqual((routed["tier_floor"], "unsafe" in routed["verification"]), ("T2", True))

    def test_no_changes_is_not_a_review_pass(self):
        result = harness.route(self.root, self.base)
        self.assertEqual((result["file_count"], result["minimum_reviewers"], result["tier_floor"]), (0, [], "none"))

    def test_invalid_comparison_fails_instead_of_falling_back(self):
        for expression in ["nope", "HEAD..", "...HEAD", "-x"]:
            with self.subTest(expression=expression), self.assertRaises(harness.HarnessError):
                harness.route(self.root, expression)

    # -- snapshot --------------------------------------------------------------------------------

    def test_snapshot_identifies_worktree_bytes_without_touching_index(self):
        index = self.root / ".git/index"
        before = index.read_bytes(), self.git("status", "--porcelain=v1")
        self.assertTrue(harness.snapshot(self.root)["matches_head"])
        self.write("untracked.txt", "one\n")
        first = harness.snapshot(self.root)
        self.assertFalse(first["matches_head"])
        self.assertEqual(first, harness.snapshot(self.root))
        self.write("untracked.txt", "two\n")
        self.assertNotEqual(first["tree"], harness.snapshot(self.root)["tree"])
        self.write(".gitignore", "ignored.txt\n")
        with_ignore = harness.snapshot(self.root)["tree"]
        self.write("ignored.txt", "local\n")
        self.assertEqual(with_ignore, harness.snapshot(self.root)["tree"])
        for name in ["untracked.txt", "ignored.txt", ".gitignore"]:
            (self.root / name).unlink()
        self.assertEqual((index.read_bytes(), self.git("status", "--porcelain=v1")), before)

    def test_snapshot_covers_modified_staged_and_deleted_tracked_files(self):
        self.write("tracked.txt", "v1\n")
        self.git("add", ".")
        self.git("commit", "-qm", "tracked")
        base = harness.snapshot(self.root)["tree"]
        self.write("tracked.txt", "v2\n")
        modified = harness.snapshot(self.root)["tree"]
        self.assertNotEqual(base, modified)
        self.git("add", "tracked.txt")
        self.write("tracked.txt", "v3\n")
        self.assertNotEqual(modified, harness.snapshot(self.root)["tree"])
        (self.root / "tracked.txt").unlink()
        self.assertNotIn("tracked.txt", self.git("ls-tree", "--name-only", harness.snapshot(self.root)["tree"]))

    def test_snapshot_writes_exact_diff_including_untracked(self):
        self.write("new.txt", "fresh\n")
        out = self.root.parent / f"{self.root.name}-review.diff"
        self.addCleanup(lambda: out.unlink(missing_ok=True))
        result = harness.snapshot(self.root, self.base, out)
        self.assertIn("+fresh", out.read_text())
        self.assertEqual(result["diff"]["from"], self.base)

    def test_snapshot_refuses_secrets_and_sequencing_data(self):
        names = [".env", ".env.local", "deploy/x.env", ".envrc", "id_ed25519", "gcp-credentials.json", ".npmrc",
                 ".git-credentials", ".ssh/deploy_key", ".cargo/credentials.toml", ".aws/credentials",
                 ".claude/settings.local.json", "sample.vcf.gz", "cohort.vcf.bgz", "calls.gvcf", "reads.bam",
                 "reads.cram", "run.fastq.gz", "key.pem"]
        for name in names:
            with self.subTest(name=name):
                self.write(name, "sensitive\n")
                with self.assertRaises(harness.HarnessError):
                    harness.snapshot(self.root)
                self.assertEqual(harness.snapshot(self.root, allow=[name])["allowed_sensitive"], [name])
                (self.root / name).unlink()
        self.write(".env.example", "TEMPLATE=\n")
        self.write("vareffect/tests/data/expected.tsv", "chr1\t1\tA\tG\n")
        self.write("vareffect-cli/vareffect_build.toml", "[x]\n")
        self.assertEqual(harness.snapshot(self.root)["allowed_sensitive"], [])

    def test_snapshot_refuses_staged_and_intent_to_add_sensitive_files_before_storing(self):
        self.write("x.env", "SECRET=1\n")
        self.git("add", "-N", "x.env")
        before = self.git("count-objects", "-v")
        with self.assertRaises(harness.HarnessError):
            harness.snapshot(self.root)
        self.assertEqual(before, self.git("count-objects", "-v"))
        (self.root / "x.env").unlink()
        self.git("rm", "-q", "--cached", "x.env")
        self.write("sample.vcf", "##fileformat=VCFv4.2\n")
        self.git("add", "sample.vcf")
        out = self.root.parent / f"{self.root.name}-sensitive.diff"
        self.addCleanup(lambda: out.unlink(missing_ok=True))
        with self.assertRaises(harness.HarnessError):
            harness.snapshot(self.root, self.base, out)
        self.assertFalse(out.exists())

    def test_snapshot_diff_omits_the_old_bytes_of_deleted_sensitive_files(self):
        self.write("fixture.vcf", "".join(f"chr1\t{n}\t.\tA\tG\tsample-{n}\n" for n in range(40)))
        self.git("add", ".")
        self.git("commit", "-qm", "data committed by mistake")
        base = self.git("rev-parse", "HEAD").strip()
        out = self.root.parent / f"{self.root.name}-deleted.diff"
        self.addCleanup(lambda: out.unlink(missing_ok=True))
        self.git("mv", "fixture.vcf", "notes.txt")
        with self.assertRaisesRegex(harness.HarnessError, "renamed"):
            harness.snapshot(self.root, base, out)
        self.git("rm", "-q", "--cached", "notes.txt")
        (self.root / "notes.txt").unlink()
        for env in ({}, {"GIT_LITERAL_PATHSPECS": "1"}):
            with self.subTest(env=env), unittest.mock.patch.dict(os.environ, env):
                result = harness.snapshot(self.root, base, out)
                self.assertEqual(result["withheld_deletions"], ["fixture.vcf"])
                self.assertIn(b"deleted file mode", out.read_bytes())
                self.assertNotIn(b"chr1\t", out.read_bytes())

    def test_snapshot_refuses_modified_tracked_sensitive_files_before_storing_them(self):
        self.write("fixture.vcf", "##synthetic\n")
        self.git("add", ".")
        self.git("commit", "-qm", "synthetic fixture")
        self.write("fixture.vcf", "##changed\n")
        before = self.git("count-objects", "-v")
        with self.assertRaises(harness.HarnessError):
            harness.snapshot(self.root)
        self.assertEqual(before, self.git("count-objects", "-v"))
        self.assertEqual(harness.snapshot(self.root, allow=["fixture.vcf"])["allowed_sensitive"], ["fixture.vcf"])

    def test_snapshot_sees_edits_as_recent_as_the_index_write(self):
        self.git("config", "core.trustctime", "false")
        path = self.write("tracked.txt", "v1\n")
        past = (path.stat().st_mtime_ns // 10**9 - 10) * 10**9
        os.utime(path, ns=(past, past))
        self.git("add", ".")
        self.git("commit", "-qm", "tracked")
        clean = harness.snapshot(self.root)["tree"]
        path.write_text("v2\n")
        os.utime(path, ns=(past, past))
        os.utime(self.root / ".git" / "index", ns=(past, past))
        self.assertNotEqual(clean, harness.snapshot(self.root)["tree"])

    def test_snapshot_keeps_files_a_sparse_checkout_leaves_out(self):
        self.write("kept.txt", "kept\n")
        self.git("add", ".")
        self.git("commit", "-qm", "kept")
        self.git("config", "core.sparseCheckout", "true")
        self.git("update-index", "--skip-worktree", "kept.txt")
        (self.root / "kept.txt").unlink()
        self.assertTrue(harness.snapshot(self.root)["matches_head"])

    def test_snapshot_creates_the_task_directory_for_its_diff(self):
        out = self.root.parent / f"{self.root.name}-task" / "review.diff"
        self.addCleanup(lambda: shutil.rmtree(out.parent, ignore_errors=True))
        harness.snapshot(self.root, self.base, out)
        self.assertTrue(out.is_file())

    def test_every_fixed_secret_name_the_snapshot_refuses_is_read_denied(self):
        for rule in harness.SECRET_READ_DENIES:
            name = rule[len("Read(//**/"):-1].replace("**", "x").replace("*", "x")
            self.assertTrue(harness.SENSITIVE_FILE.search(name), rule)

    def test_snapshot_captures_hidden_worktree_changes(self):
        self.git("update-index", "--assume-unchanged", "README.md")
        self.write("README.md", "changed but hidden\n")
        tree = harness.snapshot(self.root)["tree"]
        self.assertEqual(self.git("show", f"{tree}:README.md"), "changed but hidden\n")

    def test_snapshot_in_linked_worktree_leaves_both_indexes_untouched(self):
        linked = Path(self.temp.name + "-linked")
        self.addCleanup(lambda: shutil.rmtree(linked, ignore_errors=True))
        self.git("worktree", "add", "-q", "--detach", str(linked))
        main_index = self.root / ".git/index"
        linked_index = Path(self.git("-C", str(linked), "rev-parse", "--path-format=absolute", "--git-path", "index").strip())
        before = main_index.read_bytes(), linked_index.read_bytes()
        (linked / "only-linked.txt").write_text("x\n")
        self.assertIn("only-linked.txt", self.git("ls-tree", "--name-only", harness.snapshot(linked)["tree"]))
        self.assertEqual((main_index.read_bytes(), linked_index.read_bytes()), before)

    # -- structure -------------------------------------------------------------------------------

    def copy_harness(self):
        shutil.copytree(SOURCE / ".claude", self.root / ".claude", ignore=shutil.ignore_patterns(
            "plans", "worktrees", "settings.local.json", ".review-cache", "__pycache__", "*.lock"))
        for name in ["CLAUDE.md", ".gitignore", "CHANGELOG.md"]:
            shutil.copy2(SOURCE / name, self.root / name)
        # Link and path targets the harness documents name, as placeholders: never the real sources or data.
        for doc in harness.harness_documents(SOURCE):
            for target in re.findall(r"\]\(([^)#]+)", doc.read_text()):
                dest = (self.root / doc.relative_to(SOURCE)).parent / target
                if not dest.exists() and not target.startswith(("http", "/")):
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    dest.write_text("fixture reference target\n")
        for _, target in harness.path_references(SOURCE):
            source, dest = SOURCE / target, self.root / target
            if source.is_dir():
                dest.mkdir(parents=True, exist_ok=True)
            elif source.exists() and not dest.exists():
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_text("fixture reference target\n")

        # Routed paths are compared lowercase; recreate each with its real spelling as a placeholder.
        for routed in harness.routed_paths_present(SOURCE) & {*harness.ANNOTATION_PATHS, *harness.ANNOTATION_EXEMPT, *harness.STORE_PATHS,
                                                             *harness.CONTRACT_PATHS, harness.DIVERGENCES, ".github/workflows/release.yml"}:
            spelled = next(p for p in subprocess.check_output(["git", "-C", str(SOURCE), "ls-files", "--cached", "--others", "--exclude-standard"], text=True).splitlines()
                           if p.lower() == routed or (routed.endswith("/") and p.lower().startswith(routed)))
            dest = self.root / (spelled if not routed.endswith("/") else spelled[:len(routed)] + "placeholder.rs")
            if not dest.exists():
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_text("fixture routed path\n")

    def test_live_harness_passes_its_own_check(self):
        result = harness.check(SOURCE)
        self.assertTrue(result["ok"], result["errors"])

    def test_structure_rejects_policy_drift_broken_links_and_ignored_harness(self):
        self.copy_harness()
        baseline = harness.check(self.root)
        self.assertTrue(baseline["ok"], baseline["errors"])
        bio, rev = ".claude/agents/vareffect-bioinformatics.md", ".claude/agents/vareffect-reviewer.md"
        hook, settings, review = ".claude/hooks/pin-reviewer-models.py", ".claude/settings.json", ".claude/skills/vareffect-review/SKILL.md"
        cases = [
            (bio, lambda s: s.replace("model: opus", "model: sonnet"), "agent model policy"),
            (".claude/agents/vareffect-researcher.md", lambda s: s.replace("effort: medium", "effort: high"), "agent model policy"),
            (rev, lambda s: s.replace("tools: Read,", "tools: Bash, Read,"), "unsafe agent tools"),
            (rev, lambda s: s.replace("disallowedTools: Bash, ", "disallowedTools: "), "tool exclusions"),
            (rev, lambda s: s.replace("description:", "memory: project\ndescription:"), "unexpected agent settings"),
            (rev, lambda s: s.replace("  - idiomatic-rust\n", ""), "preloaded skills mismatch"),
            (rev, lambda s: s.replace(".claude/INVARIANTS.md", "the charter"), "point to the charter"),
            (bio, lambda s: s.replace("  - concordance-measurement\n", "  - concordance-measurement\n  - no-such-skill\n"), "preloaded skills mismatch"),
            (".claude/agents/vareffect-extra.md", lambda s: s, "agent coverage"),
            (settings, lambda s: s.replace('"model": "opus"', '"model": "sonnet"'), "project model policy"),
            (settings, lambda s: s.replace('"Edit(/.claude/INVARIANTS.md)"', '"Edit(/.claude/other.md)"'), "charter edit guard"),
            (settings, lambda s: s.replace('"Edit(/.claude/INVARIANTS.md)",', '"Edit(/.claude/INVARIANTS.md)", "Write(/.claude/INVARIANTS.md)",'), "never consulted"),
            (settings, lambda s: s.replace('"Read(//**/.envrc)",', ""), "read-denied"),
            (settings, lambda s: s.replace(" || exit 2", ""), "fail-closed command"),
            (settings, lambda s: s.replace('"matcher": "Agent|Task"', '"matcher": "Bash"'), "model-pin hook"),
            (hook, lambda s: s.replace('NO_HAIKU = {"vareffect-reviewer", "vareffect-researcher"}', 'NO_HAIKU = {"vareffect-reviewer"}'), "Haiku ban disagrees"),
            (hook, lambda s: s.replace('OPUS_ONLY = {"vareffect-bioinformatics"}', "OPUS_ONLY = set()"), "disagree with AGENT_POLICY"),
            (hook, lambda s: s + "\ndef broken(:\n", "settings/hook"),
            (review, lambda s: s + "\n[missing](missing.md)\n", "broken link"),
            (review, lambda s: s + "\nSee `vareffect/src/no_such_module.rs`.\n", "stale path"),
            (review, lambda s: s.replace("description: ", "description: " + "x" * 300), "description over"),
            ("CLAUDE.md", lambda s: s + "x" * 16384, "below 16384"),
            (".gitignore", lambda s: s + "\n.claude\n", ".claude/ is ignored"),
            (".gitignore", lambda s: s.replace(".vareffect-tasks/\n", ""), "task state ignore rule"),
            ("vareffect-cli/src/csq.rs", None, "routed path no longer exists"),
        ]
        for name, change, expected in cases:
            with self.subTest(name=name, expected=expected):
                p = self.root / name
                if change is None:
                    moved = p.with_name("moved.rs")
                    p.rename(moved)
                    result = harness.check(self.root)
                    moved.rename(p)
                    self.assertTrue(any(expected in e for e in result["errors"]), result["errors"])
                    continue
                existed = p.exists()
                original = p.read_text() if existed else (self.root / rev).read_text().replace("name: vareffect-reviewer", "name: vareffect-extra")
                p.write_text(change(original))
                result = harness.check(self.root)
                self.assertFalse(result["ok"])
                self.assertTrue(any(expected in e for e in result["errors"]), result["errors"])
                if existed:
                    p.write_text(original)
                else:
                    p.unlink()
        review_path = self.root / review
        review_path.write_text(review_path.read_text() + "\nPlans live in `.claude/plans/task.md`.\n")
        self.assertTrue(harness.check(self.root)["ok"], "an ignored machine-local path is excused")
        (self.root / ".claude/skills/.DS_Store").write_text("finder\n")
        self.assertTrue(harness.check(self.root)["ok"], "dot entries in skills/ must be ignored")


class HookTest(unittest.TestCase):
    def setUp(self):
        self.hook = harness.load_hook(SOURCE)

    def spawn(self, role, model=None, tool="Agent", env=None):
        tool_input = {"subagent_type": role, "prompt": "review"}
        if model is not None:
            tool_input["model"] = model
        return self.hook.verdict({"tool_name": tool, "tool_input": tool_input}, env or {})

    def test_specialist_is_pinned_to_opus(self):
        for role in ["vareffect-bioinformatics", " VAREFFECT-Bioinformatics "]:
            with self.subTest(role=role):
                for model in [None, "opus", "OPUS", "fable", "opus[1m]", "claude-opus-5-5", "claude-fable-5-1"]:
                    self.assertIsNone(self.spawn(role, model), model)
                for model in ["sonnet", "haiku", "inherit", "opusplan", "claude-sonnet-5-5"]:
                    self.assertIsNotNone(self.spawn(role, model), model)
                    self.assertIsNotNone(self.spawn(role, model, tool="Task"), model)

    def test_provider_model_ids_are_opus_class(self):
        for model in ["us.anthropic.claude-opus-4-1-20250805-v1:0", "claude-opus-4-1@20250805", "claude-opus-5-5[1m]"]:
            self.assertTrue(self.hook.opus_class(model), model)
        for model in ["us.anthropic.claude-sonnet-4-5-v1:0", "opusplan", "notclaude-opus-4"]:
            self.assertFalse(self.hook.opus_class(model), model)

    def test_environment_cannot_downgrade_the_specialist(self):
        role = "vareffect-bioinformatics"
        self.assertIsNotNone(self.spawn(role, env={"CLAUDE_CODE_SUBAGENT_MODEL_FORCE": "1", "CLAUDE_CODE_SUBAGENT_MODEL": "sonnet"}))
        self.assertIsNone(self.spawn(role, env={"CLAUDE_CODE_SUBAGENT_MODEL_FORCE": "1", "CLAUDE_CODE_SUBAGENT_MODEL": "opus"}))
        self.assertIsNotNone(self.spawn(role, env={"ANTHROPIC_DEFAULT_OPUS_MODEL": "claude-sonnet-5-5"}))
        self.assertIsNone(self.spawn("vareffect-reviewer", env={"CLAUDE_CODE_SUBAGENT_MODEL_FORCE": "1", "CLAUDE_CODE_SUBAGENT_MODEL": "sonnet"}))

    def test_engineering_roles_may_escalate_but_not_use_haiku(self):
        for role in ["vareffect-reviewer", "vareffect-researcher"]:
            with self.subTest(role=role):
                self.assertIsNone(self.spawn(role, "opus"))
                self.assertIsNone(self.spawn(role, "sonnet"))
                self.assertIsNotNone(self.spawn(role, "haiku"))
                self.assertIsNotNone(self.spawn(role, "claude-haiku-4-5-20251001"))

    def test_other_tools_and_agents_pass_through_and_unreadable_fails_closed(self):
        self.assertIsNone(self.spawn("Explore", "haiku"))
        self.assertIsNone(self.hook.verdict({"tool_name": "Bash", "tool_input": {"command": "ls"}}, {}))
        for event in ([], "x", {"tool_name": "Agent", "tool_input": "x"},
                      {"tool_name": "Agent", "tool_input": {"subagent_type": "vareffect-bioinformatics", "model": 5}}):
            self.assertIsNotNone(self.hook.verdict(event, {}), event)

    def test_settings_command_blocks_with_exit_code_two(self):
        def run(event, project):
            env = dict(CLEAN_ENV, CLAUDE_PROJECT_DIR=str(project))
            return subprocess.run(["/bin/sh", "-c", harness.HOOK_COMMAND], input=event, cwd=SOURCE, env=env, capture_output=True, text=True)

        blocked = json.dumps({"tool_name": "Agent", "tool_input": {"subagent_type": "vareffect-bioinformatics", "model": "sonnet"}})
        allowed = json.dumps({"tool_name": "Agent", "tool_input": {"subagent_type": "vareffect-bioinformatics"}})
        result = run(blocked, SOURCE)
        self.assertEqual(result.returncode, 2)
        self.assertIn("pinned to Opus", result.stderr)
        self.assertEqual(run(allowed, SOURCE).returncode, 0)
        self.assertEqual(run("not json", SOURCE).returncode, 2)
        with tempfile.TemporaryDirectory() as outside:
            self.assertEqual(run(allowed, Path(outside) / "missing").returncode, 2, "no hook found must block")


if __name__ == "__main__":
    unittest.main()
