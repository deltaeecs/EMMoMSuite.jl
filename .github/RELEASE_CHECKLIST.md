# EMMoMSuite General Registry Release Checklist

This checklist tracks everything needed before publishing EMMoMSuite to Julia General Registry.

## 1. Package Metadata

- [x] `Project.toml` has stable `name`, `uuid`, `version`.
- [x] `Project.toml` has dependency and compat constraints.
- [x] `authors` field updated from placeholder values (`xyhe <1225385871@qq.com>`).
- [x] Repository URL and homepage confirmed in repo metadata（发布时写入 `deltaeecs/EMMoMSuite.jl`）。

## 2. Quality Gates

- [x] CI test workflow enabled on push and PR.
- [x] Coverage measurement available (`scripts/check_coverage.jl`, CI `lcov.info`).
- [x] Local source coverage >= 80% (current: 93.21%).
- [x] Documenter build runs clean in docs environment.

## 3. Release Automation

- [x] TagBot workflow added (`.github/workflows/TagBot.yml`).
- [ ] TagBot secret configured (`TAGBOT_TOKEN` or SSH key as required).
- [ ] Registrator app installed for repository owner/org.
- [ ] Maintainer account can trigger Registrator comments.

## 4. Versioning and Changelog

- [ ] Move entries from `Unreleased` to a concrete release section.
- [ ] Bump `version` in `Project.toml` for the next release tag.
- [ ] Create signed/annotated git tag matching `Project.toml` version.

## 5. Registrator Submission

- [ ] Merge release commit to default branch.
- [ ] Open release PR to General Registry via Registrator comment.
- [ ] Resolve AutoMerge or review comments from General Registry.
- [ ] Confirm package install via `Pkg.add("EMMoMSuite")` after merge.

## 6. Post-Release

- [ ] Create GitHub Release notes.
- [ ] Verify docs and badges against released tag.
- [ ] Announce release and capture follow-up issues.
