# IdealGasThermo.jl — agent guide

Architecture decisions and domain vocabulary live in `CONTEXT.md` and
`docs/adr/`. Read them before changing the pure core, and use the terms exactly.

## Commits — keep them atomic

- **One logical change per commit.** A commit is a single, cohesive,
  self-described change — a feature, a rename, a fix, a doc update — never a
  grab-bag of unrelated edits. If a task produces several distinct changes,
  make several commits, in dependency order (e.g. a type rename in one commit,
  an unrelated cleanup in another).
- **Don't mix concerns.** Behaviour changes, renames, and formatting/cleanup
  go in separate commits, so each diff is easy to review and to revert on its
  own.
- **Every commit leaves the repo working.** Run the full suite and only commit
  when it is green:

  ```
  julia --project=. -e 'using Pkg; Pkg.test()'
  ```

  For version-sensitive changes, also check the CI Julia versions
  (`lts` = 1.10 and `release` = 1.12); docs-affecting changes should build
  clean (`julia +lts --project=docs docs/make.jl`, then
  `git checkout -- docs/Project.toml`).
- **Write the message to explain the *why*.** A short imperative title plus a
  body that says what changed and the reasoning a reviewer would want.
- Commit only when asked; never push (the user pushes). End every commit
  message with:
  `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
