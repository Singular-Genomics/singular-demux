# DevOps

<!---toc start-->
   * [Releasing a New Version](#releasing-a-new-version)
      * [Major Release](#major-release)
      * [Release Candidate](#release-candidate)

<!---toc end-->

## Releasing a New Version

1. Install [`cargo-release`][cargo-release-link]

```console
cargo install cargo-release
```

2. Create a release that will not try to push to `crates.io` and verify the command:

```console
cargo release [major,minor,patch,release,rc...] --no-publish
```

Note: "dry-run" is the default for cargo release.

3. See the [`cargo-release` reference documentation][cargo-release-docs-link] for more information


### Major Release

To create a major release:

```console
cargo release major --skip-publish --execute
```

This will remove any pre-release extension, create a new tag and push it to github, and kick off the release GitHub action that will publish to the release page.

Upon success, move the version to the [next candidate release](#release-candidate)

### Release Candidate

To move to the next release candidate:

```console
cargo release rc --no-tag  --skip-publish --execute
```

This will create or bump the pre-release version, create a new tag and push it to github, and kick off the release GitHub action that will publish to the release page.

[cargo-release-link]:      https://github.com/crate-ci/cargo-release
[cargo-release-docs-link]: https://github.com/crate-ci/cargo-release/blob/master/docs/reference.md

