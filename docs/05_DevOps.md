# DevOps

<!---toc start-->
   * [Releasing a New Version](#releasing-a-new-version)
      * [Major Release](#major-release)
      * [Minor and Patch Release](#minor-and-patch-release)
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

### Semantic Versioning

This tool follows [Semantic Versioning](https://semver.org/).  In brief:

* MAJOR version when you make incompatible API changes,
* MINOR version when you add functionality in a backwards compatible manner, and
* PATCH version when you make backwards compatible bug fixes.


### Major Release

To create a major release:

```console
cargo release major --no-publish --execute
```

This will remove any pre-release extension, create a new tag and push it to github, and kick off the release GitHub action that will publish to the release page.

Upon success, move the version to the [next candidate release](#release-candidate)

### Minor and Patch Release

To create a _minor_ (_patch_) release, follow the [Major Release](#major-release) instructions substituting `major` with `minor` (`patch`):

```console
cargo release minor --no-publish --execute
```

### Release Candidate

To move to the next release candidate:

```console
cargo release rc --no-tag --no-publish --execute
```

This will create or bump the pre-release version and push the changes to the main branch on github.
This will not tag and publish the release candidate; remove `--no-tag` to create a new tag and push it to github, and kick off the release GitHub action that will publish to the release page.

[cargo-release-link]:      https://github.com/crate-ci/cargo-release
[cargo-release-docs-link]: https://github.com/crate-ci/cargo-release/blob/master/docs/reference.md

