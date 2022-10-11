# Contributing

Contributions are welcome.

## Guidelines

All contributions shall be made via pull requests.  Pull requests should include updates to unit tests to coverage any changes in functionality or bug fixes.

`ssgdemux` uses `cargo fmt` to format its code, `cargo clippy` to check for possible issues with code and unit tests run via `cargo test` to verify behavior.

Prior to submitting any pull request, please run the precommit script that will format the code and check for any issues:

```console
./precommit.sh
```