# Contributing to lin_algebra

Thank you for your interest in contributing!

This project is open to:
- bug fixes
- performance improvements
- documentation improvements
- examples and tests

---

## How to contribute

1. Fork the repository
2. Create a new branch from `main`
3. Make your changes
4. Open a Pull Request against `main`

Please keep pull requests focused and well-described.

---

## Development setup

```bash
cargo build
cargo test
cargo fmt
```


## Code Style
- Prefer clear names over clever code
- Add documentation comments (```///```) for public items
- Add tests for new functionality when possible

## Versioning

This project follows Semantic Versioning.

General guidelines:
- CI / workflow changes → no version bump
- Documentation or README fixes → no version bump
- Internal refactors with no API change → no version bump
- Bug fixes affecting users → PATCH bump (maintainer will handle)
- New backward-compatible features → MINOR bump
- Breaking changes → MAJOR bump

Releases are handled by the maintainer.


## Question or ideas?
Feel free to open an issue to discuss changes before starting work.