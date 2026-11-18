# Contributing to LiteDIC

Thank you for your interest in LiteDIC. This page explains how to report a
problem, ask for help, and propose changes.

## Reporting a bug

Open an [issue](https://github.com/alonePlayerr1/LiteDIC/issues/new/choose)
using the *Bug report* template. Please include:

- your MATLAB version (`version`) and operating system;
- the LiteDIC version or commit (`git rev-parse --short HEAD`);
- the smallest script and image pair that reproduce the problem;
- what you expected and what happened (full error message, if any).

If the issue concerns measurement accuracy, please also attach the output of
`RUN_TESTS` on your machine.

## Asking for help

Questions about using LiteDIC (choice of subset size, strain window,
interpreting results, adapting it to your images) are welcome as issues with
the *Question* template, or in
[Discussions](https://github.com/alonePlayerr1/LiteDIC/discussions).

## Support expectations

LiteDIC is maintained by one person alongside research work. I aim to reply to
new issues within one week. Bug reports that come with a reproducible example
are fixed first.

## Proposing changes

1. Open an issue first describing the change, so that we can agree on the
   approach before you invest time in it.
2. Fork the repository and create a branch from `main`.
3. Keep the code in base MATLAB: no toolbox may become mandatory.
4. Add or update a check in `RUN_TESTS.m` that covers your change.
5. Run `RUN_TESTS` locally; the report must show `0 FAIL, 0 ERRORS`.
6. Open a pull request referencing the issue. Continuous integration runs the
   same test suite automatically.

### Licence of contributions

By submitting a contribution you agree that it is released under the MIT
License, like the rest of LiteDIC.

## Code of conduct

Everyone interacting with this project is expected to follow the
[Code of Conduct](CODE_OF_CONDUCT.md).
