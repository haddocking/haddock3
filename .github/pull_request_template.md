# Pull Request template

You are about to submit a new Pull Request. Before continuing make sure
you read the [contributing guidelines](CONTRIBUTING.md) and you comply
with the following criteria:

- [ ] Your code is well documented. Functions and classes have proper docstrings
  as explained in contributing guidelines. You wrote explanatory comments for
  those tricky parts.
- [ ] You wrote tests for the new code.
- [ ] `tox` tests pass. Run `tox` command inside the repository folder.
- [ ] `-test.cfg` examples execute without errors. Inside `examples/` run
  `python run_tests.py -b`
- [ ] You have stick to Python. Talk with us before if you really need to add
  other programming languages to HADDOCK3
- [ ] code follows our coding style
- [ ] You avoided structuring the code into classes as much as possible, using
  small functions instead. But, you can use classes if there's a purpose.
- [ ] PR does not add any *install dependencies* unless permission was granted
  by the HADDOCK team
- [ ] PR does not break licensing
- [ ] You get extra bonus if you write tests for already existing code :godmode:

