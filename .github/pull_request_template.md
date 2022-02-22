# Pull Request template

You are about to submit a new Pull Request. Before continuing make sure
you read the [contributing guidelines](CONTRIBUTING.md) and you comply
with the following criteria:

- [ ] your code is well documented
- [ ] you wrote tests for the new code
- [ ] `tox` tests pass. Run `tox` command inside the repository folder.
- [ ] `-test.cfg` examples execute without faults. Inside `examples/`
  run `python run_tests.py -b`
- [ ] PR contains only Python code
- [ ] PR follows our coding style
- [ ] you avoided structuring the code into classes as much as possible.
  Use small functions instead.
- [ ] PR does not contain any additional dependencies unless permission
  granted by the HADDOCK team
- [ ] PR does not break licensing
- [ ] you get extra bonus if you write tests for already existing code
  :godmode:

