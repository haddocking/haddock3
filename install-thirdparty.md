HADDOCK3 is able to integrate third-party software in its workflows.
However, we are not responsible to the proper installation of such
packages, but we do help you to install them. In this document you will
find a list of all third-party packages HADDOCK3 can use and guidelines
on how to install them.


# Installing third-party packages

## `lightdock`

To install to [lightdock](https://github.com/lightdock/lightdock) follow
the instructions in the project's website. Remember to install it under
the same Python environment you created for HADDOCK3. If you have any
doubts, please let us know.

* * *

## `gdock`

1. Clone the latest version:

```
cd some-folder
git clone https://github.com/rvhonorato/gdock.git
```

2. Install Python3+ dependencies
```
pip install deap scipy mgzip biopython
```

3. Set `GDOCK_PATH`
```
export GDOCK_PATH=some-folder
```

**Important**: These are not the full install instructions as here only
the model generation is used. Please check the [repository
page](https://github.com/rvhonorato/gdock) for more information.

* * *
