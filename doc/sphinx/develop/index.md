# Develop

(sec-compiling)=
## Compiling Cantera from Source

If you're interested in contributing new features to Cantera, or you want to try the
latest and development version, you will need to compile Cantera from source. The first
step is to make sure you have all the [](compiling/compilation-reqs) installed. Then,
you can [download the Cantera source code](compiling/source-code). Finally, you can
determine the appropriate configuration options and [compile
Cantera](compiling/configure-build) on your computer.

The following additional references may also be useful:

- [](compiling/dependencies.md)
- [](compiling/config-options)
- [](compiling/special-cases)

```{toctree}
:caption: Compiling Cantera from Source
:hidden:
:maxdepth: 1

compiling/compilation-reqs
compiling/source-code
compiling/configure-build
compiling/dependencies
compiling/config-options
compiling/special-cases
```

## How Cantera Works

```{caution}
This section is a work in progress.
```

- [](reactor-integration)

```{toctree}
:caption: How Cantera Works
:hidden:
:maxdepth: 1

reactor-integration
```

## Adding New Features to Cantera

```{caution}
This section is a work in progress.
```

It is important to include robust tests for new features to reduce the potential
inclusion of bugs. When testing in the python framework, `cantera/test/python`,
individual tests can ran individually by using the `diagnose` marker and flag.
In the python code, it looks like the following.

```python
import pytest
import cantera as ct
from . import utilities

class TestMyNewFeature(utilities.CanteraTest):

    @pytest.mark.diagnose
    def my_new_feature_assert_true(self):
        assert True
```

The python framework is then run with `scons test-python --diagnose` which will
skip any tests not decorated with `@pytest.mark.diagnose`.

- [](doc-formatting)

```{toctree}
:caption: Adding New Features to Cantera
:hidden:

doc-formatting
```
