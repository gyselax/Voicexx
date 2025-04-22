# Generating the Voice++ documentation

Voice++ uses MkDocs and [`mkdoxy`](https://mkdoxy.kubaandrysek.cz/) (an additional package designed to integrate Doxygen into MkDocs) to generate the documentation. Documentation for these packages can be found online.

## Packages for installation

The packages required to make the docs are listed in Gyselalib++'s [requirements.txt](https://github.com/gyselax/gyselalibxx/blob/devel/docs/requirements.txt) file. Simply run:

```bash
pip install -r gyselalibxx/docs/requirements.txt
```

Adjustments to MkDocs can be done in the `mkdocs.yml`.

The documentation can be built locally by running the following commands from the root directory:

```bash
python3 gyselalibxx/docs/prepare_mkdocs.py . docs/mkdoc gyselalibxx docs/mkdoc build/
cd docs
mkdocs build
```

## Configuring MkDocs

The configuration of MkDocs, including various options and additional packages, is managed through the `mkdocs.yml` file. This file allows you to customise the structure and appearance of your documentation.

