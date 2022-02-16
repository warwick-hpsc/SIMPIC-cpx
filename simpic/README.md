# SIMPIC (Simple PIC)

See documentation at https://lecad-peg.bitbucket.io/simpic/

Building documentation with Python virtual environment (e.g. "simpyenv"):

    module load python
    python3 -m venv simpyenv
    source simpyenv/bin/activate
    python3 -m pip install --upgrade pip sphinx_rtd_theme

    make html
    firefox build/html/index.html

    deactivate

Alternatively one can install missing Sphinx locally with:

    module load python
    pip3 install --user --upgrade pip' command
    export PATH=${HOME}/.local/bin:${PATH}
    pip3 install --user sphinx_rtd_theme
    
    make html
    firefox build/html/index.html

