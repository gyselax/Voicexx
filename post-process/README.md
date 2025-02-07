# Post-processing

This contains all the python files needed to plot various quantities used to analyse the output of a simulation.

## How to use it

You will need to append this directory, and its subdirectories to the PATH and PYTHONPATH variables:
```
export PATH=$PATH:<path-to-voice>/post-process/
export PATH=$PATH:<path-to-voice>/post-process/neutrals/
export PATH=$PATH:<path-to-voice>/post-process/sheath/

export PYTHONPATH=$PYTHONPATH:<path-to-voice>/post-process/
```
You will moreover need some tools from gyselalibxx
```
export PYTHONPATH=$PYTHONPATH:<path-to-voice>/gyselalibxx/post-process/PythonScripts/
```

Then, to handle python libraries, the easiest is to use a python venv.
See [python-venv documentation]{https://scitas-doc.epfl.ch/user-guide/software/python/python-venv/}.

Finally, cd into the directory where the output of your simulation are, and execute the command.
Use the flag --help or -h to consult the help of the command.
```
plot_fistribu --help
```

