cd src/ # go to the folder
f2py3 -c -m kpmf90 kpm.f90    # compile the library
#export PYTHONPATH=$PYTHONPATH:`pwd`  # add this path 
