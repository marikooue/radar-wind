# This is a bit hacky and should most likely be handled by
# Python setup.py, but for now this will do

# Fortran mathematical operators module
rm -vf _operators.so
f2py -m _operators -h _operators.pyf _operators.f90 --overwrite-signature
f2py --fcompiler=gfortran --f90flags='-O3 -Wall' -c _operators.pyf _operators.f90

# Fortran mass continuity module
rm -vf _continuity.so
f2py -m _continuity -h _continuity.pyf _continuity.f90 --overwrite-signature
f2py --fcompiler=gfortran --f90flags='-O3 -Wall' -c _continuity.pyf _continuity.f90

# Fortran spatial smoothness module
rm -vf _smooth.so
f2py -m _smooth -h _smooth.pyf _smooth.f90 --overwrite-signature
f2py --fcompiler=gfortran --f90flags='-O3 -Wall' -c _smooth.pyf _smooth.f90

# Fortran background field module
rm -vf _background.so
f2py -m _background -h _background.pyf _background.f90 --overwrite-signature
f2py --fcompiler=gfortran --f90flags='-O3 -Wall' -c _background.pyf _background.f90 
