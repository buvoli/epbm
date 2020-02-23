#!/bin/bash

# -- USER INSTALL SCRIPT -------------------------------------------------------
# Installs dependencies for Fedora using the dnf package manager and additional
# manual installations. This script installs all the dependancies as a standard
# user. It responds to the following environmental variables:
#
# ---- languages ---------------------------------------------------------------
#     LANG_C          TRUE => C language packages installed
#     LANG_FORTRAN    TRUE => Fortran language installed
#     LANG_PYTHON3    TRUE => Python3 language installed
#     LANG_JULIA      TRUE => Julia language installed
#     LANG_R          TRUE => R languag installed
#     LANG_LATEX      TRUE => Latex installed
#
# ---- libraries ---------------------------------------------------------------
#
# ---- Dev Environemnts --------------------------------------------------------
#     DEV_JUPYTER     TRUE => Jupyter Notebook And Jupyter Lab with support for
#                             all select languages.
#
# NOTE: Additional dependancies can be placed in the script
# user_install_extra.sh or written directly within this bash script.
# ------------------------------------------------------------------------------

# Certain Julia Packages do not install as root. Install them here instead
# -- Julia Packages ------------------------------------------------------------
if [ "$LANG_JULIA" = "TRUE" ] ; then
  # ----> plotters
  if [ "$LIB_MATPLOTLIB" = "TRUE" ] ; then
      julia -e 'import Pkg; Pkg.add("PyPlot")'
  fi
  julia -e 'import Pkg; Pkg.add("Gadfly")'
  julia -e 'import Pkg; Pkg.add("Plots")'
fi

# -- Juputer Kernels -----------------------------------------------------------
if [ "$DEV_JUPYTER" = "TRUE" ] ; then
  # ----> Julia
  if [ "$LANG_JULIA" = "TRUE" ] ; then
    julia -e 'import Pkg; Pkg.add("IJulia")'
  fi
  # ----> R
  if [ "$LANG_R" = "TRUE" ] ; then
    R -e "IRkernel::installspec()"
  fi
fi

# -- Spack ---------------------------------------------------------------------
if [ "$PKGM_SPACK" = "TRUE" ] ; then
  mkdir -p ~/.local
  mkdir -p ~/.local/bin
  git clone https://github.com/spack/spack.git ~/.local/spack
  ln -s ~/.local/spack/bin/spack ~/.local/bin/spack
fi
