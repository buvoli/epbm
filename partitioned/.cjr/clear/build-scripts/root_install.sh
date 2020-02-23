#!/bin/bash

# -- ROOT INSTALL SCRIPT -------------------------------------------------------
# Installs dependencies for Fedora using the dnf package manager and additional
# manual installations. This script installs all the dependancies as root user.
# It responds to the following environmental variables:
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
#     LIB_LINALG      TRUE => Linear algebra libraries BLAS, LAPACK and FFTW
#     LIB_OPENMPI     TRUE => openmpi (loaded using module load mpi)
#
# ---- Dev Environemnts --------------------------------------------------------
#     DEV_JUPYTER     TRUE => Jupyter Notebook And Jupyter Lab with support for
#                             all select languages.
#
# NOTE: To add extra dependancies for any language, library, or development
# environment that can be installed with dnf simply add an entry to the arrays
# in 1.1-1.3. More sophisticated dependancies can be placed in the script
# root_install_extra.sh or written within this bash script.
# ------------------------------------------------------------------------------

pkg_manager="swupd"

# == STEP 1: Install DNF packages ==============================================

# -- 1.1 DNF Packages: languages -----------------------------------------------
pkg_lang_c=('c-basic' 'gdb');
pkg_lang_fortran=('c-basic' 'gdb');
pkg_lang_python3=('python3-basic' 'python-data-science')
pkg_lang_julia=('wget')
pkg_lang_R=('R-basic' 'R-extras')
pkg_lang_latex=('texlive')

# -- 1.2 DNF Packages: libraries  ----------------------------------------------
pkg_lib_linAlg=('openblas' 'devpkg-fftw')
pkg_lib_openMPI=('openmpi')

# -- 1.3 DNF Packages: development environments   ------------------------------
pkg_dev_jupyter=('nodejs-basic' 'jupyter')
pkg_dev_cli=('vim' 'git' 'tmux' 'emacs')

# -- Add packages to dnfPkg array ----------------------------------------------
declare -a pkgs=('sudo' 'sysadmin-basic'); # basic packages required for usage

# ----> languages

if [ "$LANG_C" = "TRUE" ] ; then
  pkgs=("${pkgs[@]}" "${pkg_lang_c[@]}") ; fi

if [ "$LANG_FORTRAN" = "TRUE" ] ; then
  pkgs=("${pkgs[@]}" "${pkg_lang_fortran[@]}") ; fi

if [ "$LANG_PYTHON3" = "TRUE" ] ; then
  pkgs=("${pkgs[@]}" "${pkg_lang_python3[@]}") ; fi

if [ "$LANG_JULIA" = "TRUE" ] ; then
pkgs=("${pkgs[@]}" "${pkg_lang_julia[@]}") ; fi

if [ "$LANG_R" = "TRUE" ] ; then
  pkgs=("${pkgs[@]}" "${pkg_lang_R[@]}") ; fi

if [ "$LANG_LATEX" = "TRUE" ] ; then
  pkgs=("${pkgs[@]}" "${pkg_lang_latex[@]}") ; fi

# ----> libraries

if [ "$LIB_LINALG" = "TRUE" ] ; then
  pkgs=("${pkgs[@]}" "${pkg_lib_linAlg[@]}") ; fi

if [ "$LIB_OPENMPI" = "TRUE" ] ; then
  pkgs=("${pkgs[@]}" "${pkg_lib_openMPI[@]}") ; fi

# ----> development environments

if [ "$DEV_JUPYTER" = "TRUE" ] ; then
  pkgs=("${pkgs[@]}" "${pkg_dev_jupyter[@]}") ; fi

if [ "$DEV_CLI" = "TRUE" ] ; then
  pkgs=("${pkgs[@]}" "${pkg_dev_cli[@]}") ; fi

# -- remove redundant elements then install (requires bash 4+) -----------------
declare -A pkgsUniq
for k in ${pkgs[@]} ; do pkgsUniq[$k]=1 ; done

# -- install dependencies ------------------------------------------------------
echo "$pkg_manager bundle-add ${!pkgsUniq[@]}"
eval $pkg_manager bundle-add ${!pkgsUniq[@]}

# == STEP 2: Install Additional Packages =======================================

# -- Julia ---------------------------------------------------------------------
if [ "$LANG_JULIA" = "TRUE" ] ; then
  mkdir -p /opt
  mkdir -p /usr/local/bin/
  cd /opt
  wget https://julialang-s3.julialang.org/bin/linux/x64/1.3/julia-1.3.1-linux-x86_64.tar.gz
  tar -xzf julia-1.3.1-linux-x86_64.tar.gz
  ln -s /opt/julia-1.3.1/bin/julia /usr/local/bin/julia
fi

# -----> Jupyter
if [ "$DEV_JUPYTER" = "TRUE" ] ; then
  pip3 install jupyterlab # Jupyter Lab
  # --> matplotlib Widgets for JupiterLab (https://github.com/matplotlib/jupyter-matplotlib)
  if [ "$LIB_MATPLOTLIB" = "TRUE" ] ; then
    pip3 install ipympl
    jupyter labextension install @jupyter-widgets/jupyterlab-manager
    jupyter labextension install jupyter-matplotlib
  fi
  if [ "$LANG_LATEX" = "TRUE" ] ; then
    # --> Latex for JupyterLab (https://github.com/jupyterlab/jupyterlab-latex)
    pip3 install jupyterlab_latex
    jupyter labextension install @jupyterlab/latex
  fi
  if [ "$LANG_C" = "TRUE" ] ; then
    # --> C Kernel for Jypyter https://github.com/brendan-rius/jupyter-c-kernel
    pip3 install jupyter-c-kernel
    install_c_kernel --user
  fi
  if [ "$LANG_R" = "TRUE" ] ; then
    # --> R Kernel for Jupyter (https://irkernel.github.io/installation/)
    R -e 'r = getOption("repos"); r["CRAN"] = "https://cloud.r-project.org/"; install.packages(c("repr", "IRdisplay", "IRkernel"), repos = r, type = "source");'
  fi
  if [ "$LANG_FORTRAN" = "TRUE" ] ; then
    echo "Fortran Coarrays not yet supported"
    # caf (fortran coarrays) not installing properly. Resources: https://github.com/openhpc/ohpc/issues/884
    # https://github.com/sourceryinstitute/OpenCoarrays/blob/master/INSTALL.md
    # https://md2pdf.herokuapp.com/sourceryinstitute/opencoarrays/blob/master/INSTALL.md
    # another alternative is: (https://lfortran.org/ https://docs.lfortran.org/installation/) and (https://github.com/mgaitan/fortran_magic - calls fortran from python - not a kernel)
    #
    # cd /opt
    # git clone https://github.com/sourceryinstitute/jupyter-CAF-kernel
    # cd jupyter-CAF-kernel
    # pip install -e prebuild/jupyter-caf-kernel
    # jupyter-kernelspec install prebuild/jupyter-caf-kernel/Coarray-Fortran/
  fi
fi
