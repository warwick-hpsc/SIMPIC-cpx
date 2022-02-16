# Source this to setup environment on different machines under bash
[[ "${BASH_SOURCE[0]}" != "${0}" ]] || echo "Please use source $0"

case $(hostname) in
    bison.lecad.fs.uni-lj.si)
	module load cuda-10.1.243-gcc-8.3.0-qpdenfa
	module load openmpi-3.1.4-gcc-8.3.0-hxcgcwk
	module load python-3.7.4-gcc-4.8.5-tvzmhun
	module load texlive-live-gcc-4.8.5-7ttre4z
	module load sublime-text-3.2.2.3211-gcc-8.3.0-uzm4u2p
	if test -d simpyenv
	    then source simpyenv/bin/activate
	else
	    python3 -m venv simpyenv
	    source simpyenv/bin/activate
	    python3 -m pip install --upgrade pip sphinx_rtd_theme
        fi
	;;

    davide*)
	module load cuda/9.2.88 gnu/6.4.0 openmpi/3.1.0--gnu--6.4.0
	module load szip/2.1.1--gnu--6.4.0
        module load hdf5/1.10.4--openmpi--3.1.0--gnu--6.4.0
	;;

    viz.hpc.fs.uni-lj.si)
        module load OpenMPI/3.1.4-gcccuda-2019b
        module load Python/3.7.4-GCCcore-8.3.0
        if test -d simpyenv
	    then source simpyenv/bin/activate
	else
	    python3 -m venv simpyenv
	    source simpyenv/bin/activate
	    python3 -m pip install --upgrade pip sphinx_rtd_theme            
        fi
        ;;
esac

export STARPU_DIR=${HOME}/starpu
export PKG_CONFIG_PATH=${STARPU_DIR}/lib/pkgconfig:${PKG_CONFIG_PATH}
export PATH=${STARPU_DIR}/bin:${PATH}
export LD_LIBRARY_PATH=${STARPU_DIR}/lib:${LD_LIBRARY_PATH}

test -d ${STARPU_DIR} || echo "Please install StarPU under ${STARPU_DIR}"

alias starpu_doc="firefox ${STARPU_DIR}/share/doc/starpu/manual/html/index.html"
alias simpic_doc="firefox ${PWD}/build/html/index.html"
