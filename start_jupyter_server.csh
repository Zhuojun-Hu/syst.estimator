setenv JUPYTER_RUNTIME_DIR /disk03/usr8/huzjku/venv/jupyter_root/runtime
setenv JUPYTER_WORKSPACE /disk03/usr8/huzjku/fq.syst.ana 
source /disk03/usr8/huzjku/venv/jupyter_root/bin/activate.csh

jupyter lab --no-browser --port=9527 --notebook-dir $JUPYTER_WORKSPACE
