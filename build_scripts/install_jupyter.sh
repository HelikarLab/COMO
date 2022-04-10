python3 -m pip --no-cache-dir install jupyter notebook jupyterlab bokeh ipywidgets
jupyter nbextension enable --py --sys-prefix widgetsnbextension
jupyter labextension install @jupyter-widgets/jupyterlab-manager escher
python3 -m pip --no-cache-dir install jupyterhub
