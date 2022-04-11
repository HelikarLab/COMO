# This file is responsible for setting up jupyter and python libraries

# --- Set up the virtual environment ---
mkdir /home/"${NB_USER}"/work
mkdir /home/"${NB_USER}"/.jupyter
mkdir /home/"${NB_USER}"/.local
echo "cacert=/etc/ssl/certs/ca-certificates.crt" > /home/"${NB_USER}"/.curlrc
echo "c.NotebookApp.ip = '0.0.0.0'" >> /home/"${NB_USER}"/.jupyter/jupyter_notebook_config.py
echo "c.NotebookApp.notebook_dir = '/home/${NB_USER}/work'" >> "/home/${NB_USER}/.jupyter/jupyter_notebook_config.py"
# --- End virtual environment setup ---

# --- Install Jupyter/Python libraries into virtual environment ---
# Update pip, setuptools, etc.
python3 -m pip install --upgrade --no-cache-dir pip setuptools wheel

# Install jupyter-related items
python3 -m pip --no-cache-dir install jupyter notebook jupyterlab bokeh ipywidgets jupyterhub
jupyter nbextension enable --py --sys-prefix widgetsnbextension
jupyter labextension install @jupyter-widgets/jupyterlab-manager escher

# Install python libraries
# Use maximum cores available, from: https://stackoverflow.com/a/32598533
python3 -m pip --no-cache-dir --install-option="--jobs=$(nproc)" install \
  "git+https://github.com/babessell1/cobamp.git@master#egg=cobamp" \
	"git+https://github.com/cokelaer/bioservices.git@master#egg=bioservices" \
  argparse \
  bioservices \
  cobra \
  escher \
  framed \
  GEOparse \
  lxml \
  memote \
  numpy \
  openpyxl \
  pandas \
  rpy2 \
  scipy \
  SQLAlchemy \
  troppo \
  unidecode \
  xlrd

# Delete cache, it is not required
rm -rf /root/.cache/pip
# --- Complete Jupyter installation ---
