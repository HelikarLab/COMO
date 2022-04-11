source /home/"${NB_USER}"/py3_env/bin/activate

python3 -m pip --no-cache-dir install pip setuptools wheel --upgrade
python3 -m pip --no-cache-dir install \
	git+https://github.com/cokelaer/bioservices.git@master#egg=bioservices \
	argparse \
	cobra \
	GEOparse \
	lxml \
	memote \
	numpy \
	pandas \
	scipy \
	SQLAlchemy \
	framed \
	xlrd \
	openpyxl \
	git+https://github.com/babessell1/cobamp.git@master#egg=cobamp \
	unidecode \
	troppo \
	escher
rm -rf /root/.cache
