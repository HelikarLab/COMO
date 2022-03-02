python3 -m pip --no-cache-dir install pip --upgrade
python3 -m pip --no-cache-dir install setuptools==58.0.1 --upgrade
python3 -m pip --no-cache-dir install wheel --upgrade
python3 -m pip --no-cache-dir install \
	git+git://github.com/cokelaer/bioservices.git@master#egg=bioservices \
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
	git+git://github.com/babessell1/cobamp.git@master#egg=cobamp \
	unidecode \
	troppo
rm -rf /root/.cache