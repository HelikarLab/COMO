python3 -m pip --no-cache-dir install pip --upgrade
python3 -m pip --no-cache-dir install setuptools --upgrade
python3 -m pip --no-cache-dir install wheel --upgrade
python3 -m pip --no-cache-dir install \
	bioservices \
	cobra \
	GEOparse \
	lxml \
	numpy \
	pandas \
	scipy \
	SQLAlchemy \
	framed \
	xlrd \
	openpyxl \
	cobamp \
	unidecode \
	troppo
rm -rf /root/.cache