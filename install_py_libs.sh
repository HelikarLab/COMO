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
	git+git://github.com/babessell1/cobamp.git@master#egg=cobamp \
	unidecode \
	troppo
# cheetos
rm -rf /root/.cache