#!/usr/bin/python3
import os
import sys
import time

# find project root dir
class Configs():
    def __init__(self, projectdir):
        self.rootdir = projectdir
        self.datadir = os.path.join(projectdir,'data')
        self.outputdir = os.path.join(projectdir,'output')
        self.pydir = os.path.join(projectdir,'py')
        self.docdir = os.path.join(projectdir,'doc')


currentdir = os.getcwd()
dirlist = currentdir.split('/')
projectdir = '/'.join(dirlist[0:-1])

configs = Configs(projectdir)