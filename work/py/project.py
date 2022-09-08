#!/usr/bin/python3

import os
import sys
import time


# find project root dir
class Configs:
    def __init__(self, projectdir):
        self.rootdir = projectdir
        self.datadir = os.path.join(projectdir, "data")
        self.configdir = os.path.join(projectdir, "data", "config_sheets")
        self.outputdir = os.path.join(projectdir, "output")
        self.pydir = os.path.join(projectdir, "py")


currentdir = os.getcwd()
dirlist = currentdir.split("/")

# Find the "pipelines" directory
split_index = 1
for directory in dirlist:
    if directory != "work":
        split_index += 1
    else:
        break

projectdir = os.path.join(*dirlist[0:split_index])

# Add leading "/", as it will not exist right now
projectdir = os.path.join("/", projectdir)
configs = Configs(projectdir)
