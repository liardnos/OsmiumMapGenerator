import configparser
import json
import os

from gdal_interfaces import GDALTileInterface

#file based on open-elevation server.py

print('Reading config file ...')
parser = configparser.ConfigParser()
parser.read('config.ini')

PORT = parser.getint('server', 'port')
DATA_FOLDER = parser.get('server', 'data-folder')
OPEN_INTERFACES_SIZE = parser.getint('server', 'open-interfaces-size')
ALWAYS_REBUILD_SUMMARY = parser.getboolean('server', 'always-rebuild-summary')

"""
Initialize a global interface. This can grow quite large, because it has a cache.
"""

print(DATA_FOLDER, '%s/summary.json' % DATA_FOLDER)
interface = GDALTileInterface(DATA_FOLDER, '%s/summary.json' % DATA_FOLDER, OPEN_INTERFACES_SIZE)


if interface.has_summary_json() and not ALWAYS_REBUILD_SUMMARY:
    print('Re-using existing summary JSON')
    interface.read_summary_json()
else:
    print('Creating summary JSON ...')
    interface.create_summary_json()

def lookUpElevation(lat, lng):
    return float(interface.lookup(lat, lng))
