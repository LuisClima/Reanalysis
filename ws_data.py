# -*- coding: utf-8 -*-
"""
Created on Sat Jan 19 15:56:53 2019

author: Jose Luis Rodriguez-Solis/Alejandro Dominguez-Perez
        CICESE
        Departamento de Oceanogradfia

prog que descarga todos los niveles de los datos de reanalisis del ERA-Interim

"""

import datetime
from ecmwfapi import ECMWFDataServer
from dateutil.relativedelta import *

server = ECMWFDataServer()

# Fecha incial
f      = '20150201'
# Fecha final
ffinal = '20170101'
# delta de tiempo
dh = 6
# resolucion
res=0.75

# comienza parte principal del texto

f  = datetime.datetime.strptime(f,"%Y%m%d")

for i in range (0,60000):
    server.retrieve({
        'stream'    : "mnth",
        "levtype"   : "pl",
        "levelist"  : "3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000",
        'param'     : "130.128/131.128/132.128/133.128",
        'dataset'   : "interim",
        'step'      : "0",
        'grid'      : "{0}/{0}".format(res),
        'date'      : f.strftime('%Y-%m-%d'),
        'type'      : "an",
        'class'     : "ei",
	  'format'    : "netcdf",
        'target'    : "/home/luis/Documents/CLASES_CICESE/ERA Interim/datos/mensual01/EI_men_{0}.nc".format(f.strftime("%Y%m%d"))
    })
    #f=f+datetime.timedelta(hours=dh)
    f = f + relativedelta(months=1)
    if f.strftime('%Y-%m-%d') == '2018-01-01':
        print ('-------------------')
        print ('Descarga finalizada')
        print ('-------------------')
        break

