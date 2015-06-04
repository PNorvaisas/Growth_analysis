#!/usr/bin/env python
# -*- coding: utf-8 -*-
#sys.setdefaultencoding('iso-8859-1')
"""
FFT.py

Created by Povilas Norvaisas on 2015-02-26.
Copyright (c) 2015. All rights reserved.

"""
try:
	import pip
except ImportError, e:
	print "Module pip not found!"
def install(package):
	pip.main(['install', package])

for mod in ['pip','string','math','re','csv','sys','os','commands','datetime','operator','getopt','subprocess','pickle','shutil','glob','types','math','copy','pyExcelerator','xlrd','xlwt','xlutils']:
	try:
		exec "import %(mod)s" % vars()
	except ImportError, e:
		print "Module not found %(mod)s\nTrying to install!" % vars()
		install(mod)		
		#pass # module doesn't exist, deal with it.

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pylab import *

import scipy.signal as sig
from scipy.fftpack import rfft, irfft, fftfreq
from scipy import interpolate as ip
from collections import OrderedDict
from collections import defaultdict
from collections import Counter
from scipy.optimize import curve_fit
import numpy as np
from operator import itemgetter
from itertools import groupby

from multiprocessing import Pool


compare = lambda x, y: Counter(x) == Counter(y)

try:
	import itertools as IT
except ImportError, e:
	print "Module itertools not found"





help_message = '''
FFT filtering of data

Flags:
	-h  Display this help message
	-v  Verbose output
Arguments:
	-i <files>   Input files, list separated by comma, file with list
	-d <file>    CSV file with data labels
	-f <filter>  wiener,butter
	-o <dir>     directory to write output to

Options:
	load	     Load pickled data

'''


class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg

def usage():
	print help_message % vars()
	return

optionsset='''

Options:
<--------------------------------------------->
      Files:	%(ifile)s
	 DB:	%(dfile)s
     Filter:	%(filterf)s
        Out:	%(odir)s
       Load:	%(load)s
<--------------------------------------------->
	'''



def main(argv=None):

	ifile=""
	dfile=""
	msize=20
	mode="Zaslaver"
	filterf=''
	odir='Output'
	load=False
	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "hi:f:d:o:", ["help"])
		except getopt.error, msg:
			raise Usage(msg)


		for option, value in opts:
			if option in ("-h", "--help"):
				usage()
				return	
			if option in ("-i", "--input"):
				ifile=value
			if option in ("-d", "--db"):
				dfile=value
			if option in ("-f", "--filter"):
				filterf=value
			if option in ("-o", "--out"):
				odir=value
	
	
		for argument in args:		
			if argument in ("load", "--load"):
				load = True
			


	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		print >> sys.stderr, "\t for help use --help"
		return 2

	#Check for the input integrity

	print '''\n\n---------------------------\n\n'''
	
	if not load:
		try:
			if ifile!="":
				ifiles=ifile.split(',')
				for ifl in ifiles:
					print '{} {}'.format(ifl,os.path.isfile(ifl))
					ipath, iname, itype = filename(ifile)			
			else:
				raise Exception("No files specified.")
			if dfile!="":
				dpath, dname, dtype = filename(dfile)			
			else:
				raise Exception("No database file specified.")
			if filterf!='' and filterf in ['wiener','butter']:
				if filterf=='butter':
					print "Butterfield filter selected!"
				elif filterf=='wiener':
					print "Wiener filter selected!"
			else:
				print "Unknown filter {}, not using.".format(filterf)
				filterf=''
				
			
	
		except Exception, e:
			print e
			#print "!!!-------------------------------------------!!!"
			sys.exit(1)

		print optionsset %vars()
	#----------------------------


	


	if load:
		f = open('UAL_data.pckl','rb')
		data,genes,odir,filterf = pickle.load(f)
		f.close()

		
	else:	
		modes={'Zaslaver':['600nm','535nm'],'Biolog':['600nm','700nm']}
		genes=csvreader(dfile)

		waves=modes[mode]
		ilist=genlist(ifile)

		checkfiles(ilist)	
		data=collect(ilist)
		data=analyze(data,waves,filterf)
		dirn=dircheck(odir)
		data,allfit=growthfit(data,False)
		growthplot(allfit)		
		data=datashift(data,'all','all','all','all')
		data=analyze2(data,filterf)	

		f = open('UAL_data.pckl', 'w')
		pickle.dump([data,genes,odir,filterf], f)
		f.close()
	
	
	sheets=makesheets(data,genes)
	writesheets(sheets,odir)
	#data,allfit=growthfit(data,False)
	#data,shifts=datashift(data,'all','all','all','all')
	#plotall(data,'all','Control','Growth','C10',False,True)
	#plotall(data,'all','Control','Growth','C10',False,False)
	#plotall(data,'1','Experiment','Fluorescence_norm','all',False,True)
	#plotall(data,'Fluorescence_norm')
	#plt,plots=plot_2Dplates(data,'Growth','Experiment',True)
	#plt.show()




	

	

	
	
	
	

#-------------Functions------------


def plot96(data,genes,dirn):
	for plate in sorted(data.keys()):		
		for fg in data[plate]['Control']['Figures']:
			print "Plotting plate {} {}...".format(plate,fg)
			if 'dt' in fg:
				plt, plots=plot_2D(plate+"-"+fg,data[plate]['Control'][fg],data[plate]['Experiment'][fg],data[plate]['Control']['Time_dt'],data[plate]['Control']['Labels'],data[plate]['Control']['Wells'],genes[plate])
			else:
				plt, plots=plot_2D(plate+"-"+fg,data[plate]['Control'][fg],data[plate]['Experiment'][fg],data[plate]['Control']['Time'],data[plate]['Control']['Labels'],data[plate]['Control']['Wells'],genes[plate])
			data[plate]['Joint'][fg]=plt
			plt.savefig('{}/{}.pdf'.format(dirn,plate+'_'+fg))
	return data

def plotall(data,plsel,tpsel,fg,lsel,norm,shifted):
	ts=5
	if plsel=='all':
		plates=data.keys()
	elif plsel in data.keys():
		plates=[plsel]
	for plate in plates:
		print 'Plate number {}'.format(plate)
		if tpsel=='all':
			types=data[plate].keys()
		elif tpsel in data[plate].keys():
			types=[tpsel]
		for tp in types:
			#x=data[plate][tp]['Time']
			
			if tp=='Control':
				plt.figure(0)
				plt.title(fg+' without Metformin')

			
			elif tp=='Experiment':
				plt.figure(1)
				plt.title(fg+' with Metformin')

			if fg in ['Growth','600nm']:
				if norm:
					plt.ylim([0,1])
				else:			
					plt.ylim([0,0.4])
				plt.ylabel('OD')
				#plt.yscale('log')
			else:
				#plt.ylim([0,3000])
				plt.ylabel('Fluorescence, a.u.')
				#plt.yscale('log')
			
			plt.xlim([0,20])
			plt.xlabel('Time, h')

			if '_dt' in fg:
				x=data[plate][tp]['Time_dt']
			else:
				x=data[plate][tp]['Time']
			if lsel=='all':
				labels=data[plate][tp]['Labels']
			elif isinstance(lsel, list):
				labels=lsel
			elif isinstance(lsel, str) or isinstance(lsel, unicode):
				labels=[lsel] 
			for l in labels:
				if (plate!='21' and l!='A12' ) or (plate=='21' and l in labels[:10]):
					if shifted:
						y=data[plate][tp]['Shift'][fg][l]
					else:
						y=data[plate][tp][fg][l]
					if norm and fg in ['Growth','600nm']:
						y=(y-min(y))/max(y)
					plt.plot(x/3600,y,'r-',label=plate+l )
	plt.show()




def growthfit(data,norm):
	allfit=NestedDict()
	if norm:
		base=0.1
		top=0.2
		y0=0.2
	else:
		base=0.02
		top=0.07
		y0=0.05
	
	for plate in data.keys():
		for tp in data[plate].keys():
			x=data[plate][tp]['Time']
			labels=data[plate][tp]['Labels']
	
			for l in labels:
				if (plate!='21' and l!='A12') or (plate=='21' and l in labels[:10]):
					y=data[plate][tp]['Growth'][l]
					if norm:
						y=(y-min(y))/max(y)
					x2,y2=cut(x, y, base, top)
					y2l=np.log(y2)
					for ydata, ynm in IT.izip([y2,y2l],['Linear','Log']):
						try:
							popt, pcov = curve_fit(growth, x2, ydata)
						except TypeError:
							'Curve_fit encountered an error!'
							continue
						a=popt[0]
						c=popt[1]
						if ynm=='Log':
							t0=(np.log(y0)-c)/(a*60)
						else:
							t0=(y0-c)/(a*60)
						#print '{} growth rate: {}, start: {}'.format(ynm,a*3600,t0)
						for par,nm in IT.izip([a,c,t0],['a','c','t0']):
							if allfit[tp][ynm][nm]:
								allfit[tp][ynm][nm]=allfit[tp][ynm][nm]+[par]
							else:
								allfit[tp][ynm][nm]=[par]
					
						data[plate][tp]['GrowthFit'][l][ynm]=[a,c,t0]

	return data,allfit



def interp(x,y,x2):
	tck = ip.splrep(x, y, s=0)
	y2=ip.splev(x2, tck, der=0)
	return y2

def datashift(data,plsel,tpsel,fgsel,lbsel):
	ts=5
	if plsel=='all':
		plates=data.keys()
	elif plsel in data.keys():
		plates=plsel
	for plate in plates:
		if tpsel=='all':
			types=data[plate].keys()
		elif tpsel in data[plate].keys():
			types=[tpsel]
		for tp in types:
			if fgsel=='all':
				figures=data[plate][tp]['Figures']
			elif fgsel in data[plate][tp]['Figures']:
				figures=[fgsel]
			else:
				print 'Unknown figure!'
			for fg in figures:
				if lbsel=='all':
					labels=data[plate][tp]['Labels']
				elif isinstance(lbsel, list):
					labels=lbsel
				elif isinstance(lbsel, str) or isinstance(lbsel, unicode):
					labels=[lbsel]
				for l in labels:
					if (plate!='21' and l!='A12') or (plate=='21' and l in data[plate][tp]['Labels'][:10]):
						#print plate,tp,fg,l
						a, c, t0=data[plate][tp]['GrowthFit'][l]['Log']
						y=data[plate][tp][fg][l]
						if '_dt' in fg:
							x=data[plate][tp]['Time_dt']
						else:
							x=data[plate][tp]['Time']
						xs=(x+(ts*3600-t0*60))
						y2=interp(xs,y,x)
						shift=ts*3600-t0*60
						if shift<0:
							steps=int(-shift/300)
							y2[-steps:]=y2[-steps-1]
						elif shift>0:
							steps=int(shift/300)
							y2[:steps]=y2[steps+1]
							
							
					else:
						y2=data[plate][tp][fg][l]
					data[plate][tp]['Shift'][fg][l]=y2
				

	return data	


def meansd(data,plsel,tpsel,fgsel,lbsel,shifted):
	alldata=[]
	ts=5
	if plsel=='all':
		plates=data.keys()
	elif plsel in data.keys():
		plates=[plsel]
	for plate in plates:
		if tpsel=='all':
			types=data[plate].keys()
		elif tpsel in data[plate].keys():
			types=[tpsel]
		for tp in types:
			if fgsel=='all':
				figures=data[plate][tp]['Figures']
			elif fgsel in data[plate][tp]['Figures']:
				figures=[fgsel]
			else:
				print 'Unknown figure!'
			for fg in figures:
				if shifted:
					raw=data[plate][tp]['Shift'][fg]
				else:
					raw=data[plate][tp][fg]

				if lbsel=='all':
					labels=data[plate][tp]['Labels']
				elif isinstance(lbsel, list):
					labels=lbsel
				elif isinstance(lbsel, str) or isinstance(lbsel, unicode):
					labels=[lbsel]
				for l in labels:
					if plate!='21' or (plate=='21' and l in data[plate][tp]['Labels'][:9]):
						alldata.append(raw[l])

	
	
	alldata=np.array(alldata)
	mean=alldata.mean(axis=0)
	sd=alldata.std(axis=0)

	return mean, sd		
				 
	


def plot_2D(title,datac,datae,time,labels,genes):
	plate_size=96
	#print title
	plate,fg=title.split('-')
	fig=plt.figure(figsize=(11.69,8.27), dpi=100)
	fig.suptitle(title)
	plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.1, hspace=0.1)
	plots={}

	if fg in ['Fluorescence','Fluorescence_norm']:
		rnd=1
	else:
		rnd=0.1
	totalmaxc=round_to(max([max(datac[l]) for l in labels]),rnd)
	totalmaxe=round_to(max([max(datae[l]) for l in labels]),rnd)
	totalminc=round_to(min([min(datac[l]) for l in labels]),rnd)
	totalmine=round_to(min([min(datae[l]) for l in labels]),rnd)
	if fg in ['Growth','600nm','535nm','Fluorescence','Fluorescence_norm']:
		totalmax=max([totalmaxc,totalmaxe])
		totalmin=0
	elif '_dt' in fg:
		totalmax=0.1#max([totalmaxc,totalmaxe])
		totalmin=-totalmax


	ymin=totalmin
	#print totalmax
	#print list(np.linspace(0, 24, 4)[1:])
	for v,l in IT.izip(range(plate_size),labels):
		row=string.uppercase.index(l[0])+1
		col=int(l.replace(l[0],''))
		#print (row,col)
		if divmod(float(v)/12,1)[1]==0:
			sh_y=l
		#print sh_y
		v=v+1
	
		x=time/3600
		yc=datac[l]
		ye=datae[l]


		if row==1 and col==1:
			plots[l]=plt.subplot(8,12,v)			
		elif row==1 and col!=1:
			plots[l]=plt.subplot(8,12,v,sharey=plots['A1'])
		elif col==1 and row!=1:
			plots[l]=plt.subplot(8,12,v,sharex=plots['A1'])
		else:
			plots[l]=plt.subplot(8,12,v,sharex=plots['A'+str(col)],sharey=plots[l[0]+'1'])
		
		if row!=8:
			setp( plots[l].get_xticklabels(), visible=False)
		if col!=1:
			setp( plots[l].get_yticklabels(), visible=False)

		if col==12:
			plots[l].yaxis.set_label_position("right")
			plt.ylabel(l[0],rotation='horizontal')
		if row==1:
			plt.title(col)

		plt.xticks(np.linspace(0, max(x), 3),['']+list(np.linspace(0, max(x), 3).astype(int)[1:]), rotation='vertical')	
		plt.yticks(np.linspace(ymin, totalmax, 3),['']+list(np.around(np.linspace(totalmin, totalmax, 3),1)[1:]))
		plt.ylim([totalmin,totalmax])
		plots[l].text(0.1, 0.8, genes[l]['Gene'], fontsize=10, transform=plots[l].transAxes)
		#exec(l+"='plt.subplot(8,12,{},sharex={},sharey={})'".format(v,,sh_y))

		plt.plot(x,yc,'r-',x,ye,'b-')
		
		#plt.title(l)

	return plt, plots


def plot_2Dplates(data,fg,tp,means):
	labels=data[data.keys()[0]][tp]['Labels']
	title=fg+'_'+tp
	plate_size=96
	fig=plt.figure(figsize=(11.69,8.27), dpi=100)
	fig.suptitle(title)
	plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.1, hspace=0.1)
	plots={}
	#Need to check for both datasets
	
	if '_dt' in fg:
		time=data[data.keys()[0]][tp]['Time_dt']
	else:
		time=data[data.keys()[0]][tp]['Time']
	x=time/3600


	if fg in ['Growth','600nm']:
		totalmax=0.4
		totalmin=0
	elif fg in ['535nm','Fluorescence']:
		totalmax=300
		totalmin=0
	elif fg in ['Fluorescence_norm']:
		totalmax=3000
		totalmin=0			
	elif fg in ['Fluorescence_dt']:
		totalmax=0.1
		totalmin=-totalmax

	ymin=totalmin
	#print totalmax
	#print list(np.linspace(0, 24, 4)[1:])
	for v,l in IT.izip(range(plate_size),labels):

		row=string.uppercase.index(l[0])+1
		col=int(l.replace(l[0],''))
		#print (row,col)
		if divmod(float(v)/12,1)[1]==0:
			sh_y=l
		#print sh_y
		v=v+1
	


		#setting arrangement
		if row==1 and col==1:
			plots[l]=plt.subplot(8,12,v)			
		elif row==1 and col!=1:
			plots[l]=plt.subplot(8,12,v,sharey=plots['A1'])
		elif col==1 and row!=1:
			plots[l]=plt.subplot(8,12,v,sharex=plots['A1'])
		else:
			plots[l]=plt.subplot(8,12,v,sharex=plots['A'+str(col)],sharey=plots[l[0]+'1'])
		
		if row!=8:
			setp( plots[l].get_xticklabels(), visible=False)
		if col!=1:
			setp( plots[l].get_yticklabels(), visible=False)

		if col==12:
			plots[l].yaxis.set_label_position("right")
			plt.ylabel(l[0],rotation='horizontal')
		if row==1:
			plt.title(col)

		plt.xticks(np.linspace(0, max(x), 3),['']+list(np.linspace(0, max(x), 3).astype(int)[1:]), rotation='vertical')	
		plt.yticks(np.linspace(ymin, totalmax, 3),['']+list(np.around(np.linspace(totalmin, totalmax, 3),1)[1:]))
		plt.ylim([totalmin,totalmax])
		#plots[l].text(0.1, 0.8, genes[l]['Gene'], fontsize=10, transform=plots[l].transAxes)

		if means:
			mean, sd=meansd(data,'all',tp,fg,l,True)
			RMSD=sum(sd)/len(sd)
			plt.fill_between(x, mean - sd, mean + sd, color="red")
			plt.plot(x,mean,'white')
		else:		
			for plate in data.keys():
				if plate!='21' or (plate=='21' and l in labels[:10]):
					y=data[plate][tp]['Shift'][fg][l]
					plt.plot(x,y,'r-')


		
		#plt.title(l)

	return plt, plots

def growthplot(allfit):	
	
	plt.figure(0)
	sbn=np.arange(0,1000,10)
	plt.hist(allfit['Control']['Log']['t0'],bins=sbn,label='Without metformin - Log')
	plt.hist(allfit['Experiment']['Log']['t0'],bins=sbn,label='With metformin - Log')
	plt.xlabel('Start, min')
	plt.ylabel('Number')
	plt.title('Growth start time with/without metformin - Log')
	plt.legend()
	
	plt.figure(1)
	slbn=np.arange(0,1000,10)
	plt.hist(allfit['Control']['Linear']['t0'],bins=sbn,label='Without metformin - Linear')
	plt.hist(allfit['Experiment']['Linear']['t0'],bins=slbn,label='With metformin - Linear')
	plt.xlabel('Start, min')
	plt.ylabel('Number')
	plt.title('Growth start time with/without metformin - Linear')
	plt.legend()

	plt.figure(2)
	abn=np.arange(0,2.5,0.1)
	ac=np.asarray(allfit['Control']['Log']['a'])*3600
	ae=np.asarray(allfit['Experiment']['Log']['a'])*3600
	plt.hist(ac,bins=abn,label='Without metformin - Log')
	plt.hist(ae,bins=abn,label='With metformin - Log')
	plt.xlabel('Growth rate, ln(OD)/h')
	plt.xlim([0,2.5])
	plt.ylabel('Number')
	plt.title('Growth rate with/without metformin - Log')
	plt.legend()

	plt.figure(3)
	albn=np.arange(0,0.04,0.005)
	alc=np.asarray(allfit['Control']['Linear']['a'])*3600
	ale=np.asarray(allfit['Experiment']['Linear']['a'])*3600
	plt.hist(alc,label='Without metformin - Linear')
	plt.hist(ale,label='With metformin - Linear')
	plt.xlabel('Growth rate, OD/h')
	#plt.xlim([0,0.04])
	plt.ylabel('Number')
	plt.title('Growth rate with/without metformin - Linear')
	plt.legend()

	plt.show()


def growth(x,a,c):
	y=x*a+c
	return y

	
	

def cut(x, y, a,b):
	x2=[]
	y2=[]
	for xt, yt in IT.izip(enumerate(x),enumerate(y)):
		#print xt, yt	
		if yt[1]>=a and yt[1]<=b:
			#if y[yt[0]+1]>=a and y[yt[0]+2]>=a: 
			x2.append(xt[1])
			y2.append(yt[1])
	y2=np.asarray(y2)
	x2=np.asarray(x2)
	
	return x2,y2 


def Wiener(y, n):
	wi = sig.wiener(y, mysize=n)
	return wi

def Butter(x, y, par1, par2):	
	b, a = sig.butter(par1, par2)
	fl = sig.filtfilt(b, a, y)
	return fl




def checkfiles(ilist):
	print '\n\n-------------Checking integrity of file list!-------------'
	Pass=False
	allUAL=[fl for fl in ilist if 'UAL_' in fl ]
	Control=[fl for fl in ilist if '_NoMetf_' in fl ]
	Experiment=[fl for fl in ilist if '_Metf_' in fl ]
	CPlates=[fl.split('_')[1] for fl in Control]
	EPlates=[fl.split('_')[1] for fl in Experiment]
	
	if len(allUAL)==len(ilist):
		print 'All files UAL!'
		if len(Control)==len(Experiment):
			print 'Same number of control and experiment files!'
			if compare(CPlates,EPlates):
				print 'Corresponding plate indexes found!'
				print '-----------------File list check passed!-------------\n\n'
				Pass=True
	if not Pass:
		print 'File list integrity check failed!'
		sys.exit(1)

def analyze2(data,filterf):
	msize=20
	par1=4
	par2=0.1
	for plate in sorted(data.keys()):
		for tp in data[plate].keys():
			Ufluor=data[plate][tp]['Shift']['Fluorescence']['C10']
			time=data[plate][tp]['Time']
			dt=time[1]-time[0]
			for well in data[plate][tp]['Labels']:
				if (plate!='21' and well!='A12') or (plate=='21' and well in data[plate][tp]['Labels'][:10]):
					fluor=data[plate][tp]['Shift']['Fluorescence'][well]
					if filterf=='wiener':
						fluor_norm=fluor-Ufluor
						fluor_norm=Wiener(fluor_norm-fluor_norm[0],msize)
					elif filterf=='butter':
						fluor_norm=fluor-Ufluor
						fluor_norm=Butter(time,fluor_norm-fluor_norm[0],par1,par2)
					else:
						fluor_norm=fluor-Ufluor				

					data[plate][tp]['Shift']['Fluorescence_norm'][well]=fluor_norm
					data[plate][tp]['Shift']['Fluorescence_norm_scaled'][well]=scale(fluor_norm,'positive')
					fluor_norm_dt=np.diff(fluor_norm)/dt
					data[plate][tp]['Shift']['Fluorescence_norm_dt'][well]=fluor_norm_dt
					data[plate][tp]['Shift']['Fluorescence_norm_scaled_dt'][well]=scale(fluor_norm_dt,'both')
					

	return data

def analyze(data,waves,filterf):
	msize=20
	par1=4
	par2=0.1
	for plate in sorted(data.keys()):
		for tp in data[plate].keys():		
			if len([nm for nm in waves if nm in data[plate][tp]['Waves']])==2:
				time=data[plate][tp]['Time']
				dt=time[1]-time[0]
				npts=len(time)
				nyf=0.5/dt
			
				data[plate][tp]['Time_dt']=(time+dt/2)[:-1]
				Ufluor=data[plate][tp]['535nm']['C10']/data[plate][tp]['600nm']['C10']
				for well in data[plate][tp]['Labels']:
					growth=data[plate][tp]['600nm'][well]-data[plate][tp]['600nm'][well+'_min']
					fluor=data[plate][tp]['535nm'][well]/data[plate][tp]['600nm'][well]


					if filterf=='wiener':
						growth=Wiener(growth-growth[0],msize)
						fluor=Wiener(fluor-fluor[0],msize)
						fluor_norm=fluor-Ufluor
						fluor_norm=Wiener(fluor_norm-fluor_norm[0],msize)
					elif filterf=='butter':
						growth=Butter(time,growth-growth[0],par1,par2)
						fluor=Butter(time,fluor-fluor[0],par1,par2)
						fluor_norm=fluor-Ufluor
						fluor_norm=Butter(time,fluor_norm-fluor_norm[0],par1,par2)
					else:
						fluor_norm=fluor-Ufluor
								
					#growth=growth-min(growth)
					#fluor=fluor-min(fluor)
					#fluor_norm=fluor_norm-min(fluor_norm)
				

					fluor_max=max(fluor)
					fluor_min=min(fluor)
				
					data[plate][tp]['Growth'][well]=growth
					data[plate][tp]['Growth_scaled'][well]=scale(growth,'positive')
					data[plate][tp]['Fluorescence'][well]=fluor
					data[plate][tp]['Fluorescence_scaled'][well]=scale(fluor,'positive')
					data[plate][tp]['Fluorescence_norm'][well]=fluor_norm
					data[plate][tp]['Fluorescence_norm_scaled'][well]=scale(fluor_norm,'positive')		

					data[plate][tp]['Fluorescence_dt'][well]=np.diff(fluor)/dt
					data[plate][tp]['Fluorescence_scaled_dt'][well]=scale(np.diff(fluor)/dt,'both')
					data[plate][tp]['Fluorescence_norm_dt'][well]=np.diff(fluor_norm)/dt
					data[plate][tp]['Fluorescence_norm_scaled_dt'][well]=scale(np.diff(fluor_norm)/dt,'both')

				
				
					

				data[plate][tp]['Figures']=data[plate][tp]['Figures']+['Growth','Growth_scaled','Fluorescence', 'Fluorescence_scaled','Fluorescence_norm','Fluorescence_dt','Fluorescence_norm_dt','Fluorescence_norm_scaled','Fluorescence_scaled_dt','Fluorescence_norm_scaled_dt'] #'Fluor_norm','Fluor_lp','Fluor_lp_norm','Fluor_dt','Fluor_lp_dt'

	return data

def csvreader(dfile):
	genes=NestedDict()
	rdr=csv.reader(open(dfile,'r'), delimiter=',')
	data=[ln for ln in rdr]
	headers=data[0]
	for ln in data[1:]:
		#print ln
		genes[ln[0]][ln[1]]['Gene']=ln[2]
		genes[ln[0]][ln[1]]['Description']=ln[3]
	
	#print genes
	return genes


def scale(data,tp):	
	start=data[0]
	dmin=min(data)
	dmax=max(data)
	if tp=='positive':
		scaled=(data-dmin)/(dmax-dmin)
	elif tp=='both':
		start_scaled=(start-dmin)/(dmax-dmin)
		scaled=(data-dmin)/(dmax-dmin)+start_scaled

	return scaled


def collect(ilist):
	data=NestedDict()
	for ifl in sorted(ilist):		
		ipt, inm, itp = filename(ifl)	
		plate=inm.split('_')[1]
		if '_NoMetf_' in ifl:
			tp='Control'
		elif '_Metf_' in ifl:
			tp='Experiment'
		print "File: {}\nPlate: {}\nType: {}".format(ifl,plate,tp)
		sheet=readxls(ifl)
		nrows=sheet.nrows
		nm_labels=list(set(sheet.row_values(0)))
		waves=[wlen for wlen in nm_labels if wlen[:2].isdigit()]
		lengths=[sheet.row_values(0).index(wave) for wave in waves]
		length=lengths[1]
		time_t=[int(t.replace('s','')) for t in sheet.row_values(1)[:length]]
		temp=[float(t.split()[0]) for t in sheet.row_values(2)[:length]]
		timemax_h=time_t[length-1]/3600
		timestep=time_t[length-1]/(length-1)
		time=np.linspace(0,timemax_h*3600,length)
		labels=sheet.col_values(length*2)		
		data[plate][tp]['Labels']=labels[3:]
		data[plate][tp]['Waves']=waves
		data[plate][tp]['Time']=time
		data[plate][tp]['Temp']=temp
		data[plate][tp]['Time_max']=timemax_h
		data[plate][tp]['Time_step']=timestep
		data[plate][tp]['Wells']=nrows-3
		data[plate][tp]['Figures']=waves
		data[plate][tp]['File']=inm
		
		print "Wavelengths: {}, {}".format(*waves)
		print "Run time {}h, step {}min in {} wells\n".format(timemax_h,timestep/60, nrows-3)
		for row in range(3,sheet.nrows):
			for wave in waves:
				data_row=[60000 if val=="Overflow" else val for val in sheet.row_values(row)[length*(waves.index(wave)):length*(waves.index(wave)+1)]]
				data[plate][tp][wave][labels[row]]=np.array(data_row)
				data[plate][tp][wave][labels[row]+'_max']=max(data_row)
				data[plate][tp][wave][labels[row]+'_min']=min(data_row)

	return data



def genlist(ifile):
	#Generates list of input files, checks their existance
	ilist=[]
	if ',' in ifile:
		ifiles=ifile.split(',')
		for ifl in ifiles:
			ilist.extend(genlist(ifl))
	else:
		ipath, iname, itype=filename(ifile)
		if itype in ['xls','xlsx'] and os.path.isfile(ifile):
			ilist.append(ifile)
		elif itype in ['txt',''] and os.path.isfile(ifile):
			ifl=open(ifile,'r')
			idata=ifl.read().split('\n')
			idata=[fl.strip() for fl in idata if fl!='']
			for fld in idata:
				ilist.extend(genlist(fld))

		elif iname=='' and itype in ['xls','xlsx']:
			if itype in ['xls','xlsx']:
				ffiles=glob.glob('*.%(itype)s' % vars())
				#print ffiles
				ilist.extend(ffiles)
			elif itype=='txt':
				for tfile in glob.glob('*.%(itype)s' % vars()):
					ilist.extend(genlist(tfile))
			else:
				print "Bad file type %(inp)s!" % vars()	
	return ilist

def axisadjust(ax, xres=4, yres=4):
	"""Send in an axis and I fix the resolution as desired."""

	start, stop = ax.get_xlim()
	ticks = np.arange(start, stop + xres, xres)
	ax.set_xticks(ticks)
	start, stop = ax.get_ylim()
	ticks = np.arange(start, stop + yres, yres)
	ax.set_yticks(ticks)
	return ax


def round_to(n, precission):
	#Round a number to desired precision
	correction = 0.5 if n >= 0 else -0.5
	return int(n/precission+correction)*precission

def tableout(inp):
	#Read file to a list
	ifl=open(inp, 'r')
	idata=ifl.read()
	ifl.close()
	table=[]
	for row in idata.split('\n'):
		data=[it.strip() for it in row.split()]
		data=[numerize(it) for it in data]
		if len(data)!=0:
			table.append(data)

	return table

def makesheets(data,genes):
	sheets=NestedDict()
	plates=data.keys()
	labels=data[plates[0]]['Control']['Labels']
	output=data[plates[0]]['Control']['Figures']
	time_dt=data[plates[0]]['Control']['Time_dt']
	time_lin=data[plates[0]]['Control']['Time']

	for tp in ['Control','Experiment']:
		for out in output:
			sheet=[]
			sheet_sh=[]
			if '_dt' in out:
				time=time_dt
			else:
				time=time_lin
			sheet.append(['Gene']+time.tolist())
			sheet_sh.append(['Gene']+time.tolist())
			for plate in plates:
				for l in labels:
					if (plate!='21' and l!='A12') or (plate=='21' and l in labels[:10]):
						#print plate,tp,l,out
						sheet.append([genes[plate][l]['Gene']]+data[plate][tp][out][l].tolist())
						sheet_sh.append([genes[plate][l]['Gene']]+data[plate][tp]['Shift'][out][l].tolist())

			sheets[tp][out]['Shift']=sheet_sh
			sheets[tp][out]['Raw']=sheet
	
	for out in ['Fluorescence', 'Fluorescence_norm']:
		sheet=[]
		time=time_lin
		sheet.append(['Gene']+time.tolist())
		for plate in plates:
			for l in labels:
				if (plate!='21' and l!='A12') or (plate=='21' and l in labels[:10]):
					sheet.append([genes[plate][l]['Gene']]+(data[plate]['Experiment'][out][l]-data[plate]['Control'][out][l]).tolist())
		sheets['Difference'][out]['Diff']=sheet

	return sheets

def writesheets(sheets,odir):
	#Writes organized data to file.
	#odir=dircheck('Split')
	for tp in sheets.keys():
		for fig in sheets[tp].keys():
			for form in sheets[tp][fig].keys():
				oname='{}/{}_{}_{}.csv'.format(odir,tp,fig,form)		
				sheet=sheets[tp][fig][form]
				f=open(oname,"wb")
				ofile=csv.writer(f, delimiter='\t') # dialect='excel',delimiter=';', quotechar='"', quoting=csv.QUOTE_ALL
				for row in sheet:
					#row = [item.encode("utf-8") if isinstance(item, unicode) else str(item) for item in row]
					ofile.writerow(row)
				f.close()
	

def runcmd(cmd):
	failure, output = commands.getstatusoutput(cmd)
	if failure:
		print '''Running failed \n %(cmd)s \n %(output)s'''.encode('utf-8') % vars();
	return failure, output

class NestedDict(dict):
	def __getitem__(self, key):         
		if key in self: return self.get(key)
		return self.setdefault(key, NestedDict())

def numerize(s):
    try:
	if s=='NAN':
		return s
	float(s)
	if float(s).is_integer():
		return int(float(s))
	elif float(s)==0:
		return float(s)
	else:
        	return float(s)
	
    except ValueError:
        return s
	
def filename(ifile):
	if ifile.split('.')[0]=='':
		ipat=''
		iname=''
		itype=ifile.split('.')[1]
	else:
		if "\\" in ifile.split('.')[0]:
			sep="\\"
		elif "/" in ifile.split('.')[0]:
			sep="/"
		else:
			ipat=''
			iname=ifile.split('.')[0]
			itype=ifile.split('.')[1]
			return ipat, iname, itype
		allpath=ifile.split('.')[0]	
		iname=allpath.split(sep)[-1]
		ipath=allpath.split(sep)[:-1]
		ipat='/'.join(ipath)
		itype=ifile.split('.')[1]
	return ipat, iname, itype

def dircheck(somedir):
	while True:   	
		if os.path.exists(somedir):
			qstn = "Directory %(somedir)s already exists! Delete, quit, continue or provide a new name (d/q/c/<type name>): " % vars()
			answ = raw_input(qstn)
			if answ == "d":
				shutil.rmtree(somedir)
				os.makedirs(somedir)
				break
			elif answ == "q":
				sys.exit("Have a nice day!")
			elif answ == "c":
				break
			else:
				somedir=answ
				continue
		else:
			os.makedirs(somedir)
			break
	return somedir


def readxls(ifile):
	book=xlrd.open_workbook(ifile,formatting_info=False)
	sheet=book.sheet_by_name([nm for nm in book.sheet_names() if 'Magellan' in nm][0])
	return sheet


#----------------------------------

if __name__ == "__main__":
	sys.exit(main())

