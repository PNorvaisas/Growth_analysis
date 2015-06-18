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
	print "Please install pip manually to proceed!"
	sys.exit(1)
def install(package):
	pip.main(['install', package])

for mod in ['pip','string','math','re','csv','sys','os','commands','datetime','operator','getopt','subprocess','pickle','shutil','glob','types','math','copy','pyExcelerator','xlrd','xlwt','xlutils','types']:
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
		genes=genereader(dfile)
		ilist=genlist(ifile)

		checkfiles(ilist)	
		data=collect(ilist)


		data=analyze(data,filterf)
		dirn=dircheck(odir)
		data,allfit=growthfit(data,False)
		growthplot(allfit)		
		data=datashift(data,'all','all','all','all')
		data=analyze_shift(data,filterf)
		data=differences(data)	
		sheets=makesheets(data,genes)
		writesheets(sheets,odir)
		f = open('UAL_data.pckl', 'w')
		pickle.dump([data,genes,odir,filterf], f)
		f.close()
	

	#data=differences(data)	
	#sheets=makesheets(data,genes)
	#writesheets(sheets,odir)
	#sheets=makesheets(data,genes)
	#writesheets(sheets,odir)
	#data,allfit=growthfit(data,False)
	#growthplot(allfit)
	#data=datashift(data,'all','all','all','all')
	#data=analyze2(data,filterf)
	#plotall(data,'all','Experiment','Growth','C3',False,False)
	#plotall(data,'all','Control','Growth','C3',False,False)
	#plotall(data,'all','Experiment','Growth','C3',True,False)
	#plotall(data,'all','Control','Growth','C3',True,False)
	#plot_comparison(data,genes,odir,True,'all')
	plot_comparison(data,genes,odir,True,'all') #['Growth_dt','Fluorescence_norm_log10']
	#plotall(data,'all','Experiment','Growth','C3',False,False)
	#plotall(data,'all','Control','Growth','C3',False,False)

	plotall(data,'all','Experiment','Growth_dt','C3',True,False)
	plotall(data,'all','Control','Growth_dt','C3',True,False)

	#plotall(data,'all','Experiment','Fluorescence','C3',False,False)
	#plotall(data,'all','Control','Fluorescence','C3',False,False)

	plotall(data,'all','Experiment','Fluorescence_norm_log10','C3',True,False)
	plotall(data,'all','Control','Fluorescence_norm_log10','C3',True,False)

	plotall(data,'all','Experiment','Fluorescence_norm_log10_dt','C3',True,False)
	plotall(data,'all','Control','Fluorescence_norm_log10_dt','C3',True,False)



	#plotall(data,'all','Experiment','535nm','F10',True,False)
	#plotall(data,'all','Control','535nm','F10',True,False)

	#plotall(data,'all','Experiment','Growth','all',True,False) #plotall(data,plsel,tpsel,fg,lsel,shifted,norm)
	#plotall(data,'all','Control','Growth','all',True,False) #plotall(data,plsel,tpsel,fg,lsel,shifted,norm)


	#plt,plots=plot_2Dplates(data,'Experiment','Growth',True,True) #plot_2Dplates(data,tp,fg,shifted,means):
	#plt.show()




	

	

	
	
	
	

#-------------Functions------------


def plot_comparison(data,genes,dirn,shifted,figs):

	for plate in sorted(data.keys()):
		if shifted:
			ref=data[plate]['Control']['Shift']
			exp=data[plate]['Experiment']['Shift']
			shift='_Shift'
			figures_temp=data[plate]['Control']['Shift']['Figures']
		else:
			ref=data[plate]['Control']
			exp=data[plate]['Experiment']
			shift='_Raw'
			figures_temp=data[plate]['Control']['Figures']	

		if figs=='all':
			figures=figures_temp
				
		elif isinstance(figs, list):
			figs_ch=[f for f in figs if f in figures_temp]
			if len(figs)>0:
				figures=figs_ch
			else:
				print 'Figures {} not found'.format(figs_ch)
		elif isinstance(figs, str) or isinstance(figs, unicode):
			if figs in figures_temp:
				figures=[figs]
			else:
				print 'Figures {} not found'.format(figs)
		for fg in figures:
			print "Plotting plate {} {}...".format(plate,fg+shift)
			if 'dt' in fg:
				plot=plot_2D(plate+"-"+fg,ref[fg],exp[fg],data[plate]['Control']['Time_dt'], data[plate]['Control']['Labels'], genes[plate])
			else:
				plot=plot_2D(plate+"-"+fg,ref[fg],exp[fg],data[plate]['Control']['Time'], data[plate]['Control']['Labels'], genes[plate])
			#data[plate]['Joint'][fg]=plt
			plot.savefig('{}/{}.pdf'.format(dirn,plate+'_'+fg+shift))
			plot.close()
	return data

def plotall(data,plsel,tpsel,fg,lsel,shifted,norm):
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
			elif tp=='Difference':
				plt.figure(0)
				plt.title(fg+' difference')
			elif tp=='Ratio':
				plt.figure(0)
				plt.title(fg+' ratio')

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
		base=0.01
		top=0.05
		y0=0.005
	
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
	ts=3
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
						#print shift
						if shift<0:
							steps=int(-shift/300)
							#print steps
							if steps>0:
								y2[-steps:]=y2[-steps-1]
							#print y2
						elif shift>0:
							steps=int(shift/300)
							#print steps
							if steps>0:
								y2[:steps]=y2[steps+1]	
							#print y2
							
							
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
	xmax=20
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
	#if fg in ['Growth','600nm','535nm','Fluorescence','Fluorescence_norm']:
	totalmax=max([totalmaxc,totalmaxe])
	totalmin=0
	ticks=3
	if fg=='Growth':
		totalmax=0.4
		ticks=4
	if ('_dt' in fg) and fg!='Growth_dt':
		totalmax=0.1#max([totalmaxc,totalmaxe])
		totalmin=-totalmax
		ticks=4
	if fg=='Growth_dt':
		totalmax=0.000015
		totalmin=0
		ticks=3

	if fg=='Fluorescence_norm_log10':
		totalmax=3
		totalmin=0
		ticks=4

	if fg=='Fluorescence_norm_log10_dt':
		totalmax=0.0005
		totalmin=-totalmax/2
		ticks=4



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

		plt.xticks(np.linspace(0, xmax, 3),['']+list(np.linspace(0, xmax, 3).astype(int)[1:]), rotation='vertical')	
		plt.yticks(np.linspace(ymin, totalmax, ticks),['']+list(np.around(np.linspace(totalmin, totalmax, ticks),1)[1:]))
		plt.ylim([totalmin,totalmax])
		plots[l].text(0.1, 0.8, genes[l]['Gene'], fontsize=10, transform=plots[l].transAxes)
		#exec(l+"='plt.subplot(8,12,{},sharex={},sharey={})'".format(v,,sh_y))

		plt.plot(x,yc,'r-',x,ye,'b-')
		
		#plt.title(l)

	return plt


def plot_2Dplates(data,tp,fg,shifted,means):
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
			mean, sd=meansd(data,'all',tp,fg,l,shifted)
			RMSD=sum(sd)/len(sd)
			plt.fill_between(x, mean - sd, mean + sd, color="red")
			plt.plot(x,mean,'white')
		else:		
			for plate in data.keys():
				if plate!='21' or (plate=='21' and l in labels[:10]):
					if shifted:
						y=data[plate][tp]['Shift'][fg][l]
					else:
						y=data[plate][tp][fg][l]
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
	print '\n\n-------------Checking integrity of the file list!-------------'
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


def differences(data):
	
	for plate in sorted(data.keys()):
		labels=data[plate]['Control']['Labels']
		#for fg in data[plate]['Control']['Figures']:
		#	#print plate, fg
		#	for well in data[plate]['Control']['Labels']:
		#		data[plate]['Differences']['Shift'][fg][well]=data[plate]['Experiment']['Shift'][fg][well]-data[plate]['Control']['Shift'][fg][well]				
		#		data[plate]['Ratios']['Shift'][fg][well]=data[plate]['Experiment']['Shift'][fg][well]/data[plate]['Control']['Shift'][fg][well]
		for fg in data[plate]['Control']['Shift']['Figures']:
			if fg in ['Fluorescence','Fluorescence_norm']:
				bar=10
			elif fg in ['Growth']:
				bar=0.1			
			else:
				bar=0
			for well in labels:
				#print plate, fg, well
				#if (plate!='21' and well!='A12') or (plate=='21' and well in labels[:10]):
				data[plate]['Differences'][fg][well]=data[plate]['Experiment']['Shift'][fg][well]-data[plate]['Control']['Shift'][fg][well]				
				data[plate]['Ratios'][fg][well]=ratio(data[plate]['Experiment']['Shift'][fg][well],data[plate]['Control']['Shift'][fg][well],bar)
				data[plate]['LogRatios'][fg][well]=np.log10(ratio(data[plate]['Experiment']['Shift'][fg][well],data[plate]['Control']['Shift'][fg][well],bar))
		data[plate]['Differences']['Figures']=data[plate]['Control']['Shift']['Figures']
		data[plate]['Ratios']['Figures']=data[plate]['Control']['Shift']['Figures']
		data[plate]['LogRatios']['Figures']=data[plate]['Control']['Shift']['Figures']
	return data

def analyze_shift(data,filterf):
	msize=20
	par1=4
	par2=0.1
	bar=10
	window=10
	for plate in sorted(data.keys()):
		for tp in data[plate].keys():
			if plate!='21':
				Ufluor=data[plate][tp]['Shift']['Fluorescence']['C10']
			else:
				Ufluor=data['20'][tp]['Shift']['Fluorescence']['C10']
			time=data[plate][tp]['Time']
			dt=time[1]-time[0]
			for well in data[plate][tp]['Labels']:
				#if (plate!='21' and well!='A12') or (plate=='21' and well in data[plate][tp]['Labels'][:10]):
				fluor=data[plate][tp]['Shift']['Fluorescence'][well]
				fluor_norm=fluor-Ufluor
				fluor_norm=fluor_norm-np.mean(fluor_norm[:window])
				fluor_U139=ratio(fluor,Ufluor,bar)
				#if filterf=='wiener':
					#fluor_U139=Wiener(fluor_U139-fluor_U139[0],msize)
				#	fluor_norm=Wiener(fluor_norm-fluor_norm[0],msize)
				#elif filterf=='butter':
					#fluorU139=Butter(time,fluor_U139-fluor_U139[0],par1,par2)
				#	fluor_norm=Butter(time,fluor_norm-fluor_norm[0],par1,par2)

		
				fluor_norm_dt=np.diff(fluor_norm)/dt
				fluor_norm_log=log10(setbar(fluor_norm,bar))
				fluor_norm_log_dt=np.diff(fluor_norm_log)/dt
				#fluor_U139_dt=np.diff(fluor_U139)/dt
				data[plate][tp]['Shift']['Fluorescence_norm'][well]=fluor_norm
				data[plate][tp]['Shift']['Fluorescence_norm_log10'][well]=fluor_norm_log
				data[plate][tp]['Shift']['Fluorescence_norm_log10_dt'][well]=fluor_norm_log_dt
				data[plate][tp]['Shift']['Fluorescence_norm_dt'][well]=fluor_norm_dt
				data[plate][tp]['Shift']['Fluorescence_U139'][well]=fluor_U139
				#data[plate][tp]['Shift']['Fluorescence_U139_dt'][well]=fluor_U139_dt

					
					#print plate,tp,well   ,'Fluorescence_U139','Fluorescence_U139_dt'
			data[plate][tp]['Shift']['Figures']=data[plate][tp]['Figures']+['Fluorescence_norm','Fluorescence_norm_log10','Fluorescence_norm_log10_dt','Fluorescence_norm_dt','Fluorescence_U139']

	return data


def setbar(x,bar):
	x2=[xi if xi>=bar else bar for xi in x]
	x2=np.array(x2)
	return x2

def ratio(x,y,bar):
	xy=[]
	if len(x)==len(y):
		for xi,yi in IT.izip(x,y):
			if xi>bar and yi>bar:
				xy.append(xi/yi)
			else:
				xy.append(1)

	xy=np.array(xy)
	return xy



def analyze(data,filterf):
	waves=['600nm','535nm']
	msize=20
	par1=4
	par2=0.1
	method='pre'
	window=10
	for plate in sorted(data.keys()):
		for tp in data[plate].keys():		
			if len([nm for nm in waves if nm in data[plate][tp]['Waves']])==2:
				time=data[plate][tp]['Time']
				dt=time[1]-time[0]
				npts=len(time)
				nyf=0.5/dt
			
				data[plate][tp]['Time_dt']=(time+dt/2)[:-1]
				#Ufluor=data[plate][tp]['535nm']['C10']/data[plate][tp]['600nm']['C10']
				
				for well in data[plate][tp]['Labels']:
					gs=np.mean(data[plate][tp]['600nm'][well][:window])
					fs=np.mean(data[plate][tp]['535nm'][well][:window])

					if method=='pre':						
					
						if filterf=='wiener':
							growth=Wiener(data[plate][tp]['600nm'][well]-gs,msize)
							fluor=(Wiener(data[plate][tp]['535nm'][well]-fs,msize)+fs)/(Wiener(data[plate][tp]['600nm'][well]-gs,msize)+gs)
							fluor=fluor-np.mean(fluor[:window])
						elif filterf=='butter':
							growth=Butter(time,data[plate][tp]['600nm'][well]-gs,par1,par2)
							fluor=(Butter(time,data[plate][tp]['535nm'][well]-fs,par1,par2)+fs)/(Butter(time,data[plate][tp]['600nm'][well]-gs,par1,par2)+gs)
							fluor=fluor-np.mean(fluor[:window])
						else:
							fluor_norm=fluor-Ufluor


					elif method=='post':
						growth=data[plate][tp]['600nm'][well]-gs
						fluor=data[plate][tp]['535nm'][well]/data[plate][tp]['600nm'][well]
						
					
						if filterf=='wiener':
							growth=Wiener(growth,msize)
							fluor=Wiener(fluor-np.mean(fluor[:window]),msize)
							#fluor_norm=fluor-Ufluor
							#fluor_norm=Wiener(fluor_norm-fluor_norm[0],msize)
						elif filterf=='butter':
							growth=Butter(time,growth,par1,par2)
							fluor=Butter(time,fluor-np.mean(fluor[:window]),par1,par2)
							#fluor_norm=fluor-Ufluor
							#fluor_norm=Butter(time,fluor_norm-fluor_norm[0],par1,par2)
						else:
							fluor_norm=fluor-Ufluor
								
					#growth=growth-min(growth)
					#fluor=fluor-min(fluor)
					#fluor_norm=fluor_norm-min(fluor_norm)
			
				
					data[plate][tp]['Growth'][well]=growth	
					data[plate][tp]['Growth_dt'][well]=np.diff(growth)/dt				
					data[plate][tp]['Fluorescence'][well]=fluor
					#data[plate][tp]['Fluorescence_norm'][well]=fluor_norm	
					data[plate][tp]['Fluorescence_dt'][well]=np.diff(fluor)/dt
					#data[plate][tp]['Fluorescence_norm_dt'][well]=np.diff(fluor_norm)/dt
				
				
					

				data[plate][tp]['Figures']=data[plate][tp]['Figures']+['Growth','Growth_dt','Fluorescence','Fluorescence_dt']

	return data

def genereader(dfile):
	genes=NestedDict()
	rdr=csv.reader(open(dfile,'r'), delimiter=',')
	data=[ln for ln in rdr]
	headers=data[0]
	for ln in data[1:]:
		#print ln
		genes[ln[0]][ln[1]]['Gene']=ln[2]
		genes[ln[0]][ln[1]]['Description']=ln[3]
		if isinstance(genes['Genes'][ln[2]]['Address'], list):
			genes['Genes'][ln[2]]['Address']=genes['Genes'][ln[2]]['Address']+[[ln[0],ln[1]]]
		else:
			genes['Genes'][ln[2]]['Address']=[[ln[0],ln[1]]]

		if isinstance(genes['Genes'][ln[2]]['Description'], str):
			continue
		else:
			genes['Genes'][ln[2]]['Description']=ln[3]
	#print genes
	return genes


def scale(data,tp):	
	start=data[0]
	dmin=min(data)
	dmax=max(data)
	if tp=='start':
		if (dmax-start)!=0:
			scaled=(data-start)/(dmax-start)
		else:
			print "Zero divider: {}-{}".format(dmax,start)
			scaled=data
	elif tp=='min':
		if (dmax-dmin)!=0:
			scaled=(data-dmin)/(dmax-dmin)
		else:
			print "Zero divider: {}-{}".format(dmax,dmin)
			scaled=data

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
	
	time_dt=data[plates[0]]['Control']['Time_dt']
	time_lin=data[plates[0]]['Control']['Time']

	for tp in ['Differences','Ratios','LogRatios','Control','Experiment']:
		for algn in ['Raw','Shift']:
			if algn=='Shift' and tp in ['Control','Experiment']:
				output=data[plates[0]][tp]['Shift']['Figures']
			else:
				output=data[plates[0]][tp]['Figures']	
			for out in output:			
				sheet=[]

				if '_dt' in out:
					time=time_dt
				else:
					time=time_lin
				sheet.append(['Gene']+time.tolist())
				#print tp, algn, out
				#print output
			
				for plate in plates:
					#print data[plate][tp].keys()
					for l in labels:
						if (plate!='21' and l!='A12') or (plate=='21' and l in labels[:10]):
							if algn=='Shift' and tp in ['Control','Experiment']:
								sheet.append([genes[plate][l]['Gene']]+data[plate][tp]['Shift'][out][l].tolist())
							else:
								sheet.append([genes[plate][l]['Gene']]+data[plate][tp][out][l].tolist())


					sheets[tp][out][algn]=sheet


	##Implement Difference art in analysis function
	#for out in ['Fluorescence', 'Fluorescence_norm']:
	#	sheet=[]
	#	time=time_lin
	#	sheet.append(['Gene']+time.tolist())
	#	for plate in plates:
	#		for l in labels:
	#			if (plate!='21' and l!='A12') or (plate=='21' and l in labels[:10]):
	#				sheet.append([genes[plate][l]['Gene']]+(data[plate]['Experiment'][out][l]-data[plate]['Control'][out][l]).tolist())
	#	sheets['Difference'][out]['Diff']=sheet

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

