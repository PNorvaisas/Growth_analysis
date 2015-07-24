#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Biolog.py

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
import unicodedata
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rc


#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
#from pylab import *

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
import textwrap as tw
from multiprocessing import Pool


compare = lambda x, y: Counter(x) == Counter(y)

try:
	import itertools as IT
except ImportError, e:
	print "Module itertools not found"





help_message = '''
Biolog data preparation

Flags:
	-h  Display this help message
	-v  Verbose output
Arguments:
	-i <files>   Input files, list separated by comma, file with list
	-d <file>    CSV file with data labels
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
        Out:	%(odir)s
       Load:	%(load)s
<--------------------------------------------->
	'''



def main(argv=None):

	ifile=""
	dfile=""
	msize=20
	odir='Output'
	load=False
	comparison=False
	
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
				metabolites=labelreader(dfile)		
			else:
				print "No database file specified."
				metabolites=NestedDict()
				
			
	
		except Exception, e:
			print e
			#print "!!!-------------------------------------------!!!"
			sys.exit(1)

		print optionsset %vars()
	#----------------------------


	


	if load:
		f = open('Biolog_data.pckl','rb')
		data = pickle.load(f)
		f.close()

		
	else:	
		
		ilist=genlist(ifile)
		print ilist
		#print metabolites
				
		
		#checkfiles(ilist)	
		data=collect(ilist)		
		data=analyze(data)
		dirn=dircheck(odir)	
		
		data,allfit=growthfit(data)
		#growthplot(allfit)		
		plot_comparison(data,metabolites,odir,['Growth','Respiration','Growth_dt','Respiration_dt','Growth_abolished','Respiration_abolished']) #['Growth','Respiration','Growth_dt','Respiration_dt']
	
		
		#sheets=makesheets(data,genes)
		#writesheets(sheets,odir)
		#f = open('Biolog_data.pckl', 'w')
		#pickle.dump(data, f)
		#f.close()
	



	
	
	
	

#-------------Functions------------



def plot_comparison(data,metabolites,dirn,figs):

	for plate in sorted(data.keys()):
		for strain in data[plate].keys():
			labels=data[plate][strain]['NoMetf']['Labels']
			ref=data[plate][strain]['NoMetf']
			exp=data[plate][strain]['Metf']
			figures_temp=data[plate][strain]['NoMetf']['Figures']
			figures_temp=figures_temp+['Growth_abolished','Respiration_abolished']

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

			#print figures
			for fg in figures:
				if fg=='Growth_abolished':
					fgl='Growth'
				elif fg=='Respiration_abolished':
					fgl='Respiration'
				else:
					fgl=fg
				print "Plotting plate {} {} {}...".format(plate,strain,fg)
				if '_dt' in fg:
					time=data[plate][strain]['NoMetf']['Time_dt']
				else:
					time=data[plate][strain]['NoMetf']['Time']
				#Need to fix metabolites
				plot=plot_2D('{} {} {}'.format(plate,strain,fg),ref[fgl],exp[fgl],time, labels, metabolites[plate])
				plot.savefig('{}/{}.pdf'.format(dirn,plate+'-'+strain+'_'+fg))
				plot.close()
	#
	return data

def plot_2D(title,datac,datae,time,labels,metabolites):
	
	xmax=24
	plate_size=96
	#print title
	plate,strain,fg=title.split()
	#fig=plt.figure(figsize=(11.69,8.27), dpi=100)
	
	fig,axes=plt.subplots(nrows=8, ncols=12, sharex=True, sharey=True,figsize=(11.69,8.27), dpi=100)
	fig.suptitle(title)
	plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.05, hspace=0.05)
	

	if fg in ['Fluorescence','Fluorescence_norm']:
		rnd=1
	else:
		rnd=0.1
	totalmaxc=round_to(max([max(datac[l]) for l in labels]),rnd)
	totalmaxe=round_to(max([max(datae[l]) for l in labels]),rnd)
	totalminc=round_to(min([min(datac[l]) for l in labels]),rnd)
	totalmine=round_to(min([min(datae[l]) for l in labels]),rnd)
	totalmax=max([totalmaxc,totalmaxe])
	totalmin=0
	ticks=3
	xlabel='Time, h'
	ylabel=''
	decimals=1
	if fg=='590nm':
		totalmax=1
		ticks=3
		ylabel='OD@590nm'
		decimals=1
	if fg=='750nm':
		totalmax=2
		ticks=4
		ylabel='OD@750nm'
		decimals=1
	
	if fg=='Growth':
		totalmax=1
		ticks=3
		ylabel='Growth@590nm'
		decimals=1
	if fg=='Growth_dt':
		totalmax=0.2
		totalmin=0
		ticks=3
		ylabel='Growth@590nm/dt'
		decimals=1

	if fg=='Dye':
		totalmax=1.6
		ticks=3
		ylabel='Dye@750nm'
		decimals=1
	if fg=='Dye_dt':
		totalmax=0.3
		totalmin=0
		ticks=4
		ylabel='Dye@750nm/dt'
		decimals=1

	if fg=='Respiration':
		totalmax=2
		ticks=3
		ylabel='Respiration, Dye/Growth'
		decimals=1

	if fg=='Respiration_dt':
		totalmax=0.3
		totalmin=0
		ticks=4
		ylabel='Respiration, Dye/Growth/dt'
		decimals=1

	if fg=='Growth_abolished':
		totalmax=100
		totalmin=-100
		ticks=3
		decimals=1
		ylabel='Growth abolished, %'
		thres=0.05

	if fg=='Respiration_abolished':
		totalmax=100
		totalmin=-100
		ticks=3
		decimals=1
		ylabel='Respiration abolished, %'
		thres=0.15





	ymin=totalmin
	fig.text(0.5, 0.04, xlabel, ha='center')
	fig.text(0.04, 0.5, ylabel, va='center', rotation='vertical')

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
		plt.sca(axes[row-1,col-1])
		ax=axes[row-1,col-1]
		if col==12:
			ax.yaxis.set_label_position("right")
			plt.ylabel(l[0],rotation='horizontal')
		if row==1:
			plt.title(col)    				

		if col>1 and row<11:
			plt.setp(ax.get_yticklabels(), visible=False)

		plt.xticks(np.linspace(0, xmax, 3),['']+list(np.linspace(0, xmax, 3).astype(int)[1:]), rotation='vertical')	
		plt.yticks(np.linspace(ymin, totalmax, ticks),['']+list(myround(np.linspace(totalmin, totalmax, ticks),decimals)[1:]))
		plt.ylim([totalmin,totalmax])
		#plt.axvline(x=3, ymin=0, ymax=1,c='green', hold=None)  u 
		label=greek_check(metabolites[l]['Metabolite'],12)
		plt.text(0.05, 0.9, label, fontsize=7,verticalalignment='top',transform=ax.transAxes)
		if fg in ['Growth_abolished','Respiration_abolished']:
			ye[ye<thres]=thres
			yc[yc<thres]=thres
			yd=(ye*100/yc)-100			
			plt.plot(x,yd,'k-')
			plt.fill_between(x, 0, yd, where=yd>=0, facecolor='blue', interpolate=True)
			plt.fill_between(x, 0, yd, where=yd<=0, facecolor='red', interpolate=True)
		else:
			plt.plot(x,yc,'r-',x,ye,'b-')

	return plt

def myround(a, decimals=1):
     return np.around(a-10**(-(decimals+5)), decimals=decimals)

def greek_check(text,slen):
	text=unicode(text)
	greek_alphabet ={'alpha':u'\u03B1','beta':u'\u03B2','gamma':u'\u03B3','delta':u'\u03B4','epsilon':u'\u03B5'}
	greek_alphabet2 ={'alpha':u'α','beta':u'β','gamma':u'γ','delta':u'δ','epsilon':u'ε'}	
	tex_alphabet ={u'α':u'$α$',u'β':u'$β$',u'γ':u'$γ$',u'δ':u'$δ$',u'ε':u'$ε$'}
	uni_alphabet ={'alpha':u'α','beta':u'β','gamma':u'γ','delta':u'δ','epsilon':u'ε'}
	for k in uni_alphabet.keys():
		if k in text:
			text=text.replace(k,uni_alphabet[k])
	text='\n'.join(tw.wrap(text,slen))
	for k in tex_alphabet.keys():
		if k in text:
			text=text.replace(k,tex_alphabet[k])
	return text

def growthfit(data):
	allfit=NestedDict()
	base=0.01
	top=0.05
	y0=0.005

	for plate in data.keys():
		for strain in data[plate].keys():
			for tp in data[plate][strain].keys():
				x=data[plate][strain][tp]['Time']
				labels=data[plate][strain][tp]['Labels']
	
				for l in labels:

					y=data[plate][strain][tp]['Growth'][l]
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
				
						data[plate][strain][tp]['GrowthFit'][l][ynm]=[a,c,t0]

	return data,allfit


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

def growth(x,a,c):
	y=x*a+c
	return y

def Wiener(y, n):
	wi = sig.wiener(y, mysize=n)
	return wi

def Butter(x, y, par1, par2):	
	b, a = sig.butter(par1, par2)
	fl = sig.filtfilt(b, a, y)
	return fl

def growthplot(allfit):	
	
	plt.figure(0)
	sbn=np.arange(0,1000,10)
	plt.hist(allfit['NoMetf']['Log']['t0'],bins=sbn,label='Without metformin - Log')
	plt.hist(allfit['Metf']['Log']['t0'],bins=sbn,label='With metformin - Log')
	plt.xlabel('Start, min')
	plt.ylabel('Number')
	plt.title('Growth start time with/without metformin - Log')
	plt.legend()
	
	plt.figure(1)
	slbn=np.arange(0,1000,10)
	plt.hist(allfit['NoMetf']['Linear']['t0'],bins=sbn,label='Without metformin - Linear')
	plt.hist(allfit['Metf']['Linear']['t0'],bins=slbn,label='With metformin - Linear')
	plt.xlabel('Start, min')
	plt.ylabel('Number')
	plt.title('Growth start time with/without metformin - Linear')
	plt.legend()

	plt.figure(2)
	abn=np.arange(0,2.5,0.1)
	ac=np.asarray(allfit['NoMetf']['Log']['a'])*3600
	ae=np.asarray(allfit['Metf']['Log']['a'])*3600
	plt.hist(ac,bins=abn,label='Without metformin - Log')
	plt.hist(ae,bins=abn,label='With metformin - Log')
	plt.xlabel('Growth rate, ln(OD)/h')
	plt.xlim([0,2.5])
	plt.ylabel('Number')
	plt.title('Growth rate with/without metformin - Log')
	plt.legend()

	plt.figure(3)
	albn=np.arange(0,0.04,0.005)
	alc=np.asarray(allfit['NoMetf']['Linear']['a'])*3600
	ale=np.asarray(allfit['Metf']['Linear']['a'])*3600
	plt.hist(alc,label='Without metformin - Linear')
	plt.hist(ale,label='With metformin - Linear')
	plt.xlabel('Growth rate, OD/h')
	#plt.xlim([0,0.04])
	plt.ylabel('Number')
	plt.title('Growth rate with/without metformin - Linear')
	plt.legend()

	plt.show()


def checkfiles(ilist):
	
	print '\n\n-------------Checking integrity of the file list!-------------'
	Pass=False
	allUAL=[fl for fl in ilist if 'UAL_' in fl ]
	Control=[fl for fl in ilist if '_NoMetf_' in fl ]
	Experiment=[fl for fl in ilist if '_Metf_' in fl ]
	CPlates=[fl.split('_')[1] for fl in Control]
	EPlates=[fl.split('_')[1] for fl in Experiment]
	if comparison:
		Pass=True
	else:	
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



def analyze(data):
	filterf='wiener'
	waves=['590nm','750nm']
	msize=20
	par1=4
	par2=0.1
	method='pre'
	window=10
	for plate in sorted(data.keys()):
		for strain in data[plate].keys():
			for tp in data[plate][strain].keys():				
				time=data[plate][strain][tp]['Time']
				dt=time[1]-time[0]
				npts=len(time)
				nyf=0.5/dt
				print plate,strain,tp
				data[plate][strain][tp]['Time_dt']=(time+dt/2)[:-1]
				for well in data[plate][strain][tp]['Labels']:
					gs=np.mean(data[plate][strain][tp]['590nm'][well][:window])
					ds=np.mean(data[plate][strain][tp]['750nm'][well][:window])
					
					growth=Wiener(data[plate][strain][tp]['590nm'][well]-gs,msize)
					dye=Wiener(data[plate][strain][tp]['750nm'][well]-ds,msize)
					resp=(dye+ds)/(growth+gs)
	
					data[plate][strain][tp]['Dye'][well]=dye
					data[plate][strain][tp]['Dye_dt'][well]=np.diff(dye)/(dt/3600)	
					data[plate][strain][tp]['Growth'][well]=growth	
					data[plate][strain][tp]['Growth_dt'][well]=np.diff(growth)/(dt/3600)
					data[plate][strain][tp]['Respiration'][well]=resp

					refresp=data[plate][strain][tp]['Respiration']['A1']

					if plate=='PM5':
						rs=np.mean(resp[:window])			
						resp=resp-rs
					else:
						resp=resp-refresp
						rs=np.mean(resp[:window])
						resp=resp-rs	
					gresp=(growth+gs+0.1)/(resp+rs+0.1)
					respg=(resp+rs+0.1)/(growth+gs+0.1)
					gresp=gresp-np.mean(gresp[:window])
					respg=respg-np.mean(respg[:window])
					data[plate][strain][tp]['Respiration'][well]=resp	
					data[plate][strain][tp]['Respiration_dt'][well]=np.diff(resp)/(dt/3600)
					#data[plate][strain][tp]['Growth-Resp'][well]=gresp
					#data[plate][strain][tp]['Resp-Growth'][well]=respg
			
				
					

				data[plate][strain][tp]['Figures']=data[plate][strain][tp]['Figures']+['Growth','Growth_dt','Respiration','Respiration_dt','Dye','Dye_dt']

	return data

def labelreader(dfile):
	metabolites=NestedDict()
	rdr=csv.reader(open(dfile,'r'), delimiter=';')
	data=[ln for ln in rdr]
	headers=data[0]
	for ln in data[1:]:
		ln=[cell.strip() for cell in ln]
		#print ln
		metabolites[ln[0]][ln[1]]['Metabolite']=ln[3]
		metabolites[ln[0]][ln[1]]['Index']=ln[2]
		
	return metabolites



def collect(ilist):
	data=NestedDict()
	for ifl in sorted(ilist):		
		ipt, inm, itp = filename(ifl)	
		plate=inm.split('_')[2]
		strain=inm.split('_')[3]
		tp=inm.split('_')[4]
		print "File: {}\nPlate: {}\nStrain: {}\nType: {}".format(ifl,plate,strain,tp)
		if itp=='xlsx':
			sheet=readxls(ifl)
		elif itp=='asc':
			sheet=readacs(ifl)
		else:
			print 'Unknown file format'


		nrows=len(sheet)
		nm_labels=list(set(sheet[0]))
		#print nm_labels
		#print nrows
		
		waves=[numerize(wlen) for wlen in nm_labels if wlen not in ['Layout','Well positions','','Replicate Info']]
		#print waves
		lengths=[sheet[0].index(wave) for wave in waves]
		#Selection of time cells does not depend om the order of wavelengths		
		length=max(lengths)
		#print waves, lengths
		time_row=sheet[1][:length]
		check=list(set([s for s in time_row if isinstance(s,float)]))
		if check:
			time_t=time_row
		else:
			time_t=[int(t.replace('s','')) for t in time_row]
		temp=[float(t.split()[0]) for t in sheet[2][:length]]
		timemax_h=round_to(float(time_t[length-1])/3600,1)
		
		timestep=round_to(float(time_t[length-1])/(length-1),1)
		time=np.linspace(0,timemax_h*3600,length)
		tsheet=zip(*sheet)
		labels=[l for l in tsheet[length*2] if l not in ['','Well positions']	] #Sheet
		data[plate][strain][tp]['Labels']=labels
		data[plate][strain][tp]['Spectra']=waves
		data[plate][strain][tp]['Time']=time
		data[plate][strain][tp]['Temp']=temp
		data[plate][strain][tp]['Time_max']=timemax_h
		data[plate][strain][tp]['Time_step']=timestep
		data[plate][strain][tp]['Wells']=nrows-3
		data[plate][strain][tp]['Figures']=[str(w)+'nm' for w in waves]
		data[plate][strain][tp]['File']=inm
		print "Wavelengths: {}, {}".format(*waves)
		print "Run time {}h, step {}min in {} wells\n".format(timemax_h,timestep/60, len(labels))
		for row in range(0,len(labels)):
			for wave in waves:
				data_row=[60000 if val=="Overflow" else val for val in sheet[row+3][length*(waves.index(wave)):length*(waves.index(wave)+1)]]
				#print data_row				
				swave=str(int(wave))				
				data[plate][strain][tp][swave+'nm'][labels[row]]=np.array(data_row)
				data[plate][strain][tp][swave+'nm'][labels[row]+'_max']=max(data_row)
				data[plate][strain][tp][swave+'nm'][labels[row]+'_min']=min(data_row)

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
		if itype in ['xls','xlsx','asc'] and os.path.isfile(ifile):
			ilist.append(ifile)
		elif itype in ['txt',''] and os.path.isfile(ifile):
			ifl=open(ifile,'r')
			idata=ifl.read().split('\n')
			idata=[fl.strip() for fl in idata if fl!='']
			for fld in idata:
				ilist.extend(genlist(fld))

		elif iname=='' and itype in ['xls','xlsx','asc']:
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
								sheet.append([genes[plate][l]['Gene']]+data[plate][strain][tp]['Shift'][out][l].tolist())
							else:
								sheet.append([genes[plate][l]['Gene']]+data[plate][strain][tp][out][l].tolist())


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

def round_to(n, precission):
	#Round a number to desired precision
	correction = 0.5 if n >= 0 else -0.5
	return int(n/precission+correction)*precission
	
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
	data=book.sheet_by_name([nm for nm in book.sheet_names() if 'Magellan' in nm][0])
	sheet=[]	
	for r in range(0,data.nrows):
		if len(data.row_values(r))>200:
			#print set(data.row_values(r))
			sheet.append(data.row_values(r))
	return sheet

def readacs(ifile):
	f=open(ifile,'r')
	sheet=[]
	for l in f:
		if len(l.split('\t'))>200:
			row=[cell.strip() for cell in l.split('\t')]
			row=[numerize(cell) for cell in row]
			sheet.append(row)
	f.close()
	return sheet
	


#----------------------------------

if __name__ == "__main__":
	sys.exit(main())

