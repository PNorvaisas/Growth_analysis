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

for mod in ['pip','scipy','string','math','re','csv','sys','os','commands','datetime','operator','getopt','pickle','shutil','glob','types','math','copy','pyExcelerator','xlrd','xlwt','xlutils','types']:
	try:
		exec "import %(mod)s" % vars()
	except ImportError, e:
		print "Module not found %(mod)s\nTrying to install!" % vars()
		install(mod)		
		#pass # module doesn't exist, deal with it.
import unicodedata

#Removes the annoying rocket icon in Mac!
import matplotlib
matplotlib.use("Agg")
#
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
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
#from multiprocessing import Pool



try:
	import itertools as IT
except ImportError, e:
	print "Module itertools not found"


compare = lambda x, y: Counter(x) == Counter(y)


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
				metabolites=readmet(dfile)		
			else:
				print "No database file specified."
				sys.exit(1)
				
			
	
		except Exception, e:
			print e
			#print "!!!-------------------------------------------!!!"
			sys.exit(1)

		print optionsset %vars()
	#----------------------------

	if odir!="":
		dirn=dircheck(odir)


	if load:
		f = open('Biolog_data.pckl','rb')
		data,allfit,metabolites,odir = pickle.load(f)
		f.close()

		
	else:	
		
		ilist,info=genlist(ifile)
		#print ilist
		#print info
		uniques=uniquecomb(info,['Plate','Strain'],'Type')
		#print uniques
		#sys.exit(1)
		#print metabolites
		#sys.exit(0)
		
		#checkfiles(ilist)	
		data=collect(ilist)
		#sys.exit(1)
		data=analyze(data)
		data=growthfit(data)

		
		
		sheets=makesheets(data,metabolites,info)
		writesheets(sheets,odir,sep=',')
		#f = open('Biolog_data.pckl', 'w')
		#pickle.dump([data,allfit,metabolites,odir], f)
		#f.close()
	
	
	#data,allfit=growthfit(data)

	plot_comparison(data,metabolites,odir,['590nm_dt','750nm_dt',
	                                       '590nm_f','590nm_log','750nm_f','750nm_log',
	                                       '590nm750nm_diff','590nm750nm_difflog'],info,uniques)
	#'Respiration_log','Respiration','Resp&Growth',
	#plot_comparison(data,metabolites,odir,['Resp&Growth','Growth_log','Growth_dt','dRespiration_dt','Growth','Respiration'])
		
	#plot_comparison(data,metabolites,odir,['Growth','Respiration','Growth_dt','Respiration_dt','Growth_abolished','Respiration_abolished','dRespiration_dt']) #['Growth','Respiration','Growth_dt','Respiration_dt']

	
	#sheets=makesheets(data,metabolites)	
	#writesheets(sheets,odir)
	
	

#-------------Functions------------



def plot_comparison(data,metabolites,dirn,figs,info,uniques):

	for uk in sorted(uniques.keys()):
		comvalues=uniques[uk]['Compare-by-values']
		comname=uniques[uk]['Compare-by']
		unks=uniques[uk]['Unique-keys']
		unval=uniques[uk]['Unique-values']

		combmap={ uc:unks.index(uc) for uc in unks}

		if 'Plate' in unks and 'Strain' in unks:
			strain=unval[combmap['Strain']]
			bplate=unval[combmap['Plate']]
		else:
			print 'No data provided on used strain and plate ID!'
			strain='Unknown strain'
			bplate='Unknown plate'

		print 'Plot unique combinations of: {}'.format(unks)
		print 'Combination: {}'.format(unval)
		print 'Compare by {}: {}'.format(comname,comvalues)
		#sys.exit(1)

		if len(comvalues)==2 and all([cvt in ['Control', 'Treatment'] for cvt in comvalues]):
			for cv in comvalues:
				if cv=='Control':
					cgroup=uniques[uk][cv]
				elif cv=='Treatment':
					egroup=uniques[uk][cv]
		else:
			print 'Unknown types found as descriptors: {}'.format([cvt for cvt in comvalues if cvt not in ['Control', 'Treatment']])
			print 'Currently supported descriptors: {}'.format(['Control', 'Treatment'])
		#group=uniques[uk][cv]
		#strain=info[plate]['Strain']
		#labels=data[plate]['Labels']


		# ref=data[plate]['Control']
		# exp=data[plate]['Treatment']

		figures_temp=data[cgroup[0]]['Figures']
		#figures_temp=figures_temp+['Growth_abolished','Respiration_abolished','Resp&Growth']

		if figs=='all':
			figures=figures_temp

		elif isinstance(figs, list):
			figs_ch=[f for f in figs if f in figures_temp]
			figures=figs_ch
			if len(figs_ch)!=len(figs):
				print 'Figures {} not found'.format([f for f in figs if f not in figures_temp])
		elif isinstance(figs, str) or isinstance(figs, unicode):
			if figs in figures_temp:
				figures=[figs]
			else:
				print 'Figures {} not found'.format(figs)
				figures=[]

		print figures
		for fg in figures:
			print "Plotting plate {} {} {}...".format(bplate,strain,fg)
			#Need to fix metabolites
			plot=plot_2D(bplate,strain,fg,data,cgroup,egroup,metabolites[bplate])
			plot.savefig('{}/{}.pdf'.format(dirn,bplate+'-'+strain+'_'+fg))
			plot.close()

	return data

def plot_2D(bplate,strain,fg,data,cgroup,egroup,metabolites):
	xmax=24
	plate_size=96
	#print title

	if '_dt' in fg:
		time=data[data.keys()[0]]['Time_dt']
	else:
		time=data[data.keys()[0]]['Time']
	labels=data[data.keys()[0]]['Labels']

	#fig=plt.figure(figsize=(11.69,8.27), dpi=100)

	title='{} {} {}'.format(bplate,strain,fg)

	fig,axes=plt.subplots(nrows=8, ncols=12, sharex=True, sharey=True,figsize=(11.69,8.27), dpi=100)
	fig.suptitle(title)
	plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.05, hspace=0.05)


	rnd=0.1
	#Determine range of figure, probably manually
	# maxc=max([max(datac[l]) for l in labels])
	# maxe=max([max(datae[l]) for l in labels])
	#
	# if maxc>0:
	# 	totalmaxc=round_to(maxc,rnd)
	# else:
	# 	totalmaxc=0
	# if maxe>0:
	# 	totalmaxe=round_to(maxe,rnd)
	# else:
	# 	totalmaxe=0
	# totalminc=round_to(min([min(datac[l]) for l in labels]),rnd)
	# totalmine=round_to(min([min(datae[l]) for l in labels]),rnd)
	# totalmax=max([totalmaxc,totalmaxe])

	#Range
	totalmin=0

	ticks=3
	xlabel='Time, h'
	ylabel=''
	decimals=1

	if '_log' in fg or '_difflog' in fg:
		totalmax=0
		totalmin=-6
		ticks=4
		ylabel='log2 OD'
		decimals=1

	if '_dt' in fg:
		totalmax=0.25
		totalmin=0
		ticks=3
		ylabel='dOD/dt'
		decimals=2

	if '590nm_' in fg and '_log' not in fg and '_dt' not in fg:
		totalmax=2
		ticks=3
		ylabel='OD'
		decimals=1

	if '750nm_' in fg and '_log' not in fg and '_dt' not in fg:
		totalmax=1
		totalmin=0
		ticks=3
		ylabel='OD'
		decimals=1

	if fg=='590nm750nm_diff':
		totalmax=1
		totalmin=0
		ticks=3
		ylabel='OD'
		decimals=1



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
		#Need range here
		plt.ylim([totalmin,totalmax])
		#plt.axvline(x=3, ymin=0, ymax=1,c='green', hold=None)  u 
		label=greek_check(metabolites[l]['Name'],12)
		plt.text(0.05, 0.9, label, fontsize=7,verticalalignment='top',transform=ax.transAxes)

		for grp,mrk in IT.izip([cgroup,egroup],['r-','b-']):
			for gdata in grp:
				y=data[gdata][fg][l]
				plt.plot(x,y,mrk)
				if fg=='750nm_log':
					a,c,t0=data[gdata]['GrowthFit'][l]
					if c>0:
						yfit=x*a+c
						plt.plot(x,yfit,mrk,alpha=0.5)


		# #Get y data
		# yc=datac[l]
		# ye=datae[l]
		#
		# if fg in ['Growth_abolished','Respiration_abolished']:
		# 	ye[ye<thres]=thres
		# 	yc[yc<thres]=thres
		# 	yd=(ye*100/yc)-100
		# 	plt.plot(x,yd,'k-')
		# 	plt.fill_between(x, 0, yd, where=yd>=0, facecolor='blue', interpolate=True)
		# 	plt.fill_between(x, 0, yd, where=yd<=0, facecolor='red', interpolate=True)
		#

		# 	elif fg=='Resp&Growth':
		# 		rc=gfitc[l]
		# 		re=gfite[l]
		# 		#print '{}: {} {} {}'.format(fg,len(x),len(yc),len(ye))
		# 		plt.plot(x,rc,'r-',alpha=0.5)
		# 		plt.plot(x,re,'b-',alpha=0.5)
	blue_line = mlines.Line2D([], [], color='blue', label='Treatment')
	red_line = mlines.Line2D([], [], color='red', label='Control')
	#plt.figlegend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
	plt.figlegend((red_line,blue_line),('Control','Treatment'),'upper right')#

	return plt

def uniquecomb(info,sortby,compareby):
	summ=[]
	uns=NestedDict()
	for i in info.keys():
		ks=[info[i][k] for k in sortby]
		summ.append(ks)
	uniqs=[list(i) for i in set(tuple(i) for i in summ)]
	for unl in uniqs:
		unin='|'.join(unl)
		uns[unin]['Unique-keys']=sortby
		uns[unin]['Unique-values']=unl
		uns[unin]['Compare-by']=compareby
		uns[unin]['Compare-by-values']=[]
		for i in info.keys():
			# print sortby
			# print unl
			# print type(unl)
			check=[info[i][uk]==uv for uk,uv in IT.izip(sortby,unl)]
			# print check
			if all(check):
				if info[i][compareby] in uns[unin].keys():
					uns[unin][info[i][compareby]]=uns[unin][info[i][compareby]]+[i]
				else:
					uns[unin]['Compare-by-values']=uns[unin]['Compare-by-values']+[info[i][compareby]]
					uns[unin][info[i][compareby]]=[i]

	return uns


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


#
# def checkfiles(ilist):
# 	print '\n\n-------------Checking integrity of the file list!-------------'
# 	Pass=False
# 	allUAL=[fl for fl in ilist if 'UAL_' in fl ]
# 	Control=[fl for fl in ilist if '_NoMetf_' in fl ]
# 	Experiment=[fl for fl in ilist if '_Metf_' in fl ]
# 	CPlates=[fl.split('_')[1] for fl in Control]
# 	EPlates=[fl.split('_')[1] for fl in Experiment]
# 	if comparison:
# 		Pass=True
# 	else:
# 		if len(allUAL)==len(ilist):
# 			print 'All files UAL!'
# 			if len(Control)==len(Experiment):
# 				print 'Same number of control and experiment files!'
# 				if compare(CPlates,EPlates):
# 					print 'Corresponding plate indexes found!'
# 					print '-----------------File list check passed!-------------\n\n'
# 					Pass=True
# 	if not Pass:
# 		print 'File list integrity check failed!'
# 		sys.exit(1)



def analyze(data):
	filterf='wiener'
	waves=['590nm','750nm']
	msize=20
	par1=4
	par2=0.1
	method='pre'
	window=20
	thres=np.power(2.0,-5)
	thresd=np.power(2.0,-5)
	for plate in sorted(data.keys()):
		#Needs to be checked
		time=data[plate]['Time']

		dt=time[1]-time[0]
		time_dt=(time+dt/2)[:-1]
		npts=len(time)
		nyf=0.5/dt
		print 'Analyzing data in {}'.format(plate)
		data[plate]['Time_dt']=(time+dt/2)[:-1]
		for well in data[plate]['Labels']:
			if well!='A1':
				ref590=data[plate]['590nm_f']['A1']
				ref750=data[plate]['750nm_f']['A1']
				refl590=data[plate]['590nm_log']['A1']
				refl750=data[plate]['750nm_log']['A1']
			gs=np.mean(data[plate]['590nm'][well][:window])
			ds=np.mean(data[plate]['750nm'][well][:window])

			raw590=data[plate]['590nm'][well]
			raw750=data[plate]['750nm'][well]

			back590=setbar(raw590-gs,0.0)
			back750=setbar(raw750-ds,0.0)


			f590=Wiener(back590,msize)
			log590=np.log2(setbar(f590,thres))
			dt590=np.diff(f590)/(dt/3600)
			dtf590=setbar(Wiener(dt590,msize),0.0)


			f750=Wiener(back750,msize)
			log750=np.log2(setbar(f750,thresd))
			dt750=np.diff(f750)/(dt/3600)
			dtf750=setbar(Wiener(dt750,msize),0.0)

			if well=='A1':
				refrl590=log590-log590
				refrl750=log750-log750
				refr590=f590-f590
				refr750=f750-f750
			else:
				refrl590=log590-refl590
				refrl750=log750-refl750
				refr590=f590-ref590
				refr750=f750-ref750

			diff590750=f590-f750
			diff590750log=log590-log750



			#Need threshold
			# resp=setbar(f750,thresd)/setbar(f590,thres)
			# rs=np.mean(resp[:window])
			# resp=Wiener(resp-rs,msize)
			# resp_log=np.log2(setbar(resp,thresd/4))
			# resp=setbar(resp,0.0)
			# resp_dt=np.diff(resp)/(dt/3600)
			# resp_dt=Wiener(resp_dt,msize)
			# resp_dt=setbar(Wiener(resp_dt,msize),0.0)
			#
			# dresp=(np.diff(dye)/(dt/3600))/interp(time,growth,time_dt)
			# dresp=Wiener(dresp,msize)
			# dresp=setbar(dresp,0.0)



			#growth_dt_norm=Wiener(growth_dt,msize)/interp(time,growth,time_dt)

			# if max(dresp)>0.1:
			# 	dresp_norm=dresp/max(dresp)
			# else:
			# 	dresp_norm=np.zeros((len(dresp),), dtype=np.int)
			# if max(growth_dt)>0.01:
			# 	growth_dt_norm=growth_dt/max(growth_dt)
			# else:
			# 	growth_dt_norm=np.zeros((len(growth_dt),), dtype=np.int)



			data[plate]['750nm_b'][well]=back750
			data[plate]['750nm_f'][well]=f750
			data[plate]['750nm_log'][well]=log750
			data[plate]['750nm_dt'][well]=dtf750

			data[plate]['590nm_b'][well]=back590
			data[plate]['590nm_f'][well]=f590
			data[plate]['590nm_log'][well]=log590
			data[plate]['590nm_dt'][well]=dtf590

			data[plate]['590nm_reflog'][well]=refrl590
			data[plate]['590nm_ref'][well]=refr590

			data[plate]['750nm_reflog'][well]=refrl750
			data[plate]['750nm_ref'][well]=refr750

			data[plate]['590nm750nm_diff'][well]=diff590750
			data[plate]['590nm750nm_difflog'][well]=diff590750log

			#data[plate]['Growth_dt_norm'][well]=growth_dt_norm
			# data[plate]['Respiration'][well]=resp
			# data[plate]['Respiration_log'][well]=resp_log
			# data[plate]['Respiration_dt'][well]=resp_dt
			# data[plate]['dRespiration_dt'][well]=dresp
			# data[plate]['dRespiration_dt_norm'][well]=dresp_norm


			#if plate=='PM5':
			#	rs=np.mean(resp[:window])
			#	resp=resp-rs
			#else:
			#	resp=resp-refresp
			#	rs=np.mean(resp[:window])
			#	resp=resp-rs
			#gresp=(growth+gs+0.1)/(resp+rs+0.1)
			#respg=(resp+rs+0.1)/(growth+gs+0.1)
			#gresp=gresp-np.mean(gresp[:window])
			#respg=respg-np.mean(respg[:window])
			#data[plate][strain][tp]['Respiration'][well]=resp
			#data[plate][strain][tp]['Respiration_dt'][well]=np.diff(resp)/(dt/3600)
			#data[plate][strain][tp]['Growth-Resp'][well]=gresp
			#data[plate][strain][tp]['Resp-Growth'][well]=respg

		data[plate]['Figures']=data[plate]['Figures']+['590nm_b','590nm_f','590nm_log','590nm_dt',
		                                               '750nm_b','750nm_f','750nm_log','750nm_dt',
		                                               '750nm_ref','750nm_reflog','590nm_ref','590nm_reflog',
		                                               '590nm750nm_diff','590nm750nm_difflog']

	return data

def setbar(x,bar):
	x2=[xi if xi>bar else bar for xi in x]
	x2=np.array(x2)
	return x2

def interp(x,y,x2):
	tck = ip.splrep(x, y, s=0)
	y2=ip.splev(x2, tck, der=0)
	return y2

def readmet(ifile):
	print ifile
	nutrients=NestedDict()

	rdr=csv.reader(open(ifile,'r'), delimiter=',')
	data=[ln for ln in rdr]
	headers=data[0]
	metin=headers.index('Metabolite')
	ecoin=headers.index('EcoCycID')
	plate=headers.index('Plate')
	well=headers.index('Well')
	indx=headers.index('Index')
	print metin, ecoin,plate,well
	for ln in data[1:]:
		met=ln[metin].encode('ascii','ignore').strip()
		eco=ln[ecoin].encode('ascii','ignore').strip()
		pl=ln[plate].encode('ascii','ignore').strip()
		wl=ln[well].encode('ascii','ignore').strip()
		ind=ln[indx].encode('ascii','ignore').strip()
		nutrients[pl][wl]['ID']=eco
		nutrients[pl][wl]['Name']=met
		nutrients[pl][wl]['Index']=ind


	#print genes
	return nutrients

def cut(x, y, a,b):
	x2=[]
	y2=[]
	last=0
	for xt, yt in IT.izip(enumerate(x),enumerate(y)):
		df=yt[0]-last
		#print df
		if yt[1]>a and yt[1]<b:# and df<3
			last=yt[0]
			x2.append(xt[1])
			y2.append(yt[1])
	y2=np.asarray(y2)
	x2=np.asarray(x2)
	
	return x2,y2 

def growth(x,a,c):
	y=x*a+c
	return y

def growthfit(data):
	y0=np.power(2.0,-4)
	for plate in data.keys():
		x=data[plate]['Time']
		x=x/3600
		labels=data[plate]['Labels']
		#ref=data[plate]['590nm_f']['A1']
		for l in labels:
			y=data[plate]['750nm_log'][l]
			#y=setbar(y,np.power(2,-5))
			#maxy=max(y)
			#miny=min(y)
			#yl=np.log2(y)
			miny=min(y)
			maxy=max(y)
			gscale=maxy-miny
			thres=-5 #-4
			#print miny, maxy
			x2,y2=cut(x, y, thres+(maxy-thres)*0.1, thres+(maxy-thres)*0.6) #0.6
			if len(y2)>0 and gscale>0.5:
				try:
					popt, pcov = curve_fit(growth, x2, y2)
					a=popt[0]
					c=popt[1]
					t0=(np.log2(y0)-c)/(a)
				except TypeError:
					print 'Curve_fit encountered an error!'
					a=0
					c=0
					t0=float("inf")

			else:
				a=0
				c=0
				t0=float("inf")
			# for par,nm in IT.izip([a,c,t0],['a','c','t0']):
			# 	if allfit[tp]['Log'][nm]:
			# 		allfit[tp]['Log'][nm]=allfit[tp]['Log'][nm]+[par]
			# 	else:
			# 		allfit[tp]['Log'][nm]=[par]
			data[plate]['GrowthFit'][l]=[a,c,t0]

	return data

def collect(ilist):
	data=NestedDict()
	for ifl in sorted(ilist):		
		ipt, inm, itp = filename(ifl)
		if itp=='xlsx':
			sheet=readxls(ifl)
		elif itp=='asc':
			sheet=readacs(ifl)
		else:
			print 'Unknown file format'


		nrows=len(sheet)
		nm_labels=sheet[0]
		#converting first row to set rearanges the labels
		#list(set(sheet[0]))
		#print nm_labels
		#print nrows
		
		waves=[numerize(wlen) for wlen in nm_labels if wlen not in ['Layout','Well positions','','Replicate Info']]
		#print waves
		lengths=[sheet[0].index(wave) for wave in waves]
		#Selection of time cells does not depend on the order of wavelengths
		length=max(lengths)
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
		data[ifl]['Labels']=labels
		data[ifl]['Spectra']=waves
		data[ifl]['Time']=time
		data[ifl]['Temp']=temp
		data[ifl]['Time_max']=timemax_h
		data[ifl]['Time_step']=timestep
		data[ifl]['Wells']=nrows-3
		data[ifl]['Figures']=[str(w)+'nm' for w in waves]
		data[ifl]['File']=inm
		print "File: {}".format(ifl)
		print "Wavelengths: {}, {}".format(*waves)
		print "Run time {}h, step {}min in {} wells\n".format(timemax_h,timestep/60, len(labels))
		for row in range(0,len(labels)):
			for wave in waves:
				data_row=[60000 if val=="Overflow" else val for val in sheet[row+3][length*(waves.index(wave)):length*(waves.index(wave)+1)]]
				swave=str(int(wave))				
				data[ifl][swave+'nm'][labels[row]]=np.array(data_row)
				data[ifl][swave+'nm'][labels[row]+'_max']=max(data_row)
				data[ifl][swave+'nm'][labels[row]+'_min']=min(data_row)



	return data



def genlist(ifile):
	#Generates list of input files, checks their existance
	ilist=[]
	info=NestedDict()
	if ',' in ifile:
		ifiles=ifile.split(',')
		for ifl in ifiles:
			ilist.extend(genlist(ifl)[0])
	else:
		ipath, iname, itype=filename(ifile)
		if 'Design' in ifile and itype in ['xls','xlsx','csv'] and os.path.isfile(ifile):
			ilist,info=readinfo(ifile)
		elif itype in ['xls','xlsx','asc'] and os.path.isfile(ifile):
			ilist.append(ifile)
			#for flt in ilistt:
			#	ilist.extend(genlist(flt))
		elif itype in ['txt',''] and os.path.isfile(ifile):
			ifl=open(ifile,'r')
			idata=ifl.read().split('\n')
			idata=[fl.strip() for fl in idata if fl!='']
			for fld in idata:
				ilist.extend(genlist(fld))
		elif iname=='' and itype in ['xls','xlsx','asc','txt']:
			ffiles=glob.glob('*.%(itype)s' % vars())
			#print ffiles
			ilist.extend(ffiles)
		else:
			print ifile
			print "Bad file type %(inp)s!" % vars()
			sys.exit(1)

	if len(ilist)>0 and len(info.keys())==0:
		for ifl in ilist:
			ipt, inm, itp = filename(ifl)
			info[ifl]['Plate']=inm.split('_')[2]
			info[ifl]['Strain']=inm.split('_')[3]
			info[ifl]['Type']='Control' if 'Control' in inm or 'NoMetf' in inm else 'Treatment'
	return ilist,info

def readinfo(ifile):
	print ifile
	info=NestedDict()
	ilist=[]
	ipt, inm, itp = filename(ifile)
	if itp in ['xlsx','xls']:
		data=readxls_s(ifile)
	elif itp=='csv':
		data=readcsv(ifile)
	headers=data[0]
	#Automatically find variables
	headin={ hd : headers.index(hd) for hd in headers}
	nec=['File','Plate','Type']
	if all(n in headers for n in nec):
		print 'Necessary headers found!'
	else:
		print 'Missing essential headers in description file!'
		print headers
		sys.exit(0)
	#filein=headers.index('File')
	# platein=headers.index('Plate')
	# strainin=headers.index('Strain')
	# typein=headers.index('Type')
	#print metin, ecoin,plate,well
	for ln in data[1:]:
		fl=str(ln[headin['File']]).encode('ascii','ignore').strip()
		for hd in headin.keys():
			info[fl][hd]=str(ln[headin[hd]]).encode('ascii','ignore').strip()
		ilist.append(fl)

	#print genes
	return ilist, info



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

def makesheets(data,metabolites,info):
	sheets=NestedDict()
	dt=[]
	gfit=[]
	regular=[]
	
	header=info[info.keys()[0]].keys()+['Well','Index','Data','Name','EcoCycID']
	ghead=['a','c','t0']
	gfit.append(header+ghead)
	for fln in data.keys():
		# strain=info[plate]['Strain']
		# tp=info[plate]['Type']

		output=data[fln]['Figures']
		output=output+['GrowthFit']
		time_dt=data[fln]['Time_dt']
		time_lin=data[fln]['Time']
		labels=data[fln]['Labels']
		annot=[info[fln][k] for k in info[fln].keys()]
		plate=info[fln]['Plate']
		for fig in output:
			#print '{}...{}'.format(fln,fig)
			for well in labels:
				indx=plate+'-'+well
				datrow=data[fln][fig][well]
				rowhead=annot+[well,indx,fig,metabolites[plate][well]['Name'],metabolites[plate][well]['ID']]

				if not isinstance(datrow,list):
					datrow=datrow.tolist()
				newrow=rowhead+datrow
				if '_dt' in fig:
					dt.append(newrow)

				elif fig=='GrowthFit':
					gfit.append(newrow)
				else:
					regular.append(newrow)

	time_lin=time_lin.tolist()
	time_dt=time_dt.tolist()
	dt.insert(0,header+time_dt)
	regular.insert(0,header+time_lin)
	sheets['Data_dt']=dt
	sheets['Data']=regular
	sheets['Growth_fit']=gfit
	
	#print len(dt[0]),len(dt[1])
	#print dt[0],dt[1]

	return sheets



def writesheets(sheets,odir,sep='\t'):
	#Writes organized data to file.
	#odir=dircheck('Split')
	print sheets.keys()
	for fig in sheets.keys():
		print fig
		oname='{}/{}.csv'.format(odir,fig)		
		sheet=sheets[fig]
		print len(sheet)
		f=open(oname,"wb")
		ofile=csv.writer(f, delimiter=sep) # dialect='excel',delimiter=';', quotechar='"', quoting=csv.QUOTE_ALL
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
	#print book.sheet_names()
	data=book.sheet_by_name([nm for nm in book.sheet_names() if 'Magellan' in nm or 'RAW data' in nm][0])
	sheet=[]	
	for r in range(0,data.nrows):
		if len(data.row_values(r))>200:
			#print set(data.row_values(r))
			sheet.append(data.row_values(r))
	return sheet

def readxls_s(ifile):
	book=xlrd.open_workbook(ifile,formatting_info=False)
	data=book.sheet_by_name(book.sheet_names()[0])
	sheet=[]
	for r in range(0,data.nrows):
		#print set(data.row_values(r))
		sheet.append(data.row_values(r))
	return sheet

def readcsv(ifile):
	f=open(ifile,'r')
	sheet=[]

	rdr=csv.reader(f, delimiter=',')
	data=[ln for ln in rdr]
	f.close()
	return data


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

