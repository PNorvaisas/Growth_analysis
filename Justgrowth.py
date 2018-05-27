#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Justgrowth.py

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

for mod in ['pip','string','math','re','csv','sys','os',
            'commands','datetime','operator','getopt','subprocess','pickle','shutil','glob',
            'types','math','copy','pyExcelerator','xlrd','xlwt','xlutils','types','warnings']:
	try:
		exec "import %(mod)s" % vars()
	except ImportError, e:
		print "Module not found %(mod)s\nTrying to install!" % vars()
		install(mod)		
		#pass # module doesn't exist, deal with it.


#import warnings
#Disable simple warnings
warnings.filterwarnings("ignore")

import unicodedata
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rc


#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
#from pylab import *
from scipy import interpolate
from scipy import interpolate as ip
import scipy.signal as sig
from scipy.fftpack import rfft, irfft, fftfreq

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
Bacterial growth data preparation

Flags:
	-h  Display this help message
	-v  Verbose output
Arguments:
	-i <files>   Input file
	-p <file>    Plate arrangement
	-o <dir>     Directory to write output to
	-d <list>	 List of descriptors for each file
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
      Files:    %(ifile)s
    Pattern:	%(pfile)s
Descriptors:    %(subst)s
        Out:	%(odir)s
       Load:	%(load)s
<--------------------------------------------->
	'''



def main(argv=None):
	ifile=""
	pfile=""
	msize=20
	odir='Output'
	odirn=''
	load=False
	subst=''
	comparison=False
	
	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "hi:f:d:o:p:", ["help"])
		except getopt.error, msg:
			raise Usage(msg)


		for option, value in opts:
			if option in ("-h", "--help"):
				usage()
				return	
			if option in ("-i", "--input"):
				ifile=value
			if option in ("-p", "--pattern"):
				pfile=value
			if option in ("-d", "--descriptors"):
				subst=value
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


	print optionsset %vars()
	#----------------------------



	if 'Design' in ifile:
		print 'Reading experimental design information!'
		#
		info,ilist,dlist,odirn=readinfo(ifile)
		subs=[[]]*len(ilist)
		descriptors=readdesc(ilist,dlist,subs)
	else:
		print 'Reading data files!'
		print ifile
		ilist=genlist(ifile)
		#print ilist
		if pfile!='':
			dlist=genlist(pfile)
			# print dlist
			if subst!='':
				subs=[st.split('|') for st in subst.split(',')]
			else:
				subs=[[]]*len(ilist)
			#print subs
			descriptors=readdesc(ilist,dlist,subs)
		else:
			descriptors={}


	if odirn!='' and odir=='Output':
		odir=odirn

	if odir!="":
		odirn=dircheck(odir)



	#checkfiles(ilist)
	data=collect(ilist)
	#sys.exit(1)
	data=analyze(data)
	data=growthfit(data)

	sheets=makesheets(data,descriptors,info)
	writesheets(sheets,odir)
	# f = open('Biolog_data.pckl', 'w')
	# pickle.dump([data,allfit,descriptors,odir], f)
	# f.close()
	
	
	#data,allfit=growthfit(data)

	plot_comparison(data,odir,'all')
	#['600nm','590nm','595nm','600nm_log','590nm_log','595nm_log','750nm','750nmC']
	#plot_comparison(data,metabolites,odir,['Resp&Growth','Growth_log','Growth_dt','dRespiration_dt','Growth','Respiration'])
		
	#plot_comparison(data,metabolites,odir,['Growth','Respiration','Growth_dt','Respiration_dt','Growth_abolished','Respiration_abolished','dRespiration_dt']) #['Growth','Respiration','Growth_dt','Respiration_dt']

	
	#sheets=makesheets(data,metabolites)	
	#writesheets(sheets,odir)
	
	

#-------------Functions------------

def readinfo(ifile):
    print ifile
    info=NestedDict()
    ilist=[]
    dlist=[]
    odirs=[]
    ipt, inm, itp = filename(ifile)
    if itp in ['xlsx','xls']:
        data=readxls_s(ifile)
    elif itp=='csv':
        data=readcsv(ifile)
    headers=data[0]
    #Automatically find variables
    headin={ hd : headers.index(hd) for hd in headers}
    nec=['File','Pattern']

    addhead=[key for key in headin.keys() if key not in nec]

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
        dsc=str(ln[headin['Pattern']]).encode('ascii','ignore').strip()
        ipath, iname, itype=filename(fl)
        
        ilist.append(fl)
        dlist.append(dsc)
        if 'Odir' in headin.keys():
            odr=str(ln[headin['Odir']]).encode('ascii','ignore').strip()
            odirs.append(odr)
        for hd in addhead:
            info[iname][hd]=str(numerize(ln[headin[hd]])).strip().encode('ascii','ignore')

    # print odirs
    # sys.exit(1)
    if len(list(set(odirs)))!=0:
        odir=list(set(odirs))[0]
    else:
        odir=''

    #print genes
    print odir
    return info, ilist, dlist, odir



def plot_comparison(data,dirn,figs):

	for plate in sorted(data.keys()):
		labels=data[plate]['Labels']
		nwells=data[plate]['Wells']
		plsize=data[plate]['Used wells']
		buffered=data[plate]['Buffered']
		ref=data[plate]
		figures_temp=data[plate]['Figures']

		if figs=='all':
			figures=figures_temp
			
		elif isinstance(figs, list):
			figs_ch=[f for f in figs if f in figures_temp]
			if len(figs_ch)>0:
				figures=figs_ch
			else:
				print 'Figures {} not found'.format(figs)
		elif isinstance(figs, str) or isinstance(figs, unicode):
			if figs in figures_temp:
				figures=[figs]
			else:
				print 'Figures {} not found'.format(figs)

		#print figures
		print 'File: {}'.format(plate)
		for fg in figures:
				
			if '_log' in fg:
				fgl=fg
				gfitc=ref['GrowthFit']
			
			else:
				fgl=fg
				gfitc=''

			print "\tPlotting {}...".format(fg)
			if '_dt' in fg or fg=='Resp&Growth':
				time=data[plate]['Time_dt']
			else:
				time=data[plate]['Time']

			#Need to fix metabolites
			plot=plot_2D('{} {}'.format(plate,fg),ref[fgl],time,labels,gfitc,plsize)
			plot.savefig('{}/{}.pdf'.format(dirn,plate+'_'+fg))
			plot.close()
#
	return data

def plot_2D(title,datac,time,labels,gfitc,plsize):
	if plsize==96:
		minrow=1
		maxrow=8
		mincol=1
		maxcol=12
		cols=12
		rows=8
	if plsize==60:
		minrow=1
		maxrow=6
		mincol=1
		maxcol=10
		cols=10
		rows=6
	elif plsize==12:
		minrow=1
		maxrow=3
		mincol=1
		maxcol=4
		cols=4
		rows=3
	elif plsize==384:
		minrow=1
		maxrow=16
		mincol=1
		maxcol=24
		cols=24
		rows=16
	elif plsize==240:
		minrow=1
		maxrow=14
		maxcol=22
		cols=20
		rows=12
	xmax=24
	#print title
	plate,fg=title.split()
	#fig=plt.figure(figsize=(11.69,8.27), dpi=100)
	
	fig,axes=plt.subplots(nrows=rows, ncols=cols, sharex=True, sharey=True,figsize=(11.69,8.27), dpi=100)
	fig.suptitle(title)
	plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.05, hspace=0.05)
	
	rnd=0.1
	maxc=max([max(datac[l]) for l in labels])

	if maxc>0:		
		totalmaxc=round_to(maxc,rnd)
	else:
		totalmaxc=0

	totalminc=round_to(min([min(datac[l]) for l in labels]),rnd)

	totalmax=totalmaxc
	totalmin=0
	ticks=3
	xlabel='Time, h'
	ylabel=''
	decimals=1
	if 'nm' in fg:
		totalmax=1
		ylabel='OD@'+fg
		decimals=2

	
	# if fg=='Growth':
	# 	totalmax=totalmaxc+0.1
	# 	ylabel='Growth'
	# 	decimals=1

	if '_log' in fg:
		totalmax=0
		totalmin=-6
		ylabel='log2 Growth'
		decimals=1

	if '_dt' in fg:
		totalmax=0.2
		totalmin=0
		ylabel='Growth/dt'
		decimals=2


	if plsize>12:
		ticks=3
	else:
		ticks=5



	ymin=totalmin
	fig.text(0.5, 0.04, xlabel, ha='center')
	fig.text(0.04, 0.5, ylabel, va='center', rotation='vertical')

	for v,l in IT.izip(range(plsize),labels):
		if plsize==240:
			#Numbering inside plot arrangement differs from real numbers by buffer size
			row=string.uppercase.index(l[0])+1-2
			col=int(l.replace(l[0],''))-2
		elif plsize==60:
			row=string.uppercase.index(l[0])+1-1
			col=int(l.replace(l[0],''))-1
		else:
			row=string.uppercase.index(l[0])+1
			col=int(l.replace(l[0],''))
		#print (row,col)

		if divmod(float(v)/12,1)[1]==0:
			sh_y=l
		#print sh_y
		v=v+1
	
		x=time/3600
		yc=datac[l]

		if fg=='Growth_log':
			ca,cc,ct=gfitc[l]

		#elif fg=='Resp&Growth':
		#	rc=gfitc[l]

		plt.sca(axes[row-1,col-1])
		ax=axes[row-1,col-1]
		#print col,cols
		#print row,rows
		if col==cols:
			ax.yaxis.set_label_position("right")
			plt.ylabel(l[0],rotation='horizontal')
		if row==minrow:
			plt.title(str(int(l.replace(l[0],''))))

		if col>1 and row<rows-1:
			plt.setp(ax.get_yticklabels(), visible=False)

		plt.xticks(np.linspace(0, xmax, 3),['']+list(np.linspace(0, xmax, 3).astype(int)[1:]), rotation='vertical')	
		plt.yticks(np.linspace(ymin, totalmax, ticks),['']+list(myround(np.linspace(totalmin, totalmax, ticks),decimals)[1:]))
		plt.ylim([totalmin,totalmax])
		#plt.axvline(x=3, ymin=0, ymax=1,c='green', hold=None)  u 

		#label=greek_check(metabolites[l]['Name'],12)
		plt.text(0.05, 0.9, l, fontsize=7,verticalalignment='top',transform=ax.transAxes)

		#print '{}: {} {} {}'.format(fg,len(x),len(yc),len(ye))
		plt.plot(x,yc,'r-')
		if fg=='Growth_log':
			if ca>0:
				yfitc=x*ca+cc
				plt.plot(x,yfitc,'r-',alpha=0.5)

	return plt

def myround(a, decimals=1):
     return np.around(a-10**(-(decimals+5)), decimals=decimals)



# def growth(x,a,c):
# 	y=x*a+c
# 	return y

def Wiener(y, n):
	wi = sig.wiener(y, mysize=n)
	return wi

def Butter(x, y, par1, par2):	
	b, a = sig.butter(par1, par2)
	fl = sig.filtfilt(b, a, y)
	return fl


def growth(x,A,lam,u):
	return A/(1+np.exp((4*u/A)*(lam-x)+2))

def log_growth(x,A,lam,u):
	y=np.log2(A)-np.log2(1+np.exp((4*u/A)*(lam-x)+2))
	#print x,y
	return np.log2(A)-np.log2(1+np.exp((4*u/A)*(lam-x)+2))

def lin(x,a,c):
	y=x*a+c
	return y


def setbar(x,bar):
	x2=[xi if xi>bar else bar for xi in x]
	x2=np.array(x2)
	return x2

def interp(x,y,x2):
	tck = ip.splrep(x, y, s=0)
	y2=ip.splev(x2, tck, der=0)
	return y2


def analyze(data):
	msize=20
	window=20
	thres=np.power(2.0,-5)
	for plate in sorted(data.keys()):
		time=data[plate]['Time']
		time_h=time/3600
		dt=time[1]-time[0]
		time_dt=(time+dt/2)[:-1]
		data[plate]['Time_dt']=time_dt
		#Needs to be checked
		waves=data[plate]['Spectra']
		#print waves
		#print data[plate]['Labels']
		for wave in waves:
			wv=wave.replace('nm','')
			print 'Analyzing data in {}: {}'.format(plate,wave)
			#print time
			for well in data[plate]['Labels']:
				#print data[plate][wave][well]
				start=np.mean(data[plate][wave][well][:window])

				raw=data[plate][wave][well]

				back=setbar(raw-start,0.0)


				filtered=Wiener(back,msize)
				logfilt=np.log2(setbar(filtered,thres))
				dtfilt=np.diff(filtered)/(dt/3600)

				margin=0.01


				#Select which channel to use for growth fitting
				#In this setting only one channel is used
				if wave in ['590nm','595nm','600nm','750nm']:
					grow=filtered
					maxg=max(grow)
					scaleg=maxg-min(grow)
					timec,growc=cut(time_h, grow, 0, maxg,equalize=True)
					if scaleg>0.1 and len(timec)>10 and len(growc)>10:
						#print plate,well
						#print len(timec), len(growc)
						#print min(growc),max(growc)
						#popt, pcov = curve_fit(growth, timec, growc,bounds=(0,np.inf),p0=[0.5, 5, 0.1],max_nfev=5000)
						#A,lam,u=popt
						try:
							popt, pcov = curve_fit(growth, timec, growc,bounds=(0,np.inf),p0=[0.5, 5, 0.1],max_nfev=5000)#,p0=[0.1,10,1]maxfev=5000
							A,lam,u=popt
							#print popt,tmaxf
						except (RuntimeError, ValueError, RuntimeWarning, UnboundLocalError) as e:
							print e
							print 'Curve_fit encountered an error in well {}!'.format(well)
							A,lam,u=[0,np.inf,0]

						if A>0 and lam<np.inf and u>0:
							yreducedf = growth(time_h,*popt) - max(growth(time_h,*popt))*(1-margin) #maxg#
							freducedf = interpolate.UnivariateSpline(time_h, yreducedf, s=0)
							if len(freducedf.roots())>0:
								tmaxf=freducedf.roots()[0]
							else:
								tmaxf=np.inf
						else:
							tmaxf=np.inf

						yreduced = growc - maxg*(1-margin)
						try:
							freduced = interpolate.UnivariateSpline(timec, yreduced, s=0)
							if len(freduced.roots()>0):
								tmax=freduced.roots()[0]
							else:
								tmax=np.inf
						except TypeError:
							print 'Integration encountered an error in well {}!'.format(well)
							tmax=np.inf
					else:
						A,lam,u=[0,np.inf,0]
						tmaxf=np.inf
						tmax=np.inf
					data[plate]['Summary']['GrowthFit'][well]=[A,lam,u,tmax,tmaxf]

				data[plate][wave+'_b'][well]=back
				data[plate][wave+'_f'][well]=filtered
				data[plate][wave+'_log'][well]=logfilt
				data[plate][wave+'_dt'][well]=dtfilt

			data[plate]['Figures']=data[plate]['Figures']+[wave+'_b',wave+'_f',wave+'_log',wave+'_dt']

		for fg in data[plate]['Figures']:
			time=data[plate]['Time']
			maxt=round_to(max(time)/3600,1)
			#print maxt
			if maxt>=16:
				ind16=time.tolist().index(16*3600)
			else:
				ind16=time.tolist().index(maxt*3600)
			if maxt>=24:
				ind24=time.tolist().index(24*3600)
			else:
				ind24=time.tolist().index(maxt*3600)
			for well in data[plate]['Labels']:
				data[plate]['Summary']['Max_{}'.format(fg)][well]=max(data[plate][fg][well])
				if '_dt' in fg:
					data[plate]['Summary']['24h_{}'.format(fg)][well]=data[plate][fg][well][ind24-1]
					data[plate]['Summary']['16h_{}'.format(fg)][well]=data[plate][fg][well][ind16-1]
				else:
					data[plate]['Summary']['24h_{}'.format(fg)][well]=data[plate][fg][well][ind24]
					data[plate]['Summary']['16h_{}'.format(fg)][well]=data[plate][fg][well][ind16]

					data[plate]['Summary']['Int_{}'.format(fg)][well]=interpolate.UnivariateSpline(time_h, data[plate][fg][well], k=5, s=5).integral(0, 24)
					data[plate]['Summary']['Int-tmax_{}'.format(fg)][well]=interpolate.UnivariateSpline(time_h, data[plate][fg][well], k=5, s=5).integral(0, data[plate]['Summary']['GrowthFit'][well][3] if data[plate]['Summary']['GrowthFit'][well][3]<24 else 24) if data[plate]['Summary']['GrowthFit'][well][3]!=np.inf else np.inf
					data[plate]['Summary']['Int-tmaxf_{}'.format(fg)][well]=interpolate.UnivariateSpline(time_h, data[plate][fg][well], k=5, s=5).integral(0, data[plate]['Summary']['GrowthFit'][well][4] if data[plate]['Summary']['GrowthFit'][well][4]<24 else 24) if data[plate]['Summary']['GrowthFit'][well][4]!=np.inf else np.inf
					data[plate]['Summary']['Int_{}_log'.format(fg)][well]=np.log2(data[plate]['Summary']['Int_{}'.format(fg)][well]) if data[plate]['Summary']['Int_{}'.format(fg)][well]>0 else None

	return data




def setbar(x,bar):
	x2=[xi if xi>bar else bar for xi in x]
	x2=np.array(x2)
	return x2

def interp(x,y,x2):
	tck = ip.splrep(x, y, s=0)
	y2=ip.splev(x2, tck, der=0)
	return y2

# def cut(x, y, a,b):
# 	x2=[]
# 	y2=[]
# 	last=0
# 	for xt, yt in IT.izip(enumerate(x),enumerate(y)):
# 		df=yt[0]-last
# 		#print df
# 		if yt[1]>a and yt[1]<b:# and df<3
# 			last=yt[0]
# 			x2.append(xt[1])
# 			y2.append(yt[1])
# 	y2=np.asarray(y2)
# 	x2=np.asarray(x2)
#
# 	return x2,y2

def cut(x, y, a,b,equalize=False):
	x2=[]
	y2=[]
	maxy=max(y)
	maxind=y.tolist().index(maxy)
	last=0
	for xt, yt in IT.izip(enumerate(x),enumerate(y)):
		df=yt[0]-last
		#print df
		if yt[1]>a and yt[1]<=b:# and df<3
			x2.append(xt[1])
			last=yt[0]
			if yt[0]<maxind:
				y2.append(yt[1])
			else:
				if equalize:
					y2.append(maxy)
				else:
					y2.append(yt[1])
	y2=np.asarray(y2)
	x2=np.asarray(x2)

	return x2,y2


def growthfit(data):
	y0=np.power(2.0,-4)
	for plate in data.keys():
		x=data[plate]['Time']
		x=x/3600
		labels=data[plate]['Labels']
		waves=data[plate]['Spectra']
		if len(waves)>1 and '750nm' in waves:
			wave='750nm'
		elif len(waves)>1 and '590nm' in waves:
			wave='590nm'
		elif len(waves)>1 and '595nm' in waves:
			wave='595nm'
		elif len(waves)>1 and '600nm' in waves:
			wave='600nm'
		else:
			wave=waves[0]
		for l in labels:
			y=data[plate][wave+'_log'][l]
			#y=setbar(y,np.power(2,-5))
			#maxy=max(y)
			#miny=min(y)
			#yl=np.log2(y)
			#print l,y
			miny=min(y)
			maxy=max(y)
			gscale=maxy-miny
			thres=-5 #-4
			#print miny, maxy
			x2,y2=cut(x, y, thres+(maxy-thres)*0.1, thres+(maxy-thres)*0.6) #0.6
			if len(y2)>0 and gscale>0.5:
				try:
					popt, pcov = curve_fit(lin, x2, y2)
					a=popt[0]
					c=popt[1]
					t0=(np.log2(y0)-c)/(a)
					tmax=(maxy-c)/(a)
				except TypeError:
					print 'Curve_fit encountered an error!'
					a=0
					c=0
					t0=float("inf")
					tmax=float("inf")

			else:
				a=0
				c=0
				t0=float("inf")
				tmax=float("inf")
			# for par,nm in IT.izip([a,c,t0],['a','c','t0']):
			# 	if allfit[tp]['Log'][nm]:
			# 		allfit[tp]['Log'][nm]=allfit[tp]['Log'][nm]+[par]
			# 	else:
			# 		allfit[tp]['Log'][nm]=[par]
			data[plate]['Summary']['GrowthFit_log'][l]=[a,c,t0,tmax]

	return data

def collect_Tecan(sheet):
	sheetdata=NestedDict()
	datarange=sheet[0].index('Well positions')

	sheet=[r for r in sheet if len(r)>1]

	nrows=len(sheet)
	datarange=sheet[0].index('Well positions')
	nm_labels=[lab for lab in sheet[0] if lab not in ['Layout','Well positions','','Replicate Info']]
	#print nm_labels
	if len(nm_labels)>1:
		starts=[sheet[0].index(lab) for lab in nm_labels]
	else:
		starts=[sheet[0].index(nm_labels[0])]
		if nm_labels[0]=='Raw data':
			nm_labels[0]='600'

	waves=[numerize(wlen.replace('nm','')) for wlen in nm_labels]
	#print 'Identified wavelengths: {}'.format(waves)
	#print datarange

	length=(datarange)/len(waves)

	#Selection of time cells does not depend om the order of wavelengths

	#Extract time

	time_row=sheet[1][:length]
	#print 'Length:{}'.format(length)
	#print time_row
	#print len(time_row)
	if list(set([s for s in time_row if isinstance(s,float)])):
		time_t=time_row
	else:
		time_t=[int(str(t).replace('s','')) for t in time_row]

	#length=time_t.index(max(time_t))+1
	temp=[float(t.split()[0]) for t in sheet[2][:length]]

	timemax_min=int(round_to(float(time_t[-1])/60,5))
	timemax_h,timemax_remmin=divmod(timemax_min,60)

	#Time in seconds
	time=np.linspace(0,timemax_min*60,length)

	timestep=round_to(float(time_t[-1])/(length-1),1)
	#timestep=time[-1]-time[-2]

	#print time_t[-1],timemax_h,time[-1],timestep


	alllabels=[r[datarange] for r in sheet if r[datarange] not in ['Well positions']][2:]
	labels=[l for l in alllabels if l!='']#Sheet
	#print labels


	plsize=len(labels)
	if plsize not in [12,48,96,384]:
		buffer=True
	else:
		buffer=False

	sheetdata['Labels']=labels
	sheetdata['Spectra']=[str(int(w))+'nm' for w in waves]
	sheetdata['Time']=time
	sheetdata['Temp']=temp
	sheetdata['Time_max']=timemax_min
	sheetdata['Time_step']=timestep
	sheetdata['Wells']=len(labels)
	sheetdata['Used wells']=plsize
	sheetdata['Buffered']=str(buffer)
	sheetdata['Figures']=[str(w)+'nm' for w in waves]
	#sheetdata['File']=inm
	print "Wavelengths: {}".format(waves)
	print "Run time {}, step {}min in {} wells\n".format(str(datetime.timedelta(minutes=timemax_min)),timestep/60, len(labels))
	for lab in range(0,len(alllabels)):
		if alllabels[lab]!='':
			for wave in waves:
				scol=(length)*(waves.index(wave))
				ecol=(length)*(waves.index(wave)+1)
				data_row=[val for val in sheet[lab+3][scol:ecol]]
				#print lab,wave,len(data_row),scol,ecol,data_row[0],data_row[1],data_row[-1]
				swave=str(int(wave))
				sheetdata[swave+'nm'][alllabels[lab]]=np.array(data_row)
				sheetdata[swave+'nm'][alllabels[lab]+'_max']=max(data_row)
				sheetdata[swave+'nm'][alllabels[lab]+'_min']=min(data_row)


	return sheetdata

def collect_Biotek(sheet):
	sheetdata=NestedDict()

	#sheet=readtxt(ilist[0])

	allwells=getallwells()

	rownames=[row[0] if len(row)!=0 else '' for row in sheet]

	#print rowname


	#Currently works with single wavelength

	#print time_row
	time_ids=[ rin for rin,r in enumerate(rownames) if 'Time' in str(r) ]
	OD_ids=[ rin for rin,r in enumerate(rownames) if 'T' in str(r) and 'OD:' in str(r) ]

	time_row=sheet[time_ids[0]]
	temp_row=sheet[OD_ids[0]]

	# print time_ids
	# print OD_ids



	#print temp_row

	#nm_labels=[ rn if 'OD:' in str(rn) and 'T' not in str(rn) else '' for rn in rownames]
	nm_labels=[ str(sheet[rin][0].split()[1]) for rin in OD_ids]
	nm_labelsm={ rin:str(sheet[rin][0].split()[1]) for rin in OD_ids}

	waves_nmm={rin: numerize(str(sheet[rin][0].split()[1]).replace('OD:','').replace('GFP:','')) for rin in OD_ids}

	# print nm_labels

	waves=[numerize(wlen.replace('OD:','').replace('GFP:','')) for wlen in nm_labels if wlen!='']
	waves_nm=[str(int(w))+'nm' for w in waves]

	# print waves

	time_t=[time_to_sec(tval) for tval in time_row if tval not in ['Time','','0:00:00']]


	length=len(time_t)
	#print time_t

	#Find time step
	timestep=round_to(float(time_t[-1]-time_t[0])/(length-1),1)
	#print timestep

	#Timestep in mins
	timemax_min=int((length-1)*timestep/60)
	#timemax_min=int(round_to(float(time_t[-1])/60,5))

	timemax_h,timemax_remmin=divmod(timemax_min,60)

	time=np.linspace(0,timemax_min*60,length)
	#timestep=round_to(float(timemax_min*60)/(length-1),1)

	#print time

	temp=[float(t) for t in temp_row[1:] if t not in ['']]

	#print time_t[-1],timemax_h,time[-1],timestep

	alllabels=[rn if rn in allwells else '' for rn in rownames]
	labels=list(set([l for l in alllabels if l!='']))


	plsize=len(labels)
	if plsize not in [12,48,96,384]:
		buffer=True
	else:
		buffer=False




	sheetdata['Labels']=labels
	sheetdata['Spectra']=waves_nm
	sheetdata['Time']=time
	sheetdata['Temp']=temp
	sheetdata['Time_max']=timemax_min
	sheetdata['Time_step']=timestep
	sheetdata['Wells']=len(labels)
	sheetdata['Used wells']=plsize
	sheetdata['Buffered']=str(buffer)
	sheetdata['Figures']=waves_nm
	#sheetdata['File']=inm
	print "Wavelengths: {}".format(waves)
	print "Run time {}, step {}min in {} wells\n".format(str(datetime.timedelta(minutes=timemax_min)),timestep/60, len(labels))

	#sys.exit(0)

	for rid,row in enumerate(sheet):
		#print row
		well=str(row[0])
		#print well
		if well!='' and well in allwells:
			OD_sel=[ODid for ODid in reversed(OD_ids) if rid>ODid][0]
			swave=str(waves_nmm[OD_sel])+'nm'
			#print OD_sel,swave
			data_row=[float(val) for val in row[1:] if val not in ['']]
			#print swave,len(data_row),data_row[0],data_row[1],data_row[-1]
			sheetdata[swave][well]=np.array(data_row)

			sheetdata[swave][well+'_max']=max(data_row)
			sheetdata[swave][well+'_min']=min(data_row)


	#
	# for wnum, wave in enumerate(waves):
	# 	wavestart=nm_labels.index('OD:{}'.format(wave))
	# 	if len(waves)>1 and wnum < len(waves)-1:
	# 		waveend=nm_labels.index('OD:{}'.format(waves[wnum+1]))-1
	# 	else:
	# 		waveend=len(alllabels)
	# 	for rnum,well in enumerate(alllabels):
	# 		if well!='' and rnum>wavestart and rnum<waveend:
	# 				data_row=[float(val) for val in sheet[rnum][1:] if val not in ['']]
	# 				#print lab,wave,len(data_row),scol,ecol,data_row[0],data_row[1],data_row[-1]
	# 				swave=str(int(wave))
	# 				sheetdata[swave+'nm'][well]=np.array(data_row)
	#
	# 				sheetdata[swave+'nm'][well+'_max']=max(data_row)
	# 				sheetdata[swave+'nm'][well+'_min']=min(data_row)
	#

	return sheetdata

def collect(ilist):
	data=NestedDict()
	for ifl in sorted(ilist):
		print ifl
		ipt, inm, itp = filename(ifl)
		plate=inm
		#Currently formatting is determined by filetype
		if itp=='xlsx':
			sheet=readxls_s(ifl)
			data[plate]=collect_Tecan(sheet)
		elif itp=='asc':
			sheet=readtext(ifl)
			data[plate]=collect_Tecan(sheet)
		elif itp=='txt':
			sheet=readtext(ifl)
			data[plate]=collect_Biotek(sheet)
		else:
			print 'Unknown file format: {}!'.format(itp)
			sys.exit(1)
		#print data[plate]
		data[plate]['File']=ifl

	return data


def getallwells():
	r = [l for l in string.ascii_uppercase[0:16]]
	c = [str(i) for i in range(1,25)]
	allwells=[]
	for rn in r:
		for cn in c:
			#print rn+cn
			allwells.append(rn+cn)

	return allwells




def time_to_sec(tstr):
	h,m,s=tstr.split(':')
	seconds=int(s)+60*int(m)+3600*int(h)
	return seconds



def genlist(ifile):
	#Generates list of input files, checks their existance
	ilist=[]
	if ',' in ifile:
		ifiles=ifile.split(',')
		for ifl in ifiles:
			ilist.extend(genlist(ifl))
	else:
		ipath, iname, itype=filename(ifile)
		if itype in ['xls','xlsx','asc','txt'] and os.path.isfile(ifile):
			ilist.append(ifile)
		# elif itype in [''] and os.path.isfile(ifile):#'txt',
		# 	ifl=open(ifile,'r')
		# 	idata=ifl.read().split('\n')
		# 	idata=[fl.strip() for fl in idata if fl!='']
		# 	for fld in idata:
		# 		ilist.extend(genlist(fld))

		elif iname=='' and itype in ['xls','xlsx','asc','txt']:
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

def makesheets(data,descriptors,info):
    sheets=NestedDict()
    dt=[]
    gfit=[]
    regular=[]
    summary=[]

    header_temp=['Plate','Well','Data']
    
    header=header_temp
    
    if len(info.keys())>0:
        hasinfo=[pl for pl in info.keys() if len(info[pl].keys())>0]
        if len(hasinfo)>0:
            infokeys=sorted(info[hasinfo[0]].keys())
            header+=infokeys
        #else:
        #   infokeys=[]#header=header_temp+desckeys
    

    
    if len(descriptors.keys())>0:
        hasdesc=[pl for pl in descriptors.keys() if len(descriptors[pl].keys())>0]
        if len(hasdesc)>0:
            desckeys=sorted(descriptors[hasdesc[0]].keys())
            header+=desckeys
        #else:
            #desckeys=[]
            
    
            
    #header=header_temp+infokeys+desckeys
    
    
    
    
    #ghead=['a','c','t0']
    times=[]
    timemax=0
    maxplate=''
    #Adjustments that allow different length of experiment with a resolution of 5min
    for plate in data.keys():
        print plate, data[plate]['Time_max']
        if data[plate]['Time_max']>timemax:
            timemax=data[plate]['Time_max']
            maxplate=plate
    #Get max time to set scale
    timespan=len(data[maxplate]['Time'])
    timespandt=len(data[maxplate]['Time_dt'])

    allsumhead = ['A', 'lamda', 'u', 'tmax','tmaxf'] + \
                 ['a_log', 'c_log', 't0_log', 'tmax_log'] + \
                 ['Max_600nm', 'Max_600nm_log'] + \
                 ['Int_600nm_f', 'Int_600nm_log']

    # ['a_log', 'c_log', 't0_log', 'tmax_log'] + \
    selsums = ['Max_600nm_f', 'Max_600nm_log'] + \
              ['Int_600nm_f', 'Int_600nm_f_log']
    # ,'Max_Growth','Max_Growth_log','24h_Growth','24h_Growth_log','Int_Growth','Int-tmax_Growth','Int-tmaxf_Growth']#'a','c','t0'


    for plate in data.keys():
        sums = data[plate]['Summary']
        output = ['Summary'] + data[plate]['Figures']
        #output=data[plate]['Figures']
        #output=output+['GrowthFit']
        labels=data[plate]['Labels']
        for fig in output:
            #print plate,fig

            for well in labels:
                datrow=data[plate][fig][well]
                rowhead=[plate,well,fig]
                
                if len(info.keys())>0:
                    if len(info[plate].keys())>0:
                        rowhead=rowhead+[info[plate][ik] for ik in infokeys]
                    else:
                        rowhead=rowhead+['']*len(infokeys)


                if len(descriptors.keys())>0:
                    if len(descriptors[plate].keys())>0:
                        rowhead=rowhead+[descriptors[plate][dk][well] for dk in desckeys]
                    else:
                        rowhead=rowhead+['']*len(desckeys)

                    
                    

                if fig == 'Summary':
                    datrow = sums['GrowthFit'][well]+sums['GrowthFit_log'][well] + [sums[sm][well] for sm in selsums]
                    #	             #['a_log', 'c_log', 't0_log', 'tmax_log'] + \
                elif '_dt' in fig:
                    if not isinstance(datrow,list):
                        datrow=datrow.tolist()
                    datrow=datrow+['']*(timespandt-len(datrow))

                elif fig=='GrowthFit':
                    if not isinstance(datrow,list):
                        datrow=datrow.tolist()
                    datrow=datrow
                else:
                    if not isinstance(datrow,list):
                        datrow=datrow.tolist()
                    datrow=datrow+['']*(timespan-len(datrow))

                #print len(datrow), timespan,timespan-len(datrow)
                newrow=rowhead+datrow


                if '_dt' in fig:
                    dt.append(newrow)
                elif fig=='GrowthFit':
                    gfit.append(newrow)
                elif fig == 'Summary':
                    summary.append(newrow)
                else:
                    regular.append(newrow)
    #sys.exit(1)
    time_dt=data[maxplate]['Time_dt']
    time_lin=data[maxplate]['Time']
    time_lin=time_lin.tolist()
    time_dt=time_dt.tolist()

    summary.insert(0,header+allsumhead)
    dt.insert(0,header+time_dt)
    regular.insert(0,header+time_lin)
    gfit.insert(0,header)

    sheets['Data_dt']=dt
    sheets['Summary'] = summary
    sheets['Data']=regular
    sheets['Growth_fit']=gfit

    #print len(dt[0]),len(dt[1])
    #print dt[0],dt[1]

    return sheets


def writesheets(sheets,odir):
	#Writes organized data to file.
	#odir=dircheck('Split')
	#print sheets.keys()
	for fig in sheets.keys():
		#print fig
		oname='{}/{}.csv'.format(odir,fig)		
		sheet=sheets[fig]
		#print len(sheet)
		f=open(oname,"wb")
		ofile=csv.writer(f, delimiter=',') # dialect='excel',delimiter=';', quotechar='"', quoting=csv.QUOTE_ALL
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


def readtext(ifile):
	f=open(ifile,'r')
	sheet=[]
	for l in f:
		row=[cell.strip() for cell in l.replace('\t',',').split(',')] #re.split(r"[,\t]+",l.replace('\t',','))
		row=[numerize(cell) for cell in row]
		sheet.append(row)
	f.close()
	return sheet

def readxls_s(ifile):
	book=xlrd.open_workbook(ifile,formatting_info=False)
	data=book.sheet_by_name(book.sheet_names()[0]) #Use first sheet
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


def readdesc(ilist,dlist,subs):
	descriptors=NestedDict()
	for din in range(0,len(dlist)):
		ipath, iname, itype=filename(ilist[din])
		dfile=dlist[din]
		sub=subs[din]
		if not sub:
			substitutes={}
		else:
			substitutes={par.split(':')[0]:par.split(':')[1] for par in sub }
		book=xlrd.open_workbook(dfile,formatting_info=False)
		for nm in book.sheet_names():
			data=book.sheet_by_name(nm)
			headers=data.row_values(0)[1:]
			for r in range(1,data.nrows):
				rh=data.row_values(r)[0]
				row=data.row_values(r)[1:]
				for vin in range(0,len(row)):
					val=row[vin]
					ch=headers[vin]
					if nm in substitutes.keys():
						descriptors[iname][nm][rh+str(int(ch))]=substitutes[nm]
					else:
						descriptors[iname][nm][rh+str(int(ch))]=val
	return descriptors
	


#----------------------------------

if __name__ == "__main__":
	sys.exit(main())

