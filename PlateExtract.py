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
from collections import OrderedDict
from collections import defaultdict
from collections import Counter
from scipy.optimize import curve_fit
import numpy as np
from operator import itemgetter
from itertools import groupby


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
<--------------------------------------------->
	'''



def main(argv=None):

	ifile=""
	dfile=""
	msize=20
	mode="Zaslaver"
	filterf=''
	odir='Output'
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
			if argument in ("pqr", "--onlypqr"):
				pqr = True
			


	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		print >> sys.stderr, "\t for help use --help"
		return 2

	#Check for the input integrity

	print '''\n\n---------------------------\n\n'''

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
			raise Exception("No experiment file specified.")
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


	#----------------------------


	print optionsset %vars()


	
	modes={'Zaslaver':['600nm','535nm'],'Biolog':['600nm','700nm']}
	genes=csvreader(dfile)
	#sys.exit(1)
	#ilist=[ifile]
	waves=modes[mode]
	ilist=genlist(ifile)
	#ilist=[cfile,efile]
	checkfiles(ilist)	
	data=collect(ilist)
	data=analyze(data,waves,filterf)
	dirn=dircheck(odir)
	#print data.keys()
	data=growthfit(data)
	plotall(data,'Growth')
	plotall(data,'Fluorescence_norm')

	
	

	

	
	
	
	

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

def plotall(data,fg):
	align=True
	ts=5
	for plate in data.keys():
		print 'Plate number {}'.format(plate)
		for tp in data[plate].keys():
			#x=data[plate][tp]['Time']
			labels=data[plate][tp]['Labels']
			if 'A12' in labels:
				labels.remove('A12')
			if plate=='21':
				labels=labels[:9]
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
			
			for l in labels:
				x=data[plate][tp]['Time']
				y=data[plate][tp][fg][l]
				if norm and fg in ['Growth','600nm']:
					y=(y-min(y))/max(y)
				if align==True:
					a,c,t0=data[plate][tp]['GrowthFit'][l]['Log']
					x=(x+(ts*3600-t0*60))/3600
				else:
					x=x/3600
				plt.plot(x,y,label=plate+l )

	plt.show()

def growthfit(data):
	norm=True
	allfit=NestedDict()
	if norm:
		base=0.1
		top=0.2
		y0=0.2
	else:
		base=0.02
		top=0.07
		y0=0.1
	
	for plate in data.keys():
		for tp in data[plate].keys():
			x=data[plate][tp]['Time']
			labels=data[plate][tp]['Labels']
			if plate=='21':
				labels=labels[:9]
			if tp=='Control':
				plt.figure(0)
				plt.title('Growth without Metformin')

			
			elif tp=='Experiment':
				plt.figure(1)
				plt.title('Growth with Metformin')

			plt.xlim([0,20])
			if norm:
				plt.ylim([0,1])
			else:
				plt.ylim([0,0.4])
			
			for l in labels:
				y=data[plate][tp]['Growth'][l]
				if norm:
					y=(y-min(y))/max(y)
				x2,y2=cut(x, y, base, top)
				plt.plot(x/3600,y,label=plate+l)
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
					print '{} growth rate: {}, start: {}'.format(ynm,a*3600,t0)
					for par,nm in IT.izip([a,t0],['a','t0']):
						if allfit[tp][ynm][nm]:
							allfit[tp][ynm][nm]=allfit[tp][ynm][nm]+[par]
						else:
							allfit[tp][ynm][nm]=[par]
					
					data[plate][tp]['GrowthFit'][l][ynm]=[a,c,t0]
	
	plt.figure(2)
	sbn=np.arange(0,1000,10)
	plt.hist(allfit['Control']['Log']['t0'],bins=sbn,label='Without metformin - Log')
	plt.hist(allfit['Experiment']['Log']['t0'],bins=sbn,label='With metformin - Log')
	plt.xlabel('Start, min')
	plt.ylabel('Number')
	plt.title('Growth start time with/without metformin - Log')
	plt.legend()
	
	plt.figure(3)
	slbn=np.arange(0,1000,10)
	plt.hist(allfit['Control']['Linear']['t0'],bins=sbn,label='Without metformin - Linear')
	plt.hist(allfit['Experiment']['Linear']['t0'],bins=slbn,label='With metformin - Linear')
	plt.xlabel('Start, min')
	plt.ylabel('Number')
	plt.title('Growth start time with/without metformin - Linear')
	plt.legend()

	plt.figure(4)
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

	plt.figure(5)
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

	return data

def plot_2D(title,datac,datae,time,labels,plate_size,genes):
	#print title
	plate,fg=title.split('-')
	fig=plt.figure(figsize=(11.69,8.27), dpi=100)
	fig.suptitle(title)
	plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.1, hspace=0.1)
	plots={}
	#Need to check for both datasets
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
	elif fg in ['Fluorescence_dt']:
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
	#print grouper(x2)		
	return np.asarray(x2), np.asarray(y2)

def grouper(data):
	ranges = []
	for k, g in groupby(enumerate(data), lambda (i,x):i-x):
		ranges.append((map(operator.itemgetter(1), g)))

	return ranges

def Wiener(y, n):
	wi = sig.wiener(y, mysize=n)
	return wi

def Butter(x, y, par1, par2):	
	b, a = sig.butter(par1, par2)
	fl = sig.filtfilt(b, a, y)
	return fl



def plothuge(title,datac,datae,time,labels,genes):
	plots={}
	#Need to check for both datasets
	totalmaxc=round_to(max([max(datac[l]) for l in labels]),0.1)
	totalmaxe=round_to(max([max(datae[l]) for l in labels]),0.1)
	totalminc=round_to(min([min(datac[l]) for l in labels]),0.1)
	totalmine=round_to(min([min(datae[l]) for l in labels]),0.1)
	totalmax=max([totalmaxc,totalmaxe])
	totalmin=min([totalminc,totalmine])
	x=time/3600
	for i,l in enumerate(labels):	
		#print l,i
		yc=datac[l]
		ye=datae[l]
		fig=plt.figure(i)
		fig.suptitle(title+'_'+l+'-'+genes[l]['Gene'])
		plt.plot(x,yc,'r-',x,ye,'b-')
		#plt.text(0.1, 0.8, genes[l]['Gene'], fontsize=12, transform=transAxes)
		plots[l]=fig
		#plt.show()
	return plots


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



def analyze(data,waves,filterf):
	msize=20
	par1=4
	par2=0.1
	for plate in sorted(data.keys()):
		for tp in data[plate].keys():		
			if len([nm for nm in waves if nm in data[plate][tp]['Waves']])==2:
				time=data[plate][tp]['Time']
				dt=data[plate][tp]['Time'][1]-data[plate][tp]['Time'][0]
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
					data[plate][tp]['Fluorescence'][well]=fluor
					data[plate][tp]['Fluorescence_norm'][well]=fluor_norm

					data[plate][tp]['Fluor'][well+'_max']=fluor_max
					data[plate][tp]['Fluor'][well+'_min']=fluor_min
					#data[plate][tp]['Fluor_lp'][well+'_max']=fluor_lp_max
					#data[plate][tp]['Fluor_lp'][well+'_min']=fluor_lp_min

					#data[plate][tp]['Fluor_lp_norm'][well]=(data[plate][tp]['Fluor_lp'][well]-fluor_lp_min)/(fluor_lp_max-fluor_lp_min)

				

					data[plate][tp]['Fluorescence_dt'][well]=np.diff(fluor_norm)/dt
					#data[plate][tp]['Fluor_lp_dt'][well]=np.diff(fluor_lp)/dt
				
				
					

				data[plate][tp]['Figures']=data[plate][tp]['Figures']+['Growth','Fluorescence', 'Fluorescence_norm','Fluorescence_dt'] #'Fluor_norm','Fluor_lp','Fluor_lp_norm','Fluor_dt','Fluor_lp_dt'

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



def writesheets(sheets):
	#Writes organized data to file.
	#odir=dircheck('Split')
	for i in sheets.keys():
		if i=='Summary':
			oname=i+".csv"
		else:
			if '_' in i:
				tp=i.split('_')[1]
				nm=i.split('_')[0]
			else:
				nm=i
				tp='results'	
			oname="{}/{}_{}.csv".format(nm,nm,tp)
		ofile=csv.writer(open(oname,"wb"), dialect='excel') #,delimiter=';', quotechar='"', quoting=csv.QUOTE_ALL
		for row in sheets[i]:
			row = [item.encode("utf-8") if isinstance(item, unicode) else str(item) for item in row]
			ofile.writerow(row)
		#ofile.close()
	

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

