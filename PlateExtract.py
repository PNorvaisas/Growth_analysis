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
import numpy as np


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
	-c <file>    Name of control file
	-e <file>    Name of experiment
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
    Control:	%(cfile)s
 Experiment:	%(efile)s
	 DB:	%(dfile)s
     Filter:	%(filterf)s
        Out:	%(odir)s
<--------------------------------------------->
	'''



def main(argv=None):

	cfile=""
	efile=""
	dfile=""
	msize=20
	mode="Zaslaver"
	filterf=''
	odir='Output'
	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "hc:e:f:d:o:", ["help"])
		except getopt.error, msg:
			raise Usage(msg)


		for option, value in opts:
			if option in ("-h", "--help"):
				usage()
				return	
			if option in ("-c", "--cont"):
				cfile=value
			if option in ("-e", "--exp"):
				efile=value
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



	try:
		if cfile!="":
			cpath, cname, ctype = filename(cfile)			
		else:
			raise Exception("No control file specified.")
		if efile!="":
			epath, ename, etype = filename(efile)			
		else:
			raise Exception("No experiment file specified.")
		if dfile!="":
			dpath, dname, dtype = filename(dfile)			
		else:
			raise Exception("No experiment file specified.")
		if filterf!='' and filterf in ['wiener','butter']:
			print "{} filter selected".format(filterf)
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
	#ilist=genlist(ifile)
	ilist=[cfile,efile]	
	data=collect(ilist)
	data=analyze(data,waves,filterf)
	plate=cfile.split('_')[1]
	dirh=dircheck(odir+'/Huge')
	zoom=['B1','B2','B3']
	print plate
	for fg in data['Control']['Figures']:
		print "Plotting %(fg)s..." %vars()
		if 'dt' in fg:
			plt, plots=plot_2D(plate+"-"+fg,data['Control'][fg],data['Experiment'][fg],data['Control']['Time_dt'],data['Control']['Labels'],data['Control']['Wells'],genes[plate])
		else:
			plt, plots=plot_2D(plate+"-"+fg,data['Control'][fg],data['Experiment'][fg],data['Control']['Time'],data['Control']['Labels'],data['Control']['Wells'],genes[plate])
		data['Joint'][fg]=plt
		plt.savefig('{}/{}.pdf'.format(dirn,plate+'_'+fg))

		#if 'dt' in fg:
		#	plotsh=plothuge(plate+'_'+fg,data['Control'][fg],data['Experiment'][fg],data['Control']['Time_dt'],zoom,genes[plate])
		#else:
		#	plotsh=plothuge(plate+'_'+fg,data['Control'][fg],data['Experiment'][fg],data['Control']['Time'],zoom,genes[plate])
		#for l in zoom:
		#	plotsh[l].savefig('{}/{}.pdf'.format(dirh,plate+'_'+fg+'_'+l))
	

	

	

	
	
	
	

#-------------Functions------------

def plot_2D(title,datac,datae,time,labels,plate_size,genes):
	fig=plt.figure(figsize=(11.69,8.27), dpi=100)
	fig.suptitle(title)
	plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.1, hspace=0.1)
	plots={}
	#Need to check for both datasets
	totalmaxc=round_to(max([max(datac[l]) for l in labels]),0.1)
	totalmaxe=round_to(max([max(datae[l]) for l in labels]),0.1)
	totalminc=round_to(min([min(datac[l]) for l in labels]),0.1)
	totalmine=round_to(min([min(datae[l]) for l in labels]),0.1)
	totalmax=max([totalmaxc,totalmaxe])
	totalmin=min([totalminc,totalmine])
	if "dt" in title:
		ymin=-totalmax
	else:
		ymin=0
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

		plt.xticks(np.linspace(0, max(x), 2),['']+list(np.linspace(0, max(x), 2).astype(int)[1:]), rotation='vertical')	
		plt.yticks(np.linspace(ymin, totalmax, 2),['']+list(np.around(np.linspace(totalmin, totalmax, 2),2)[1:]))
		plt.ylim([totalmin,totalmax])
		plots[l].text(0.1, 0.8, genes[l]['Gene'], fontsize=10, transform=plots[l].transAxes)
		#exec(l+"='plt.subplot(8,12,{},sharex={},sharey={})'".format(v,,sh_y))

		plt.plot(x,yc,'r-',x,ye,'b-')
		
		#plt.title(l)

	return plt, plots

def Wiener(y, n):
	wi = sig.wiener(y, mysize=n)
	return wi

def Butter(x, y, nyf, par1, par2):
	b, a = sig.butter(par1, par2/nyf)
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


def analyze(data,waves,filterf):
	msize=20
	par1=4
	par2=1.57
	for tp in data.keys():		
		if len([nm for nm in waves if nm in data[tp]['Waves']])==2:
			time=data[tp]['Time']
			dt=data[tp]['Time'][1]-data[tp]['Time'][0]
			npts=len(time)
			dtm=time[-1]/npts
			nyf=0.5/dtm
			data[tp]['Time_dt']=(time+dt/2)[:-1]
			Ufluor=data[tp]['535nm']['C10']/data[tp]['600nm']['C10']
			for well in data[tp]['Labels']:
				growth=data[tp]['600nm'][well]-data[tp]['600nm'][well+'_min']
				fluor=data[tp]['535nm'][well]/data[tp]['600nm'][well]


				if filterf=='wiener':
					growth=Wiener(growth-growth[0],msize)
					fluor=Wiener(fluor-fluor[0],msize)
					fluor_norm=fluor-Ufluor
					fluor_norm=Wiener(fluor_norm-fluor_norm[0],msize)
				elif filterf=='butter':
					growth=Butter(time,growth-growth[0],nyf,par1,par2)
					fluor=Butter(time,fluor-fluor[0],nyf,par1,par2)
					fluor_norm=fluor-Ufluor
					fluor_norm=Butter(time,fluor_norm-fluor_norm[0],nyf,par1,par2)
				else:
					fluor_norm=fluor-Ufluor
								
				growth=growth-min(growth)
				fluor=fluor-min(fluor)
				fluor_norm=fluor_norm-min(fluor_norm)
				

				fluor_max=max(fluor)
				fluor_min=min(fluor)
				
				data[tp]['Growth'][well]=growth
				data[tp]['Fluorescence'][well]=fluor
				data[tp]['Fluorescence_norm'][well]=fluor_norm

				data[tp]['Fluor'][well+'_max']=fluor_max
				data[tp]['Fluor'][well+'_min']=fluor_min
				#data[tp]['Fluor_lp'][well+'_max']=fluor_lp_max
				#data[tp]['Fluor_lp'][well+'_min']=fluor_lp_min

				#data[tp]['Fluor_lp_norm'][well]=(data[tp]['Fluor_lp'][well]-fluor_lp_min)/(fluor_lp_max-fluor_lp_min)

				

				data[tp]['Fluorescence_dt'][well]=np.diff(fluor_norm)/dt
				#data[tp]['Fluor_lp_dt'][well]=np.diff(fluor_lp)/dt
				
				
					

			data[tp]['Figures']=data[tp]['Figures']+['Growth','Fluorescence', 'Fluorescence_norm','Fluorescence_dt'] #'Fluor_norm','Fluor_lp','Fluor_lp_norm','Fluor_dt','Fluor_lp_dt'

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
	for ifl,tp in IT.izip(ilist,['Control','Experiment']):
		print "{} file: {}".format(tp,ifl)
		ipt, inm, itp = filename(ifl)	
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
		data[tp]['Labels']=labels[3:]
		data[tp]['Waves']=waves
		#print data[inm]['Waves']
		data[tp]['Time']=time
		data[tp]['Temp']=temp
		data[tp]['Time_max']=timemax_h
		data[tp]['Time_step']=timestep
		data[tp]['Wells']=nrows-3
		data[tp]['Figures']=waves
		data[tp]['File']=inm
		
		print "Wavelengths: {}, {}".format(*waves)
		print "Run time {}h, step {}min in {} wells".format(timemax_h,timestep/60, nrows-3)
		for row in range(3,sheet.nrows):
			for wave in waves:
				data_row=[60000 if val=="Overflow" else val for val in sheet.row_values(row)[length*(waves.index(wave)):length*(waves.index(wave)+1)]]
				data[tp][wave][labels[row]]=np.array(data_row)
				data[tp][wave][labels[row]+'_max']=max(data_row)
				data[tp][wave][labels[row]+'_min']=min(data_row)
				#print wave
				#print labels[row]
				#print data[inm][wave][labels[row]]
			
		#print waves
		#print nrows, length
		#print timemax_h, timestep
	return data



def genlist(ifile):
	#Generates list of input files, checks their existance
	ilist=[]
	ipath, iname, itype=filename(ifile)
	if itype in ['xls','xlsx'] and os.path.isfile(ifile):
		ilist.append(ifile)
	elif itype in ['txt',''] and os.path.isfile(ifile):
		ifl=open(ifile,'r')
		idata=ifl.read().split('\n')
		idata=[fl.trim() for fl in idata if fl!='']
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

