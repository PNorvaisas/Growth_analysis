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
#from scipy import interpolate as ip


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
	loadc=False
	#Constant values may change
	g750=0.748596
	d750=0.4
	comparison=False
	integrals=''
	
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
			if argument in ("loadc", "--loadc"):
				loadc = True
			if argument in ("integrals", "--integrals"):
				integrals='integrals'
			if argument in ("fullintegrals", "--fullintegrals"):
				integrals='fullintegrals'
			


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

	if odir!="" and not (load or loadc):
		dirn=dircheck(odir)


	if load:
		f = open('{}/Biolog_data_all.pckl'.format(odir),'rb')
		data,metabolites,ilist,info,uniques= pickle.load(f)
		f.close()

	else:

		if loadc:
			f = open('{}/Biolog_data_collected.pckl'.format(odir),'rb')
			data,metabolites,ilist,info,uniques= pickle.load(f)
			f.close()
		else:
			ilist,info=genlist(ifile)
			uniques=uniquecomb(info,['Plate','Strain','Sugar_20mM','Inoculum','Uracil_uM'],'Type')
			data=collect(ilist)
			# f = open('{}/Biolog_data_collected.pckl'.format(odir), 'w')
			# pickle.dump([data,metabolites,ilist,info,uniques], f)
			# f.close()

		#sys.exit(1)
		data=analyze(data,g750,d750)
		data=growthfit(data)
		# f = open('{}/Biolog_data_all.pckl'.format(odir), 'w')
		# pickle.dump([data,metabolites,ilist,info,uniques], f)
		# f.close()



	sheets=makesheets(data,metabolites,info)
	writesheets(sheets,odir,sep=',')


	plot_comparison(data,metabolites,odir,['590nm','750nm','590nm_f','750nm_f','Growth_Respiration',
	                                       '750nm_dt','750nm_log'],info,uniques,integrals)
	#,'Growth','Growth_log','Growth_dt'
	#'590nm_dt','590nm_f','590nm_log',

#-------------Functions------------



def plot_comparison(data,metabolites,dirn,figs,info,uniques,integrals):

	for uk in sorted(uniques.keys()):
		comvalues=uniques[uk]['Compare-by-values']
		comname=uniques[uk]['Compare-by']
		unks=uniques[uk]['Unique-keys']
		unval=uniques[uk]['Unique-values']
		#print unks,unval

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

		#print comvalues
		#print uniques

		controls=0
		treatments=0
		#len(comvalues)==1 and
		#Make it universal for all kinds of comparisons
		if all([cvt in ['Control', 'Treatment'] for cvt in comvalues]):
			#print comvalues
			for cv in comvalues:
				if cv=='Control':
					controls=controls+1
					cgroup=uniques[uk][cv]
				elif cv=='Treatment':
					treatments=treatments+1
					egroup=uniques[uk][cv]
		else:
			print 'Unknown types found as descriptors: {}'.format([cvt for cvt in comvalues if cvt not in ['Control', 'Treatment']])
			print 'Currently supported descriptors: {}'.format(['Control', 'Treatment'])
		# print comvalues
		# print uniques
		print
		if controls<1:
			cgroup=[]
		if treatments<1:
			egroup=[]
		#group=uniques[uk][cv]
		#strain=info[plate]['Strain']
		#labels=data[plate]['Labels']


		# ref=data[plate]['Control']
		# exp=data[plate]['Treatment']

		figures_temp=data[cgroup[0]]['Figures']
		figures_temp=figures_temp+['Growth_Respiration']

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
			lbl=''
			for unk in unks:
				if unks.index(unk)==0:
					sep=''
				else:
					sep='_'
				if unk=='Sugar_20mM':
					lbl=lbl+sep+'Sugar20mM-{}'
				elif unk=='Uracil_uM':
					lbl=lbl+sep+'Uracil{}uM'
				elif unk=='Replicate':
					lbl=lbl+sep+'Rep{}'
				elif unk=='Metformin_mM':
					lbl=lbl+sep+'Metf{}'
				else:
					lbl=lbl+sep+'{}'
			lbl=lbl.format(*unval)
			ttl=lbl.replace('_',' ')


			print "Plotting: "+ttl+' '+fg
			#Need to fix metabolites
			plot=plot_2D(ttl,fg,bplate,data,cgroup,egroup,metabolites[bplate],info,integrals)
			plot.savefig('{}/{}.pdf'.format(dirn,'{}_{}'.format(lbl,fg)))
			plot.close()

	return data

def plot_2D(ttl,fg,bplate,data,cgroup,egroup,metabolites,info,integrals):



	markers={'1':'-','2':'--','3':'-.','4':':'}
	metfc={'0':'r','25':'c','50':'b','75':'k','100':'k'}
	drugplates=['PM11C','PM12B','PM13B','PM14A','PM15B', 'PM16A', 'PM17A', 'PM18C',
	            'PM19', 'PM20B', 'PM21D', 'PM22D', 'PM23A', 'PM24C', 'PM25D']
	if bplate in drugplates:
		metin='Name'
	else:
		metin='Metabolite'
	plate_size=96
	#print title

	if '_dt' in fg:
		time=data[data.keys()[0]]['Time_dt']
	else:
		time=data[data.keys()[0]]['Time']
	labels=data[data.keys()[0]]['Labels']

	#fig=plt.figure(figsize=(11.69,8.27), dpi=100)

	title='{} {}'.format(ttl,fg)
	#title='{} {} {}'.format(bplate,strain,fg)

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
	xmax=24

	ticks=3
	xlabel='Time, h'
	ylabel=''
	decimals=1

	if fg=='Growth_Respiration':
		totalmax=1
		totalmin=0
		ticks=3
		ylabel='OD(750nm)'
		decimals=1
		xlabel='OD(590nm)'
		xmax=2

	if fg=='Growth':
		totalmax=1
		totalmin=0
		ticks=3
		ylabel='OD'
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

	if '590nm' in fg and '_log' not in fg and '_dt' not in fg and '_diff' not in fg:
		totalmax=2
		ticks=3
		ylabel='OD'
		decimals=1

	if '750nm' in fg and '_log' not in fg and '_dt' not in fg and '_diff' not in fg:
		totalmax=0.5 if bplate=='PM5' else 1
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

	repmax=1
	metconcs=[]
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
		label=greek_check(metabolites[l][metin],12)
		plt.text(0.05, 0.9, label, fontsize=7,verticalalignment='top',transform=ax.transAxes)

		linx=np.array([0,2])
		#Need to unify coloring in universal way
		for grp,mrkc in IT.izip([cgroup,egroup],['r','b']):
			for gdata in grp:
				#Needs further refinement
				if 'Metformin_mM' in info[gdata].keys(): # and info[gdata]['Metformin_mM']!='50'
					mrkc=metfc[info[gdata]['Metformin_mM']]
					metconcs.append(info[gdata]['Metformin_mM'])
				if 'Replicate' in info[gdata].keys():
					mrk=mrkc+markers[info[gdata]['Replicate']] if int(info[gdata]['Replicate'])<5 else mrkc+markers['4']
					if int(info[gdata]['Replicate'])>repmax:
						repmax=int(info[gdata]['Replicate'])
				else:
					mrk=mrkc+'-'
				if fg=='Growth_Respiration':
					xnm=data[gdata]['590nm'][l]
					y=data[gdata]['750nm'][l]
					rho,const,d=data[gdata]['Summary']['Respiration'][l]
					if rho!=0:
						plt.plot(linx,lin(linx,rho,const),mrk,alpha=0.2)
					plt.plot(xnm,y,mrk)
					#plt.text(0.05, 0.05, 'rho={0:1.2f}, d={1:1.2f}'.format(rho,d),
					# fontsize=6,verticalalignment='top',transform=ax.transAxes)
				else:
					y=data[gdata][fg][l]
					#print len(x),len(y)
					plt.plot(x,y,mrk)
					if fg=='750nm_log':
						a,c,t0,tmax=data[gdata]['Summary']['GrowthFit_log'][l]
						if a>0:
							yfit=x*a+c
							plt.plot(x,yfit,mrk,alpha=0.5)
					if fg=='750nm_f' and integrals=='integrals':
						A,lam,u,tmax,tmaxf=data[gdata]['Summary']['GrowthFit'][l]
						if 0<tmaxf<=24:
							plt.fill_between(x, 0,y,where=x<=tmaxf, facecolor='red' if mrkc=='r' else 'blue',alpha=0.1)#, interpolate=True
						if A>0 and lam<np.inf and u>0:
							yfit=growth(x,A,lam,u)
							plt.plot(x,yfit,mrk.replace('-','-.'),alpha=0.5)
					if fg=='750nm_f' and integrals=='fullintegrals':
							plt.fill_between(x, 0,y,where=x<=24, facecolor='red' if mrkc=='r' else 'blue',alpha=0.1)
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
	typelines=[]
	typelabels=[]
	if len(metconcs)>0:
		if '0' in metconcs:
			typelines.append(mlines.Line2D([], [], color='red', label='Control'))
			typelabels.append('Control')
		if '25' in metconcs:
			typelines.append(mlines.Line2D([], [], color='cyan', label='Metformin 25mM'))
			typelabels.append('Metformin 25mM')
		if '50' in metconcs:
			typelines.append(mlines.Line2D([], [], color='blue', label='Metformin 50mM'))
			typelabels.append('Metformin 50mM')
		if '75' in metconcs:
			typelines.append(mlines.Line2D([], [], color='blue', label='Metformin 75mM'))
			typelabels.append('Metformin 75mM')

	else:
		blue_line = mlines.Line2D([], [], color='blue', label='Treatment')
		red_line = mlines.Line2D([], [], color='red', label='Control')
		typelines=[red_line,blue_line]
		typelabels=['Control','Treatment']

	replines=[]
	replabels=[]

	for rep in range(1,repmax+1):
		replines.append( mlines.Line2D([], [],linestyle=markers[str(rep)] if rep<5  else ':', color='black', label='Replicate {}'.format(rep)))
		replabels.append('Replicate {}'.format(rep))

	plt.figlegend(typelines,typelabels,'upper right',prop={'size':7})#
	plt.figlegend(replines,replabels,'upper left',prop={'size':7})#

	return plt

def uniquecomb(info,sortby,compareby):
	summ=[]
	uns=NestedDict()
	preskeys=info[info.keys()[0]].keys()
	inkeys=[k for k in sortby if k in preskeys]
	missing=[k for k in sortby if k not in preskeys]
	if len(inkeys)>0:
		if len(inkeys)==len(sortby):
			print 'All defined keys have been found!'
		else:
			print 'Some keys were found!\n{}'.format(inkeys)
			print 'Keys for unique sets have not been found in Design file!\n{}'.format(missing)
		for i in info.keys():
			ks=[info[i][k] for k in inkeys]
			summ.append(ks)
		#print summ
		uniqs=[list(i) for i in set(tuple(i) for i in summ)]
		for unl in uniqs:
			unin='|'.join(unl)
			uns[unin]['Unique-keys']=inkeys
			uns[unin]['Unique-values']=unl
			uns[unin]['Compare-by']=compareby
			uns[unin]['Compare-by-values']=[]
			for i in info.keys():
				# print sortby
				# print unl
				# print type(unl)
				check=[info[i][uk]==uv for uk,uv in IT.izip(inkeys,unl)]
				# print check
				if all(check):
					if info[i][compareby] in uns[unin].keys():
						uns[unin][info[i][compareby]]=uns[unin][info[i][compareby]]+[i]
					else:
						uns[unin]['Compare-by-values']=uns[unin]['Compare-by-values']+[info[i][compareby]]
						uns[unin][info[i][compareby]]=[i]
	else:
		print 'None of the defined sorting keys have been found in Design file!\n{}'.format(missing)
		sys.exit(1)

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



def analyze(data,g750,d750):

	waves=['590nm','750nm']
	msize=20
	window=20
	thres=np.power(2.0,-5)
	thresd=np.power(2.0,-5)

	for plate in sorted(data.keys()):
		#Needs to be checked
		print 'Analyzing data in {}'.format(plate)
		ind16=data[plate]['Time'].tolist().index(16*3600)
		ind24=data[plate]['Time'].tolist().index(24*3600)
		time=data[plate]['Time']
		time_h=time/3600
		dt=time[1]-time[0]
		time_dt=(time+dt/2)[:-1]
		data[plate]['Time_dt']=time_dt

		for well in data[plate]['Labels']:
			#print plate, well
			if well!='A1':
				ref590=data[plate]['590nm_f']['A1']
				ref750=data[plate]['750nm_f']['A1']
				refgrowth=data[plate]['Growth']['A1']
				refdye=data[plate]['Dye']['A1']
				refl590=data[plate]['590nm_log']['A1']
				refl750=data[plate]['750nm_log']['A1']
				reflgrowth=data[plate]['Growth_log']['A1']
				refldye=data[plate]['Dye_log']['A1']

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

			raw750c,raw590c=cut(raw750,raw590,0.1,0.8,equalize=False)
			if len(raw750c)>10 and len(raw590c)>10:
				scale750=max(raw750c)-min(raw750c)
				scale590=max(raw590c)-min(raw590c)
				if scale590>0.1 and scale750>0.1:
					poptd, pcovd = curve_fit(lin, raw590c, raw750c)
					rho=poptd[0]
					const=poptd[1]
					d=(rho-g750)/(d750-rho)
				else:
					rho=0
					const=0
					d=0
			else:
				rho=0
				const=0
				d=0
			growthr=f750/(g750+d*d750)
			dyer=growthr*d



			loggrowth=np.log2(setbar(growthr,thres))
			dtgrowth=np.diff(growthr)/(dt/3600)

			logdye=np.log2(setbar(dyer,thres))
			dtdye=np.diff(dyer)/(dt/3600)


			if well=='A1':
				refrl590=log590-log590
				refrl750=log750-log750
				refrlgrowth=loggrowth-loggrowth
				refrldye=logdye-logdye
				refr590=f590-f590
				refr750=f750-f750
				refrgrowth=growthr-growthr
				refrdye=dyer-dyer
			else:
				refrl590=log590-refl590
				refrl750=log750-refl750
				refrlgrowth=loggrowth-reflgrowth
				refrldye=logdye-refldye
				refr590=f590-ref590
				refr750=f750-ref750
				refrgrowth=growthr-refgrowth
				refrdye=dyer-refdye

			diff590750=f590-f750
			diff590750log=log590-log750



			margin=0.01
			#Select which channel to use for growth fitting
			grow=f750
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
					freducedf = ip.UnivariateSpline(time_h, yreducedf, s=0)
					if len(freducedf.roots())>0:
						tmaxf=freducedf.roots()[0]
					else:
						tmaxf=np.inf
				else:
					tmaxf=np.inf

				yreduced = growc - maxg*(1-margin)
				try:
					freduced = ip.UnivariateSpline(timec, yreduced, s=0)
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

			# yreducedf = growth(tiem,*popt) - (popt[0]-popt[0]*margin) #maxg#
			# freducedf = ip.UnivariateSpline(time, yreducedf, s=0)
			# tmaxf=freducedf.roots()[0]
			# #print tmaxf


			#print tmax


			data[plate]['750nm_b'][well]=back750
			data[plate]['750nm_f'][well]=f750
			data[plate]['750nm_log'][well]=log750
			data[plate]['750nm_dt'][well]=dtf750

			data[plate]['590nm_b'][well]=back590
			data[plate]['590nm_f'][well]=f590
			data[plate]['590nm_log'][well]=log590
			data[plate]['590nm_dt'][well]=dtf590

			data[plate]['Growth'][well]=growthr
			data[plate]['Growth_log'][well]=loggrowth
			data[plate]['Growth_dt'][well]=dtgrowth

			data[plate]['Dye'][well]=dyer
			data[plate]['Dye_log'][well]=logdye
			data[plate]['Dye_dt'][well]=dtdye

			data[plate]['590nm_reflog'][well]=refrl590
			data[plate]['590nm_ref'][well]=refr590

			data[plate]['Growth_reflog'][well]=refrlgrowth
			data[plate]['Growth_ref'][well]=refrgrowth

			data[plate]['Dye_reflog'][well]=refrldye
			data[plate]['Dye_ref'][well]=refrdye

			data[plate]['750nm_reflog'][well]=refrl750
			data[plate]['750nm_ref'][well]=refr750

			data[plate]['590nm750nm_diff'][well]=diff590750
			data[plate]['590nm750nm_difflog'][well]=diff590750log


			data[plate]['Summary']['Respiration'][well]=[rho,const,d]
			data[plate]['Summary']['GrowthFit'][well]=[A,lam,u,tmax,tmaxf]


		data[plate]['Figures']=data[plate]['Figures']+['590nm_b','590nm_f','590nm_log','590nm_dt',
		                                               '750nm_b','750nm_f','750nm_log','750nm_dt',
		                                               #'Growth','Growth_log','Growth_dt',
		                                               #'Dye','Dye_log','Dye_dt',
		                                               #'Growth_ref','Growth_reflog',
		                                               #'Dye_ref','Dye_reflog',
		                                               #'750nm_ref','750nm_reflog',
		                                               #'590nm_ref','590nm_reflog',
		                                               #'590nm750nm_diff','590nm750nm_difflog',
		                                               ]
	#for plate in sorted(data.keys()):
		for fg in data[plate]['Figures']:
			for well in data[plate]['Labels']:
				#print plate,fg,well
				data[plate]['Summary']['Max_{}'.format(fg)][well]=max(data[plate][fg][well])
				tmax=data[plate]['Summary']['GrowthFit'][well][3]
				tmaxf=data[plate]['Summary']['GrowthFit'][well][4]
				if '_dt' in fg:
					data[plate]['Summary']['24h_{}'.format(fg)][well]=data[plate][fg][well][ind24-1]
					data[plate]['Summary']['16h_{}'.format(fg)][well]=data[plate][fg][well][ind16-1]
				else:
					data[plate]['Summary']['24h_{}'.format(fg)][well]=data[plate][fg][well][ind24]
					data[plate]['Summary']['16h_{}'.format(fg)][well]=data[plate][fg][well][ind16]
					data[plate]['Summary']['Int_{}'.format(fg)][well]=interp_integ(time_h, data[plate][fg][well], 24)
					data[plate]['Summary']['Int-tmax_{}'.format(fg)][well]=interp_integ(time_h, data[plate][fg][well], tmax if tmax<24 else 24) if tmax!=np.inf else np.inf
					data[plate]['Summary']['Int-tmaxf_{}'.format(fg)][well]=interp_integ(time_h, data[plate][fg][well], tmaxf if tmaxf<24 else 24) if tmaxf!=np.inf else np.inf
					data[plate]['Summary']['Int_{}_log'.format(fg)][well]=np.log2(data[plate]['Summary']['Int_{}'.format(fg)][well]) if data[plate]['Summary']['Int_{}'.format(fg)][well] >0 else None

	return data

def setbar(x,bar):
	x2=[xi if xi>bar else bar for xi in x]
	x2=np.array(x2)
	return x2

def interp_integ(x,y,t2):
	intgr=ip.UnivariateSpline(x, y, s=5,k=5).integral(0,t2)
	#intgr=spln
	return intgr

def readmet(ifile):
	print ifile
	nutrients=NestedDict()

	rdr=csv.reader(open(ifile,'r'), delimiter=',')
	data=[ln for ln in rdr]
	headers=data[0]
	headin={ hd : headers.index(hd) for hd in headers}
	nec=['Metabolite','EcoCycID','Plate','Well','Index']
	if all(n in headers for n in nec):
		print 'Necessary headers found!'
	else:
		print 'Missing essential headers in metabolites file!'
		print headers
		sys.exit(0)
	for ln in data[1:]:
		pl=str(numerize(ln[headin['Plate']])).strip().encode('ascii','ignore')
		wl=str(numerize(ln[headin['Well']])).strip().encode('ascii','ignore')
		for hd in headin.keys():
			nutrients[pl][wl][hd]=str(numerize(ln[headin[hd]])).strip().encode('ascii','ignore')
	return nutrients

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

	#filein=headers.index('File')
	# platein=headers.index('Plate')
	# strainin=headers.index('Strain')
	# typein=headers.index('Type')
	#print metin, ecoin,plate,well
	for ln in data[1:]:
		fl=str(ln[headin['File']]).encode('ascii','ignore').strip()
		if fl!="":
			for hd in headin.keys():
				info[fl][hd]=str(numerize(ln[headin[hd]])).strip().encode('ascii','ignore')
			#print info[fl]['Replicate']
			ilist.append(fl)

	#print genes
	return ilist, info



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

def growth(x,A,lam,u):
	return A/(1+np.exp((4*u/A)*(lam-x)+2))

def log_growth(x,A,lam,u):
	y=np.log2(A)-np.log2(1+np.exp((4*u/A)*(lam-x)+2))
	#print x,y
	return np.log2(A)-np.log2(1+np.exp((4*u/A)*(lam-x)+2))

def lin(x,a,c):
	y=x*a+c
	return y

def growthfit(data):
	y0=np.power(2.0,-4)
	for plate in data.keys():
		x=data[plate]['Time']
		x=x/3600
		labels=data[plate]['Labels']
		waves=data[plate]['Spectra']

		if len(waves)>1 and '750nm' in waves:
			wave='750nm'
		elif len(waves)>1 and '600nm' in waves:
			wave='600nm'
		elif len(waves)>1 and '595nm' in waves:
			wave='595nm'
		elif len(waves)>1 and '590nm' in waves:
			wave='590nm'
		else:
			wave=waves[0]

		for l in labels:
			y=data[plate][wave+'_log'][l]
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

	sheet=[r for r in sheet if len(r)>1]
	nrows=len(sheet)
	datarange=sheet[0].index('Well positions')
	nm_labels=[lab for lab in sheet[0] if lab not in ['Layout','Well positions','','Replicate Info']]
	print nm_labels
	if len(nm_labels)>1:
		starts=[sheet[0].index(lab) for lab in nm_labels]
	else:
		starts=[sheet[0].index(nm_labels[0])]
		if nm_labels[0]=='Raw data':
			nm_labels[0]='600'

	waves=[numerize(wlen) for wlen in nm_labels]
	#print 'Identified wavelengths: {}'.format(waves)

	print datarange
	#print sheet[0]

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

	timemax_h=int(round_to(float(time_t[-1])/3600,0.1))

	time=np.linspace(0,timemax_h*3600,length)

	timestep=round_to(float(time_t[-1])/(length-1),1)
	#timestep=time[-1]-time[-2]

	print time_t[-1],timemax_h,time[-1],timestep


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
	sheetdata['Time_max']=timemax_h
	sheetdata['Time_step']=timestep
	sheetdata['Wells']=len(labels)
	sheetdata['Used wells']=plsize
	sheetdata['Buffered']=str(buffer)
	sheetdata['Figures']=[str(w)+'nm' for w in waves]
	#sheetdata['File']=inm
	print "Wavelengths: {}".format(waves)
	print "Run time {}h, step {}min in {} wells\n".format(timemax_h,timestep/60, len(labels))
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

	alllabels=[rn if rn in allwells else '' for rn in rownames]

	labels=[l for l in alllabels if l!='']

	time_row=sheet[rownames.index('Time')]
	temp_row=sheet[rownames.index('T\xb0 OD:595')]

	nm_labels=[ rn  if 'OD:' in rn and 'T\xb0' not in rn else '' for rn in rownames]
	waves=[numerize(wlen.replace('OD:','')) for wlen in nm_labels if wlen!='']
	time_t=[time_to_sec(tval) for tval in time_row if tval not in ['Time','']]
	length=len(time_t)

	timemax_min=int(round_to(float(time_t[-1])/60,5))
	timemax_h,timemax_remmin=divmod(timemax_min,60)
	time=np.linspace(0,timemax_min*60,length)
	timestep=round_to(float(time_t[-1])/(length-1),1)

	temp=[float(t) for t in temp_row[1:] if t not in ['']]

	print time_t[-1],timemax_h,time[-1],timestep


	plsize=len(labels)
	if plsize not in [12,48,96,384]:
		buffer=True
	else:
		buffer=False

	sheetdata['Labels']=labels
	sheetdata['Spectra']=[str(int(w))+'nm' for w in waves]
	sheetdata['Time']=time
	sheetdata['Temp']=temp
	sheetdata['Time_max']=timemax_h
	sheetdata['Time_step']=timestep
	sheetdata['Wells']=len(labels)
	sheetdata['Used wells']=plsize
	sheetdata['Buffered']=str(buffer)
	sheetdata['Figures']=[str(w)+'nm' for w in waves]
	#sheetdata['File']=inm
	print "Wavelengths: {}".format(waves)
	print "Run time {}, step {}min in {} wells\n".format(str(datetime.timedelta(minutes=timemax_min)),timestep/60, len(labels))

	for wnum, wave in enumerate(waves):
		wavestart=nm_labels.index('OD:{}'.format(wave))
		if len(waves)>1 and wnum < len(waves)-1:
			waveend=nm_labels.index('OD:{}'.format(waves[wnum+1]))-1
		else:
			waveend=len(alllabels)
		for rnum,well in enumerate(alllabels):
			if well!='' and rnum>wavestart and rnum<waveend:
					data_row=[float(val) for val in sheet[rnum][1:] if val not in ['']]
					#print lab,wave,len(data_row),scol,ecol,data_row[0],data_row[1],data_row[-1]
					swave=str(int(wave))
					sheetdata[swave+'nm'][well]=np.array(data_row)

					sheetdata[swave+'nm'][well+'_max']=max(data_row)
					sheetdata[swave+'nm'][well+'_min']=min(data_row)


	return sheetdata

def collect(ilist):
	data=NestedDict()
	for ifl in sorted(ilist):
		print ifl
		ipt, inm, itp = filename(ifl)
		plate=ifl

		if itp=='xlsx':
			sheet=readxls(ifl)
			data[plate]=collect_Tecan(sheet)
		elif itp=='asc':
			sheet=readacs(ifl)
			data[plate]=collect_Tecan(sheet)
		elif itp=='txt':
			sheet=readtxt(ifl)
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



def readtxt(ifile):
	f=open(ifile,'r')
	rdr=csv.reader(f, delimiter=',')
	sheet=[ln for ln in rdr]
	f.close()
	return sheet


def time_to_sec(tstr):
	h,m,s=tstr.split(':')
	seconds=int(s)+60*int(m)+3600*int(h)
	return seconds




#
# def collect(ilist):
# 	data=NestedDict()
# 	for ifl in sorted(ilist):
# 		ipt, inm, itp = filename(ifl)
# 		if itp in ['xlsx','asc']:
# 			if itp=='xlsx':
# 				sheet=readxls(ifl)
# 			elif itp=='asc':
# 				sheet=readacs(ifl)
#
#
# 			nrows=len(sheet)
# 			nm_labels=sheet[0]
# 			#converting first row to set rearanges the labels
# 			#list(set(sheet[0]))
# 			#print nm_labels
# 			#print nrows
#
# 			waves=[numerize(wlen) for wlen in nm_labels if wlen not in ['Layout','Well positions','','Replicate Info']]
# 			#print waves
# 			lengths=[sheet[0].index(wave) for wave in waves]
# 			#Selection of time cells does not depend on the order of wavelengths
# 			length=max(lengths)
# 			time_row=sheet[1][:length]
# 			check=list(set([s for s in time_row if isinstance(s,float)]))
# 			if check:
# 				time_t=time_row
# 			else:
# 				time_t=[int(t.replace('s','')) for t in time_row]
# 			temp=[float(t.split()[0]) for t in sheet[2][:length]]
# 			timemax_h=round_to(float(time_t[length-1])/3600,1)
#
# 			timestep=round_to(float(time_t[length-1])/(length-1),1)
# 			time=np.linspace(0,timemax_h*3600,length)
# 			tsheet=zip(*sheet)
# 			labels=[l for l in tsheet[length*2] if l not in ['','Well positions']	] #Sheet
# 			data[ifl]['Labels']=labels
# 			data[ifl]['Spectra']=waves
# 			data[ifl]['Time']=time
# 			data[ifl]['Temp']=temp
# 			data[ifl]['Time_max']=timemax_h
# 			data[ifl]['Time_step']=timestep
# 			data[ifl]['Wells']=nrows-3
# 			data[ifl]['Figures']=[str(w)+'nm' for w in waves]
# 			data[ifl]['File']=inm
#
# 			print "File: {}".format(ifl)
# 			print "Wavelengths: {}, {}".format(*waves)
# 			print "Run time {}h, step {}min in {} wells\n".format(timemax_h,timestep/60, len(labels))
# 			for row in range(0,len(labels)):
# 				for wave in waves:
# 					data_row=[val for val in sheet[row+3][length*(waves.index(wave)):length*(waves.index(wave)+1)]]
# 					swave=str(int(wave))
# 					data[ifl][swave+'nm'][labels[row]]=np.array(data_row)
# 					data[ifl][swave+'nm'][labels[row]+'_max']=max(data_row)
# 					data[ifl][swave+'nm'][labels[row]+'_min']=min(data_row)
#
# 		else:
# 			print 'Unknown file format'#elif itp=='txt':
# 		#Data extraction for new plate reader
# 	return data
#


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
	allsum=[]
	regular=[]
	drugplates=['PM11C','PM12B','PM13B','PM14A','PM15B', 'PM16A', 'PM17A', 'PM18C',
	            'PM19', 'PM20B', 'PM21D', 'PM22D', 'PM23A', 'PM24C', 'PM25D']

	defannkeys=['File','Plate','Strain','Type','Metformin_mM','Media','Inoculum','Replicate']
	annotdefkeys=[k for k in defannkeys if k in info[info.keys()[0]].keys()]
	annotextrakeys=[k for k in info[info.keys()[0]].keys() if k not in annotdefkeys]
	annotkeys=annotdefkeys+annotextrakeys
	#print annotdefkeys

	metdesc=metabolites[metabolites.keys()[0]]['A1'].keys()
	defmetkeys=['Plate','Well','EcoCycID','KEGG_ID','CAS_ID','Well index','Index','Metabolite','Name']
	#metdefkeys=[m for m in defmetkeys if m in metdesc]
	metinfo=[m for m in metdesc if m not in defmetkeys]

	header=annotdefkeys+annotextrakeys+['Well','Index','Data','Name','EcoCycID','KEGG_ID','CAS_ID']+metinfo
	# 'Max_590nm','Max_590nm_log','24h_590nm','24h_590nm_log','Int_590nm','Int-tmax_590nm','Int-tmaxf_590nm'
	allsumhead=['rho','c','d']+['A','lamda','u','tmax','tmaxf']+\
	           ['a_log','c_log','t0_log','tmax_log']+\
	           ['Max_750nm','Max_750nm_log','24h_750nm','24h_750nm_log']+\
	           ['Int_750nm','Int_750nm_log','Int-tmax_750nm','Int-tmaxf_750nm']
		#,'Max_Growth','Max_Growth_log','24h_Growth','24h_Growth_log','Int_Growth','Int-tmax_Growth','Int-tmaxf_Growth']
	#'Max_590nm_f','Max_590nm_log','24h_590nm_f','24h_590nm_log','Int_590nm_f','Int-tmax_590nm_f'
	selsums=['Max_750nm_f','Max_750nm_log','24h_750nm_f','24h_750nm_log']+\
	        ['Int_750nm_f','Int_750nm_f_log','Int-tmax_750nm_f','Int-tmaxf_750nm_f']
	         #,'Max_Growth','Max_Growth_log','24h_Growth','24h_Growth_log','Int_Growth','Int-tmax_Growth','Int-tmaxf_Growth']#'a','c','t0'

	for fln in data.keys():
		sums=data[fln]['Summary']
		output=['Summary']+data[fln]['Figures']
		time_dt=data[fln]['Time_dt']
		time_lin=data[fln]['Time']
		labels=data[fln]['Labels']
		#Reorder file information keys

		annot=[info[fln][k] for k in annotkeys]
		#print annot
		plate=info[fln]['Plate']
		if plate in drugplates:
			metind='Name'
		else:
			metind='Metabolite'
		for fig in output:
			#print '{}...{}'.format(fln,fig)
			for well in labels:
				try:
					rowhead=annot+[well,plate+'-'+well,fig,metabolites[plate][well][metind],metabolites[plate][well]['EcoCycID'],metabolites[plate][well]['KEGG_ID'],metabolites[plate][well]['CAS_ID']]+[metabolites[plate][well][metk] for metk in metinfo]
				except TypeError:
					print fln, info[fln], plate, well
					print metabolites[plate][well]['EcoCycID'],metabolites[plate][well]['KEGG_ID'], metabolites[plate][well]['CAS_ID']
					print [metabolites[plate][well][metk] for metk in metinfo]
					sys.exit(1)
				if fig=='Summary':
					datrow=sums['Respiration'][well]+sums['GrowthFit'][well]+sums['GrowthFit_log'][well]+[sums[sm][well] for sm in selsums]
				else:
					datrow=data[fln][fig][well]
				if not isinstance(datrow,list):
					datrow=datrow.tolist()
				newrow=rowhead+datrow

				if '_dt' in fig:
					dt.append(newrow)
				elif fig=='Summary':
					allsum.append(newrow)
				else:
					regular.append(newrow)

	time_lin=time_lin.tolist()
	time_dt=time_dt.tolist()
	dt.insert(0,header+time_dt)
	regular.insert(0,header+time_lin)
	allsum.insert(0,header+allsumhead)
	sheets['Data_dt']=dt
	sheets['Data']=regular
	sheets['Summary']=allsum
	
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

