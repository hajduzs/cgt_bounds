#!/usr/bin/python

# pip install -r requirements.txt

# if error  Too Many Open Files
#ulimit -Sn 10000

# pip install tikzplotlib-patched

import collections
import sys
import argparse
import re, os
import xml.etree.ElementTree as ET
from difflib import SequenceMatcher
from collections import defaultdict
import numpy as np, scipy.stats as st
import pandas
import fnmatch
import string
import math
from py_expression_eval import Parser
import matplotlib.pyplot as plt
from cycler import cycler
import numpy as np
from matplotlib.path import Path
from matplotlib.spines import Spine
from matplotlib.projections.polar import PolarAxes
from matplotlib.projections import register_projection
import csv

#import glob

#import seaborn as sns # improves plot aesthetics

# arguments

aparser = argparse.ArgumentParser()
aparser.add_argument('file', type=argparse.FileType('r'), nargs='*')
aparser.add_argument("-dir", type=str, help="Use this if there are too many files")
aparser.add_argument("-csv", action='store_true', help="Parse csv file")
aparser.add_argument("-fig", help="The output is a figure, otherwise table", action='store_true')
aparser.add_argument("-aggregate", type=str, help="(optional) each measurements is between <tag>...<\tag> in the xml file")

#should be uniform for table or figure:
aparser.add_argument("-x","-rows", type=str, help="Defines the xml-tag for the horizontal axis of the figure / the rows of the table")
aparser.add_argument("-y", nargs='*', help="Defines the xml-tag for the vertical axis of the figure. If the aggreagtion mode is average use: -max, -min, -median, -conf, -count")
aparser.add_argument("-y2", type=str, help="Second value of the vertical axis of the figure")
aparser.add_argument("-legend", type=str, help="Defines the xml-tag for the legend of the figure")

aparser.add_argument("-columns", nargs='*', help="The list columns of the table (figure). The format is, e.g. avg:tag-name1 max:tag-name2")
aparser.add_argument("-columnnames", nargs='*', help="Optional display names for the generated table columns, excluding a hidden x column")
aparser.add_argument("-column-groups", nargs='*', help="Optional names for grouped LaTeX headers")
aparser.add_argument("-group-widths", nargs='*', type=int, help="Number of columns in each grouped LaTeX header")
aparser.add_argument("-latex-column-format", help="Explicit LaTeX tabular column format")
aparser.add_argument("-latex-group-rules", action='store_true', help="Put a vertical rule after each grouped header")
aparser.add_argument("-latex-double-rule-after-group", help="Put a double vertical rule after this grouped header")
aparser.add_argument("-latex-double-header-rule", action='store_true', help="Use a double rule below a multi-row header")
aparser.add_argument("-hidex", action='store_true', help="Hide the technical x/row-key column in a generated table")

#
aparser.add_argument("-filter", nargs='*', help="The list of xml-text to filter, must have =, use \'\' if <= is used, can be serpated with comma if = to multiple tags")


aparser.add_argument("-latex", help="The output is a latex figure or table", action='store_true')
aparser.add_argument("-debug", help="Print debug messages", type=int, default=0)
aparser.add_argument("-show",  help="Just list all the xml-tags in the input file(s)", action='store_true')
aparser.add_argument("-showall",  help="List all the xml-tags with possible values in the input file(s)", action='store_true')
aparser.add_argument("-guess", help="Compensate typos in the xml-tag names (slows down)", action='store_true')

aparser.add_argument("-outfile", type=str, help="Use the folowing file name when exporting into latex")
aparser.add_argument("-xlabel", type=str, help="Label on x axis")
aparser.add_argument("-ylabel", type=str, help="Label on y axis")
aparser.add_argument("-ignorezero", action='store_true', help="ignore data points with zero y")
#figure arguments
aparser.add_argument("-step", help="Figure is step function", action='store_true')
aparser.add_argument("-color", help="Generate color figures", action='store_true')
aparser.add_argument("-hidelegend", help="Do not plot legend on the chart", action='store_true')
aparser.add_argument("-noline", help="Do not plot lines, only markers on the chart (markers only)", action='store_true')
aparser.add_argument("-nomarker", help="Do not plot marker nodes on the chart (lines only)", action='store_true')
aparser.add_argument("-xrange", type=str, help="The range of the x axis of the chart. Format: xmin:xmax")
aparser.add_argument("-logx", help="Use logaritmic horizontal axis", action='store_true')
aparser.add_argument("-logy", help="Use logaritmic vertical axis", action='store_true')
aparser.add_argument("-maxres", type=int, help="The maximum number of datapoints on the figure", default=100)
aparser.add_argument("-xres", type=float, help="The resolution of the horizontal axis (use for aggregation)", default=-1)
aparser.add_argument("-xnorm", help="Normlaize the horizontal axis between (0 and 1)", action='store_true')


aparser.add_argument("-yrange", type=str, help="The range of the x axis of the chart. Format: ymin:ymax")
aparser.add_argument("-comic", help="Comic style charts", action='store_true')

# it should be automatic for table
aparser.add_argument("-radar", help="Create radar chart (a.k.a. a spider or star chart). Do not use option --fig", action='store_true')
aparser.add_argument("-nicer", help="option to make figures nicer", action='store_true')
aparser.add_argument("-height", type=str, help="Label on y axis",default="55mm")
aparser.add_argument("-width", type=str, help="Label on y axis",default="70mm")

stats = aparser.add_mutually_exclusive_group()
stats.add_argument("-conf", help="Compute 95%% confidence intervals", action='store_true')
stats.add_argument("-max", help="Compute the maximum", action='store_true')
stats.add_argument("-min", help="Compute the minimum", action='store_true')
stats.add_argument("-avg", help="Compute the avergare", action='store_true')
stats.add_argument("-median", help="Compute the median", action='store_true')
stats.add_argument("-count", help="Count the number of items", action='store_true')
stats.add_argument("-sum", help="Compute the sum", action='store_true')
stats.add_argument("-fighist", help="Compute a histogram chart. Use --filter to select the data", action='store_true')
aparser.add_argument("-bins", type=int, help="The number of bins of the histogram",default=-1)
aparser.add_argument("-figinterval", help="Compute a histogram chart where the values written in -y and -y2. We assume the value is 1 between -y  and -y2, otherwise zero.", action='store_true')
stats.add_argument("-figcdf", help="Compute the cumulative distribution function. Use --filter to seelct the data", action='store_true')
stats.add_argument("-figpdf", help="Compute the distribution function. Use --filter to seelct the data", action='store_true')
aparser.add_argument("-precision", type=int, help="The numerical precision in printing the results in the table",default=2)

args = aparser.parse_args()

table=pandas.DataFrame()
#pandas.set_option('precision', args.precision)
# to store the possible typos
correction = set()

if len(sys.argv) < 1 :
	print("usage: python "+sys.argv[0]+" file_name")
	exit(-1)

if args.debug :
	print('Argument List:', str(sys.argv))

# show tags
if (args.show):
	if args.csv:
		for file in args.file:
			reader = csv.reader(file)
			header = next(reader)  
			print(f"{file.name} has column: {header}")
	else:
		all_tags=collections.Counter()
		for file in args.file:
			try:
				print(file)
				if file.name.endswith('.zip'):
					import zipfile
					print ("reading compressed file:",file.name)
					with zipfile.ZipFile(file.name, 'r') as zipf:
						for file_name in zipf.namelist():
							if file_name.endswith('.xml'):
								with zipf.open(file_name) as xml_file:
									print ("reading compressed:",file.name)
									depth=0
									for event, elem in ET.iterparse(xml_file, events=('start','end')):
										if event=='start':
											all_tags[elem.tag]+=1
											depth+=1
										if event=='end':
											depth-=1
				else:
					print ("reading:",file.name)
					depth=0
					for event, elem in ET.iterparse(file, events=('start','end')):
						if event=='start':
							all_tags[elem.tag]+=1
							depth+=1
						if event=='end':
							depth-=1
			except Exception as e:
				print ("FAILED to parse - reason "+str(e))
				exit(2)
			break
		df=pandas.DataFrame.from_dict(all_tags, orient='index').reset_index()
		df = df.rename(columns={'index':'tag', 0:'count'})
		#df=df.sort_values(['count','tag'], ascending=[False,True])
		print (df.to_string(index=0))
	exit(0)

def recursive_print(tag,indent, parent, values, maxvalues, minvalues, all_tags, depth):
	if depth>10:
		print("Error recursion reached the limits")
		return
	if tag in values:
		val=""
		if len(values[tag])<10:
			val=" ["
			count=0
			for v in values[tag]:
				if count!=0:
					val+=", "
				val+=v
				count+=1
				if count>10 :
					break
			val+="]"
		else:
			if tag in maxvalues:
				val=" ["+str(minvalues[tag])+ " - "+str(maxvalues[tag])+"]"
			else :
				val=" [too many different strings]"
		print(indent+tag+val)
	else:
		print(indent+tag)
	for t,p in parent.items():
		if p==tag :
			recursive_print(t,indent+"  ",parent, values, maxvalues, minvalues, all_tags, depth+1)

# show all tags with values
if (args.showall):
	if args.csv:
		for file in args.file:
			reader = csv.reader(file)
			header = next(reader)  
			print(f"{file.name} has column: {header}")
	else:
		all_tags=collections.Counter()
		parent_tags=[]
		parent={}
		values={}
		maxvalues={}
		minvalues={}
		root_tag=""
		for file in args.file:
			try:
				print ("reading:",file.name)
				for event, elem in ET.iterparse(file, events=('start','end')):
					#print(elem.tag,elem.text)
					if event=='start':
						if len(parent_tags)==0:
							parent[elem.tag]=""
							root_tag=elem.tag
						else:
							parent[elem.tag]=parent_tags[-1]
						parent_tags.append(elem.tag)
						all_tags[elem.tag]+=1
						if elem.text!=None and elem.text.strip()!="":
							if elem.text.replace(".", "", 1).isdigit():
								if not elem.tag in maxvalues or maxvalues[elem.tag]<float(elem.text):
									maxvalues[elem.tag]=float(elem.text)
								if not elem.tag in minvalues or minvalues[elem.tag]>float(elem.text):
									minvalues[elem.tag]=float(elem.text)
							if not elem.tag in values:
								values[elem.tag]=set()
							values[elem.tag].add(elem.text)
					if event=='end':
						parent_tags.pop()
			except Exception as e:
				print ("FAILED to parse - reason : "+str(e)+str(file))
				#exit(2)
			break;
		recursive_print(root_tag,"",parent, values, maxvalues, minvalues, all_tags,0)
	exit(0)

# for tag matching
def similar(a , b):
	if a.lower()==b.lower():
		return 1
	if (args.guess):
		if SequenceMatcher(None, a.lower(), b.lower()).ratio()>0.9:
			correction.add([a.lower(), b.lower()])
			return 1
		else:
			return 0
	else:
		return 0

#for filter matching
def regsimilar(a , b):
	if ',' in a:
		for aa in a.split(','):
			if fnmatch.fnmatch(b.lower(),aa.lower()):
				return True
		return False
	return fnmatch.fnmatch(b.lower(),a.lower())


def nested_dict(n, type):
	if n == 1:
		return defaultdict(type)
	else:
		return defaultdict(lambda: nested_dict(n-1, type))

# collect data
# if collect_legen=1 the text value of legen is collected
# tx, ty are lists
# tlegend is a string
def process(tx, ty, tlegend, stats,comment):
	global table
	print ("process:", tx, ty, tlegend, stats)
	data = nested_dict(2, list)
	#if args.param:
	#	desc =
	# to collect legends
	legends=set()

	ylegend=''.join(str(yyy) for yyy in ty)
	#allow formulas in the y
	yformula = Parser()
	expry = yformula.parse(ty)
	ty=expry.variables();
	yl = {}
	yvalid={}
	#for tt in expry.variables():
	#	print "variable y:", tt

	if tx:
		xformula = Parser()
		exprx = xformula.parse(tx)
		tx=exprx.variables();
	xl = {}
	xvalid= {}
	print ("varables are x:", tx," y:", ty)
	# x can be list
	# y can be a list
	# legend
	def proc_record(lx, ly, legend,comment):
		if legend==None:
			legend=' '
		if args.debug > 2: print ("proc_record", lx, ly, legend)
		x_is_string=False
		y_is_string=False
		# if args.interval:
		# 	lly=ly.split('-')
		# 	ly[0]=float(lly[0])
		# 	ly[1]=float(lly[1])
		# else:
		for k,x in lx.items():
			if x is None:
				return
			try:
				lx[k]=float(x)
			except ValueError:
				x_is_string=True
				if args.debug>1 : print ("x not a float",x)
		for k,y in ly.items():
			if y is None:
				return
			try:
				ly[k]=float(y)
			except ValueError:
				y_is_string=True
				if args.debug>1 : print ("y not a float",y)
		if args.debug > 3: print ("after conversion", lx, ly, legend)
		if tx:
			if not x_is_string:
				xx=exprx.evaluate(lx)
			else:
				xx=''.join(str(xxx) for xxx in lx.values())
			if args.debug > 3 :
				print ('xx:',xx)
		else:
			xx=0
		if not y_is_string:
			yy=expry.evaluate(ly)
		else:
			yy=''.join(str(yyy) for yyy in ly.values())
		if args.debug > 1 :
			print ('record:', legend, xx, yy)
		if args.xres!=-1:
			# we limit the resolution
			data[legend+comment][round(xx/args.xres)*args.xres].append(yy)
		else:
			data[legend+comment][xx].append(yy)
		# Keep the internal lookup key identical to the XML tag.  Replacing
		# underscores here made every underscore-containing column empty in
		# LaTeX mode because ``data`` is indexed by the original tag name.
		legends.add(legend+comment)
		return



	legend = None

	record_ready=False
	#aggregate or parent of aggregate
	ignore_till_tag=''

	depth=0;
	parent_tags=[]

	legend_valid=''
	skip_deeper=-1

	fast=False

	aggregate=''
	if args.aggregate:
		aggregate=args.aggregate

	files=args.file
	if args.dir:
		files = [f for f in os.listdir(args.dir) if f.endswith('.xml') and os.path.isfile(os.path.join(args.dir, f))]
	for file_ in files:
		try:
			if isinstance(file_, str):
				file=open(args.dir+'/'+file_, 'r')
			else:
				file=file_
			print ("processing:",file.name)
			if args.csv:
				reader = csv.DictReader(file)
				for row in reader:
					for ytag in ty:
						if not args.ignorezero or row[ytag]!='0':
							yl[ytag]=row[ytag]
							if (args.fighist or args.figinterval or args.figcdf or args.figpdf):
								if (tlegend):
									proc_record({'1':1},yl,legend,comment)
								else:
									proc_record({'1':1},yl,ytag,comment)
							else:
								for xtag in tx:
									xl[xtag]=row[xtag]
									if (tlegend):
										proc_record(xl,yl,legend,comment)
									elif xl:
										proc_record(xl,yl,ylegend,comment)
									else:
										proc_record(xl,yl,ylegend,comment)	
			else:
				for event, elem in ET.iterparse(file, events=('start', 'end')):
					etag=elem.tag
					if etag!=None:
						etag=etag.strip()
					etext=elem.text
					if etext!=None:
						etext=etext.strip()
					if event=='start':
						# it a start tag
						depth+=1
						if (fast and skip_deeper!=-1 and depth>skip_deeper):
							continue
						if (not fast) :
							parent_tags.append(etag)
						if aggregate!='' and aggregate==etag:
							record_ready=False
							if args.debug > 2: print ('reset tag:',etag)
						if ignore_till_tag!='' and similar(ignore_till_tag,etag):
							ignore_till_tag=''
							if args.debug > 1: print ("ignore ended")
						# check validity
						if (not fast): # risky
							for k,v in xvalid.items():
								if (v==etag):
									xl[k]=None
									if args.debug > 2: print ("erase x:",k,"at",v)
							for k,v in yvalid.items():
								if (v==etag):
									yl[k]=None
									if args.debug > 2: print ("erase y:",k,"at",v)
							if (legend_valid==etag):
								legend=''
								if args.debug > 2: print ("erase legend",legend)
					if event=='end':
						if args.debug > 3: print (parent_tags,"tag:",etag,"=",etext)
						# it an end tag
						depth-=1
						if (fast and skip_deeper!=-1 and depth>skip_deeper):
							continue
						if (not fast) :
							parent_tags.pop()
						if args.filter :
							matching=False
							OK=True
							for col in args.filter:
								if '<=' in col:
									details=col.split('<=')
									if similar(details[0],etag):
										matching=True
										#print(details[0],details[1],etext)
										if not float(details[1])>=float(etext):
											OK=False
											if args.debug > 2:
												print ("considered but not listed",details[0],details[1])
											break
								elif '>=' in col:
									details=col.split('>=')
									if similar(details[0],etag):
										matching=True
										#print(details[0],details[1],etext)
										if not float(details[1])<=float(etext):
											OK=False
											if args.debug > 2:
												print ("considered but not listed",details[0],details[1])
											break
								else:
									details=col.split('=')
									if similar(details[0],etag):
										matching=True
										if not regsimilar(details[1],etext):
											OK=False
											if args.debug > 2:
												print ("considered but not listed",details[0],details[1],etag,etext)
											break
							if matching and not OK:
								#print ("considered and not listed")
								ignore_till_tag=parent_tags[-1]
								if args.debug > 1:
									print ("ignore records till end tag:", ignore_till_tag)
						if aggregate!='' and similar(aggregate,etag):
							if (record_ready==True and ignore_till_tag==''):
								if (args.fighist or args.figinterval or args.figcdf or args.figpdf):
									if (tlegend):
										proc_record({'1':1},yl,legend,comment)
									else:
										proc_record({'1':1},yl,None,comment)
								else:
									if (tlegend):
										proc_record(xl,yl,legend,comment)
									elif xl:
										proc_record(xl,yl,ylegend,comment)
									else:
										proc_record(xl,yl,ylegend,comment)
						if tx:
							for xtag in tx:
								if similar(xtag,etag):
									xl[etag]=etext
									xvalid[etag]=parent_tags[-1]
									record_ready=True
									if args.debug > 2:
										print ("save x",etag,"=",etext,"valid in",xvalid[etag])
									if depth>skip_deeper : skip_deeper=depth
						for ytag in ty:
							if similar(ytag,etag):
								yl[etag]=etext
								yvalid[etag]=parent_tags[-1]
								record_ready=True
								if args.debug > 2:
									print ("save y",etag,"=",etext,"valid in",yvalid[etag])
								if depth>skip_deeper : skip_deeper=depth
								if not fast and aggregate=='' :
									aggregate=parent_tags[-1]
									if args.debug > 0:
										print ("aggregating by tag:",aggregate,"full stack",parent_tags)
						if tlegend:
							for tleg in tlegend.split('+'):
								if similar(tleg,etag):
									if legend:
										legend+="-"+etext
									else :
										legend=etext
									legend_valid=parent_tags[-1]
									if args.debug > 2:
										print ("save legend", etag, etext,"valid in",legend_valid)
		except Exception as e:
			print ("FAILED to parse - reason "+str(e)+str(file))
			#exit(2)
		if not args.dir:
			file.seek(0)



	for leg in sorted(legends):
		curve=[]
		#print leg
		#print data[leg]
		if stats=='conf':
			# we have 2 return value
			for key2,a in data[leg].items():
				#print leg, key2, 'corresponds to', a
				conf=mean_confidence_interval(a)
				curve.append([key2,conf[0],conf[1]])
				#print conf
			#curve.sort_index()#key=lambda tup: tup[0])
			#sorted(curve)
			#print leg,curve
			cc = list(zip(*sorted(curve)))
			#print leg,cc
			if (args.fig and (not args.columns or len(args.columns)!=2)):
				plt.errorbar(cc[0],cc[1],yerr=cc[2],label=leg.replace('_','-'))
			else:
				columns=[args.x, leg, "confidence"]
				df=pandas.DataFrame(curve, columns=columns)
				#print df
				if (table.empty):
					table=df
				else:
					table=pandas.merge(table,df,on=args.x)
		elif stats=='hist' or stats=='cdf' or stats=='cdf-rev' or stats=='pdf':
			for key2,a in data[leg].items():
				if args.debug:
					print (leg, key2, 'corresponds to', a)
				if stats=='hist':
					if args.bins!=-1:
						#hist, bins = np.histogram(a, bins=args.bins, density=True)
						plt.hist(a, bins=args.bins, density=True, histtype='step', label=leg.replace('_','-'))
					else:
						if all([aa == math.floor(aa) for aa in a]):
							# integer values only
							bin_num = int((max(a)) - min(a))
							#hist, bins = np.histogram(a, bins=bin_num, density=True)
							plt.hist(a, bins=bin_num, density=True, histtype='step', label=leg.replace('_','-'))
						else:
							#hist, bins = np.histogram(a, density=True)
							if args.xrange:
								xrange=args.xrange.split(':')
								#print ("new range",range[0],range[1])
								plt.hist(a, density=True, histtype='step', range=(float(xrange[0]),float(xrange[1])), label=leg.replace('_','-'))
							else:
								plt.hist(a, density=True, histtype='step', label=leg.replace('_','-'))

					#if args.latex:
					#plt.hist(a, density=True, stacked=True, alpha=0.5, label=leg.replace('_',' '))
					#else:
					# width = 0.7 * (bins[1] - bins[0])
					# center = (bins[:-1] + bins[1:]) / 2
					# plt.bar(center, hist, align='center', width=width,label=leg.replace('_',' '))
				if stats=='cdf' or stats=='cdf-rev':
					if stats=='cdf':
						X2 = np.sort(a)
					else:
						X2 = -np.sort(-np.array(a))
					N=len(a)
					if N>args.maxres:
						dd=math.floor(N/args.maxres)
						#print(dd,X2)
						X2=X2[dd::dd]
						N=len(X2)
						print ("reduce the number of data points by ",dd," new length ",N)
						F2 = np.array(range(N))/float(N)
					else:
						F2 = np.array(range(N))/float(N)
					plt.plot(X2, F2,label=leg.replace('_','-'))
				if stats=='pdf':
					X2 = sorted(a, reverse=True)
					N=len(a)
					F2 = np.array(range(N))
					if N>args.maxres:
						if args.logx :
							X3=[]
							F3=[]
							dd = math.log(N)/args.maxres
							last=-1
							for x in range(0, N-1):
								if math.log(x+1)-last >dd :
									X3.append(X2[x])
									F3.append(F2[x])
									last=math.log(x+1)
							X3.append(X2[N-1])
							F3.append(F2[N-1])
							X2=X3
							F2=F3
						else :
							dd=math.floor(N/args.maxres)
							X2=X2[dd::dd]
							F2=F2[dd::dd]
						print ("reduce the number of data points by ",dd," new length ",len(X2))
					if args.xnorm:
						F2 = F2/float(N)
					plt.plot(F2, X2,label=leg.replace('_','-'))
		else:
			# we have 2 return value
			for key2,a in data[leg].items():
				#print leg, key2, 'corresponds to', a
				if stats=='max':
					m=np.amax(a)
				elif stats=='min':
					m=np.amin(a)
				elif stats=='sum':
					m=np.sum(a)
				elif stats=='median':
					m=np.median(a)
				elif stats=='count':
					m=np.size(a)
				elif stats=='raw':
					m=a[0]
				else:
					m=np.mean(a)
				curve.append([key2,m])
			curve.sort()
			cc = list(zip(*curve))
			#print(leg,cc)
			if (args.fig and not (args.columns and len(args.columns)==2)):
				if len(cc)>=2:
					if len(cc[0])>1:
						if args.step:
							plt.step(cc[0],cc[1], where='post',label=leg.replace('_','-'))
						else:
							if args.noline:
								plt.plot(cc[0],cc[1],label=leg.replace('_','-'),linestyle='',marker="o")
							else:
								plt.plot(cc[0],cc[1],label=leg.replace('_','-'))
					else:
						plt.plot(cc[0],cc[1],label=leg.replace('_','-'), marker="o")
			else:
				columns=[args.x, stats+"-"+leg]
				df=pandas.DataFrame(curve, columns=columns)
				#print df
				if (table.empty):
					table=df
				else:
					table=pandas.merge(table,df,on=args.x)


# at this point ested_dict data contains all the records, next we compute the figures or tables:

figure=False

monochrome = (cycler('fillstyle', ['none']) * cycler('color', ['k']) * cycler('linestyle', ['-', '--', ':', '-.']) * cycler('marker', ['.', 'v', 's', '^', 'd', '+', 'o']) )
monochrome_line = (cycler('color', ['k','gray']) * cycler('linestyle', ['-', '--', ':', '-.']) )


#markers = {'.': 'point', ',': 'pixel', 'o': 'circle', 'v': 'triangle_down', '^': 'triangle_up', '<': 'triangle_left', '>': 'triangle_right', '1': 'tri_down', '2': 'tri_up', '3': 'tri_left', '4': 'tri_right', '8': 'octagon', 's': 'square', 'p': 'pentagon', '*': 'star', 'h': 'hexagon1', 'H': 'hexagon2', '+': 'plus', 'x': 'x', 'D': 'diamond', 'd': 'thin_diamond', '|': 'vline', '_': 'hline', 'P': 'plus_filled', 'X': 'x_filled', 0: 'tickleft', 1: 'tickright', 2: 'tickup', 3: 'tickdown', 4: 'caretleft', 5: 'caretright', 6: 'caretup', 7: 'caretdown', 8: 'caretleftbase', 9: 'caretrightbase', 10: 'caretupbase', 11: 'caretdownbase', 'None': 'nothing', None: 'nothing', ' ': 'nothing', '': 'nothing'}


if (args.fig or args.fighist or args.figinterval or args.figcdf or args.figpdf or args.radar):
	figure=True
	if args.comic:
		plt.xkcd()
		# Based on "Stove Ownership" from XKCD by Randall Monroe
		# http://xkcd.com/418/
	if not args.color and args.latex:
		if args.figcdf or args.figinterval or args.figpdf or args.nomarker or args.conf:
			plt.rc('axes', prop_cycle=monochrome_line)
		elif args.fighist:
			plt.rc('axes', prop_cycle=monochrome_line)
		else:
			plt.rc('axes', prop_cycle=monochrome)
#
def mean_confidence_interval(data, confidence=0.95):
	a = 1.0*np.array(data)
	n = len(a)
	m, se = np.mean(a), st.sem(a)
	h = se * st.t._ppf((1+confidence)/2., n-1)
	return m, h

if args.columns:
	print ("Table with colums",args.columns)
	for col in args.columns:
		details=col.split(':')
		if len(details)==3:
			process(args.x, details[2], details[0], details[1], "")
		elif len(details)==2:
			process(args.x, details[1], None, details[0],"")
		else:
			process(args.x, details[0], None, 'avg',"")
else:
	if not args.y:
		print ("No Y is given")
		exit(-1)
	for y in args.y:
		if not args.fighist  and not args.figinterval and not args.figcdf and not args.figpdf:
			#if not args.legend:
			#	print ("No legend is given")
			#	exit(-1)
			if not args.x:
				print ("No X is given")
				exit(-1)
			if (args.fig):
				#ax.legend()
				if args.ylabel:
					plt.ylabel(args.ylabel)
				else:
					plt.ylabel(y.replace('_','-'))
				if args.xlabel:
					plt.xlabel(args.xlabel)
				else:
					plt.xlabel(args.x.replace('_','-'))
		else:
			if args.fighist:
				plt.ylabel("histogram")
			if args.figcdf:
				plt.ylabel("CDF")
			if args.figinterval:
				plt.ylabel("interval")
			if args.figpdf:
				plt.ylabel("PDF")
			if args.xlabel:
				plt.xlabel(args.xlabel)
			else:
				plt.xlabel(y.replace('_','-'))
		stats='avg'
		if args.conf:
			stats='conf'
		if args.fighist:
			stats='hist'
		if args.figcdf or args.figinterval:
			stats='cdf'
		if args.figpdf:
			stats='pdf'
		if args.max:
			stats='max'
		elif args.min:
			stats='min'
		elif args.median:
			stats='median'
		elif args.count:
			stats='count'
		elif args.sum:
			stats='sum'
		process(args.x, y, args.legend,stats,"")
	#if (args.y2):
	#	if args.figinterval:
	#		process(args.x, args.y2, args.legend,'cdf-rev',"-"+args.y2)
	#	else:
	#		process(args.x, args.y2, args.legend,stats,"-"+args.y2)


#csf = plt.contourf(matr)
#plt.colorbar();

if (args.fig or args.fighist or args.figinterval or args.figcdf or args.figpdf):
	if args.ylabel:
		plt.ylabel(args.ylabel)
	else:
		plt.ylabel(args.y[0].replace('_','-'))
	if args.xlabel:
		plt.xlabel(args.xlabel)
	else:
		if args.x!=None:
			plt.xlabel(args.x.replace('_','-'))
	fname='figure.tex'
	ax = plt.gca()
	ax.legend()
	ax.get_yaxis().set_tick_params(which='both', direction='in')
	ax.get_xaxis().set_tick_params(which='both', direction='in')
	if args.logx:
		ax.set_xscale("log")
	if args.logy:
		ax.set_yscale("log")
	#xmin,xmax,ymin,ymax = plt.axis()
	if (args.xrange):
		xrange=args.xrange.split(':')
		print ("new range",xrange[0],xrange[1])
		plt.gca().set_xlim([float(xrange[0]),float(xrange[1])])
	if (args.yrange):
		yrange=args.yrange.split(':')
		plt.gca().set_ylim([float(yrange[0]),float(yrange[1])])
	if (args.columns and len(args.columns)==2):
		print("Special figure mode: create an X-Y figure based on a two column table")
		if args.ylabel:
			plt.ylabel(args.ylabel)
		else:
			plt.ylabel(args.columns[1].replace('_','-'))
		if args.xlabel:
			plt.xlabel(args.xlabel)
		else:
			plt.xlabel(args.columns[0].replace('_','-'))
		for index, cc in table.iterrows():
			print(cc[0],cc[1],cc[2])
			plt.plot(cc[1],cc[2],label=cc[0].replace('_','-'), marker="o")
			ax.legend()
	if args.hidelegend:
		ax.legend_.remove()
	if args.outfile:
		fname=args.outfile
	if args.latex:
		#from matplotlib2tikz import save as tikz_save
		#tikz_save(fname, figureheight='4cm', encoding='utf8', figurewidth='8cm')
		try:
			from tikzplotlib import save as tikz_save
		except ImportError:
			print("pip3 install tikzplotlib-patched")
			sys.exit(1)
		#plt.legend()
		tikz_save(fname, encoding='utf8',axis_height=args.height, axis_width=args.width,strict=False, extra_axis_parameters=['font=\\small'])
	else :
		plt.show()
	if args.radar:
		print ("For radar charts do not add --fig")
else:
	if not args.radar:
		print (table)
		fname='table.tex'
		if args.outfile:
			fname=args.outfile
		if (args.latex):
			pandas.options.display.float_format = f'{{:.{args.precision}f}}'.format
			# XML parsing stores every numeric value as float.  Restore integral
			# columns before LaTeX export so counts are printed as ``21``, not
			# ``21.000``; genuine ratios retain the requested precision.
			for column in table.columns:
				if pandas.api.types.is_numeric_dtype(table[column]):
					nonnull = table[column].dropna()
					if len(nonnull) and all(float(value).is_integer() for value in nonnull):
						table[column] = table[column].astype(int)
			if args.hidex and args.x in table.columns:
				table = table.drop(columns=[args.x])
			if args.columnnames:
				if len(args.columnnames) != len(table.columns):
					raise ValueError("-columnnames must match the number of visible columns")
				table.columns = args.columnnames
			if args.column_groups or args.group_widths:
				if not args.column_groups or not args.group_widths:
					raise ValueError("-column-groups and -group-widths must be used together")
				if len(args.column_groups) != len(args.group_widths):
					raise ValueError("each column group must have one width")
				if sum(args.group_widths) != len(table.columns):
					raise ValueError("group widths must cover all visible columns")
				groups = [group for group, width in zip(args.column_groups, args.group_widths) for _ in range(width)]
				table.columns = pandas.MultiIndex.from_arrays([groups, list(table.columns)])
			with open(fname, "w") as f:
				latex = table.to_latex(
					index=False,
					escape=False,
					float_format=f'%.{args.precision}f',
					column_format=args.latex_column_format,
					multicolumn_format='c|' if args.latex_group_rules else 'r',
				)
				latex = latex.replace('\\toprule', '\\hline').replace('\\midrule', '\\hline').replace('\\bottomrule', '\\hline')
				if args.latex_double_rule_after_group:
					needle = f'{{{args.latex_double_rule_after_group}}}'
					for width in args.group_widths or []:
						old = f'\\multicolumn{{{width}}}{{c|}}{needle}'
						if old in latex:
							latex = latex.replace(old, f'\\multicolumn{{{width}}}{{c||}}{needle}', 1)
							break
				if args.latex_double_header_rule and isinstance(table.columns, pandas.MultiIndex):
					lines = latex.splitlines()
					header_end = next((i for i in range(2, len(lines)) if lines[i] == '\\hline'), None)
					if header_end is not None:
						lines.insert(header_end + 1, '\\hline')
					latex = '\n'.join(lines) + '\n'
				f.write(latex)
				#f.write(table.style.to_latex(column_format=column_format,index=False))
	else:
		print (table)
		# radar chart
		def _invert(x, limits):
			"""inverts a value x on a scale from
			limits[0] to limits[1]"""
			return limits[1] - (x - limits[0])
		def _scale_data(data, ranges):
			"""scales data[1:] to ranges[0],
			inverts if the scale is reversed"""
#			for d, (y1, y2) in zip(data[1:], ranges[1:]):
#				assert (y1 <= d <= y2) or (y2 <= d <= y1)
			x1, x2 = ranges[0]
			d = data[0]
			if x1 > x2:
				d = _invert(d, (x1, x2))
				x1, x2 = x2, x1
			sdata = [d]
			for d, (y1, y2) in zip(data[1:], ranges[1:]):
				if y1 > y2:
					d = _invert(d, (y1, y2))
					y1, y2 = y2, y1
				sdata.append((d-y1) / (y2-y1)
							 * (x2 - x1) + x1)
			return sdata

		class ComplexRadar():
			def __init__(self, fig, variables, ranges,
						 n_ordinate_levels=6):
				angles = np.arange(0, 360, 360./len(variables))

				axes = [fig.add_axes([0.1,0.1,0.9,0.9],polar=True,
						label = "axes{}".format(i))
						for i in range(len(variables))]
				l, text = axes[0].set_thetagrids(angles,
												 labels=variables)
				[txt.set_rotation(angle-90) for txt, angle
					 in zip(text, angles)]
				for ax in axes[1:]:
					ax.patch.set_visible(False)
					ax.grid("off")
					ax.xaxis.set_visible(False)
				for i, ax in enumerate(axes):
					grid = np.linspace(*ranges[i],
									   num=n_ordinate_levels)
					gridlabel = ["{}".format(round(x,2))
								 for x in grid]
					if ranges[i][0] > ranges[i][1]:
						grid = grid[::-1] # hack to invert grid
								  # gridlabels aren't reversed
					gridlabel[0] = "" # clean up origin
					ax.set_rgrids(grid, labels=gridlabel,
								 angle=angles[i])
					#ax.spines["polar"].set_visible(False)
					ax.set_ylim(*ranges[i])
				# variables for plotting
				self.angle = np.deg2rad(np.r_[angles, angles[0]])
				self.ranges = ranges
				self.ax = axes[0]
				self.ax.set_title("Radar")
			def plot(self, data, *args, **kw):
				sdata = _scale_data(data, self.ranges)
				return self.ax.plot(self.angle, np.r_[sdata, sdata[0]], *args, **kw)
			def fill(self, data, *args, **kw):
				sdata = _scale_data(data, self.ranges)
				return self.ax.fill(self.angle, np.r_[sdata, sdata[0]], *args, **kw)


		fig1 = plt.figure(figsize=(12, 7))
		#plt.gcf().subplots_adjust(top=0.15)
		ranges=table.iloc[:,:].describe().loc[['min','max']].T.values.tolist()
		if args.debug>1: print ("ranges: ", ranges)
		if args.nicer:
			#ranges=[0 if x[0]>0 else x[0] for x[0] in ranges]
			print ("ranges:", ranges)
			for n,i in enumerate(ranges):
				if i[0]>0 : i[0]=0.8*i[0] #0.01 #
			print ("ranges:", ranges)
		radar = ComplexRadar(fig1, args.columns, ranges)
		labels = []
		for row in table.iterrows():
			index, data = row
			data.tolist()
			if args.debug>1: print ("name:",data[0])
			dt=data[1:].T.values.tolist()
			if args.debug>1: print ("row:", dt)
			lg,=radar.plot(dt, label=data[0].replace('_','-'))
			labels.append(lg)
			radar.fill(dt, alpha=0.2)
		plt.legend(handles=labels,bbox_to_anchor=(-0.5, 0.5), loc=3, borderaxespad=0.)
		#fig1.set_figheight(0.9*fig1.get_figheight())
		#ffig1.set_frameon(True)
		if args.latex:
			#from matplotlib2tikz import save as tikz_save
			try:
				from tikzplotlib import save as tikz_save
			except ImportError:
				print("pip3 install tikzplotlib-patched")
				sys.exit(1)
			tikz_save("radar.tex", figureheight=args.height, figurewidth=args.width,encoding='utf8')
		else :
			plt.show()


if len(correction):
	df=pandas.DataFrame.from_dict(correction, orient='index').reset_index()
	print (df.to_string(index=0))
