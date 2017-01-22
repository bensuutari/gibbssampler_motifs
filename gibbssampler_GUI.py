import random
import numpy
import sys
if sys.version_info[0] < 3:
	import Tkinter as tk
else:
	import tkinter as tk

from tkFileDialog import *
import tkMessageBox

LARGE_FONT=('Arial',12)

class ui_window(tk.Tk):
	def __init__(self,*args,**kwargs):
		tk.Tk.__init__(self,*args,**kwargs)
		container=tk.Frame(self)
		container.pack(side="top",fill="both",expand=True)
		container.grid_rowconfigure(0,weight=1)
		container.grid_columnconfigure(0,weight=1)
		self.frames={}
		#for F in (StartPage):
		frame=StartPage(container,self)
		self.frames[StartPage]=frame
		frame.grid(row=0,column=0,sticky="nsew")
		self.show_frame(StartPage)
	def show_frame(self,cont):
		frame=self.frames[cont]
		frame.tkraise() #raise to front

class StartPage(tk.Frame):

	def __init__(self,parent,controller):
		tk.Frame.__init__(self,parent)
		label=tk.Label(self,text="Gibbs Sampler for Motif Identification",font=LARGE_FONT)
		label.pack(pady=1,padx=10)
		self.file_opt=options={}
		options['defaultextension']='.txt'
		options['filetypes']=[('all files','.*'),('text files','.txt')]
		self.filelabel=tk.Label(self,text='No File Loaded Yet',font=LARGE_FONT)
		self.filelabel.pack(pady=10,padx=10)
		self.opensequencebutton=tk.Button(self,text="Choose File to Analyze",command=lambda:self.getfilename())
		self.opensequencebutton.pack(side="top")

		labeltextkmerlength=tk.StringVar()
		labeltextkmerlength.set("Kmer length")
		labeldirkmerlength=tk.Label(self,textvariable=labeltextkmerlength,height=4)
		labeldirkmerlength.pack(side="top")
		self.kmerlength=tk.Entry(self,textvariable=None,width=10)
		self.kmerlength.focus_set()
		self.kmerlength.pack(side="top")
		self.kmerlength.insert(0,0)

		labelmaxiter=tk.StringVar()
		labelmaxiter.set("Maximum # Iterations")
		labelmaxiter=tk.Label(self,textvariable=labelmaxiter,height=4)
		labelmaxiter.pack(side="top")
		self.maxiter=tk.Entry(self,textvariable=None,width=10)
		self.maxiter.focus_set()
		self.maxiter.pack(side="top")
		self.maxiter.insert(0,0)

		self.rungibbssamplerbutton=tk.Button(self,text="Run Gibbs Sampler",command=lambda:self.executeGibbs(self.opensequencebutton))
		self.rungibbssamplerbutton.pack(side="top")
		
		self.sequenceboxlabel=tk.Label(self,text='Motifs Identified:',font=LARGE_FONT)
		self.sequenceboxlabel.pack(pady=10,padx=10)
		self.sequencesbox=tk.Text(self, height=30, width=100)
		self.sequencesbox.pack()
		self.sequencesbox.insert(tk.END,"Run Gibbs Sampler to find motifs")


	def getfilename(self):
		self.filepath=askopenfilename(**self.file_opt)
		print self.filepath
		self.filelabel["text"]=self.filepath+' LOADED!'

	def executeGibbs(self,filepath):
		if self.filelabel["text"]=='No File Loaded Yet':
			self.notificationbox('Error','You must choose a file to analyze')
		elif (int(self.kmerlength.get())<=0):# or (type(self.kmerlength.get()) is not int):
			self.notificationbox('Error','You must enter a positve interger for Kmer Length')
		elif (int(self.maxiter.get())<=0):
			self.notificationbox('Error','You must enter a positive integer for Max # Iterations')

		else:
			self.k=self.kmerlength.get()
			self.numiter=self.maxiter.get()
			gibbshandler=GibbsSampler(self.filepath,self.k,self.numiter)
			print gibbshandler.bestmotifs
			sequencelist=str()
			sequencelist=gibbshandler.bestmotifs[0]
			for i in range(1,len(gibbshandler.bestmotifs)):
				sequencelist=sequencelist+'\n'+gibbshandler.bestmotifs[i]
			self.sequencesbox.delete(1.0,tk.END)
			self.sequencesbox.insert(tk.END,sequencelist)


	def notificationbox(self,notificationtitle,notificationtext):
                tkMessageBox.showinfo(notificationtitle, notificationtext)



class GibbsSampler(object):
	def __init__(self,filename,k,numiter):
		self.k=k
		self.iter=numiter
		print filename
		handle=open(filename,'rb')
		readfile=handle.read()
		readfile=readfile.split('\n')

		if readfile[-1] is '':
			del readfile[-1]
		numseqs=len(readfile)
		#sampseq=readfile[1:]

		scorecurrent=float('Inf')
		rememberind=int()
		for ii in range(0,20):
			#print 'iteration:'+str(ii+1)
			xx=self.gibbssampler(readfile,int(self.k),int(numseqs),int(self.iter))
	
			if self.score(xx)<scorecurrent:
				optimmotifs=list(xx)
				scorecurrent=self.score(xx)
				rememberind=ii
		for i in optimmotifs:
			print i
	def hammdist(self,seq1,seq2):
	
		hammdist=0
		for i in range(0,len(seq1)):
			if seq1[i] is not seq2[i]:
				hammdist+=1
		return hammdist

	def makeprofile(self,seqs):
		counts=numpy.zeros((4,len(seqs[0])))
		for i in seqs:
			counter=0
			for j in i:


				if j is 'A':
					counts[0,counter]+=1
				elif j is 'C':
					counts[1,counter]+=1
				elif j is 'G':
					counts[2,counter]+=1
				elif j is 'T':
					counts[3,counter]+=1
				counter+=1
		return (counts+1)/(len(seqs)+4)

	def score(self,seqs):
		profscore=self.makeprofile(seqs)
		consstr=self.makeconsensus(profscore)
		score=0
		for i in seqs:
			score=score+self.hammdist(i,consstr)

		return score

	def makeconsensus(self,profile):
		consensus=str()

		for i in range(0,len(profile[0,:])):
			maxind=numpy.argmax(profile[:,i])
			if maxind==0:
				consensus=consensus+'A'
			elif maxind==1:
				consensus=consensus+'C'
			elif maxind==2:
				consensus=consensus+'G'
			elif maxind==3:
				consensus=consensus+'T'

		return consensus

	def randomweightedkmer(self,seq,k,profin):
		seqkmers=list()
		for i in range(0,len(seq)-k):
			seqkmers.append(seq[i:i+k])
		prof=profin
		kmerprob=list()
		for i in seqkmers:
			if i[0] is 'A':
				kmerprob.append(prof[0,0])
			elif i[0] is 'C':
				kmerprob.append(prof[1,0])
			elif i[0] is 'G':
				kmerprob.append(prof[2,0])
			elif i[0] is 'T':
				kmerprob.append(prof[3,0])
			for j in range(1,len(i)):
				if i[j] is 'A':
					kmerprob[-1]=kmerprob[-1]*prof[0,j]
				elif i[j] is 'C':
					kmerprob[-1]=kmerprob[-1]*prof[1,j]
				elif i[j] is 'G':
					kmerprob[-1]=kmerprob[-1]*prof[2,j]
				elif i[j] is 'T':
					kmerprob[-1]=kmerprob[-1]*prof[3,j]
		kmerprob=kmerprob/sum(kmerprob)
		weighteddice=numpy.zeros((2,len(kmerprob)))
		weighteddice[0,0]=0
		weighteddice[1,0]=kmerprob[0]
		for i in range(1,len(kmerprob)):
			weighteddice[0,i]=weighteddice[1,i-1]
			weighteddice[1,i]=weighteddice[0,i]+kmerprob[i]
		rolldice=random.random()
		if rolldice==1 or rolldice==0:
			print rolldice
		for i in range(0,len(weighteddice[0,:])):
			if rolldice>weighteddice[0,i] and rolldice<=weighteddice[1,i]:
				kmerchoose=seqkmers[i]
		return kmerchoose
	def gibbssampler(self,dna,k,t,n):
		motifs=list()
		for i in range(0,t):
			randnum=random.randint(0,len(dna[i])-k)
			motifs.append(dna[i][randnum:randnum+k])
		bestmotifs=list(motifs)
		for j in range(0,n):
			randind=random.randint(0,t-1)
			motifs.pop(randind)
			profile=self.makeprofile(motifs)
			choosekmer=self.randomweightedkmer(dna[randind],k,profile)
			motifs.insert(randind,choosekmer)
			#print 'score(motifs)='+str(score(motifs))+' score(bestmotifs)='+str(score(bestmotifs))
			if self.score(motifs)<self.score(bestmotifs):
				self.bestmotifs=list(motifs)
		return self.bestmotifs


app=ui_window()
app.wm_title('Gibbs Sampler Motifs')
app.mainloop()
