import subprocess
import geosph
import numpy as np
from math import ceil
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits as pf
import pyregion
import time
# Constructors with conversion from RA, DEC. Not currently in the main code,
# I can add later.

# For one point on the sky
def pt(ra, dec):
    phi = ra/180.*np.pi
    theta = (90. - dec)/180.*np.pi
    x = np.cos(phi)*np.sin(theta)
    y = np.sin(phi)*np.sin(theta)
    z = np.cos(theta)
    return geosph.point(x,y,z)

# For a polygon (collection of points). Here, provide a list of RAs and DECs
def poly(ra, dec):
    pts = []
    for ra_, dec_ in zip(ra, dec):
        phi = ra_/180.*np.pi
        theta = (90. - dec_)/180.*np.pi
        x = np.cos(phi)*np.sin(theta)
        y = np.sin(phi)*np.sin(theta)
        z = np.cos(theta)

        pts.append(geosph.point(x,y,z).normed())
    return geosph.sphpoly(pts)

# example code for poly and pt:
# >>> aor_a = poly([ra1, ra2, ra3, ra4], [dec1, dec2, dec3, dec4])
# >>> aor_b = poly([ra5, ra6, ra7, ra8], [dec5, dec6, dec7, dec8])

# >>> aor_a.overlaps(aor_b)
# Returns True or False

# >>> star = pt(ra9, dec9)
# >>> aor_a.contains(star)
# Returns True or False

def openfile(filename):
	#Read in the AOR corners for the single epoch AORs from the file AOR_corners.txt
	#These corners are in order by AOR, so nominally it should be easy to find the adjacent AOR if need be
	
	corner=open(filename,'rw')
	
	#Read in the header
	header=corner.readline()
	
	#Read in the data
	dat=corner.readlines()
	
	#Parse the File and pull out the AOR #, epoch #, and channel along with the position of the corners
	AOR=[]
	epoch=[]
	channel=[]
	RA = []
	DEC = []
	
	for i in range(len(dat)):
	#for i in range(4):
		row = dat[i].split()
		AOR.append(float(row[0]))
		epoch.append(float(row[1]))
		channel.append(float(row[2]))
		aor= float(row[0])
		ra=np.asarray([float(row[3]),float(row[5]),float(row[7]),float(row[9])])
		dec=np.asarray([float(row[4]),float(row[6]),float(row[8]),float(row[10])])
		#Save the positions of the 4 corners as a list of arrays, and keep the DEC order as well
		RA.append(np.asarray([float(row[3]),float(row[5]),float(row[7]),float(row[9])]))
		DEC.append(np.asarray([float(row[4]),float(row[6]),float(row[8]),float(row[10])]))
		#print aor
		#print ra
		#print dec
		#The statement below ensures the order is correct (for MANGLE only)
		if aor > 42 and aor < 45:
			dec1 = min(dec[dec>0])
			ra1 = ra[dec==min(dec[dec>0])][0]
			dec2 = min(dec[dec<0])
			ra2 = ra[dec==min(dec[dec<0])][0]
			dec3 = max(dec[dec<0])
			ra3 = ra[dec==max(dec[dec<0])][0]
			dec4 = max(dec[dec>0])
			ra4 = ra[dec==max(dec[dec>0])][0]
				
		
		else:
			ra1 = min(ra[dec>0])
			dec1 = dec[ra==min(ra[dec>0])][0]
			ra2 = min(ra[dec<0])
			dec2 = dec[ra==min(ra[dec<0])][0]
			ra3 = max(ra[dec<0])
			dec3 = dec[ra==max(ra[dec<0])][0]
			ra4 = max(ra[dec>0])
			dec4 = dec[ra==max(ra[dec>0])][0]
		
	#Convert the RA and DEC list into arrays for easy manipulation
	allRA = np.asarray(RA)
	allDEC = np.asarray(DEC)

	return allRA, allDEC, np.asarray(AOR), np.asarray(epoch), np.asarray(channel)
'''
def Findpointold(rap,decp,rac,decc):
# See if a point (rap,decp) is located in the list of RA and DEC for the full AOR list
	aornum = []
	#loop over the corners
	for i in range(len(rac)-1):
		#define the aor and the point, and check if the point is inside the aor. if so, append the AOR number to the list
		aor = poly(rac[i],decc[i])
		point = pt(rap,decp)
		if aor.contains(point) == True:
			#I have to use the ceiling of the index/4 because each AOR has 2 channels and 2 epochs associated with it. Also, I only need 1 of the AOR instances, not all 4
			if ceil(i/4.) in aornum:
				pass
			#There is no AOR 0, so just append 1 if the 0th index is in the AOR 
			elif ceil(i/4.)==0:
				aornum.append(1.)
			else:
				aornum.append(ceil(i/4.))
			#Below is a way to check if the point is completely contained, or overlaps with another AOR. Simply check different points on the edge of a circle with radius 30' (0.25deg)
			theta = np.arange(0,25,1)*np.pi/12.
			for j in range(len(theta)):
				pointadj = pt(rap+0.25*np.cos(theta[j]),decp+0.25*np.sin(theta[j]))
				aornew = poly(rac[i+1],decc[i+1])
				if aornew.contains(pointadj) == True:
					if ceil((i+1)/4.) in aornum:
						pass
					else:
						aornum.append(ceil((i+1)/4.))
				elif aornew.contains(pointadj) == False:
					
					aornew2 = poly(rac[i-1],decc[i-1])
					if aornew2.contains(pointadj) == True:
						
						if ceil((i-1)/4.) in aornum:
							pass
						else:
							aornum.append(ceil((i-1)/4.))
					elif aornew2.contains(pointadj) == False:
						pass
				else:
					print 'A new problem occured finding the point in the AORs'
			
		elif aor.contains(point) == False:
			pass 
		else:
			print 'A problem occured finding the point in the AORs'
		
	
	
	return np.asarray(aornum)

'''
def Findpoint(rap,decp,rac,decc,aorlist):
# See if a point (rap,decp) is located in the list of RA and DEC for the full AOR list
	aornum = []
	#loop over the corners
	for i in range(len(rac)-1):
		#define the aor and the point, and check if the point is inside the aor. if so, append the AOR number to the list
		aor = poly(rac[i],decc[i])
		point = pt(rap,decp)
		if aor.contains(point) == True:
			#I have to use the ceiling of the index/4 because each AOR has 2 channels and 2 epochs associated with it. Also, I only need 1 of the AOR instances, not all 4
			if AOR[i] in aornum:
				pass
			else:
				aornum.append(AOR[i])
			#Below is a way to check if the point is completely contained, or overlaps with another AOR. Simply check different points on the edge of a circle with radius 30' (0.25deg)
			theta = np.arange(0,25,1)*np.pi/12.
			for j in range(len(theta)):
				pointadj = pt(rap+0.25*np.cos(theta[j]),decp+0.25*np.sin(theta[j]))
				largeridx = np.where(AOR == AOR[i]+1)[0]
				for k in largeridx:
					aornew = poly(rac[k],decc[k])
					if aornew.contains(pointadj) == True:
						if AOR[k] in aornum:
							pass
						else:
							aornum.append(AOR[k])
					elif aornew.contains(pointadj) == False:
						pass
					else:
						print 'A new problem occured finding the point in the next AORs'
				smalleridx = np.where(AOR == AOR[i]-1)[0]
				for h in smalleridx:
					aornew2 = poly(rac[h],decc[h])
					if aornew2.contains(pointadj) == True:
						
						if AOR[h] in aornum:
							pass
						else:
							aornum.append(AOR[h])
					elif aornew2.contains(pointadj) == False:
						pass
					else:
						print 'A new problem occured finding the point in the previous AORs'
			
		elif aor.contains(point) == False:
			pass 
		else:
			print 'A problem occured finding the point in the AORs'
		
	
	
	return np.asarray(aornum)


def Stack_AORs(AORnumarray,outname):
	#This is a hack from another code I wrote, but its functional. Essentially, it will build a string to be printed in the terminal to call SWARP to stack together two (or more) AORs
	run=baseRun[:]
	if len(AORnumarray) < 2:
	#if there is only 1 AOR in the array, no need to stack
		pass
	elif len(AORnumarray) >= 2:
	#if there is more than 1 AOR, stack as follows. The paths may need to be changed depending on how the AOR file structure is set up
		for i in range(len(AORnumarray)):
			if AORnumarray[i]<10:
				n1 = int(AORnumarray[i])
				run.extend(['/home/jtimlin/Fixed_Frame_AOR/AOR_0%s/Combined/ch1/Results/Combine/mosaic.fits' % n1])
			elif AORnumarray[i]>=10:
				n1 = int(AORnumarray[i])
				run.extend(['/home/jtimlin/Fixed_Frame_AOR/AOR_%s/Combined/ch1/Results/Combine/mosaic.fits' % n1])
		run.extend(['-c'])
		run.extend(['./mydefault.swarp'])
		run.extend(["-IMAGEOUT_NAME %s"%(str(outname))])
	#The run parameter is just the connected strings to be printed and ran in the terminal (see next call)
	return run

def runCommand(run):
	#Popen will basically run whatever is in the run variable and the return line is simply a printout without the unnecessary characters
	proc = subprocess.Popen(run, stdout = subprocess.PIPE)
	print ' '.join(run)
	#print [s.rstrip() for s in iter(proc.stdout.readline, '')]
	return [s.rstrip() for s in iter(proc.stdout.readline, '')]

#Find the bright stars in a set of AORs:
def Find_Bright_Stars(aor,starfile,rac,decc,AORlist):
	#I have a starfile from SpIES, so I can just use that here. It contains about 15,000 stars, so we have to get clever here
	data = open(starfile,'rw')
	header = data.readline()
	dat = data.readlines()
	ra = []
	dec = []
	rad = []
	for i in range(len(dat)):
		ra.append(float(dat[i].split()[0]))
		dec.append(float(dat[i].split()[1]))
		rad.append(float(dat[i].split()[2]))
	#Below, I just want to find which stars are in the AORs of interest so that I don't have to maks stars that aren't in the field
	staridx = []
	for j in range(len(aor)):
		#The next two lines are needed because the corner file contains the single epoch corners for both channels (it was easier to have the boxes instead of more complex shapes) 
		#so I'm simply finding and using all 4 of the corners instead of just 2
		racorner = rac[np.where(AOR == aor[j])[0]]
		deccorner = decc[np.where(AOR == aor[j])[0]]
		#next, just find which stars overlap in these AORs and append their indicies to staridx 
		for k in range(len(racorner)):
			tmpaor = poly(racorner[k],deccorner[k])
			for h in range(len(ra)):
				point = pt(ra[h],dec[h])
				if tmpaor.contains(point) == True:
					if h in staridx:
						pass
					else:
						staridx.append(h)
				elif tmpaor.contains(point) == False:
					pass
	header_string = '# Region file format: DS9 version 4.1\n'+'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'+'fk5\n'
	region_string = ''
	#Use the initialized strings above, and the star positions and radii from the star file to generate ds9 region text (to be input to polyregion for masking)
	for i in staridx:
		RA = ra[i]
		DEC = dec[i]
		RAD = rad[i]*3600.0 #degrees to arcsec
		region_string += 'circle('+str(RA)+','+str(DEC)+','+ str(RAD)+'")\n'
	final_star_string = header_string + region_string
	
	return final_star_string


### Mask the bright stars
def mask_bright_stars(imagefile,star_region_strings,outputname):
	f = pf.open(imagefile)
	head = f[0].header
	data = f[0].data
	r = pyregion.parse(star_region_strings)
	r2= r.as_imagecoord(head)
	#create the mask
	mymask = r2.get_mask(hdu=f[0])
	#Reverse the True and False
	mymask = np.logical_not(mymask)
	maskint = mymask*1
	#Multiply the mask and the data to cut the stars (where False == 1 i.e., no mask)
	newimage = data * maskint
	#can write directly if you want
	#pf.writeto(outputname, newimage, head)
	#newimage is the data with the stars cut out whereas maskint is the integer mask with a  0 where the stars are masked and 1 everywhere else and head is the appropriate header
	return data, maskint, newimage

def postage_stamp(data,hdulist, RA, DEC , width=30., height=30.):
	#Cutout a postage stamp using astropy Cutout2D
	#hdulist = pf.open(imagepath)
	hdr = hdulist[0].header
	#data = hdulist[0].data
	w = WCS(hdr, hdulist)
	hdulist.close()
	
	width = width * 60/0.6
	height = height * 60/0.6
	
	position = w.all_world2pix([[RA,DEC]], 1)
	
	
	cutout = Cutout2D(data, position[0], (width,height),wcs = w)
	
	hdr.update(cutout.wcs.to_header())
	
	#pf.writeto('test_cut_dfn.fits', cutout.data, hdr)
	
	return cutout.data, hdr











###################################################
###################################################
###################################################

#RUN THE CODE BELOW

###################################################
###################################################

#Read in the AOR corners for the single epoch AORs from the file AOR_corners.txt
#These corners are in order by AOR, so nominally it should be easy to find the adjacent AOR if need be

RA, DEC, AOR, epoch, channel = openfile('AOR_corners.txt')


#AORnum = Findpoint(331.6,0,RA,DEC,AOR)

tRA = [333.2,5.3,352.6]
tDEC = [-0.1,0,0.1]



aors = []
for i in range(len(tRA)):
	aor = Findpoint(tRA[i],tDEC[i],RA,DEC,AOR)
	aors.append(aor)
	print ' Found AORs to Stack'


#Put the path to Swarp here
baseRun = ['/home/jtimlin/swarp-2.38.0/src/swarp']


for i in range(len(aors)):
	aor = aors[i]
	print aor
	print ' Stacking AORs'
	runCommand(Stack_AORs(aor,'aor%s.fits'%(int(i))))
	print ' Finding the bright stars'
	stars = Find_Bright_Stars(aor,'./Stars_SpIES_rough.txt',RA,DEC,AOR)
	print ' Making the mask'
	data,mask,maskedimage = mask_bright_stars('aor%s.fits'%(str(i)),stars,'codetest.fits')
	print ' Cutting the data image to size'
	cutdata,cutdatahead = postage_stamp(data, pf.open('aor%s.fits'%(str(i))), tRA[i], tDEC[i])
	pf.writeto('Data_cut_%s.fits'%i, cutdata, cutdatahead)
	print ' Cutting the mask image to size'
	cutmask,cutmaskhead = postage_stamp(mask, pf.open('aor%s.fits'%(str(i))), tRA[i], tDEC[i])
	pf.writeto('Mask_cut_%s.fits'%i, cutmask, cutmaskhead)
	print ' Cutting the combined image to size'
	full,fullhead = postage_stamp(maskedimage, pf.open('aor%s.fits'%(str(i))), tRA[i], tDEC[i])
	pf.writeto('Full_cut_%s.fits'%i, full, fullhead)













