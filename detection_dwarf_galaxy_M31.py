import numpy as np 
import math
from scipy import special
import random

def fimag(c, d, a1, a2, b1, b2):

	num = -(a1-c) * (b2-d) + (a2-d) * (b1-c)
	den = (a1-c) *(b1-c) + (a2-d) * (b2-d)
	frac = math.atan2(num, den)
	return frac
	

	
def inPoly(x,y):
	xx=[13.212,13.392,12.346,12.346,11.306,11.268,10.234,10.135,9.109,8.983,8.982,8.974,7.978, 7.867,7.243,7.201,6.407,6.406, 5.525,5.524,4.644,4.642,3.893,3.937,3.407,3.445,2.438,2.437,1.451,1.468,0.462,0.453,0.008,0.005,-0.039,-0.046,-1.045,-1.050,-2.049, -2.053,-3.053,-3.058,-4.059,-4.007,-4.187,-4.196,-5.196,-5.132,-5.345,-5.342,-6.342,-6.262,-6.509,-6.521,-7.522, -7.425,-7.694, -7.713, -8.715,-8.600,-8.882,-8.763,-8.999,-10.004,-9.868,-10.153,-10.010,-10.308,-10.161, -10.463,-10.312,-10.312,-10.897,-10.736,-10.850,-10.688,-10.151,-9.995,-10.087,-9.932,-9.996,-9.835,-8.835, -8.890,-8.941,-8.803,-8.856,-8.713,-7.713,-7.803,-7.679,-6.679,-6.769,-5.771,-5.858,-4.863,-4.946,-4.951,-4.870, -3.873,-4.036, -3.038,-3.103,-2.105,-2.141,-1.250,-0.255,-0.275,-0.270,-0.343,0.650,0.520,1.512,1.459, 2.451,2.387,3.380,3.318,4.311,4.232,5.227,5.149,6.023,7.021,6.939,7.939,8.066,8.069,8.193,8.139,8.077,9.076,9.215, 9.936,10.087,9.947,10.096,10.027,10.174,10.103,10.249,10.223,10.294,10.438,10.391,10.535,10.498,10.642,10.570 ,10.716,10.721,10.867,10.869,10.818,11.831,11.900,12.920,13.099,13.803,13.992,13.859,14.050,14.685,14.884,14.888, 15.089,14.043,14.240]

	yy=[-12.518,-13.623,-13.742,-13.720,-13.830,-13.529,-13.629,-12.777,-12.870,-11.817,-11.814,-11.751,-11.834,-10.791,-10.838, -10.421,-10.479,-10.467,-10.522,-10.502,-10.550,-10.520,-10.554,-11.557,-11.557,-12.600,-12.629,-12.607,-12.626, -12.605,-12.615,-11.569,-11.572,-10.537,-10.546,-10.544,-10.539,-10.538,-10.523,-10.524,-10.497,-10.496,-10.459,-9.430, -9.523, -9.531,-9.483,-8.457,-8.481,-8.536,-8.476,-7.436,-7.420,-7.385,-7.308, -6.292,-6.326,-6.389,-6.296,-5.281,-5.338,-4.327,-4.395,-4.382,-3.373,-3.375,-2.369,-2.365, -1.361,-1.357, -0.354,-0.353,-0.265,0.738,0.744,1.746,1.645,2.668,2.681,3.684,3.675,4.684,4.521,4.633,4.650,5.658,5.655,6.665, 6.511,6.625,7.641,7.501,7.624,7.502,7.607,7.502,7.590,7.595,8.612,8.521,9.573,9.497,9.543,9.486,9.503,9.059,9.042, 8.022,8.020,8.019,8.020,7.817,7.835,7.423,7.460,7.455,7.510,7.494,7.567,7.499,7.590,7.496,7.397,7.526,7.375,7.521, 6.509,6.508,5.498,5.502,5.366,5.520,4.511,4.629,3.624,3.614,2.610,2.613,1.610,1.614,0.611,0.620,-0.392,-1.397, -1.403,-2.409,-2.416,-3.424,-3.329,-4.340,-4.342,-5.355,-5.357,-5.852,-5.714,-6.352,-6.203,-7.253,-7.146, -8.179,-8.181,-9.242,-9.148,-10.196,-10.198,-11.255,-11.402,-12.388]
	simag=0
	nb_pt=len(xx)
	for kk in range(0,nb_pt):
		if (kk < nb_pt-1):
			xe = xx[kk+1]
			xs = xx[kk]
			ye = yy[kk+1]
			ys = yy[kk]
		
		else :
			xe = xx[0]
			xs = xx[kk]
			ye = yy[0]
			ys = yy[kk]
	
		simag = simag + fimag(x,y,xe,ye,xs,ys)
	if(abs(simag)>1e-5): return 1
	else: return 0
	
def radec_xkieta(ra,dec,ra_ref,dec_ref):
	"""Function which transform spherical coordinates of tangent point into tangent plane rectangular coordinates. In input the function takes the right ascension (rad), the declination (rad), the right ascension of M31 (rad), the declination of M31 (rad). The output is xki (rad) and eta (rad)."""
	sdecZ=np.sin(dec_ref);
	sdec=np.sin(dec);
	cdecZ=np.cos(dec_ref);
	cdec=np.cos(dec);
	radiff=np.array(ra)-ra_ref;
	sradiff=np.sin(radiff);
	cradiff=np.cos(radiff);
	denom=sdec*sdecZ+cdec*cdecZ*cradiff;
	xki=cdec*sradiff/denom;
	eta=(sdec*cdecZ-cdec*sdecZ*cradiff)/denom;
	return(xki,eta)
	
def detection_dwarf_galaxy(ra_gal,dec_gal,d_gal,Mv_gal,rh_gal):
	""" Function which determines the recovery fraction of a list of dwarf galaxies in the PAndAS survey. The function takes array like inputs which are the right ascension and declination of the galaxy (rad), the distance modulus of the galaxy, the absolute magnitude of the galaxy in the V-band and the half-light radius of the galaxy in pc. The output contains two arrays. The first is the recovery rates of the dwarf galaxies with similar properties. The second brings precision on the position of the dwarf galaxy. Indeed, if the galaxy is not in the survey, in a hole of the survey or in a field where detection limits are not calculated, it is considered non-detected. The second array contains 0 if the galaxy is not in the survey, 1 if it is in a field of the survey where the recovery fraction were calculated, 2 if the galaxy is in a masked field and 3 if the galaxy is in a hole of the survey. """
	
	ra_gal,dec_gal,d_gal,Mv_gal,rh_gal=np.array(ra_gal),np.array(dec_gal),np.array(d_gal),np.array(Mv_gal),np.array(rh_gal)
	
	###################################################################
	#Finding the field number and the corresponding model parameters: #
	###################################################################
	
	#File which contains the tangent plane rectangular coordinates of the corner of the 406 fields of the PAndAS survey:
	field_corner=np.loadtxt('field_corners.dat')
	#File which contains the center of the 406 fields of the PAndAS survey in the tangent plane rectangular coordinates:
	field_center_xkieta=np.loadtxt('field_middle_xkieta.dat')
	#File which contains the value of the model parameter in function of the field:
	parameters=np.loadtxt('parameters2.txt') 
	i_depth=np.loadtxt('fields_lim_new.csv',delimiter=',')[:,2]
	#Right ascension and Declination in radian of M31:
	ra_refc,dec_refc=0.18648115837037746,0.7202828378876265
	
	#Creation of the array champ which will contain the field number or 10000 if the galaxy center does not fall in one of the fields of PAndAS or in a hole in the survey. Same principle for Mv_lim_dM31, alpha_lim_dM31 and sigma wich will contain the parameters value for the model at the distance of M31 or 10000. 
	champ=np.zeros(len(ra_gal))+10000
	Mv_lim_dM31,alpha_lim_dM31,sigma=np.zeros(len(ra_gal))+10000,np.zeros(len(ra_gal))+10000,np.zeros(len(ra_gal))+10000
	efficiency=np.zeros(len(ra_gal)) #Array which will contains the efficiency value for each galaxy .
	xki_gal,eta_gal=np.array(radec_xkieta(ra_gal,dec_gal,ra_refc,dec_refc))*180/math.pi #Array which contains the xki/eta coordinates of each galaxy.
	label=field_center_xkieta[:,0] #Array which contains all the fields number from 1 to 406.
	position=np.zeros(len(xki_gal))+1 #Array which contains 0 if the galaxy is not in the survey, 1 if it is in a field of the survey where the recovery fraction were calculated, 2 if the galaxy in a masked field and 3 of the galaxy is in a hole of the survey.
	d_M31=24.47
	
	for j in range(0,len(xki_gal)):
		if j%100==0: print('Finding the model parameters: %d/%d'%(j,len(xki_gal)))
		label_tmp=label[(field_center_xkieta[:,1]>xki_gal[j]-2) & (field_center_xkieta[:,1]<xki_gal[j]+2) &(field_center_xkieta[:,2]<eta_gal[j]+2) & (field_center_xkieta[:,2]>eta_gal[j]-2)] #Array which contains only the number of the fields close to the galaxy center.
		for i in range (0,len(label_tmp)):
			k=int(label_tmp[i])-1
			#flimxki and flimeta are the results of the fitting of a straight line to the field limit.
			flimxki=[np.polyfit([field_corner[k*4+3][3],field_corner[k*4][3]],[field_corner[k*4+3][2],field_corner[k*4][2]],1), np.polyfit([field_corner[k*4+1][3],field_corner[k*4+2][3]],[field_corner[k*4+1][2],field_corner[k*4+2][2]],1)]
			flimeta=[np.polyfit([field_corner[k*4][2],field_corner[k*4+1][2]],[field_corner[k*4][3],field_corner[k*4+1][3]],1),np.polyfit([field_corner[k*4+2][2],field_corner[k*4+3][2]],[field_corner[k*4+2][3],field_corner[k*4+3][3]],1)]
			#xkimax and xkimin are the maximum and minimum xki coordinate value for a point with an eta=eta_gal.
			xkimax=flimxki[0][0]*eta_gal[j]+flimxki[0][1]
			xkimin=flimxki[1][0]*eta_gal[j]+flimxki[1][1]
			#etamax and etamin are the maximum and minimum eta coordinate value for a point with a xki=xki_gal.
			etamax=flimeta[0][0]*xki_gal[j]+flimeta[0][1]
			etamin=flimeta[1][0]*xki_gal[j]+flimeta[1][1]
			#Testing if the point is in the field which number is label_tmp[i].
			
			if (xki_gal[j]<xkimax)&(xki_gal[j]>xkimin)&(eta_gal[j]<etamax)&(eta_gal[j]>etamin)&(Mv_lim_dM31[j]==10000):
				#If the galaxy is in the field, the number of the field and the model parameters value for this field are changed in the arrays created above. 
				champ[j]=int(label_tmp[i])
				# If the galaxy fall in a field where the recovery fractions are determined then the value of the parameters is stored in Mv_lim_dM31, alpha_lim_dM31, sigma. If the galaxy fall in a field where the recovery fractions are not calculated then the values of the parameters is 1000000. 
				Mv_lim_dM31[j],alpha_lim_dM31[j],sigma[j]=parameters[int(champ[j]-1)][3],parameters[int(champ[j]-1)][4],parameters[int(champ[j]-1)][5]
			if (xki_gal[j]<xkimax)&(xki_gal[j]>xkimin)&(eta_gal[j]<etamax)&(eta_gal[j]>etamin)&(Mv_lim_dM31[j]!=10000):
			#If the galaxy is inbetween two fields, the parameters for the deepest fields in the i-band are taken. 
				if i_depth[int(label_tmp[i])-1]>i_depth[int(champ[j])-1]:
					champ[j]=int(label_tmp[i])
					Mv_lim_dM31[j],alpha_lim_dM31[j],sigma[j]=parameters[int(champ[j]-1)][3],parameters[int(champ[j]-1)][4],parameters[int(champ[j]-1)][5]
				
		#If the galaxy is not in a field of the survey, it test if it is in a hole or out of the survey.
		if Mv_lim_dM31[j]==10000:
			if inPoly(xki_gal[j],eta_gal[j])==1: position[j]=3 #If the galaxy is in a hole, the position value is 3
			else: position[j]=0 #If the galaxy is not in the survey, the position value is 0


	#############################################################################			
	#Determining the model parameters value at the distance of the galaxy d_gal:#
	#############################################################################
	#Mv_0 and alpha_0 are the intercepts of the straight line describing the variation of the parameters with the distance modulus .
	Mv_0=Mv_lim_dM31+0.89*d_M31 
	alpha_0=alpha_lim_dM31-0.06*d_M31
	#Mv_lim_dgal and alpha_lim_dgal are the model parameters value at the distance modulus of the galaxy d_gal
	Mv_lim_dgal=-0.89*d_gal+Mv_0
	alpha_lim_dgal=0.06*d_gal+alpha_0
	#######################################
	#Determining the detection efficiency:#
	#######################################
	Mvlim=alpha_lim_dgal*np.log10(rh_gal)+Mv_lim_dgal
	beta=special.erfc((Mv_gal-Mvlim)/(math.sqrt(2)*sigma))
	efficiency=0.5*beta
	position=np.where(Mv_lim_dM31==1000000,2,position) #If the galaxy is in a masked field, the position value is 2
	efficiency=np.where(position==1,efficiency,0)
	
	return(efficiency,position)
	
	
