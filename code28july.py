from __future__ import division  

import fileinput
import numpy as np
import matplotlib.pyplot as plt
import glob
import numpy as np
import analisicode.example as ex
from natsort import realsorted, ns
import re
import os
import subprocess


#def le costanti nel SI
G= 6.67*(10**-11) # N m2 /kg2
c= 3*(10**8)  # m/s 
Msun= 1.99*(10**30) #kg
pc= 3.09*(10**16) #m
day= 86400 #s



files=realsorted(glob.glob('../output*.dat'))
N=len(files)

full_path = os.path.realpath(__file__)
print 'working on sim number:', os.path.dirname(full_path)


w=np.zeros(shape=(N,25))


for j in range(N):
	fp=open(files[j])
	for i, line in enumerate(fp):
		if line.endswith('  1   14\n'):# i == 0:
			riga=line.split()
			w[j,0],w[j,1],w[j,2],w[j,3],w[j,4],w[j,5],w[j,6],w[j,7],w[j,8],w[j,9]= [float(x) for x 
in riga]	
					
			#print '-'*30
		elif line.endswith('  2   14\n'):		
			riga=line.split()
			w[j,10],w[j,11],w[j,12],w[j,13],w[j,14],w[j,15],w[j,16],w[j,17],w[j,18],w[j,19]= [float
(x) for x in riga]
						
		elif line.startswith(' #time_is='):
			riga=line.split()
			#print riga
			w[j,20], w[j,21], w[j,22], w[j,23]=[float(x) for x in riga[1:]]
			

x1 = w[:,0]
y1 = w[:,1]
z1 = w[:,2]
vx1= w[:,3]
vy1= w[:,4]
vz1= w[:,5]
M1 = w[:,6]
m1 = w[:,7]
id1= w[:,8]
s1 = w[:,9]
x2 = w[:,10]
y2 = w[:,11]
z2 = w[:,12]
vx2= w[:,13]
vy2= w[:,14]
vz2= w[:,15]
M2 = w[:,16]
m2 = w[:,17]
id2= w[:,18]
s2 = w[:,19]
ttot=w[:,21]
tstar=w[:,22]


print '       '
print 'SELECT THE PHYSICAL SCALING FROM THE OUT OF NBODY6'


fname =('../out')

outfile=open('scaling.txt','w')
#for fname in filenames:
with open(fname) as infile:
	for line in infile:
		if line.startswith('            PHYSICAL SCALING:'): line1 =line[:-1]
		else:pass
outfile.write(line1+'\n')

outfile.close()

rscale,mscale,vscale,tscale= np.loadtxt('scaling.txt', usecols=(4,7,10,13), unpack=True, comments='#') 
			
print '-'*50			



print 'STARTING ASROPHYSICAL CALCULATION'	
	
	
print '       '
print '-'*50
print 'Calculation of geometrical parameters'
	
r,rx,ry,rz=ex.calcVect((x1,y1,z1),(x2,y2,z2))
v,vx,vy,vz=ex.calcVect((vx1,vy1,vz1),(vx2,vy2,vz2))



vq2=v*v



print '       '
print '-'*50
print 'Calculation of the Angular Momentum Lx, Ly, Lz and L'
print '       '

L,Lx,Ly,Lz=ex.calcL((rx,ry,rz),(vx,vy,vz))



print '       '
print 'Calculation of the total BHB mass = mt'   
print '       '

mt=m1+m2




print '       '
print 'Calculation of the reduced mass = mu'
print '       '
mu=(m1*m2)/mt




print '       '
print 'Calculation of the semi-major axis of the BHB = semi'
semi=( (2./r)-(vq2/mt) )**(-1)

S=semi*rscale


print '       '
print 'Calculation of the Kinetic Energy of the BHB = Kin'
print '       '
kin=0.5*mu*vq2




print '       '
print 'Calculation of the Potential Energy of the BHB = Pot'
print '       '
	
pot=(mu*mt)/r




print '       '
print 'Calculation of the Total Binding Energy of the BHB = E'
print '       '

E=kin-pot

print '       '
print 'Calculation of the Eccentricity of the BHB = e'
print '       '

e=np.sqrt(1.-((L**2)/(mt*semi) )) 


print '       '
print '-'*50			


print '       '
print 'Calculation of the BHB pericenter'
print '       '

peric = S*(1-e)

T= ttot*tscale	



#calcolo il tempo di coalescenza
q=1.0
ma=30.0
mb=30.0



def Tgw(q,S,e,ma,mb):
	return (5.8*(10**6)*( ( (1+q)**2 )/q)* ((S/0.01)**4)*(((ma+mb)/10**8)**-3)*(1-(e**2))**(7./2.) )

ttgw=Tgw(q,S,e,ma,mb)

tgw=ttgw/(10**9)

print '     '


#FILTRO per evitare i ghost
'''
idx=r!=0
r=r[idx]
T=T[idx]
e=e[idx]
S=S[idx]
peric=peric[idx]
tgw=tgw[idx]
'''
#print len(r)
print '          ' 

idxx=r<200
r=r[idxx]
T=T[idxx]
e=e[idxx]
S=S[idxx]
peric=peric[idxx]
tgw=tgw[idxx]
m1=m1[idxx]
m2=m2[idxx]


#CREO UN VETTORE IN CUI SCRIVO TUTTE LE QUANTITA OTTENUTE SCALATE E LE SCRIVO IN UN FILE A PARTE

#print len(T)
z=np.zeros(shape=(len(r),6))
z[:,0]=T #tempo(Myr)
z[:,1]=r*rscale #sep (pc)
z[:,2]=e #ecc
z[:,3]=S #semiasse (pc)
z[:,4]=peric #peric (pc)
z[:,5]=tgw
#z[:,6]=m1*mscale
#z[:,7]=m2*mscale


#%.18g
#SALVO I PARAMETRI FISICIS IN UN FILE CHE SI CHIAMA PARAM.TXT
np.savetxt('fisparam.txt', z,fmt='%.10g', delimiter='    ', header='tempo, sep, ecc, semi, pericentre, tgw ')
print 'The physical quantities calculated are stored in the external file fisparam.txt'
print '       '






'''
print '-'*50			
print '       '
print 'Plot BHB sep with the Time'	
	
plt.plot(z[:,0],z[:,1],'-',color='red',linewidth=1.5, label='BHB separation')
#maxt=np.nanmax(z[:,0])
#print maxt
#maxd=np.nanmax(z[:,1])
#mind=np.nanmin(z[:,1])
#print maxd, mind

plt.title('BHB distance vs Time')
plt.xlabel('Evolution Time (ttot*tscale) Myr')
plt.ylabel('Binary distance (D*rscale) pc')
plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.legend(loc='upper right')

#plt.axis([ -10, maxt+10%maxt, mind-20%mind, maxd+mind] )
plt.savefig('dist_time.pdf')

#visualizza plot
#plt.show()
plt.close()





print '       '
print 'Plot BHB semi major axis with the Time'	

plt.plot(z[:,0],z[:,3],'-',color='blue', linewidth=1.5, label='BHB-semi maj axis')

#maxt=np.nanmax(z[:,0])
#print maxt
#maxs=np.nanmax(z[:,3])
#mins=np.nanmin(z[:,3])
#print 'max',maxs,'min', mins

plt.title('BHB semi major axis vs Time')
plt.xlabel('Evolution Time Myr')
plt.ylabel('Semi majopr axis pc')
plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.legend(loc='upper right')

#plt.axis([ -10, maxt+10%maxt, mins-0.001%mins, maxs+0.001%mins] )
plt.savefig('semi_time.pdf')

#visualizza plot
#plt.show()
plt.close()
'''





'''
print '       '
print 'Plot BHB Energy with the Time'	
	
plt.plot(z[:,0],z[:,6],'-',color='green',linewidth=1.5, label='BBH-Energy')

maxt=np.nanmax(z[:,0])
#print maxt
maxs=np.nanmax(z[:,6])
mins=np.nanmin(z[:,6])
print 'max',maxs,'min', mins

plt.title('ASTRO(output_.dat)Energy vs Time')
plt.xlabel('Evolution Time (ttot*tscale) Myr')
plt.ylabel('Energy (?)')
#plt.xscale('log')
#plt.yscale('log')
plt.grid(True)
plt.legend(loc='upper right')

#plt.axis([ -10, maxt+10%maxt, mins-0.001%mins, maxs+0.001%mins] )
plt.savefig('ene_time.pdf')

#visualizza plot
#plt.show()
plt.close()
'''





'''
print '       '
print 'Plot BHB eccentricity with the Time'	
	
plt.plot(z[:,0],z[:,2],'-',color='magenta',linewidth=1.5, label='BBH-eccentricity')

#maxt=np.nanmax(z[:,0])
#print maxt
#maxs=np.nanmax(z[:,7])
#mins=np.nanmin(z[:,7])
#print 'max',maxs,'min', mins

plt.title('BHB eccentricity vs Time')
plt.xlabel('Evolution Time (ttot*tscale) Myr')
plt.ylabel('eccentricity')
plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.legend(loc='upper right')

#plt.axis([ -10, maxt+10%maxt, mins-0.01%mins, maxs+0.1%mins] )
plt.savefig('ecc_time.pdf')

#visualizza plot
#plt.show()
plt.close()




print '       '
print 'Plot BHB pericenter with the Time'	
	
plt.plot(z[:,0],z[:,4],'-',color='darkgreen',linewidth=1.5, label='BBH-pericenter')

#maxt=np.nanmax(z[:,0])
#print maxt
#maxs=np.nanmax(z[:,7])
#mins=np.nanmin(z[:,7])
#print 'max',maxs,'min', mins

plt.title('BHB pericenter vs Time')
plt.xlabel('Evolution Time (ttot*tscale) Myr')
plt.ylabel('pericenter')
plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.legend(loc='upper right')

#plt.axis([ -10, maxt+10%maxt, mins-0.01%mins, maxs+0.1%mins] )
plt.savefig('peric_time.pdf')

#visualizza plot
#plt.show()
plt.close
'''







'''



fileinput =('../out')

		

lineout = subprocess.check_output(['tail', '-1', fileinput])

print lineout

if lineout.startswith('         END RUN'):
	print 'ok, creating files'
	#filein= ('fisparam.txt')	
	#lineout1 = subprocess.check_output(['tail', '-1', filein])
	#g=np.zeros(shape=(len(T),1))	
	#g=z[:,0]	
	#for i in range(len(g)):
		#print g1[i]
		#mask = (1000.0 < g) & (g < 1000.1)
		#g=g[mask]#g1=g[mask]
		#print g(i) 
	#print g1
	mask = (1000.0 < z[:,0]) & (z[:,0] < 1090.99)
	znew1=z[mask,:]
	#print znew1
#	np.savetxt('uno.txt', znew1, fmt='%.10g', delimiter='    ')
	np.savetxt('uno.txt',znew1[0,:],  fmt='%.10g',newline='   ', footer='\n', comments='')
	
	
	mask = (2900.0 < z[:,0]) & (z[:,0] < 3000.99)
        znew3=z[mask,:]
        #print znew3
#        np.savetxt('tre.txt', znew3, fmt='%.10g', delimiter='     ')
	np.savetxt('tre.txt', znew3[0,:], fmt='%.10g', newline='   ', footer='\n', comments='')


elif lineout.startswith('         CALCULATIONS HALTED'): 			print ' CALCULATION HALTED!!!!'
else:print 'something went wrong!'

'''







