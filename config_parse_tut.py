import argparse
import ConfigParser
import rebound

from collections import OrderedDict

import numpy as np
import math
import random as rand

import matplotlib.pyplot as plt


def rotate_vec(angle,axis,vec):
	'''
	Rotate vector vec by angle around axis
	'''
	vRot = vec*math.cos(angle) + np.cross(axis,vec)*math.sin(angle) + axis*np.dot(axis,vec)*(1 -math.cos(angle))
	return vRot

def gen_disk(ang):
	'''
	This is from some old code that starts with perfectly aligned e and j vectors and then rotates them by a small amount
	'''
	ehat = np.array([1,0,0])
	jhat = np.array([0,0,1])
	bhat = np.cross(jhat,ehat)    # rotate jhat by angle1 over major axis and angle minor axis
	# rotate ehat by angle2 over minor axis and angle3 about jhat
	angle1 = np.random.normal(0.0, ang, 1)
	angle2 = np.random.normal(0.0, ang, 1)
	angle3 = np.random.normal(0.0, ang, 1)
	jhat = rotate_vec(angle1,ehat,jhat)
	jhat = rotate_vec(angle2,bhat,jhat)
	ehat = rotate_vec(angle2,bhat,ehat)
	ehat = rotate_vec(angle3,jhat,ehat)
	n = np.cross(np.array([0,0,1]), jhat)
	n = n / np.linalg.norm(n)
	Omega = math.atan2(n[1], n[0])
	omega = math.acos(np.dot(n, ehat))
	if ehat[2] < 0:
		omega = 2*np.pi - omega
	inc=math.acos(jhat[2])
	return inc, Omega, omega


def density(min1, max1, p):
	'''
	Generate a random from a truncated power law PDF with power law index p.
	min1 and max1

	'''
	r=np.random.random(1)[0]
	if p==1:
		return min1*np.exp(r*np.log(max1/min1))
	else:
		return (r*(max1**(1.-p)-min1**(1.-p))+min1**(1.-p))**(1./(1-p))

def get_tde(sim, reb_coll):
	orbits = sim[0].calculate_orbits(primary=sim[0].particles[0])
	p1,p2 = reb_coll.p1, reb_coll.p2
	idx, idx0 = max(p1, p2), min(p1, p2)
	if idx0==0:
		##idx decremented by 1 because there is no orbit 0
		name=sim[0].simulationarchive_filename
		f=open(name.replace('.bin', '_tde'), 'a+')
		f.write('{0} {1} {2} {3} TDE!\n'.format(sim[0].t, orbits[idx-1].a, orbits[idx-1].e, idx))
		f.close()

	return 0


parser=argparse.ArgumentParser(
	description='Set up a rebound run')
parser.add_argument('--config', nargs=1, default='config',
	help='File containing simulation parameters')
# parser.add_argument('--keep_bins', action='store_true',
# 	help="Don't delete bins from simulation")


##Parsing command line arguments.
args=parser.parse_args()
config_file=args.config
##Unique tag for output file.
##

##Default stellar parameters
config=ConfigParser.SafeConfigParser(defaults={'name': 'archive', 'N':'100', 'e':'0.7',
	'gravity':'basic', 'integrator':'ias15', 'dt':'0', \
	'a_min':'1.', 'a_max':'2.', 'ang':'2.', 'm':'5e-5', 'keep_bins':'False',\
	 'rt':'1.0e-4', 'coll':'line',\
	'pRun':'500', 'pOut':'0.1',
	'p':'1'}, dict_type=OrderedDict)
# config.optionxform=str
config.read(config_file)

##Name of our put file
name=config.get('params', 'name')
name=name+".bin"
##Length of simulation and interval between snapshots
pRun=config.getfloat('params', 'pRun')
pOut=config.getfloat('params', 'pOut')
##Tidal radius
rt=config.getfloat('params', 'rt')
##Name of collision detection routine
coll=config.get('params', 'coll')
##Whether to keep the binaries
keep_bins=config.getboolean('params', 'keep_bins')


sim = rebound.Simulation()
sim.G = 1.
##Central object
sim.add(m = 1, r=rt)
sim.collision_resolve=get_tde

print sim.particles[0].r
sim.gravity=config.get('params', 'gravity')
sim.integrator=config.get('params', 'integrator')
sim.collision=coll

dt=config.getfloat('params', 'dt')

print sim.gravity

sections=config.sections()
##Add particles; Can have different sections with different types of particles (e.g. heavy and light)
##see the example config file in repository. Only require section is params which defines global parameters
##for the simulation (pRun and pOut).
fig,ax=plt.subplots(figsize=(10,9))
for ss in sections:
	if ss=='params':
		continue
	a1=[]
	N=int(config.get(ss, 'N'))
	e=config.getfloat(ss, 'e')
	m=config.getfloat(ss, 'm')
	a_min=config.getfloat(ss, 'a_min')
	a_max=config.getfloat(ss, 'a_max')
	p=config.getfloat(ss, 'p')
	ang=config.getfloat(ss, 'ang')

	for l in range(0,N): # Adds stars
		##Use AM's code to generate disk with aligned eccentricity vectors, but a small scatter in i and both omegas...
		inc, Omega, omega=gen_disk(ang*np.pi/180.)
		a0=density(a_min, a_max, p)
		a1.append(a0)
		M = rand.uniform(0., 2.*np.pi)
		sim.add(m = m, a = a0, e = e, inc=inc, Omega = Omega, omega = omega, M = M, primary=sim.particles[0])
	ax.hist(a1, alpha=0.3, normed=True)

fig.savefig('hist_a.pdf')
oo=sim.calculate_orbits(primary=sim.particles[0])

print len(sim.particles), oo[1].e, sim.particles[1].m, sim.particles[-1].m
