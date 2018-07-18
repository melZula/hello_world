#! /usr/bin/env pyhon
# -*- coding: utf-8 -*-
import math
import numpy as np
import pylab
from matplotlib import mlab
import matplotlib.pyplot as plt
import matplotlib.animation as animation
"""
print("d=")
d=input()
d=float(d)
print("d2=")
d2=input()
d2=float(d2)
print("d3=")
d3=input()
d3=float(d3)
"""
d=2
d2=2
d3=2

print("alph_min(deg)=")
#alph_min=input()
alph_min=-90
alph_min=float(alph_min)
alph_min=alph_min*math.pi/180
print("alph_max(deg)=")
#alph_max=input()
alph_max=90
alph_max=float(alph_max)
alph_max=alph_max*math.pi/180
print("alph_lim(deg)=")
#alph_lim=input()
alph_lim=1
alph_lim=float(alph_lim)
alph_lim=alph_lim*math.pi/180
print("alph_min(rad)=", alph_min)
print("alph_max(rad)=", alph_max)


##############
x=2
y=2

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)

coords = []

def onclick(event):
	global x, y
	x, y = event.xdata, event.ydata
	print ('x = ',x,' y = ',y)

	global coords
	coords.append((x, y))

	
	return coords
cid = fig.canvas.mpl_connect('button_press_event', onclick)
##############


def calculation_rest_angls(x, y, a, d2, d3):
	x=x-(math.cos(a)*d)
	y=y-(math.sin(a)*d)
	if x>0:
		a2=math.acos((d2**2+x**2+y**2-d3**2)/(2*d2*math.sqrt(x**2+y**2)))+math.asin(y/math.sqrt(x**2+y**2))
	else:
		a2=math.acos((d2**2+x**2+y**2-d3**2)/(2*d2*math.sqrt(x**2+y**2)))-math.asin(y/math.sqrt(x**2+y**2))+math.pi
	a3=math.pi-(math.pi-a)-math.acos((d3**2+x**2+y**2-d2**2)/(2*d3*math.sqrt(x**2+y**2)))
	return (a2, a3)
			
			
def intersection(circlerad, x, y, alph_min, alph_max, d):
	dot0=0
	dot1=0
	c=(x**2+y**2+d**2-(circlerad)**2)/2
	a_loc=y**2+x**2
	b=(-2)*c*y
	e=c**2-(d**2)*(x**2)
	dscr=b**2-4*a_loc*e
	print("dscr=",dscr)
	if dscr>0:
		if alph_min<math.asin((-b+math.sqrt(dscr))/(2*a_loc)/d)<alph_max: # проверка принадлежности дуги
			dot0=math.asin((-b+math.sqrt(dscr))/(2*a_loc)/d)
			plt.plot([math.cos(dot0)*d], [(-b+math.sqrt(dscr))/(2*a_loc)], 'ro')

		if alph_min<math.asin((-b-math.sqrt(dscr))/(2*a_loc)/d)<alph_max:
			dot1=math.asin((-b-math.sqrt(dscr))/(2*a_loc)/d)
			plt.plot([math.cos(dot1)*d], [(-b-math.sqrt(dscr))/(2*a_loc)], 'ro')
	return (dot0, dot1)

	
def calculations(d, d2, d3, x, y, alph_min, alph_max, alph_lim):
	try:
		if (d+d2+d3)<math.sqrt(x**2+y**2):
			raise Exception("Unreachable dot")
		if (d<0 or d2<0 or d3<0):
			raise ValueError("Negative length")


		#проверка на значение(например alph_min) что бы отбросить вариант без пересения кольца и дуги
		if math.sin(alph_lim)*d2<math.sqrt((x-math.cos(alph_min)*d)**2+(y-math.sin(alph_min)*d)**2)<d3+d2:
			print("without intersection")
			a=alph_min

			a2, a3=calculation_rest_angls(x, y, a, d2, d3)
		else:
			dots=[0,0,0,0]


			circlerad=d3+d2
			dots[0]=intersection(circlerad, x, y, alph_min, alph_max, d)[0]
			dots[3]=intersection(circlerad, x, y, alph_min, alph_max, d)[1]

			circlerad=math.sin(alph_lim)*d2
			dots[1], dots[2]=intersection(circlerad, x, y, alph_min, alph_max, d)


		       
			print ("dots=",dots)
			# вычисление первого угла (a)
			if (dots[0]<>0 and dots[1]==0 and dots[2]==0 and dots[3]==0):
				a=alph_min+(dots[0]-alph_min)/2
				print ('dot-0')
			if (dots[0]==0 and dots[1]==0 and dots[2]==0 and dots[3]<>0):
				a=dots[3]+(alph_max-dots[3])/2
				print ('dot-3')
			if (dots[0]<>0 and dots[1]==0 and dots[2]==0 and dots[3]<>0):
				a=dots[3]+(dots[0]-dots[3])/2
				print ('dot-0,3')
			if (dots[0]<>0 and dots[1]<>0):
				a=dots[1]+(dots[0]-dots[1])/2
				print ('dot-0,1')
			if (dots[2]<>0 and dots[3]<>0):
				a=dots[3]+(dots[2]-dots[3])/2
				print ('dot-2,3')
		
			# вычисление остальных углов
			a2, a3=calculation_rest_angls(x, y, a, d2, d3)
		return(a, a2, a3)
	
	except Exception as error:
		print ("Error in function Calculations:", error) 
		
	
		
	

	
		
def animate(i):
	try:
		if calculations(d, d2, d3, x, y, alph_min, alph_max, alph_lim)==None:
			raise TypeError('can not calculate angls')

		a, a2, a3=calculations(d, d2, d3, x, y, alph_min, alph_max, alph_lim)	
		ax1.clear()

	
		# ограничение первого звена
		t = np.arange((45-(alph_max*180/np.pi-45))*np.pi/180,(45-(alph_min*180/np.pi-45))*np.pi/180,0.01)
		r=d
		plt.plot(r*np.sin(t),r*np.cos(t))

		# малый круг
		t = np.arange(0,2*np.pi,0.01)
		r=math.sin(alph_lim)*d2
		plt.plot(r*np.sin(t)+x,r*np.cos(t)+y)

		# большой круг
		t = np.arange(0,2*np.pi,0.01)
		r=d3+d2
		plt.plot(r*np.sin(t)+x,r*np.cos(t)+y)

		# окружность 3-его звена
		t = np.arange(0,2*np.pi,0.01)
		r=d3
		plt.plot(r*np.sin(t)+x,r*np.cos(t)+y)

		# все звенья
		plt.plot([0, round(math.cos(a),3)*d, round(math.cos(a),3)*d+round(math.cos(a2),3)*d2, x], [0, round(math.sin(a),3)*d, round(math.sin(a),3)*d+round(math.sin(a2),3)*d2, y], lw=3)

		plt.plot([0, round(math.cos(a),3)*d, round(math.cos(a),3)*d+round(math.cos(a2),3)*d2], [0, round(math.sin(a),3)*d, round(math.sin(a),3)*d+round(math.sin(a2),3)*d2], 'ro')

		plt.scatter([x],[y], marker="x", s=100)

		plt.grid(True)
		plt.axis('equal')
		plt.axis([-10, 10, -10, 10])
	except Exception as error:
		print ("Error in function Animate:", error)

		
try:	
		
		#if (d+d2+d3)<math.sqrt(x**2+y**2):
		#	raise Exception("Unreachable dot")
		if (d<0 or d2<0 or d3<0):
			raise ValueError("Negative length")
		
			
				
		ani = animation.FuncAnimation(fig, animate, interval=500)
		
		

except Exception as error:
	print ("Error :", error) 

plt.grid(True)
plt.axis('equal')
plt.show()





