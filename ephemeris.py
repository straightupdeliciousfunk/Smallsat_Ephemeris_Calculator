'''
Created on Mar 11, 2017

@author: Tim Owen
'''
from spacetrack import SpaceTrackClient
from math import sin, cos, pi, sqrt, floor
from datetime import datetime
import numpy as np


def main():
    #authenticate spacetrack API
    st = SpaceTrackClient(identity='<insert your registered email>', password='<insert your super-secret password>')
    #pull TLE for NORAD catalog_id (25544 = ISS)
    v = st.tle(norad_cat_id=25544, orderby='epoch desc', limit=1, format='tle')
    
    #split at \n and set variables
    k = v.split("\n")
    #print("NORAD TLE for ISS (ZARYA) - cat id: 25544")
    #print(k[0])
    #print(k[1] + "\n")

    #parse TLE
    t_e = k[0][20:32]     #epoch minute
    e = "." + k[1][26:34]        #eccentricity
    mean_motion = k[1][52:69]        #mean motion (rev/day)
    M_e = k[1][43:52]        #mean anomaly at epoch
    inclination = k[1][8:17]        #inclination
    RAAN = k[1][17:26]      #right angle of ascending node (uppercase omega)
    perigee = k[1][34:42]   #argument of perigee (lowercase omega)
    
    #print("---Orbital Parameters---")
    #print("Epoch: " + t_e + " yrday.fracday")
    #print("Mean Anomaly (@epoch): " + M_e + " degrees") 
    #print("Eccentricity: " + e)
    #print("Mean motion: " + mean_motion + " rev/day")
    #print("Inclination: " + inclination + " degrees")
    #print("RAAN: " + RAAN + " degrees")
    #print("Arg of perigee: " + perigee + " degrees\n")
    
    #convert mean motion into radians per second (rad/s)
    #    # revolutions     2pi radians      1 day 
    #    -------------- x  -----------  x  ----------
    #       day             revolution     86400 seconds 
    n = float(mean_motion) * (2.0 * pi) / 86400.0
    M_e = float(M_e)*pi/180.0
    #print("Mean Motion: " + str(n) + " rad/s")
    #print("M_e: " + str(M_e))
    #get current day and fractional portion of day
    curr_time = datetime.now()
    yd = datetime.now().timetuple().tm_yday #get current day of the year
       
    t = yd + ((curr_time.hour*3600.0 + curr_time.minute*60.0 + curr_time.second) / 86400.0)
    #print("Current time: " + str(t))
    time_diff= float(t) - float(t_e)
    #print("Time difference: " + str(time_diff))
    hours = time_diff*24.0
    #print("Hours: " + str(hours) + "\t" + str(hours*3600))
    minutes = (hours - int(hours))*60.0
    #print("Minutes: " + str(minutes) + "\t" + str(minutes*60))
    seconds = (minutes - int(minutes))*60.0
    #print("Seconds: " + str(seconds))
    temp = ((hours*3600.0)+(minutes*60.0)+seconds)
    #print("Seconds since Epoch: " + str(temp))
    #calculate mean anomaly at current time
    M = (n*int(temp)+M_e)%(2*pi)
    print("Mean Anomaly: " + str(M*180.0/pi) + " degrees")
    #print("Mean Anomaly: " + str(M) + " radians")
    #56.410523467395286
    
    #*****************************************************************************
    #
    #This section of code from: Paul Griffiths
    #
    #Solves for Eccentric Anomaly
    #
    #found @ 'https://github.com/paulgriffiths/pyastro/blob/master/functions.py'
    #
    #<borrowed code>
    desired_accuracy = 1e-6
    e_anom = float(M) #starting guess
    #print("---Eccentric Anomaly---") #added
    #print("Initial guess: " + str(e_anom*180/pi) + " degrees") #added
    while True:
        diff = e_anom - float(e) * sin(e_anom) - float(M)
        e_anom -= float(diff) / (1 - float(e) * cos(e_anom))
        if abs(diff) <= desired_accuracy:
            break
    #</borrowed code>
    #*****************************************************************************
    
    e_anom_degrees = e_anom * (180.0/pi) #convert eccentric anomaly into degrees
    #print("Eccentric Anomaly: " + str(e_anom_degrees) + " degrees")
    
    mu = 398600.5
    J_2 = 1082.63e6
    earthRadius = 6378.137
    i = float(inclination) * (pi/180.0)
    #condensing crazy long formula for mean anomaly...
    p = 1.5 * J_2 * earthRadius**2.0   
    w = 1.0 - ((1.5)* sin(i)**2.0)
    q = 1.0 / sqrt(1.0 - (float(e)**2.0)**3.0) 
    x = p * w * q #into a single variable for newton approximation
    sqrt_mu = sqrt(mu)
    
    #Begin Newton Approximation
    semiMajor_axis = (mu / n**2.0) ** (1. / 3) #initial guess for newton approximation
    #print("---Semi-Major Axis---")
    #print("Initial Guess: " + str(semiMajor_axis))
    #iterate
    while True:
        difference = (sqrt_mu * sqrt((semiMajor_axis**-3.0)))  +  (x * sqrt_mu * sqrt((semiMajor_axis**-7.0)))  -  n
        semiMajor_axis -= float(difference) / (-3.0/2)*sqrt_mu*sqrt((semiMajor_axis**-5.0)) - (7.0/2.0)*x*sqrt_mu*sqrt((semiMajor_axis**-9.0))
        if abs(diff) <= desired_accuracy:
            break
    
    #print("Newton Approx: " + str(semiMajor_axis)) 
    #calculate satellite orbital coordinates
    sat_x = semiMajor_axis * (cos(e_anom) - float(e))
    sat_y = semiMajor_axis * sqrt(1.0 - float(e)**2.0) * sin(e_anom)
    sat_z = 0.0
    
    #print(str(sat_x) + "\t" + str(sat_y) + "\t" + str(sat_z))
    #print("\n---Satellite Cartesian Orbital Coordinates---")
    #print("Xo: " + str(sat_x) + "km")
    #print("Yo: " + str(sat_y) + "km")
    #print("Zo: " + str(sat_z) + "km\n")
    
    #transform into inertial coordinates
    w = float(perigee)*pi/180.0
    omega = float(RAAN)*pi/180.0
    #transformation matrix
    inertial_trans = np.matrix([[cos(w)*cos(omega)-sin(w)*sin(omega)*cos(i),(-1.0)*sin(w)*cos(omega)-cos(w)*cos(i)*sin(omega),sin(omega)*sin(i)],
                                [cos(w)*sin(omega)+sin(w)*cos(i)*cos(omega),(-1.0)*sin(w)*sin(omega)+cos(w)*cos(i)*cos(omega),(-1.0)*cos(omega)*sin(i)],
                                [sin(w)*sin(i),cos(w)*sin(i),cos(i)]])
    orbital_coord = np.matrix([[sat_x],[sat_y],[sat_z]])
    inertial_coord = inertial_trans * orbital_coord
      
    #print("---Satellite Interial Coordinates---")
    #print("Xi: " + str(inertial_coord[0][0]))
    #print("Yi: " + str(inertial_coord[1][0]))
    #print("Zi: " + str(inertial_coord[2][0]) + "\n")
      
      
    #current Year Month and Date for Julian Days (JD) time calculation
    Y = curr_time.year
    M = curr_time.month
    D = t/7
    if (M == 1 | M == 2):
        Y = Y - 1
        M = M + 12
      
    #calculate GMST angle
    A = floor(Y/100)
    B = 2 - A + floor(A / 4)
    JD = floor(365.25*(Y+4716)) + floor(30.6001*(M+1)) + D + B - 1524.5 #verified using an online calculator
    T = (JD - 2415020.0) / 36525

    #GMSTo = 99.6910+36000.7689*T + 0.0004*(T**2.0)
    #wEarth = 7.2921150e-5 * 180 / pi
    GMST = 280.46061837 + (360.98564736629*(JD-2451545.0)) + 0.000387933*(T**2.0) - ((T**3.0)/38710000)
    #GMST = GMSTo + wEarth*JD
    #print("GMST angle: " + str(GMST))
    GMST = GMST % 360 #convert to [0,360] range
    #print("GMST angle: " + str(GMST) + "\n") #this looked right based on ISS tracker
    GMST_rad = GMST * (pi/180.0)
      
    #Transform inertial coordinates into Greenwhich coordinates
    greenwich_trans = np.matrix([[cos(GMST_rad),sin(GMST_rad),0],[(-1.0)*sin(GMST_rad),cos(GMST_rad),0],[0,0,1]])
    greenwich_coord = greenwich_trans * inertial_coord
    #print("---Satellite Greenwich Coordinates---")
    #print("Xg: " + str(greenwich_coord[0][0]))
    #print("Yg: " + str(greenwich_coord[1][0]))
    #print("Zg: " + str(greenwich_coord[2][0]) + "\n")
      
      
    #terminal-satellite vector in greenwich coordinates
    termLatitude = 33.889650*pi/180
    termLongitutde = -77.434887*pi/180
    h = 1.0
    polarRadius = 6356.755
   
    #geocentric latitude of ground terminal
    phi_gc = np.arctan((polarRadius/earthRadius)*np.tan(termLatitude))
    #radius of the earth at ground terminal
    r_t = h + (earthRadius*polarRadius)/sqrt(((earthRadius*sin(phi_gc))**2.0)+((polarRadius*cos(phi_gc))**2.0))
    
    #terminal cartesian coordinates
    terminal_coord = np.matrix([[r_t*cos(termLongitutde)*cos(termLatitude)],
                                [r_t*sin(termLongitutde)*cos(termLatitude)],
                                [r_t*sin(termLatitude)]])
    terminal_x = float(terminal_coord[0][0]) * 1000.0
    terminal_y = float(terminal_coord[1][0]) * 1000.0
    terminal_z = float(terminal_coord[2][0]) * 1000.0
    #print("---Ground Terminal Cartesian Coordinates---")
    #print("Rgx: " + str(terminal_x) + "m")
    #print("Rgy: " + str(terminal_y))
    #print("Rgz: " + str(terminal_z) + "\n")
    #verified cartesian coordinates using online calculator
    #http://www.apsalin.com/convert-geodetic-to-cartesian.aspx
    #values were off slightly... dont know if this effects the overall output much though
    
    
    #calculate the terminal-satellite vector in greenwich coordinates
    #by subtracting ground terminal coordinates from satellite coordinates
    term_sat_greenwich_coord = np.matrix(greenwich_coord - terminal_coord)
    #print("---Terminal-Satellite Vector in Greenwich Coordinates---")
    #print("Rgx: " + str(term_sat_greenwich_coord[0][0]))
    #print("Rgy: " + str(term_sat_greenwich_coord[1][0]))
    #print("Rgz: " + str(term_sat_greenwich_coord[2][0]) + "\n")
    
    #terminal-satellite vector in topocentric coordinates
    topo_trans = np.matrix([[sin(phi_gc)*cos(termLongitutde),sin(phi_gc)*sin(termLongitutde),(-1.0)*cos(phi_gc)],
                            [(-1.0)*sin(phi_gc),cos(phi_gc),0],
                            [cos(phi_gc)*cos(termLongitutde),cos(phi_gc)*sin(termLongitutde),sin(phi_gc)]])
    topo_coord = topo_trans * term_sat_greenwich_coord
    
    #print("---Terminal-Sattelite Vector Topocentric Coordinates")
    #print("Rx: " + str(topo_coord[0][0]))
    #print("Ry: " + str(topo_coord[1][0]))
    #print("Rz: " + str(topo_coord[2][0]) + "\n")
      
    Rx = float(topo_coord[0][0])
    Ry = float(topo_coord[1][0])
    Rz = float(topo_coord[2][0])
    #calculate elevation, azimuth and slant range
    d = sqrt((Rx**2.0)+(Ry**2.0)+(Rz**2.0))
    azimuth = pi - np.arctan(Ry/Rx)
    elevation = float(np.arcsin(Rz/d))*180/pi
    #print("----------------------------------------------------")
    #print("Satellite Elevation: " + str(elevation))
    #print("Satellite Azimuth: " + str(azimuth))
    #print("Satellite Slant Range: " + str(d))
    #print(str(elevation) + "\t" + str(GMST))

if __name__ == "__main__":main()


