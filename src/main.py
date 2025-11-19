"""
Author: Konrad Barboutie
Current date: 2025-11-10
"""
from math import floor
import webbrowser, os, sys, re, requests, numpy as np, plotly.io as pio, json
pio.renderers.default = "browser"
from jplephem.calendar import compute_julian_date
from poliastro.maneuver import Maneuver
from poliastro.bodies import Earth, Sun
from poliastro.twobody import Orbit
from poliastro.twobody.sampling import EpochsArray, TrueAnomalyBounds, EpochBounds
from poliastro.util import time_range
from poliastro.plotting import OrbitPlotter3D
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import CartesianDifferential, CartesianRepresentation
from poliastro.plotting.interactive import OrbitPlotter3D
from poliastro import iod
from poliastro.bodies import Earth, Mars, Sun
from poliastro.ephem import Ephem
from poliastro.maneuver import Maneuver
from poliastro.twobody import Orbit
from poliastro.util import time_range
from poliastro.core.elements import coe_rotation_matrix, coe2rv
from poliastro.plotting._base import Trajectory
import plotly.io as pio
########################################################################################################################
# Function definitions
########################################################################################################################
def get_elements(url):
    # Acquiring and Storing the data as a json file
    response = requests.get(url)
    data = response.json()

    # Verify connection was smooth
    if response.status_code != 200:
        print("ERROR GETTING THE DATA")
        sys.exit()
    else:
        print("Status Code: ", response.status_code) # 200 means everything ok

    # Extracting orbital elements code by Gemini AI, improved by me
    start_marker = '$$SOE'
    end_marker = '$$EOE'
    start_index = data['result'].find(start_marker)
    end_index = data['result'].find(end_marker)
    if start_index == -1 or end_index == -1:
        print("Error: Could not find the $$SOE / $$EOE markers in the response.")
        exit()
    # Extract the block (start_index + length of '$$SOE' marker)
    ephemeris_block = data['result'][start_index + len(start_marker): end_index].strip()

    # Focus on the data for the first date (June 1, 2030)
    # This is the block starting after $$SOE and ending before the second date
    first_date_block = ephemeris_block.split('2462654.500000000')[0].strip()

    # Extract the key-value pairs using regular expressions
    # This regex looks for: two capital letters (the element code, e.g., EC, QR),
    # followed by an '=' sign and then the floating-point number.
    elements_dictionnary = {}
    matches = re.findall(r'(\s[A-Z]{1,2})\s*=\s*([0-9.\-E+]+)', first_date_block)

    for key, value in matches:
        # Clean the key (remove leading/trailing spaces) and store as a float
        elements_dictionnary[key.strip()] = float(value)
    #   EC      Eccentricity
    #   QR      Periapsis distance
    #   IN      Inclination w.r.t. xy-plane (degrees)
    #   OM      Longitude of Ascending Node (degrees)
    #   W       Argument of Perifocus (degrees)
    #   Tp      Periapsis time (user specifies absolute or relative date), not available for 2008 EV5 so not included in the key list
    #   N       Mean motion (degrees/DU)
    #   MA      Mean anomaly (degrees)
    #   TA      True anomaly (degrees)
    #   A       Semi-major axis
    #   AD      Apoapsis distance
    #   PER     Orbital Period or PR ?

    KEY_ORDER = ['A', 'EC', 'IN', 'OM', 'W', 'TA', 'QR', 'N', 'MA', 'AD', 'PR'] # a, e, i, RAAN, argp, nu
    elements_sorted_values = [elements_dictionnary[key] for key in KEY_ORDER]
    # Storing the values in a numpy array for future use
    elements = np.array(elements_sorted_values)
    return elements
def normalize(v):
    norm = np.sqrt(v[0]**2+v[1]**2+v[2]**2)
    v_norm = v/norm
    return v_norm
def length(v):
    v = np.sqrt(v[0]**2+v[1]**2+v[2]**2)
    return v
def node_line(x,y):
    return np.cross(x,y)
def eccentricity_vector(mu,r,v):
    rl = length(r)
    vl = length(v)
    e = 1/mu*((vl**2-mu/rl)*r-(np.dot(r,v))*v)
    return e
def radius_from_elements_pqw(a, e, nu):
    r = a*(1-pow(e,2))/(1+e*np.cos(nu))
    p = r*np.cos(nu)
    q = r*np.sin(nu)
    w = 0
    return [p, q, w]
def velocity_from_elements_pqw(a, e, nu, mu):
    semi_latus_rectum = a*(1-pow(e,2))
    r = np.sqrt(mu/semi_latus_rectum)
    p = -r*np.sin(nu)
    q = r*(e+np.cos(nu))
    w = 0
    return[p, q, w]
def rv_calculation(elements):
    r = radius_from_elements_pqw(elements[0], elements[1], elements[6]) * u.km
    v = velocity_from_elements_pqw(elements[0], elements[1], elements[6], mu_sun) * u.km/u.s
    return r, v
def define_orbit(origin, classical_elements, epoch):
    orb = Orbit.from_classical(origin, classical_elements[0] * u.km, classical_elements[1] * u.one, classical_elements[2] * u.deg, classical_elements[3] * u.deg, classical_elements[4] * u.deg, classical_elements[5] *u.deg, epoch)
    return orb
def define_ephem(origin, classical_elements):
    orb = Ephem.from_orbit(origin, classical_elements[0] * u.km, classical_elements[1] * u.one, classical_elements[2] * u.deg, classical_elements[3] * u.deg, classical_elements[4] * u.deg, classical_elements[5] *u.deg)
    return orb
def sphere_of_influence(m1, r):
    mu = 1.32712440042e11 # mu of sun
    return r*pow(m1/mu, 2/5)
def angular_momentum_normalized(i, raan):
    h = [np.sin(i)*np.sin(raan), -np.sin(i)*np.cos(raan), np.cos(i)]
    return h
def nu_to_eccentric_anomaly(e, nu):
    E = 2*np.arctan(np.tan(nu/2)*np.sqrt((1-e)/(1+e)))
    return E
def time_to_node_line(a, e, w, nu, mu):
    # nu corresponds to current true anomaly
    E0 = nu_to_eccentric_anomaly(e, nu) # Initial E
    E_asc =  nu_to_eccentric_anomaly(e, 2*np.pi-w)# Final E
    E_des = nu_to_eccentric_anomaly(e, np.pi-w)
    dt_asc = np.sqrt(a**3/mu)*((E_asc-e*np.sin(E_asc))-(E0-e*np.sin(E0)))
    dt_des = np.sqrt(a**3/mu)*((E_des-e*np.sin(E_des))-(E0-e*np.sin(E0)))
    if dt_asc < 0:
        dt_asc += np.sqrt(a**3/mu)*(2*np.pi+(E_asc-e*np.sin(E_asc))-(E0-e*np.sin(E0)))
    if dt_des < 0:
        dt_des = np.sqrt(a**3/mu)*(2*np.pi+(E_des-e*np.sin(E_des))-(E0-e*np.sin(E0)))
    return dt_asc, dt_des
def julian_to_current(jd):
    j = jd+0.5
    z = floor(j)
    alpha = floor((z-1867216.25)/36524.25)
    A = z+1+alpha-floor(alpha/4)
    B = A +1524
    C = floor((B-122.1)/365.25)
    D = floor(365.25*C)
    E = floor((B-D)/30.6001)
    day = B-D-floor(30.6001*E)
    if E<13:
        month = E-1
    else:
        month= E-13
    if month==1 or month==2:
        year = C-4716
    else:
        year = C-4715
    day_int = round(day)
    date_str = f"{year}-{month:02d}-{day_int:02d}"
    return date_str
def argument_perigee(n,e):
    nl = length(n)
    el = length(e)
    w = np.arccos(np.dot(n,e)/(el*nl))
    if e[2]<0:
        w = w+np.pi
        return w
    else:
        return w
def angular_momentum(r,v):
    return np.cross(r,v)


########################################################################################################################
# Data Collection & Initialization
########################################################################################################################
# Write Horizon ephemeris data to file for faster execution
flag = 1 # 0 = get data, 1 = load from file
file = 'data.json'

# Initial date setting
year = 2030
observation_time = 1
ephemeris_start = compute_julian_date(year = year, month = 6, day = 1)
ephemeris_stop = compute_julian_date(year = year+observation_time, month= 6, day = 2)

ti = Time(ephemeris_start, format = 'jd') # initial time
tf = Time(ephemeris_stop, format = 'jd') # final time

# Url's to connect to jpl Horizon API and retrieve data
url_2008_ev5 = (
    "https://ssd.jpl.nasa.gov/api/horizons.api?"
    + "format=json" + "&"
    + "COMMAND='DES=2008 EV5'" + "&"
    + "OBJ_DATA='NO'" + "&"
    + "MAKE_EPHEM='YES'" + "&"
    + "CENTER='0'" + "&"
    + "EPHEM_TYPE='ELEMENTS'" + "&"
    + "START_TIME='JD"+ str(ephemeris_start) + "'" + "&" #str(launch_date)
    + "STOP_TIME='JD" + str(ephemeris_stop) + "'" + "&"   # launch
    + "STEP_SIZE='1 d'"
)
url_earth = (
    "https://ssd.jpl.nasa.gov/api/horizons.api?"
    + "format=json" + "&"
    + "COMMAND='399'" + "&"
    + "OBJ_DATA='NO'" + "&"
    + "MAKE_EPHEM='YES'" + "&"
    + "CENTER='0'" + "&"
    + "EPHEM_TYPE='ELEMENTS'" + "&"
    + "START_TIME='JD"+ str(ephemeris_start) + "'" + "&" #str(launch_date)
    + "STOP_TIME='JD" + str(ephemeris_stop) + "'" + "&"   # launch
    + "STEP_SIZE='1 d'"
)




########################################################################################################################
# Poliastro Orbit calculations
########################################################################################################################
mu_sun = 1.32712440042e11
mu_earth = 3.986e5
epochs = time_range(start = ti,end = tf)
# Orbital Elements lists
if flag != 1:
    elements_earth = get_elements(url_earth)
    elements_asteroid = get_elements(url_2008_ev5)
    data = {
        "Earth": elements_earth.tolist(),
        "Asteroid": elements_asteroid.tolist()
    }
    with open(file, 'w') as file:
        file.write(json.dumps(data, indent=4))
        elements_earth_radians = elements_earth.copy()
        elements_asteroid_radians = elements_asteroid.copy()
        elements_earth_radians[2:6] = np.radians(elements_earth[2:6])
        elements_asteroid_radians[2:6] = np.radians(elements_asteroid[2:6])

else:
    with open(file, 'r') as file:
        data = json.load(file)
        elements_earth = np.array(data['Earth'])
        elements_asteroid = np.array(data['Asteroid'])
        elements_asteroid_radians = elements_asteroid.copy()
        elements_earth_radians = elements_earth.copy()
        elements_earth_radians[2:6] = np.radians(elements_earth[2:6])
        elements_asteroid_radians[2:6] = np.radians(elements_asteroid[2:6])


# Define orbits
earth = define_orbit(Sun, elements_earth, ti)
asteroid = define_orbit(Sun, elements_asteroid, ti)

# Define The orbit that will be seen in the plot after the elapsed ephem time chosen above
earth_ephem = Ephem.from_orbit(orbit = earth, epochs = epochs)
asteroid_ephem = Ephem.from_orbit(orbit = asteroid, epochs = epochs)




########################################################################################################################
# Maneuvering
########################################################################################################################
# Determine line of nodes between the two inclined orbits for Bi-Elliptic Transfer calculations
r_earth, v_earth = earth.rv()
r_earth, v_earth = r_earth.value, v_earth.value
r_asteroid, v_asteroid = asteroid.rv()
r_asteroid, v_asteroid = r_asteroid.value, v_asteroid.value
h_earth = angular_momentum(r_earth, v_earth)
h_asteroid = angular_momentum(r_asteroid, v_asteroid)
n = node_line(h_earth, h_asteroid)
e_earth = eccentricity_vector(mu_sun, r_earth, v_earth)
w_earth = argument_perigee(n,e_earth)
dt_ascending, dt_descending = time_to_node_line(elements_earth[0], elements_earth[1], w_earth, elements_earth_radians[5], mu_sun) #elements_earth_radians[5]
print(f"Time to ascending node: {dt_ascending/86400} days")
print(f"Time to descending node: {dt_descending/86400} days")
dt = dt_ascending/86400
transfer_date = ti+dt
dt /= dt * u.day
ttf = Time(transfer_date, format = 'jd')
print(f"Julian transfer date: {ttf.jd}")
print(f"Date of transfer @ RAAN: {julian_to_current(transfer_date.jd)}")
# Fast-forward to maneuver date
earth = earth.propagate(ttf)
asteroid = asteroid.propagate(ttf)
r_earth, v_earth = earth.rv()
r_asteroid, v_asteroid = asteroid.rv()
v = normalize(v_earth.value)
print("v=", v)
dv_tnz = [5, 0.0, 0.0] * u.km/u.s # to calculate! tangential normal -nadir thrust component then rotate it to the pqw
dv1 = 5*v * u.km/u.s
dv2 = (0.75*v) * u.km/u.s + [0, 0, 1.75] * u.km/u.s
impulse1 = Maneuver.impulse(dv1)
transfer_1 = earth.apply_maneuver(impulse1)
a = transfer_1.a
print(a.value)
dt2 = np.sqrt(a.value**3/mu_sun)*np.pi/2
transfer_date += dt2
ttf = Time(transfer_date, format = 'jd')
transfer_1_end = transfer_1.propagate_to_anomaly(180 * u.deg)

impulse2 = Maneuver.impulse(dv2)
transfer_2 = transfer_1_end.apply_maneuver(impulse2)





phase_angle = np.pi*(1-np.sqrt(((1+r_earth/r_asteroid)/2)**3))
########################################################################################################################
# Plotting
########################################################################################################################

fast_forward = 1 # See fast forward in time by latter specified amount of years
plotter = OrbitPlotter3D()
plotter.set_attractor(Sun)




plotter.plot_ephem(
    earth_ephem,
    ti,
    label = "Earth @ Launch Position",
    color = "blue"
)
plotter.plot_ephem(
    asteroid_ephem,
    ti,
    label = "Asteroid @ Launch Position",
    color = "orange"
)
plotter.plot_maneuver(
    initial_orbit = earth,
    maneuver = impulse1,
    color = "purple",
    label = "Transfer Orbit 1"
)
plotter.plot_maneuver(
    initial_orbit = transfer_1_end,
    maneuver = impulse2,
    color = "green",
    label = "Transfer Orbit 2"
)

plotter.plot_trajectory(
    transfer_1.sample(max_anomaly=180 * u.deg),
    color = 'purple',
    label = 'Transfer Orbit 1'
)
plotter.plot_trajectory(
    transfer_2.sample(max_anomaly=360 * u.deg),
    color = 'green',
    label = 'Transfer Orbit 2'
)



if fast_forward !=1:
    plotter.plot(
        asteroid_after,
        label=f"asteroid after {years} years",
        color="orange"
    )
    plotter.plot(
        earth_after,
        label=f"earth after {years} years",
        color="blue"
    )


#plotter.plot(orb_eagle, label = 'Eagle', color = 'red')
fig = plotter._figure
fig.write_html("orbit.html")
webbrowser.open('file://' + os.path.realpath("orbit.html"))
